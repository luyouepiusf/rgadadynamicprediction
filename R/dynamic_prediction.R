dynamic_prediction<-function(
    list_yykj,
    list_ttkj,
    zzk,
    state_history_vec,
    landmark){
  
  state_history_mat<-apply(outer(1:nstates,state_history_vec,`==`),c(1,2),as.integer)
  is.na(state_history_mat)<-0
  
  at_risk_1<-case_when(state_history_vec==1~TRUE,TRUE~FALSE)
  at_risk_2<-case_when(state_history_vec==2~TRUE,TRUE~FALSE)
  at_risk_3<-case_when(state_history_vec==3~TRUE,TRUE~FALSE)
  time_12<-which(dplyr::lag(state_history_vec==1)&state_history_vec==2)[1]
  time_13<-which(dplyr::lag(state_history_vec==1)&state_history_vec==3)[1]
  time_23<-which(dplyr::lag(state_history_vec==2)&state_history_vec==3)[1]
  
  temp_mat<-solve(Sigmab)
  temp_mat<-0.5*(temp_mat+t(temp_mat))
  temp_vec<-rep(0,ndimb)
  for(kk in 1:ndimy){
    Phik<-as.matrix(F_basismatrix_long[[kk]][list_ttkj[[kk]],,drop=F])
    temp_mat[yidx2bidx[[kk]],yidx2bidx[[kk]]]<-
      temp_mat[yidx2bidx[[kk]],yidx2bidx[[kk]]]+
      (t(Phik)%*%Phik)/Sigmae2[kk]
    temp_vec[yidx2bidx[[kk]]]<-
      (t(Phik)%*%(list_yykj[[kk]]-Phik%*%mub[yidx2bidx[[kk]]]))/Sigmae2[kk]
  }
  temp_mat<-0.5*(temp_mat+t(temp_mat))
  Vbi_<-solve(temp_mat)
  Vbi_<-0.5*(Vbi_+t(Vbi_))
  Vbi_half_<-expm::sqrtm(Vbi_)
  Vbi_half_<-0.5*(Vbi_half_+t(Vbi_half_))
  Ebi_<-mub+Vbi_%*%temp_vec
  
  alleeieval<-t(randtoolbox::sobol(1000,dim=47,normal=T))
  allbbieval<-sweep(Vbi_half_%*%alleeieval,1,Ebi_,'+')
  
  all_loglik<-rep(NA,ncol(alleeieval))
  sum_exploglik<-0
  sum_ppt<-matrix(0,nstates,ntimepoints)
  for(elem in 1:ncol(alleeieval)){
    bbitemp<-allbbieval[,elem]
    yykt_temp<-matrix(0,ndimy,ntimepoints)
    zzkt_temp<-matrix(0,ndimz,ntimepoints)
    for(kk in 1:ndimy){
      yykt_temp[kk,]<-as.matrix(F_basismatrix_long[[kk]])%*%bbitemp[yidx2bidx[[kk]]]
    }
    for(kk in 1:ndimz){
      zzkt_temp[kk,]<-zzk[kk]
    }
    vvkt_temp<-rbind(yykt_temp,zzkt_temp)
    
    logh_term<-0
    if(!is.na(time_12))logh_term<-logh_term+
      log(arr_hazard0[1,2,time_12]+1e-5)+sum(vvkt_temp[,time_12]*arr_betav[1,2,])
    if(!is.na(time_13))logh_term<-logh_term+
      log(arr_hazard0[1,3,time_13]+1e-5)+sum(vvkt_temp[,time_13]*arr_betav[1,3,])
    if(!is.na(time_23))logh_term<-logh_term+
      log(arr_hazard0[2,3,time_23]+1e-5)+sum(vvkt_temp[,time_23]*arr_betav[2,3,])
    
    logF_term<-0
    for(idx_begin in 1:nstates){
      for(idx_end in 1:nstates){
        if(!possible_transition[idx_begin,idx_end])next
        logF_term<-logF_term-sum(arr_hazard0[idx_begin,idx_end,at_risk_1]*exp(c(t(vvkt_temp[,at_risk_1])%*%arr_betav[idx_begin,idx_end,])))
      }
    }

    loglik_temp<-logh_term+logF_term
    all_loglik[elem]<-loglik_temp
    sum_exploglik<-sum_exploglik+exp(loglik_temp)
    
    PPt_temp<-array(0,c(nstates,nstates,ntimepoints))
    for(idx_begin in 1:nstates){
      for(idx_end in 1:nstates){
        for(tt in 1:ntimepoints){
          if(!possible_transition[idx_begin,idx_end])next
          hexp_temp<-arr_hazard0[idx_begin,idx_end,tt]*
            exp(sum(arr_betav[idx_begin,idx_end,]*vvkt_temp[,tt]))
          PPt_temp[idx_begin,idx_end,tt]<-1-exp(-hexp_temp)
        }
      }
      PPt_temp[idx_begin,idx_begin,]<-1-colSums(PPt_temp[idx_begin,,])
    }
    ppt_temp<-matrix(0,nstates,ntimepoints)
    ppt_temp[,1:landmark]<-state_history_mat[,1:landmark]
    for(tt in (landmark+1):ntimepoints){
      ppt_temp[,tt]<-t(PPt_temp[,,tt])%*%ppt_temp[,tt-1]
    }
    sum_ppt<-sum_ppt+ppt_temp*exp(loglik_temp)
  }
  ppt<-sum_ppt/sum_exploglik
  wt<-exp(all_loglik)/sum_exploglik
  Ebbi<-c(allbbieval%*%wt)
  Vbbi<-allbbieval%*%(wt*t(allbbieval))
  
  list_result<-list(
    ppt=ppt,
    list_yykj=list_yykj,
    list_ttkj=list_ttkj,
    Ebbi=Ebbi,
    Vbbi=Vbbi,
    landmark=landmark)
  return(list_result)
}
