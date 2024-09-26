plot_prediction<-function(list_result,heights=c(1.5,1)){
  
  ppt<-list_result$ppt
  list_yykj<-list_result$list_yykj
  list_ttkj<-list_result$list_ttkj
  Ebbi<-list_result$Ebbi
  Vbbi<-list_result$Vbbi
  landmark<-list_result$landmark
  
  grid.upper_lower_time<-function(upper,lower,time,color){
    upper0<-rep(upper,each=2)
    lower0<-rep(lower,each=2)
    time0<-c(0,rep(time[-length(time)],each=2),time[length(time)])
    grid.polygon(
      x=c(time0,rev(time0)),
      y=c(upper0,rev(lower0)),
      gp=gpar(fill=color),
      default.units="native")
  }
  
  plot_range<-12*11
  ppt_range<-ppt[,1:plot_range]
  horizon_lines<-seq(landmark,plot_range-1,12)
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow=2,heights=heights)))
  
  ################
  # longitudinal #
  ################
  
  pushViewport(viewport(layout.pos.row=1))
  varnamesv<-c("Fast. glu.\n(log scale)", "2-hour glu.\n(log scale)", "Fast. Cp.\n(log scale)", "2-hour Cp.\n(log scale)", 
               "HbA1c", "GADA-z\n(log scale)", "IA-2A-z\n(log scale)", "IAA-z\n(log scale)", "ZnT8A-z\n(log scale)", 
               "FDR", "Age", "DQ8X", "DQ2X")
  pushViewport(plotViewport(margins=c(0,8,0,2)+0.1))
  pushViewport(viewport(gp=gpar(cex=0.8),layout=grid.layout(nrow=ndimy)))
  for(kk in 1:ndimy){
    pushViewport(viewport(layout.pos.row=kk));grid.text(varnamesv[kk],x=unit(0,"npc")-unit(0.5,"line"),just="right",gp=gpar(lineheight=0.9))
    pushViewport(dataViewport(xscale=c(0,plot_range),yData=list_yaxis_range[[kk]],clip=T))
    grid.lines(x=unit(landmark,"native"))
    for(a_line in horizon_lines)grid.lines(x=unit(a_line,"native"),gp=gpar(lty=2,col="grey"))
    curve_mean<-c(F_basismatrix_long[[kk]]%*%Ebbi[yidx2bidx[[kk]]])
    curve_se<-sqrt(diag(F_basismatrix_long[[kk]]%*%Vbbi[yidx2bidx[[kk]],yidx2bidx[[kk]]]%*%t(F_basismatrix_long[[kk]])))
    grid.polygon(
      y=unit(c(curve_mean[1:plot_range]-curve_se[1:plot_range]*qnorm(0.975),
               rev(curve_mean[1:plot_range]+curve_se[1:plot_range]*qnorm(0.975)))*scaley[kk]+centery[kk],"native"),
      x=unit(c(1:plot_range,rev(1:plot_range)),"native"),gp=gpar(fill="grey90",col=NA))
    grid.lines(
      y=unit(curve_mean[1:plot_range]*scaley[kk]+centery[kk],"native"),
      x=unit(1:plot_range,"native"),gp=gpar(col="cornflowerblue",lwd=2))
    if(length(list_ttkj[[kk]]>0)){
      grid.points(x=list_ttkj[[kk]],y=list_yykj[[kk]]*scaley[kk]+centery[kk],default.units="native",pch=16,size=unit(2,"mm"))
    }
    grid.rect(gp=gpar(fill=NA))
    popViewport()
    pushViewport(dataViewport(xscale=c(0,plot_range),yData=list_yaxis_range[[kk]],clip=F))
    grid.yaxis(at=list_yaxis_at[[kk]],label=list_yaxis_label[[kk]],main=FALSE,gp=gpar(cex=0.75))
    # if(kk==ndimy){
    #   grid.xaxis(at=seq(0,plot_range,12),label=seq(0,plot_range,12)/12,gp=gpar(cex=0.75))
    #   grid.text("Years Since GADA positivity",y=unit(0,"npc")-unit(3,"line"),gp=gpar(cex=0.75))
    # }
    popViewport()
    popViewport()
  }
  popViewport()
  popViewport()
  popViewport()
  
  ############
  # survival #
  ############
  
  pushViewport(viewport(layout.pos.row=2))
  statenames<-c(
    "Single GADA",
    "Multiple without IA-2A",
    "Multiple with IA-2A",
    "Type 1 diabetes")
  ppt_upper<-apply(apply(ppt_range,2,dplyr::lag,default=0),2,cumsum)
  ppt_lower<-apply(ppt_range,2,cumsum)
  num_to_str<-function(x,digits=2)sprintf(paste0("%.",digits,"f"),x)
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  state_colors<-gg_color_hue(4)
  state_colors_darker<-c(
    colorRampPalette(c(state_colors[1],"black"))(100)[20],
    colorRampPalette(c(state_colors[2],"black"))(100)[20],
    colorRampPalette(c(state_colors[3],"black"))(100)[20],
    colorRampPalette(c(state_colors[4],"black"))(100)[20])
  
  pushViewport(plotViewport(margins=c(3,8,0.5,2)+0.1))
  pushViewport(viewport(layout=grid.layout(
    nrow=5,heights=unit(c(1,1.5,1.5,1.5,1.5),c("null","line","line","line","line")))))
  pushViewport(viewport(xscale=c(0,plot_range),yscale=c(1,0),layout.pos.row=1,layout.pos.col=1))
  for(idx_state in 1:4){
    grid.upper_lower_time(
      ppt_upper[idx_state,],
      ppt_lower[idx_state,],
      1:plot_range,
      state_colors[idx_state])
  }
  for(a_line in horizon_lines)grid.lines(x=unit(a_line,"native"),gp=gpar(lty=2,col="grey"))
  grid.text("Predicted probability\nof state occupation",x=unit(0,"npc")-unit(4,"line"),gp=gpar(cex=0.75),rot=90)
  grid.lines(x=unit(landmark,"native"))
  grid.yaxis(gp=gpar(cex=0.75))
  grid.yaxis(gp=gpar(cex=0.75),main=F)
  popViewport()
  for(idx_state in 1:nstates){
    pushViewport(viewport(xscale=c(0,plot_range),layout.pos.row=idx_state+1,layout.pos.col=1))
    grid.rect()
    for(a_line in horizon_lines)grid.lines(x=unit(a_line,"native"),gp=gpar(lty=2,col="grey"))
    grid.text(num_to_str(ppt_range[idx_state,horizon_lines],3),x=unit(horizon_lines,"native"),gp=gpar(col=state_colors_darker[idx_state],cex=0.75))
    grid.text(statenames[idx_state],x=unit(0,"npc")-unit(0.5,"line"),just="right",gp=gpar(cex=0.75,lineheight=0.9))
    grid.rect(
      x=unit(0,"npc")-unit(1,"line")-grobWidth(textGrob(statenames[idx_state],just="right",gp=gpar(cex=0.75))),
      width=unit(1,"line"),height=unit(1,"line"),gp=gpar(fill=state_colors[idx_state]))
    if(idx_state==nstates){
      grid.xaxis(at=seq(0,plot_range,12),label=seq(0,plot_range,12)/12,gp=gpar(cex=0.75))
      grid.text("Years Since GADA positivity",y=unit(0,"npc")-unit(3,"line"),gp=gpar(cex=0.75))
    }
    popViewport()
  }
  popViewport()
  popViewport()
}

# plot_prediction(list_result)


