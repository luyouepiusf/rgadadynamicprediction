gui_gadadynamicprediction<-function(){
  server<-function(input,output,session){

    df_hba1c = data.frame(x1=as.character(1:5),x2=rep(NA_real_,5))
    colnames(df_hba1c)<-c("Months","HbA1c\n(%)")

    df_autoab = data.frame(
      x1=as.character(1:5),
      x2=rep(NA_real_,5),
      x3=rep(NA_real_,5),
      x4=rep(NA_real_,5),
      x5=rep(NA_real_,5))
    colnames(df_autoab)<-c(
      "Months",
      "GADA",
      "IA-2A",
      "IAA",
      "ZnT8A")

    df_ogtt = data.frame(
      x1=as.character(1:5),
      x2=rep(NA_real_,5),
      x3=rep(NA_real_,5),
      x4=rep(NA_real_,5),
      x5=rep(NA_real_,5))
    colnames(df_ogtt)<-c(
      "Months",
      "Fasting\nglucose\n(mg/dL)",
      "2-hour\nglucose\n(mg/dL)",
      "Fasting\nC-peptide\n(ng/mL)",
      "2-hour\nC-peptide\n(ng/mL)")

    update_nrow<-function(df,nn){
      if(nn>nrow(df)){
        n_add<-nn-nrow(df)
        df<-df%>%bind_rows(data.frame(matrix(ncol=0,nrow=n_add)))
      }else{
        df<-df[1:nn,]
      }
      return(df)
    }
    get_state2<-function(df_autoab){
      df_select<-df_autoab%>%
        mutate(`Months`=as.numeric(`Months`))%>%
        filter(!is.na(`Months`))%>%
        arrange(`Months`)%>%
        filter(`IA-2A`>0|`IAA`>0|`ZnT8A`>0)
      return(ifelse(nrow(df_select)>=1,df_select$`Months`[1],NA))
    }
    get_state3<-function(df_autoab){
      df_select<-df_autoab%>%
        mutate(`Months`=as.numeric(`Months`))%>%
        filter(!is.na(`Months`))%>%
        arrange(`Months`)%>%
        filter(`IA-2A`>0)
      return(ifelse(nrow(df_select)>=1,df_select$`Months`[1],NA))
    }

    ## Defining a reactivevalues object so that whenever dataset value changes it affects everywhere in the scope of every reactive function
    datavalues<-reactiveValues(
      df_autoab=df_autoab,
      df_ogtt=df_ogtt,
      df_hba1c=df_hba1c,
      state2_time=NA_real_,
      state3_time=NA_real_)

    # returns rhandsontable type object - editable excel type grid data
    output$df_autoab <- renderRHandsontable({
      rhandsontable(datavalues$df_autoab,rowHeaderWidth=20)%>%
        hot_cols(colWidths=45)%>%
        hot_validate_numeric(cols=1,choices=as.character(1:500))
    })
    output$df_ogtt <- renderRHandsontable({
      rhandsontable(datavalues$df_ogtt,rowHeaderWidth=20)%>%
        hot_cols(colWidths=60)%>%
        hot_validate_numeric(cols=1,choices=as.character(1:500))
    })
    output$df_hba1c <- renderRHandsontable({
      #req(datavalues$df_hba1c)
      rhandsontable(datavalues$df_hba1c,rowHeaderWidth=20)%>%
        hot_cols(colWidths=45)%>%
        hot_validate_numeric(cols=1,choices=as.character(1:500))%>%
        hot_validate_numeric(cols=2,min=3,max=10)
    })

    observeEvent(
      input$nrow_autoab_bn,
      {
        datavalues$df_autoab<-update_nrow(datavalues$df_autoab,input$nrow_autoab)
        datavalues$state2_time<-get_state2(datavalues$df_autoab)
        datavalues$state3_time<-get_state3(datavalues$df_autoab)
      }
    )
    observeEvent(
      input$ex1_bn,
      {
        datavalues$df_ogtt<-
          tibble(
            `Months`=as.character(c(16,30)),
            `Fasting\nglucose\n(mg/dL)`=c(80,89),
            `2-hour\nglucose\n(mg/dL)`=c(60,115),
            `Fasting\nC-peptide\n(ng/mL)`=c(0.82,0.74),
            `2-hour\nC-peptide\n(ng/mL)`=c(NA_real_,NA_real_))
        datavalues$df_autoab<-
          tibble(
            `Months`=as.character(c(1,5,10,12,14,19,24,26,28,32,35)),
            "GADA"=c(0.2,2.2,1.2,1.1,1.6,3.4,4.8,4.2,2.6,2.4,1.2),
            "IA-2A"=c(-0.08,-0.07,NA,-0.08,-0.07,-0.07,-0.05,0.01,0.02,0.02,0.01),
            "IAA"=c(-0.33,-0.32,-0.33,-0.33,NA,-0.33,-0.32,-0.15,-0.25,-0.26,-0.26),
            "ZnT8A"=c(NA, NA, NA,-0.11,-0.10,-0.11,NA,-0.11,-0.11,-0.10,-0.09))
        datavalues$df_hba1c<-tibble(`Months`=32,`HbA1c\n(%)`=5.0)
        updateNumericInput(session,"nrow_autoab",value=nrow(datavalues$df_autoab))
        updateNumericInput(session,"nrow_ogtt",value=nrow(datavalues$df_ogtt))
        updateNumericInput(session,"nrow_hba1c",value=nrow(datavalues$df_hba1c))

        updateNumericInput(session,"Age",value=26)
        updateRadioButtons(session,"Sex",selected=1)
        updateRadioButtons(session,"FDR",selected=0)
        updateRadioButtons(session,"HLA",selected=0)

        updateRadioButtons(session,"event_history",selected=0)
        updateNumericInput(session,"state2_time",value=NA_real_)
        datavalues$state2_time<-NA_real_
        updateNumericInput(session,"state3_time",value=NA_real_)
        datavalues$state3_time<-NA_real_
        updateNumericInput(session,"landmark",value=36)
      }
    )
    observeEvent(
      input$ex2_bn,
      {
        datavalues$df_ogtt<-
          tibble(
            `Months`=as.character(c(16,30,45)),
            `Fasting\nglucose\n(mg/dL)`=c(80,89,89),
            `2-hour\nglucose\n(mg/dL)`=c(60,115,120),
            `Fasting\nC-peptide\n(ng/mL)`=c(0.82,0.74,0.76),
            `2-hour\nC-peptide\n(ng/mL)`=c(NA_real_,NA_real_,NA_real_))
        datavalues$df_autoab<-
          tibble(
            `Months`=as.character(c(1,5,10,12,14,19,24,26,28,32,35,39,43,45,48,51,55)),
            "GADA"=c(0.2,2.2,1.2,1.1,1.6,3.4,4.8,4.2,2.6,2.4,1.2,0.7,0.8,0.5,0.5,0.3,0.4),
            "IA-2A"=c(-0.08,-0.07,NA,-0.08,-0.07,-0.07,-0.05,0.01,0.02,0.02,0.01,NA,0.4,0.3,0.3,0.05,0.2),
            "IAA"=c(-0.33,-0.32,-0.33,-0.33,NA,-0.33,-0.32,-0.15,-0.25,-0.26,-0.26,-0.28,-0.26,-0.24,-0.24,-0.28,-0.25),
            "ZnT8A"=c(NA, NA, NA,-0.11,-0.10,-0.11,NA,-0.11,-0.11,-0.10,-0.09,0.02,NA,0.01,NA,0.00,-0.05))
        datavalues$df_hba1c<-tibble(`Months`=c(32,45),`HbA1c\n(%)`=c(5.0,5.1))
        updateNumericInput(session,"nrow_autoab",value=nrow(datavalues$df_autoab))
        updateNumericInput(session,"nrow_ogtt",value=nrow(datavalues$df_ogtt))
        updateNumericInput(session,"nrow_hba1c",value=nrow(datavalues$df_hba1c))

        updateNumericInput(session,"Age",value=26)
        updateRadioButtons(session,"Sex",selected=1)
        updateRadioButtons(session,"FDR",selected=0)
        updateRadioButtons(session,"HLA",selected=0)

        updateRadioButtons(session,"event_history",selected=0)
        updateNumericInput(session,"state2_time",value=NA_real_)
        datavalues$state2_time<-NA_real_
        updateNumericInput(session,"state3_time",value=43)
        datavalues$state3_time<-43
        updateNumericInput(session,"landmark",value=60)
      }
    )
    observeEvent(
      input$nrow_ogtt_bn,
      {
        datavalues$df_ogtt<-update_nrow(datavalues$df_ogtt,input$nrow_ogtt)
      }
    )
    observeEvent(
      input$nrow_hba1c_bn,
      {
        datavalues$df_hba1c<-update_nrow(datavalues$df_hba1c,input$nrow_hba1c)
      }
    )
    observeEvent(
      input$event_history,
      {
        output$test1 <- renderText(
          ifelse(
            is.na(datavalues$state2_time),
            "Still not multiple autoantibody positive",
            paste("Multiple autoantibody at ", datavalues$state2_time," Months")))
        output$test2 <- renderText(
          ifelse(
            is.na(datavalues$state2_time),
            "Still not IA-2A positive",
            paste("IA-2A positive at ", datavalues$state3_time," Months")))
      }
    )

    observeEvent(
      input$df_autoab$changes$changes,
      {
        datavalues$df_autoab<-hot_to_r(input$df_autoab) # convert the rhandontable to R data frame object so manupilation / calculations could be done
        updateNumericInput(session,"nrow_autoab",value=nrow(datavalues$df_autoab))
        datavalues$state2_time<-get_state2(datavalues$df_autoab)
        datavalues$state3_time<-get_state3(datavalues$df_autoab)
      }
    )
    observeEvent(
      input$df_ogtt$changes$changes,
      {
        datavalues$df_ogtt <- hot_to_r(input$df_ogtt) # convert the rhandontable to R data frame object so manupilation / calculations could be done
        updateNumericInput(session,"nrow_ogtt",value=nrow(datavalues$df_ogtt))
      }
    )
    observeEvent(
      input$df_hba1c$changes$changes,
      {
        datavalues$df_hba1c <- hot_to_r(input$df_hba1c) # convert the rhandontable to R data frame object so manupilation / calculations could be done
        updateNumericInput(session,"nrow_hba1c",value=nrow(datavalues$df_hba1c))
      }
    )

    observeEvent(
      input$predict_bn,
      {
        df_autoab_clean<-
          datavalues$df_autoab%>%
          mutate(`Months`=as.numeric(`Months`))%>%
          filter(!is.na(`Months`))%>%
          arrange(`Months`)%>%
          mutate(across(!`Months`,function(x)log(x+1)))
        df_ogtt_clean<-
          datavalues$df_ogtt%>%
          mutate(`Months`=as.numeric(`Months`))%>%
          filter(!is.na(`Months`))%>%
          arrange(`Months`)%>%
          mutate(across(!`Months`,function(x)log(x+1)))
        df_hba1c_clean<-
          datavalues$df_hba1c%>%
          mutate(`Months`=as.numeric(`Months`))%>%
          filter(!is.na(`Months`))%>%
          arrange(`Months`)

        df_logglufast_clean<-df_ogtt_clean%>%
          filter(!is.na(`Fasting\nglucose\n(mg/dL)`))%>%
          select(`Months`,`Fasting\nglucose\n(mg/dL)`)
        df_logglu120_clean<-df_ogtt_clean%>%
          filter(!is.na(`2-hour\nglucose\n(mg/dL)`))%>%
          select(`Months`,`2-hour\nglucose\n(mg/dL)`)
        df_logcpepfast_clean<-df_ogtt_clean%>%
          filter(!is.na(`Fasting\nC-peptide\n(ng/mL)`))%>%
          select(`Months`,`Fasting\nC-peptide\n(ng/mL)`)
        df_logcpep120_clean<-df_ogtt_clean%>%
          filter(!is.na(`2-hour\nC-peptide\n(ng/mL)`))%>%
          select(`Months`,`2-hour\nC-peptide\n(ng/mL)`)
        df_hba1c_clean<-df_hba1c_clean%>%
          filter(!is.na(`HbA1c\n(%)`))%>%
          select(`Months`,`HbA1c\n(%)`)
        df_loggadaz<-df_autoab_clean%>%
          filter(!is.na(`GADA`))%>%
          select(`Months`,`GADA`)
        df_logia2az<-df_autoab_clean%>%
          filter(!is.na(`IA-2A`))%>%
          select(`Months`,`IA-2A`)
        df_logmiaaz<-df_autoab_clean%>%
          filter(!is.na(`IAA`))%>%
          select(`Months`,`IAA`)
        df_logzincz<-df_autoab_clean%>%
          filter(!is.na(`ZnT8A`))%>%
          select(`Months`,`ZnT8A`)

        list_ttkj<-list(
          t1=df_logglufast_clean$`Months`,
          t2=df_logglu120_clean$`Months`,
          t3=df_logcpepfast_clean$`Months`,
          t4=df_logcpep120_clean$`Months`,
          t5=df_hba1c_clean$`Months`,
          t6=df_loggadaz$`Months`,
          t7=df_logia2az$`Months`,
          t8=df_logmiaaz$`Months`,
          t9=df_logzincz$`Months`)

        list_yykj<-list(
          y1=(df_logglufast_clean[[2]]-centery[1])/scaley[1],
          y2=(df_logglu120_clean[[2]]-centery[2])/scaley[2],
          y3=(df_logcpepfast_clean[[2]]-centery[3])/scaley[3],
          y4=(df_logcpep120_clean[[2]]-centery[4])/scaley[4],
          y5=(df_hba1c_clean[[2]]-centery[5])/scaley[5],
          y6=(df_loggadaz[[2]]-centery[6])/scaley[6],
          y7=(df_logia2az[[2]]-centery[7])/scaley[7],
          y8=(df_logmiaaz[[2]]-centery[8])/scaley[8],
          y9=(df_logzincz[[2]]-centery[9])/scaley[9])

        zzk<-(c(
          case_when(input$FDR=="1"~1,TRUE~0),
          input$Age/12,
          case_when(input$HLA=="1"~1,TRUE~0),
          case_when(input$HLA=="2"~1,TRUE~0))-centerz)/scalez

        landmark<-input$landmark

        state_history_vec<-rep(1,ntimepoints)
        if(!is.na(datavalues$state2_time))state_history_vec[datavalues$state2_time:ntimepoints]<-2
        if(!is.na(datavalues$state3_time))state_history_vec[datavalues$state3_time:ntimepoints]<-3
        state_history_vec[(landmark+1):ntimepoints]<-NA

        output$main_plot<-renderPlot(
          {
            list_result<-dynamic_prediction(
              list_yykj,list_ttkj,zzk,
              state_history_vec,landmark)
            plot_prediction(list_result)
          },width=900,height=1150)
      }
    )
  }

  ui<-function(request){
    page_fluid(
      tags$style(HTML("
    body {
      font-size: 11px; /* Adjust this value to your desired font size */
    }
    /* Adjust font size of inputs including numeric inputs */
    input, button, select, textarea {
      font-size: 11px; /* Match the global font size */
    }
    /* Specifically target numeric inputs */
    input[type='number'] {
      font-size: 11px; /* Match the global font size */
    }
    /* Adjust font size for action buttons */
    .btn {
      font-size: 11px; /* Match the global font size */
    }
    ")),
      h1("Dynamic Prediction Tool for GADA Positive Children"),
      layout_columns(
        actionButton("ex1_bn","Ex1"),
        actionButton("ex2_bn","Ex2"),
        actionButton("ex3_bn","Ex3")),
      layout_columns(
        card(
          card_header("Age"),
          numericInput("Age",label="Age in months",value=72,min=1,step=1)),
        card(
          card_header("Sex"),
          radioButtons(
            "Sex",
            "Sex",
            choices=list("Female"=0,"Male"=1),
            selected=0)),
        card(
          card_header("Family history"),
          radioButtons(
            "FDR",
            "First-degree relatives with type 1 diabetes?",
            choices=list("No"=0,"Yes"=1),
            selected=0)),
        card(
          card_header("Human leukocyte antigen (HLA)"),
          radioButtons(
            "HLA",
            "HLA genotype",
            choices=list("DQ2/8"=0,"DQ8/X, not DQ2/8"=1,"DQ2/X, not DQ2/8"=2),
            selected=0)),
        col_widths=c(3,3,3,3)),
      layout_columns(
        card(
          card_header("Autoantibody tests",height="30px"),
          a("Input Autoantibody test results, leave blank if a value is missing."),
          numericInput("nrow_autoab","Number of rows to input:",value=5,min=1,step=1),
          actionButton("nrow_autoab_bn","Update"),
          rHandsontableOutput("df_autoab")),
        card(
          card_header("OGTT variables",height="30px"),
          a("Input OGTT test results, leave blank if a value is missing."),
          numericInput("nrow_ogtt","Number of rows to input:",value=5,min=1,step=1),
          actionButton("nrow_ogtt_bn","Update"),
          rHandsontableOutput("df_ogtt")),
        card(
          card_header("HbA1c",height="30px"),
          a("Input HbA1c test results, leave blank if a value is missing."),
          numericInput("nrow_hba1c","Number of rows to input:",value=5,min=1,step=1),
          actionButton("nrow_hba1c_bn","Update"),
          rHandsontableOutput("df_hba1c")),
        col_widths=c(4,5,3)),
      layout_columns(
        card(
          card_header("Event History",height="30px"),
          a("Enter event history in months since GADA positivity"),
          radioButtons(
            "event_history",
            NULL,
            choices=list("Manual entry"=0,"Auto-fill using autoantibody results"=1),
            selected=0),
          conditionalPanel(
            condition="input.event_history==0",
            numericInput("state2_time","Time of multiple autoantibody (in months). Leave empty if single autoantibody positive.",value=NA,min=1,max=12*9,step=1),
            numericInput("state3_time","Time of IA-2A positivity in months. Leave empty if IA-2A negative.",value=NA,min=1,max=12*9,step=1)),
          conditionalPanel(
            condition="input.event_history==1",
            textOutput("test1"),
            textOutput("test2"))),
        card(
          card_header("Most recent follow-up",height="30px"),
          numericInput(
            "landmark",
            "The origin of prediction is the most recent follow-up time in months since GADA positivity, which is also called the landmark time. Please input your land mark time (in months).",value=NA)),
        col_widths=c(6,6),
        row_heights=c(1,1,1,1)),
      card(
        actionButton(
          "predict_bn","All done. Calculate Predictions!",
          style="background-color: #f5b1a6;height: 50px; font-size: 14px"),
        card_header("Predictions"),
        plotOutput(outputId = "main_plot", height = "1200px"),
      )
    )
  }

  runGadget(shinyApp(ui, server, enableBookmarking = "url"), viewer = dialogViewer(dialogName="Dynamic Prediction Tool for GADA Positive Children", width = 1400))
}
