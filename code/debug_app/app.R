###
# shiny application to run sequential designs.
# choice sets are generated adaptively sequentially using the KL criterion
# Results are saved on dropbox account.
###

rm(list=ls())
options(shiny.reactlog=TRUE)
#setwd
#setwd("C:/Users/u0105757/Desktop/thesis_Gil/survey_app")

#libraries
library(choice)
library(shinyjs)
library(mgcv)
library(plyr)
library(DT)
library(xlsx)
library(rdrop2)


### settings DCE ######################################################################
answer_options<-c("Nog niets geselecteerd", "Bank A", "Bank B")
n_alts=2
levels<-c(3,2,3,2,2,2)
mindiff<- 4

#dropbox map to store data
outputDir <- "responses"

#level names
names<- vector(mode="list", length(levels))
names[[1]]<-c("0%","1%","2%")
names[[2]]<-c("spaarrekening","spaar- en zichtrekening")
names[[3]]<-c("bankautomaten","bankautomaten + PC-banking","bankautomaten + PC-banking + banking-app")
names[[4]]<-c("ja", "nee")
names[[5]]<-c("ja", "nee")
names[[6]]<-c("ja", "nee")

#alternatives names
alternatives<-c("BANK A","BANK B")

#attribute names
attributes<-c("Jaarlijkse rente", "Aanbod producten", "Betalingswijzen", "Klimaatschadelijke investeringen",
              "Investeringen in de productie van oorlogswapens", "Investeringen in mensenrechten schendende praktijken")

#PRIOR
p_mod<- c( 0.11, 0, 0, 0.2, 0, 0.19, 0.29, 0.34)
var<-1
p_cov<- diag(rep(var,length(p_mod)))

#design
n_total<-10 #total

#introduction text
intro_text<- tags$div(

tags$h3("Hypothetische situatie:"),

tags$h5("Eerder in de vragenlijst hebt u ingevuld welke bank uw voorkeur geniet voor het
        uitvoeren van uw dagdagelijkse betalingen." ,tags$br(),
"Beeld u in dat u momenteel enkel bent aangesloten bij deze bank. Op dit ogenblik staat er 25.000 euro op uw spaarrekening bij deze bank. Op uw zichtrekening beschikt u gemiddeld over 1.000 euro.",tags$br(),tags$br(),
"Vandaag verneemt u echter het spijtige nieuws dat uw bank in de toekomst geen spaarrekening meer zal aanbieden. U bent genoodzaakt om elders een spaarrekening te openen."),

tags$h5("In wat volgt worden telkens twee nieuwe banken voorgesteld die op zes vlakken van elkaar verschillen."),

tags$h4(tags$b("Netto rendement van uw geld (1):")),h5("de rente die u jaarlijks zal overhouden op het overgezette geld. De jaarlijkse lidmaatschapskosten zijn hierin verrekend."),
tags$h4(tags$b("Aanbod financiele producten (2):")),h5("De bank biedt enkel een spaarrekening aan, of de bank biedt zowel een zicht- als een spaarrekening aan."),
tags$h4(tags$b("Betalingswijzen (3):")),
     tags$h5(tags$em("Bankautomaat:"),("u voert verrichtingen uit via de bankautomaten die u op straat kan terugvinden.")),
     tags$h5(tags$em("PC-banking:"),("u kan online inloggen op een site van uw bank om zo uw verrichtingen uit te voeren.")),
     tags$h5(tags$em("Banking-app:"),("u kan verrichtingen uitvoeren op uw smartphone, tablet of computer via een applicatie (app) die u kan downloaden.")),

tags$h4(tags$b("De bank doet al dan niet investeringen in:")),
tags$h4(tags$b("- activiteiten die een negatieve impact hebben op de opwarming van de aarde (4);")),
     tags$h4(tags$b("- de productie van wapens die gebruikt worden in oorlogsgebieden (5);")),
     tags$h4(tags$b("- bedrijven die hun personeel laten werken in mensonwaardige en gevaarlijke omstandigheden (6).")),


tags$h5("Een bank verleent kredieten aan ondernemingen en overheden en behaalt hier rendementen op.
        In dit experiment wordt verondersteld dat banken zulke kredieten verlenen met al het zicht- en spaarrekeninggeld
        van hun klanten. Zo een krediet kan verleend worden aan eender welke organisatie.
        Als een bank een krediet verleent aan een omstreden organisatie, gebeurt dat in dit experiment aan 1 of
        meerdere van de drie voorgaande activiteiten. Als ze dat doen, doen ze dat telkens voor 10 procent van al het
        zicht- en spaarrekeninggeld van hun klanten.")
)



#end text
end_txt<-c("Druk 1 maal op opslaan en afsluiten. Indien het venster niet automatisch sluit mag je het zelf afsluiten.")
end_txt2<-c("Bedankt voor uw deelname!")
###############################################################################################

#### settings questionnaire ###################################################################
Qlist <- read.xlsx("Qlist.xlsx", 1)
Qlist <- Qlist[-14,]

Qlist2 <- read.xlsx("Qlist2.xlsx", 1)
Qlist2 <- Qlist2[,-10]

intro1_text<-tags$div("Bedankt dat u even tijd wil vrijmaken voor dit onderzoek.",tags$br(),tags$br(),
"Mijn naam is Gil Op de Beeck, student Handelsingenieur aan de KU Leuven.
In het kader van mijn masterproef doe ik onderzoek naar de bankkeuze van de Belgen.
Ik zal me hiervoor baseren op de resultaten van deze vragenlijst. Uw antwoorden zijn dus een heuse bijdrage aan mijn onderzoek.
",tags$br(),tags$br(),"In wat volgt worden u enkele vragen voorgelegd en zal u deelnemen aan een klein keuze-experiment. De vragenlijst zal maximaal vijftien minuten van uw tijd in beslag nemen.
              U kan geen verkeerde antwoorden geven en uw gegevens worden strikt anoniem verwerkt.", tags$br(), tags$br(),
"Alvast hartelijk bedankt!",tags$br(),"Gil")


###############################################################################################


jscode <- "shinyjs.closeWindow = function() { window.close(); }"

#fcomb<-full_sets(lvls = levels, n_alts = n_alts, mindiff = mindiff)
#write.xlsx(x = fcomb, file = "fcomb.xlsx",sheetName = "fcomb", row.names = FALSE)

fcomb<-read.xlsx(file = "fcomb.xlsx", 1)

### user interface
ui <- fluidPage(

 # tags$style(type="text/css", "#opttext5 { height: 400px; width: 100%; text-align:center; font-size: 30px; display: block;}"),
  useShinyjs(),
  extendShinyjs(text = jscode, functions = c("closeWindow")),
  column(width=12,

         mainPanel(align = "left",

                   uiOutput("setnr"),
                   DT::dataTableOutput('choice.set'),
                   uiOutput("MainAction"),
                   actionButton("OK", label = "OK", icon("paper-plane" )),
                   br(),br(),
                   uiOutput("sideaction")




         ))
)


### server
server<-function(input, output) {

  #initialize
  n_att<-length(levels)
  fp<-profiles(lvls = levels)
  drops <- c("set","alt")
  f <- n <- 0

  makeReactiveBinding('n')
  makeReactiveBinding('f')


  reac_vars = reactiveValues(des = matrix(),
                             y_bin=vector("numeric"),
                             sam=matrix(data = NA, ncol = length(p_mod)),
                             w= vector("numeric"),
                             resp=vector("character"),
                             choice.set=matrix(),
                             resQ=vector("character"))

  #dynamic userinterface
  output$MainAction <- renderUI({ dynamicUi() })
  dynamicUi <- reactive({

    #1 Questionnaire explanation
    if (input$OK == 0 )
      return(list(h3(intro1_text, align = "left")))

    # Questionnaire questions
    if (input$OK > 0  & input$OK < nrow(Qlist)+1)
      return(list(
        h5(textOutput("questnr")),
          h3(textOutput("question")),
          uiOutput("answeroptions"),
          uiOutput("optionaltext")
        ))

    #2 explanation DCE
    if (input$OK == nrow(Qlist)+1 )
      return(list(h3(intro_text)))

    #3 DCE
    if (input$OK > nrow(Qlist)+1  & input$OK <= nrow(Qlist)+1+n_total)
      return(
        uiOutput("radio", align="center")
       )

    #4 Questionnaire 2
    if (input$OK > nrow(Qlist)+1+n_total & input$OK <= nrow(Qlist)+1+n_total + nrow(Qlist2))
      return(list(
        h5(textOutput("questnr2")),
        h3(textOutput("question2")),
        uiOutput("answeroptions2"),
        uiOutput("optionaltext4")

      ))

     #5 thanks + SAVE
    if (input$OK > nrow(Qlist) + 1 + n_total + nrow(Qlist2))
       return(list( h4(end_txt), actionButton('Save', 'Opslaan en afsluiten'), br(), h4(end_txt2)))

  })

  output$sideaction <-renderUI({ sidefun() })

  sidefun<-reactive({

    if (input$OK > nrow(Qlist)+1  & input$OK <= nrow(Qlist)+1+n_total)
      return(list(
    uiOutput("info"),
    uiOutput("help")
      ))

  })
  #radiobuttons DCE
  output$radio<-renderUI({radioButtons("survey", "Welk van deze twee banken zou u verkiezen:",
               answer_options , inline = T)})

  #radiobuttons help file
  output$info<-renderUI({radioButtons("checkinfo", "Beschrijving hypothetische situatie:",
                                       c("laat zien","verberg") , inline = T, selected = "verberg")})

  #show text if button is selected
  output$help<-renderText({
    if (input$checkinfo == "laat zien"){ paste(intro_text) }

     })


  #hide actionbutton end of survey
  observeEvent(input$OK, {
    if(n == nrow(Qlist)+1+n_total+nrow(Qlist2)){ hide('OK') }
     n <<- n + 1
  })

  #hide datatable
  observeEvent(input$OK, {
    if(f > nrow(Qlist)+n_total){ hide('choice.set') }
    f <<- f + 1
  })

  #close app
  observeEvent(input$Save, {
    saveData(d = reac_vars$des, Q= reac_vars$resQ, Y= reac_vars$y_bin )
    js$closeWindow()
    })

  output$vragen<-renderText({reac_vars$resQ})

  #store responses DCE and questionnaire
  observeEvent(input$OK, {

    ###Questionnaire 1###
    #R
    if (input$OK %in% c(2,10,13)){
      reac_vars$resQ<-c(reac_vars$resQ, input$R_Questions)
    }
    #T
    if (input$OK %in% c(3,5)){
      reac_vars$resQ<-c(reac_vars$resQ, input$T_Questions)
    }
    #RT
    if (input$OK %in% c(4,6:8,11,12) ){
      reac_vars$resQ<-c(reac_vars$resQ, paste(c(paste(input$RT_Questions),paste(input$opttext, collapse = " ")),collapse = " "))
    }
    #C
    if (input$OK %in% c(9,14)){
      reac_vars$resQ<-c(reac_vars$resQ, paste(c(paste(input$C_Questions, collapse = " "),paste(input$opttext2, collapse = " ")),collapse = " "))
    }

    ###Questionnaire 2###
    #R
    if (input$OK %in% c(nrow(Qlist)+n_total+4, nrow(Qlist)+n_total+6)){
      reac_vars$resQ<-c(reac_vars$resQ, input$R_Questions2)
    }
    #T
    if (input$OK %in% c(nrow(Qlist)+n_total+3, nrow(Qlist)+n_total+3)){
      reac_vars$resQ<-c(reac_vars$resQ, input$foo)
    }
    #RT
    if (input$OK == (nrow(Qlist)+n_total+5)){
      reac_vars$resQ<-c(reac_vars$resQ, paste(c(paste(input$RT_Questions2, collapse = " "),paste(input$opttext4, collapse = " ")),collapse = " ") )
    }
    #RTT
    if (input$OK == (nrow(Qlist)+n_total+7)){
      reac_vars$resQ<-c(reac_vars$resQ, paste(c(paste(input$RTT_quesions2, collapse = " "),paste(input$foo, collapse = " ")),collapse = " "))
    }
})
       ###DCE###

# observe({
#   if (input$OK > nrow(Qlist)+2 & input$OK <= nrow(Qlist)+2+n_total){
#
#       observeEvent(input$OK,{
#       reac_vars$resp<-isolate(c(reac_vars$resp, input$survey))
#       print(reac_vars$resp)
#       reac_vars$y_bin<-map_resp(resp = reac_vars$resp, resp_options = answer_options, n_alts = n_alts, neutral = T)
#
#       })
#   }
# })

observeEvent(input$OK,{
  if (input$OK > nrow(Qlist)+2 & input$OK <= nrow(Qlist)+2+n_total){
  reac_vars$resp<-isolate(c(reac_vars$resp, input$survey))
  reac_vars$y_bin<-map_resp(resp = reac_vars$resp, resp_options = answer_options, n_alts = n_alts, neutral = T)
  }
})


  # The option list
  option.list <- reactive({
    qlist <- Qlist[input$OK, 4:ncol(Qlist)]
    # Remove items from the qlist if the option is empty.
    as.matrix(na.omit(qlist[qlist!=""]))
  })

  # The option list2
  option.list2 <- reactive({
    qlist2 <- Qlist2[input$OK-nrow(Qlist)-1-n_total, 4:ncol(Qlist2)]
    # Remove items from the qlist if the option is empty.
    as.matrix(na.omit(qlist2[qlist2!=""]))
  })

  # This function show the question  text.
  output$question <- renderText({
    paste0( Qlist[input$OK, 3])
  })

  # This function show the question  text 2.
  output$question2 <- renderText({
    paste0( Qlist2[input$OK-nrow(Qlist)-1-n_total, 3])
  })

  # Return answer options conditional
  output$answeroptions<- renderUI({

    if (input$OK > 0  & input$OK < (nrow(Qlist)+1)){

      if (Qlist[input$OK, 2]=="R"){
      return(radioButtons("R_Questions",label=NULL, choices = option.list(), selected = character(0)))
      }

      if (Qlist[input$OK, 2]=="C"){
      return(checkboxGroupInput("C_Questions", label = NULL, choices = option.list()))
      }

      if (Qlist[input$OK, 2]=="RT"){
      return(radioButtons("RT_Questions",label=NULL, choices = option.list(), selected = character(0)))
      }

      if (Qlist[input$OK, 2]=="T"){
      return(textInput("T_Questions","",""))
      }
    }
  })

  # Return answer options conditionals 2
  output$answeroptions2<- renderUI({

    if (input$OK > nrow(Qlist)+1+n_total & input$OK <= nrow(Qlist)+1+n_total+nrow(Qlist2)){

      if (Qlist2[input$OK-nrow(Qlist)-1-n_total, 2]=="R"){
        return(radioButtons("R_Questions2",label=NULL, choices = option.list2(),inline = T, selected = character(0)))
      }
      if (Qlist2[input$OK-nrow(Qlist)-1-n_total, 2]=="RT"){
        return(radioButtons("RT_Questions2",label=NULL, choices = option.list2(),  selected = character(0)))
      }
      if (Qlist2[input$OK-nrow(Qlist)-1-n_total, 2]=="RTT"){
          return(radioButtons("RTT_Questions2",label=NULL, choices = option.list2(), selected = character(0)))
      }
    }
  })

  # Optional textboxes
  output$optionaltext<-renderUI({
    if (input$OK > 0  & input$OK < (nrow(Qlist)+1)){

    if (Qlist[input$OK, 2]=="RT"){
      return(textInput("opttext","",value=""))}
    if (Qlist[input$OK, 2]=="C"){
      return(textInput("opttext2","",value=""))}
    }
  })

  # Optional textboxes 2
  output$optionaltext4<-renderUI({
    if (input$OK > nrow(Qlist)+1+n_total & input$OK <= nrow(Qlist)+1+n_total + nrow(Qlist2)){

      if (Qlist2[input$OK-nrow(Qlist)-1-n_total, 2]=="RT"){
        return(textInput("opttext4","",value=""))}

    if (Qlist2[input$OK-nrow(Qlist)-1-n_total, 2]=="T"||Qlist2[input$OK-nrow(Qlist)-1-n_total, 2]=="RTT"){
      tags$textarea(id="foo", rows=6, cols=80, "")}
    }


  })

  # Observe conditions
  observe({
    toggle("opttext", condition= {input$RT_Questions == "Andere:" || input$RT_Questions == "Ik volg de studierichting:" || input$RT_Questions == "Ja, ik voel me aangesproken tot:"} )

  })
  observe({
    toggle("opttext2", condition= {"Andere:" %in% input$C_Questions} )

  })
  observe({
    toggle("opttext4", condition= {input$RT_Questions2 == "Ja en ik kan een voorbeeld geven:" || input$RT_Questions2 == "Ja" || input$RT_Questions2 == "Nee"  })

  })
  observe({
   if(input$OK > nrow(Qlist)+1  & input$OK <= nrow(Qlist)+1+n_total){
     toggle("OK", condition = {input$survey == "Bank A" || input$survey == "Bank B" })

   }

  })


  # Shows question nr
  output$questnr <- renderText({
    paste0("Vraag ", input$OK,"/",nrow(Qlist))
  })

  # Shows question nr 2
  output$questnr2 <- renderText({
    paste0("Vraag ", input$OK-nrow(Qlist)-1-n_total,"/",nrow(Qlist2))
  })


  ##### Produce choice set #####
  select_set <-eventReactive(input$OK, {

    #DCE fase
    if (input$OK > nrow(Qlist)+1  & input$OK <= nrow(Qlist)+1+n_total){

      #first set (no responses observed)
      if (input$OK < nrow(Qlist)+2+1){
        N<-1000
        reac_vars$sam<-MASS::mvrnorm(n=N, mu=p_mod, Sigma= p_cov)
        reac_vars$w<-rep(1,N)/N

        reac_vars$des<-KL_info(fp = fp, fcomb = fcomb, par_samples = reac_vars$sam, weights = reac_vars$w, n_alts = n_alts )

        reac_vars$choice.set<-t(reac_vars$choice.set[, 1:n_att])
      }

      else{
        #update parameters
        samples<-imp_sampling(prior_mode = p_mod, prior_covar = p_cov, design = reac_vars$des, n_alts = n_alts, Y=reac_vars$y_bin)
        reac_vars$sam<-samples[[1]]
        reac_vars$w<-samples[[2]]

        #new set
        new_set<-KL_info(fp = fp, fcomb = fcomb, par_samples = reac_vars$sam, weights = reac_vars$w, n_alts = n_alts )

        #update
        reac_vars$des<-rbind(reac_vars$des, new_set)

        reac_vars$choice.set<-present_Gil(design = new_set, levels=levels, lvl_names = names, n_alts = n_alts)
        reac_vars$choice.set<-t(reac_vars$choice.set[, 1:n_att])

      }

      #Fill in attribute names and alternatives names
      colnames(reac_vars$choice.set) <- alternatives
      rownames(reac_vars$choice.set) <- attributes

      return(reac_vars$choice.set)

    }

  })

  #plot choice set
  output$choice.set<- DT::renderDataTable({

    if (input$OK > nrow(Qlist)+1  & input$OK <= nrow(Qlist)+1+n_total){

      set<-select_set()
      datatable(set, filter = 'none', selection="multiple", escape=FALSE,
                options = list(sDom  = '<"top"><"bottom">'))
    }
    else {hide()}
  })

  #plot setnr DCE
  output$setnr <- renderText({

    if (input$OK > nrow(Qlist)+1  & input$OK <= nrow(Qlist)+1 +n_total){paste(h5("set: ", input$OK-(nrow(Qlist)+1) ,"/",n_total))}
  })

  #Save data Questionnaire + DCE
  saveData <- function(d,Y,Q) {

    #datafiles
    dataDCE<-cbind(d,Y)
    dataQ<-Q

    #unique code
    time<-as.integer(Sys.time())

    # Create a unique file name
    fileName_DCE <- sprintf("R_%s_%s_DCE.csv", n_alts, time)
    fileName_Q <- sprintf("R_%s_%s_Q.csv", n_alts, time)

    # Write the data to a temporary file locally
    filePath_DCE <- file.path(tempdir(), fileName_DCE)
    filePath_Q <- file.path(tempdir(), fileName_Q)

    # Write
    write.csv(dataDCE, file = filePath_DCE, row.names = FALSE)
    write.csv(dataQ, file = filePath_Q, row.names = FALSE)

    # Upload the file to Dropbox
    drop_upload(filePath_DCE, dest = outputDir)
    drop_upload(filePath_Q, dest = outputDir)


  }

}

shinyApp(ui, server)


#require('devtools')
#devtools::install_github('rstudio/shinyapps')
#install_github("traets/choice/choice")
#require('shinyapps')
#rsconnect::setAccountInfo(name='projectb', token='DC6FBFE968B39240BB8CC4C7BA8395F1', secret='Fo3QCehY7FrwoT5HrwUivjpPb146WitGAF0Gp8sq')

#deployApp()



