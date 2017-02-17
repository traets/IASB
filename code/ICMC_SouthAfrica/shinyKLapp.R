###
# ICMC south africa. 
# shiny application to run sequential designs.
# choice sets are generated adaptively sequentially using the KL criterion
###

rm(list=ls())

#libraries
library(choice)
library(shinyjs) #quiting
library(mgcv) #uniquecombs function
library(plyr) #mapvalues function 
library(DT) #datatable
library(xlsx) #saving
library(rdrop2) #saving

### settings DCE ######################################################################
answer_options<-c("None","Alt A","Alt B")
n_alts=2
levels<-c(3,3,3)

#dropbox map to store data
outputDir <- "responses"

#level names
names<- vector(mode="list", length(levels))
names[[1]]<-c("$50","$75","$100")
names[[2]]<-c("2min","15min","30min")
names[[3]]<-c("bad","average","good")

#alternatives names
alternatives<-c("Alternative  A","Alternative B")

#attribute names
attributes<-c("price","travel time","comfort")

#PRIOR
p_mod<- c(-1, -1, -1, -1, 1, 1)
var<-3
p_cov<- diag(rep(var=3, length(p_mod)))

#design
n_total<-10 

#introduction text
intro_text<- c("Instructions")

#end text
end_txt<-c("Thanks for taking the survey!")

#######################################################################################

profs<-profiles(lvls = levels) #possible profiles
combs<- full_sets(lvls = levels, n_alts= n_alts, mindiff = 0) #combination of profiles
#jscode <- "shinyjs.closeWindow = function() { window.close(); }"

### user interface
ui <- fluidPage(
  
  #useShinyjs(),
  #extendShinyjs(text = jscode, functions = c("closeWindow")),
  column(12,
         
         #headerPanel("Survey KL"),
         
         mainPanel(align = "center",
                   
                   uiOutput("setnr"),
                   DT::dataTableOutput('choice.set'),
                   uiOutput("MainAction"),
                   actionButton("OK", "OK")
         ))
)


### server
server<-function(input, output) {
  
  #initialize
  n_att<-length(levels)
  drops <- c("set","alt")
  n <- 0
  makeReactiveBinding('n')
  
  reac_vars = reactiveValues(des = matrix(),
                             y_bin=vector("numeric"), 
                             sam=matrix(data = NA, ncol = length(p_mod)),
                             w= vector("numeric"),
                             resp=vector("character"),
                             choice.set=matrix())
  
  #dynamic userinterface
  output$MainAction <- renderUI({ dynamicUi() })
  dynamicUi <- reactive({
    
    #1 explanation.
    if (input$OK == 0 ) return(h3(intro_text))
    
    #2 DCE
    if (input$OK > 0 & input$OK <= n_total)
      return(list(
        radioButtons("survey", "Please select the option you prefer:",
                     answer_options , inline = T, selected = "None")))
    
    #3 Endnote
    if (input$OK > n_total)
      return(list( h4(end_txt), actionButton('Save', 'Save and Quit'), br()))
    
  })
  
  #actionbutton
  observeEvent(input$OK, { 
    if(n == n_total){ hide('OK') }
    n <<- n + 1
  })
  
  #close app
  observeEvent(input$Save, { 
    saveData(d = reac_vars$des, Y= reac_vars$y_bin )
    js$closeWindow()
  })
  
  #store responses
  observeEvent(input$OK, {
    if (input$OK > 1){
      reac_vars$resp<-c(reac_vars$resp, input$survey)
      reac_vars$y_bin<-map_resp(resp = reac_vars$resp, resp_options = answer_options, n_alts = n_alts)
    }
  })
  
  ###produce choice set###
  select_set <-eventReactive(input$OK, {
    
    #survey fase
    if (input$OK > 0 & input$OK <= n_total){
      
      #first set (no responses observed)
      if (input$OK < 2){ 
        N=1000
        reac_vars$sam<-MASS::mvrnorm(n=N, mu=p_mod, Sigma= p_cov)
        reac_vars$des<-KL_info(fp= profs, fcomb= combs, par_samples = reac_vars$sam, weights = rep(1,N)/N, n_alts = n_alts )
      
        reac_vars$choice.set<-present(design = reac_vars$des, lvl_names = names, n_alts = n_alts)
        reac_vars$choice.set<-t(reac_vars$choice.set[ , 1:n_att])
      }
      
      else{  
        #update parameters
        samples<-imp_sampling(prior_mode = p_mod, prior_covar = p_cov, design = reac_vars$des, n_alts = n_alts, Y=reac_vars$y_bin)
        reac_vars$sam<-samples[[1]]
        reac_vars$w<-samples[[2]]
        
        #new set 
        new_set<-KL_info(fp= profs, fcomb= combs, par_samples = reac_vars$sam, weights = reac_vars$w, n_alts = n_alts )
        
        #update
        reac_vars$des<-rbind(reac_vars$des, new_set)
        
        reac_vars$choice.set<-present(design = new_set, lvl_names = names, n_alts = n_alts)
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
    
    if (input$OK > 0 & input$OK <= n_total){
      
      set<-select_set()
      datatable(set, filter = 'none', selection="multiple", escape=FALSE,
                options = list(sDom  = '<"top"><"bottom">'))}
  })
  
  #plot setnr
  output$setnr <- renderText({
    
    if (input$OK > 0 & input$OK <= n_total){paste(h5("set: ", input$OK,"/",n_total))}
  })
  
  #Save data
  saveData <- function(d, Y) {
    
    data<-cbind(d, Y)
    # Create a unique file name
    fileName <- sprintf("%s_%s.csv", n_alts, input$ID)
    
    # Write the data to a temporary file locally
    filePath <- file.path(tempdir(), fileName)
    write.csv(data, filePath, row.names = FALSE, quote = TRUE)
    
    # Upload the file to Dropbox
    drop_upload(filePath, dest = outputDir)
  }
  
  
}

shinyApp(ui, server)

#require('devtools')
#require('shinyapps')
#deployApp()



