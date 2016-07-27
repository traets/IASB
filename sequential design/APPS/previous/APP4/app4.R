
#####################################################################
#   survey APP 1.0                                                  #
#   Generate sequential design based on updated parameter values    #
#   (using package choice)                                          #
#####################################################################
rm(list=ls())

#Libraries
library(shiny)
library(shinyjs)
library(choice)
library(mgcv)
library(plyr)
library(DT)

# RESEARCHER #####################################
#levels attributes
levels<-c(2,3,2)
#attribute names
attributes<-c("price","travel time","comfort")
#level names
names<-list_example <- vector(mode="list", length(levels))
names[[1]]<-c("$50","$100")
names[[2]]<-c("2min","15min","30min")
names[[3]]<-c("bad","good")
#alternatives names
alternatives<-c("Alternative  A","Alternative B")
#answer options names
answer_options<-c("None","Alt A","Alt B")

#Number of sets
initial_sets<-10
total_sets<-15
#prior
p_mod<- c(1, -2, 3, 0.4)
p_cov<- diag(c(3, 3, 3, 3))

###############################################################

# PREPARATION #################################################
n_att<-length(levels)
n_alts<-length(alternatives)
#vector to store responses
resp<-character()
#generate initial design
des<-design(lvls=levels, n_sets = initial_sets, n_alts = n_alts)
#transform coded design into presentable design
des_ok<-present(design = des, lvl_names = names, n_alts = n_alts)
#drops
drops <- c("set","alt")
###############################################################

ui <- fluidPage(


  column(12,
         #  Application title
         headerPanel("Survey 1.0"),

         mainPanel(align = "center",

                   uiOutput("Question"),
                   dataTableOutput("choice.set"),
                   uiOutput("MainAction"),
                   actionButton("Click.Counter", "Next"),
                   uiOutput("Y")

         ))
)

#################################################################################################
#################################################################################################

server<-function(input, output) {


  # Main action
  output$MainAction <- renderUI({ dynamicUi() })
  dynamicUi <- reactive({

    #1 welcome message.
    if (input$Click.Counter==0)
      return( h3("Explanation of the survey") )

    #2 Survey
    if (input$Click.Counter>0 & input$Click.Counter<=total_sets)
      return(list(
        radioButtons("survey", "Please select the option you prefer:",
                     answer_options , inline = T, selected = "None")
      ))

    #3 Endnote
    if (input$Click.Counter>total_sets)
      return(list(
        h4("Thanks for taking the survey!"),
        actionButton('Save', 'Save Results on Dropbox'),
        br()))

    if (input$Click.Counter==total_sets)
      hide("Click.Counter")

  })

  # choice sets
  output$choice.set <- renderDataTable({


    if (input$Click.Counter>0 & input$Click.Counter <= initial_sets){

      #select choice set.
      choice.set<-des_ok[des_ok$set==input$Click.Counter, !(names(des_ok) %in% drops)]
      choice.set<-t(choice.set)

      #Fill in attribute names and alternatives names
      colnames(choice.set) <- alternatives
      rownames(choice.set) <- attributes

      datatable(choice.set, filter = 'none', selection="multiple", escape=FALSE,
                options = list(sDom  = '<"top"><"bottom">'))

    }

    # sequential
    else if (input$Click.Counter > initial_sets & input$Click.Counter <= total_sets){


      #update parameters
      draws<-imp_sampling(prior_mode = p_mod, prior_covar = p_cov, design = des, n_alts = n_alts, Y=input_resp())
      w<-draws[[2]]
      sam<-draws[[1]]

      #new set based on updated parameters
      new_set<-choice::DB_seq_fed(design = des, lvls = levels, n_alts = n_alts, par_samples = sam, weights = w, prior_covar = p_cov)

      #update
      des<<- rbind(des, new_set)
      p_mod<<-modes(samples = sam, weights = w, s=15)

      choice.set<-present(design = new_set, lvl_names = names, n_alts = n_alts)
      choice.set<-t(choice.set[, 1:n_att])

      #Fill in attribute names and alternatives names
      colnames(choice.set) <- alternatives
      rownames(choice.set) <- attributes


      datatable(choice.set, filter = 'none', selection="multiple", escape=FALSE,
                options = list(sDom  = '<"top"><"bottom">'))


    }

    else{}


  })

  output$Question <- renderText({

    if (input$Click.Counter>0 & input$Click.Counter<=total_sets){paste(h5("Question: ", input$Click.Counter))}
  })


  output$Y <- renderText({
    if (input$Click.Counter>0 & input$Click.Counter <= total_sets){


      input_resp <-eventReactive(input$OK, {

        resp<<-c(resp,input$survey)
        y_bin<<-map_resp(resp = resp, resp_options = answer_options, n_alts = n_alts)

      })





    }

  })
}

shinyApp(ui , server)




