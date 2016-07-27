
shinyUI(fluidPage(

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
))


# library(rsconnect)
# rsconnect::setAccountInfo(name='projectb', token='6BA49A929820DC78E1F441A22338BBC8', secret='QZhc1xDWHOUt+23y8otnMbQFDJXdgudyJXLNRYrS')
#deployApp()
