#server

#Libraries
library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  
  output$uitleg<-renderText({  paste(input$ok) })
  
  
})
  
  
  
  
  
  
  
  
  
  
  
