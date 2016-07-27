#Example choice experiment


# Define UI for application that presents alternatives  
shinyUI(
  fluidPage(
  
    # Application title
    titlePanel("Questionnaire "),
    
    #Explanation
    h3("Choose your preferred option"),
  
      #sidebar to define a choice
      sidebarLayout( 
        
         position="right",
    
    
             sidebarPanel(
      
                 br(),br(),br(),
                 checkboxGroupInput("dep.Var", label= h3("select option"),
                 choices = list ("option A"=1, "option B"=2, "option C"=3),selected=0),
                          
                 actionButton("ok", label = "confirm"),
                 textOutput("uitleg")
                          ),
    
              mainPanel(
    
              
                tags$pre(h4("Option A \t Option B  \t Option C")),
                
                img(src = "gras.A.jpg", height = 172, width = 160),
                img(src = "gras.B.jpg", height = 172, width = 160),
                img(src = "gras.C.jpg", height = 172, width = 172),
                
                tags$pre(h5("specifications \t         specifications \t specifications"))
  
   
                      )
                    )
   
            ))