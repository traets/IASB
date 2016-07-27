#Example choice experiment


# Define UI for application that presents alternatives  
shinyUI(
  fluidPage(
  
    # Application title
    titlePanel("Questionnaire "),
    
    conditionalPanel
    (
    condition = "input.ok < 5",
    
    #Explanation
    h3("Choose your preferred option")
    
    ),
  
      #sidebar to define a choice
      sidebarLayout( 
        
         position="right",
    
    
             sidebarPanel(
                  
                 #Condition
               conditionalPanel
               (
               condition = "input.ok < 5",
                
                 br(),br(),br(),
                 checkboxGroupInput("dep.Var", label= h3("select option"),
                 choices = list ("option A"=1, "option B"=2, "option C"=3),selected=0),
                          
                 actionButton("ok", label = "confirm")
               ) 
                          ),
    
              mainPanel(
                
                #Eerste fase
                conditionalPanel
                (
                condition = "input.ok == 0",
                
                tags$pre(h4("Option A \t Option B  \t Option C")),
                
                img(src = "gras.A.jpg", height = 172, width = 160),
                img(src = "gras.B.jpg", height = 172, width = 160),
                img(src = "gras.C.jpg", height = 172, width = 172),
                
                tags$pre(h5("specifications \t         specifications \t specifications"))
                ),
   
                #Tweede fase
                conditionalPanel
                (
                condition = "input.ok == 1",
                
                tags$pre(h4("Option A \t Option B  \t Option C")),
                
                
                img(src = "gras.A.jpg", height = 172, width = 160),
                img(src = "gras.A.jpg", height = 172, width = 160),
                img(src = "gras.B.jpg", height = 172, width = 172),
                
                tags$pre(h5("specifications \t         specifications \t specifications"))
                ),
                
                #Derde fase
                conditionalPanel
                (
                condition = "input.ok == 2",
                
                tags$pre(h4("Option A \t Option B  \t Option C")),
                
                img(src = "gras.A.jpg", height = 172, width = 160),
                img(src = "gras.B.jpg", height = 172, width = 160),
                img(src = "gras.B.jpg", height = 172, width = 172),
                
                tags$pre(h5("specifications \t         specifications \t specifications"))
                ),
                
                #Vierde fase
                conditionalPanel
                (
                condition = "input.ok == 3",
                
                tags$pre(h4("Option A \t Option B  \t Option C")),
                
                img(src = "gras.C.jpg", height = 172, width = 160),
                img(src = "gras.B.jpg", height = 172, width = 160),
                img(src = "gras.C.jpg", height = 172, width = 172),
                
                tags$pre(h5("specifications \t         specifications \t specifications"))
                ),
                
                #Vijfde fase
                conditionalPanel
                (
                condition = "input.ok == 4",
                
                tags$pre(h4("Option A \t Option B  \t Option C")),
                
                img(src = "gras.A.jpg", height = 172, width = 160),
                img(src = "gras.C.jpg", height = 172, width = 160),
                img(src = "gras.A.jpg", height = 172, width = 172),
                
                tags$pre(h5("specifications \t         specifications \t specifications"))
                ),
                
                
                #Vijfde fase
                conditionalPanel
                (
                condition = "input.ok == 5",
                
                print(h2("THANKS!"))
                
                )
                
                
                
                
                
                
                
                      )
                    )
   
            ))