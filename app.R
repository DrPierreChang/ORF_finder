library(shiny)
library(shinythemes)
source("ORF_finder_genome.R")


# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "",
                  tabPanel("ORF FINDER",
                           sidebarPanel(
                             #tags$h3("Input:"),
                             tags$h3("Data"),
                             textInput("txt1", "Paste DNA seq:", ""),
                             textInput("txt2", "Enter minlen:", ""),
                             
                           ), # sidebarPanel
                           mainPanel(
                             h3("Results"),
                             
                             h4("Statistics:"),
                             dataTableOutput("table"),
                             
                             h4("ORFs:"),
                             verbatimTextOutput("table2"),
                             
                           ) # mainPanel
                           
                  ) # Navbar 1, tabPanel
                  
                  
                ) # navbarPage
) # fluidPage





# Define server function  
server <- function(input, output) {
  
  test <- reactive({
    d <- input$txt1
    ml <- as.integer(input$txt2)
    if (!is.na(d) && !is.na(ml)) {
      a <- orf_finder(d, ml)
      l <- list(df1 = a[2:6], df2 = a[c(1)])
      
    }
    #return(l)
  })
  
  
      output$table <- renderDataTable({
        data.frame(test()[[1]])
      })
      output$table2 <- renderPrint({
        print(test()[[2]])
      })
      
      
    }

  
  
 # server



# Create Shiny object
shinyApp(ui = ui, server = server)