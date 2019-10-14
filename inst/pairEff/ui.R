library(shiny)
library(shinyWidgets)

ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification {
           height: 50px;
           width: 400px;
           position:fixed;
           top: calc(50% - 25px);;
           left: calc(50% - 200px);;
           }
           "
      )
    )
  ),

  titlePanel("pairEff: pairwise efficiency calculator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("inputFile", "xls"),
      actionButton("exmplFile", "Use example file"),
      width = 2
    ),
    mainPanel(
      fluidRow(
        column(6,
               selectInput("pointsSet", "Set", ""),
               checkboxInput("showInRange", "Show in range points only"),
               plotOutput("pointsPlot")),
        column(6,
               pickerInput("densitySets", "Sets", choices = "",
                           multiple = TRUE,
                           options = list("actions-box" = TRUE)),
               plotOutput("densityPlot"))
      ),
      tableOutput("resultsTbl")
    )
  )
)
