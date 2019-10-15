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
      tags$p(),
      wellPanel(
        selectInput("regionStart", "Region start",
                    c("set", "gene", "plate", "manual"),
                    "plate"),
        selectInput("regionEnd", "Region end",
                    c("set", "gene", "plate", "manual"),
                    "plate")
      ),
      width = 2
    ),
    mainPanel(
      fluidRow(
        column(6,
               selectInput("pointsSet", "Set", ""),
               checkboxInput("showInRange", "Show in range points only"),
               plotOutput("pointsPlot", click = "pointsPlot_click",
                          brush = "pointsPlot_brush",
                          hover = hoverOpts("pointsPlot_hover",
                                            delay = 100,
                                            delayType = "debounce")),
               uiOutput("pointsPlot_hover_info")
        ),
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
