#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)

# Define UI for application that draws a histogram
ui <- fluidPage(
  fileInput("xlsFile", "xls"),
  fluidRow(
    column(6,
           selectInput("includedPointsSet", "Set", ""),
           checkboxInput("showIncluded", "Show included points only"),
           plotOutput("includedPoints")),
    column(6,
           pickerInput("densitySets", "Sets", choices = "",
                       multiple = TRUE,
                       options = list("actions-box" = TRUE)),
           plotOutput("densityPlot"))
  ),
  tableOutput("resultsTbl")
)
