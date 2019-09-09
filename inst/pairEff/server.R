#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(pairEff)
library(tidyverse)
library(qpcR)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  pEff <- reactive({
    req(input$xlsFile)
    pEff <- pairEff(input$xlsFile$datapath, l6)
    updateSelectInput(session,
                      "includedPointsSet",
                      choices = pEff$result$set,
                      selected = pEff$result$set[1])
    updatePickerInput(session,
                      "densitySets",
                      choices = pEff$result$set,
                      selected = pEff$result$set[1:4])
    pEff
  })

  output$includedPoints <- renderPlot({
    req(pEff(), input$includedPointsSet)
    peff <<- pEff()
    if (input$showIncluded) {
      ggplot(pEff()$inptl %>%
               filter(set == input$includedPointsSet & usePoint) %>%
               mutate(conc = as.factor(conc))) +
      geom_line(aes(x = Cycle, y = RFU,
                    color = conc,
                    group = Well))
    } else {
      ggplot(pEff()$inptl %>%
               filter(set == input$includedPointsSet) %>%
               mutate(conc = as.factor(conc))) +
        geom_line(aes(x = Cycle, y = RFU,
                      color = conc,
                      group = Well,
                      size = usePoint))
    }

  })

  output$densityPlot <- renderPlot({
    req(pEff(), input$densitySets)
    pEff()$pEtbl %>%
      filter(set %in% input$densitySets & pE >= 0 & pE <= 2) %>%
      # rename(E = pE) %>%
      ggplot() +
      stat_density(aes(x = pE, group = set, fill = set),
                   position = "dodge", geom = "area",
                   alpha = 0.3, color = "black")

  })

  output$resultsTbl <- renderTable({
    req(pEff())
    pEff()$result
  }, digits = 8)
})
