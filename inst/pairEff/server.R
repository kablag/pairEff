library(shiny)
library(shinyWidgets)
library(pairEff)
library(tidyverse)
library(qpcR)

# source("../../R/pairEff2.R")

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })

  inputfile <- reactiveVal()

  observeEvent(
    input$exmplFile,
    {
      inputfile(list(
        name = "Gapdh-Got1-Gtf2-Gusb -  rfu.xls",
        datapath = paste0(path.package("pairEff"), "/extdata/", "Gapdh-Got1-Gtf2-Gusb -  rfu.xls")
      ))
    }
  )
  observeEvent(
    input$inputFile,
    {
      inputfile(input$inputFile)
    }
  )

  pEff <- reactive({
    req(inputfile())
    # pEff <- pairEff(input$xlsFile$datapath, l6)
    withProgress(message = 'Calculating pairwise efficiency...', value = 0, {
      pEff <- pairEff2(inputfile()$datapath)
    })
    updateSelectInput(session,
                      "pointsSet",
                      choices = pEff$result$set,
                      selected = pEff$result$set[1])
    updatePickerInput(session,
                      "densitySets",
                      choices = as.character(pEff$result$set),
                      selected = as.character(pEff$result$set[1:4]))
    pEff
  })

  output$pointsPlot <- renderPlot({
    req(pEff())
    inptl <- pEff()$input
    if (input$showInRange)
      inptl <- inptl[RFU <= geneEndPoint & RFU >= geneStartPoint]
    ggplot(inptl[set == input$pointsSet]) +
      geom_line(aes(x = Cycle, y = RFU, color = as.factor(conc), group = well)) +
      geom_point(aes(x = takeoffC, y = takeoffF, color = as.factor(conc))) +
      geom_point(aes(x = derMaxC, y = derMaxF, color = as.factor(conc))) +
      geom_hline(aes(yintercept = geneStartPoint)) +
      geom_hline(aes(yintercept = geneEndPoint)) +
      guides(color = guide_legend(title = "Conc."))
  })

  output$densityPlot <- renderPlot({
    req(pEff(), input$densitySets)
    pEff()$pE[set %in% input$densitySets] %>%
      ggplot() +
      stat_density(aes(x = pE, group = set, fill = set),
                   position = position_dodge(width = NULL),
                   geom = "area",
                   alpha = 0.3, color = "black")
  })

  output$resultsTbl <- renderTable({
    req(pEff())
    pEff()$result
  }, digits = 8)
})
