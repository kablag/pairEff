library(shiny)
library(shinyWidgets)
library(pairEff)
library(tidyverse)
library(RDML)
# library(data.table)
# library(plotly)


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
        datapath = paste0(path.package("pairEff"),
                          "/extdata/", "Gapdh-Got1-Gtf2-Gusb -  rfu.xls")
      ))
    }
  )
  observeEvent(
    input$inputFile,
    {
      inputfile(input$inputFile)
    }
  )

  clickY <- reactiveVal()

  observeEvent(input$pointsPlot_brush$ymin,
               if (length(input$pointsPlot_brush$ymin) &&
                   length(input$pointsPlot_brush$ymax))
                 clickY(c(input$pointsPlot_brush$ymin,
                          input$pointsPlot_brush$ymax)))

  pairEffOptions <- reactive({
    list(regionStart = {
      if (input$regionStart != "manual")
        input$regionStart
      else
        clickY()[1]
    },
    regionEnd = {
      if (input$regionEnd != "manual")
        input$regionEnd
      else
        clickY()[2]
    })
  })

  rdml <- reactive({
    req(inputfile())
    withProgress(message = 'Loading file...', value = 0, {
      rdml <- {
        if (tools::file_ext(inputfile()$datapath) %in% c("xls", "xlsx")) {
          # read table
          inptw <- data.frame(readxl::read_excel(inputfile()$datapath))
          # convert to long format
          wells <- colnames(inptw)[-1]
          description <- data.frame(
            fdata.name = wells,
            quantity = rep(c(100, 50, 25, 12, 6, 3), 16),
            exp.id = c("exp1"),
            run.id = c("run1"),
            react.id = 1:96,
            sample = wells,
            sample.type = c("std"),
            target = paste(
              as.vector(sapply(1:16, function(i)
                rep(paste0(wells[i * 6 - 5], "-", wells[i * 6]), 6))),
              c(sapply(1:4, function(i) rep(paste0("Gene", i), 12)),
                sapply(1:4, function(i) rep(paste0("Gene", i), 12))),
              sep = "@"),
            target.dyeId = c("Sybr"),
            stringsAsFactors = FALSE
          )
          rdml <- RDML$new()
          rdml$SetFData(inptw, description)
          rdml
        } else {
          RDML$new(inputfile()$datapath)
        }
      }
    })
  })

  pEff <- reactive({
    req(rdml())
    withProgress(message = 'Calculating pairwise efficiency...', value = 0, {
      pEff <- pairEff(rdml(),
                      regionStart = pairEffOptions()$regionStart,
                      regionEnd = pairEffOptions()$regionEnd)
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

    inptl <- reactive({
      req(pEff())
      inptl <- pEff()$input
      if (input$showInRange)
        inptl <- inptl[RFU <= regionEnd & RFU >= regionStart]
      inptl[set == input$pointsSet]
    })

    output$pointsPlot <- renderPlot({
      req(inptl(), nrow(inptl()) != 0)
      p <- ggplot(inptl()) +
        geom_line(aes(x = cyc, y = fluor, color = as.factor(conc), group = position)) +
        geom_point(aes(x = takeoffC, y = takeoffF, color = as.factor(conc))) +
        geom_point(aes(x = derMaxC, y = derMaxF, color = as.factor(conc))) +
        geom_hline(aes(yintercept = regionStart)) +
        geom_hline(aes(yintercept = regionEnd)) +
        guides(color = guide_legend(title = "Conc."))
      p +
        scale_y_continuous(
          breaks = sort(c(ggplot_build(p)$layout$panel_params[[1]]$y.major_source,
                          inptl()$regionStart[1], inptl()$regionEnd[1])))
    })

    output$pointsPlot_hover_info <- renderUI({
      req(input$pointsPlot_hover)
      tags$p(HTML(sprintf("RFU = %f", input$pointsPlot_hover$y)))
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
