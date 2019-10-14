library(data.table)
library(readxl)
library(magrittr)
library(ggplot2)

pE <- function(dilutionRate,
               dilutionDelta,
               RFU.i, RFU.j,
               cycle.i, cycle.j) {
  (dilutionRate ^ ((log(RFU.j/RFU.i, base = dilutionRate) +
                      dilutionDelta) /
                     (cycle.j - cycle.i))) - 1
}

detectIndividualStartPoint <- function(fPoints, pval = 0.05, nsig = 3)
{
  cycles <- seq_along(fPoints)
  res <- sapply(5:length(cycles),
                function(i) {
                  mod <- lm(fPoints[1:i] ~ cycles[1:i], na.action = na.exclude)
                  1 - pt(tail(rstudent(mod), 1), df = mod$df.residual)
                })
  sig <- sapply(res, function(x) x < pval)
  sig[is.na(sig)] <- FALSE
  selTakeoff <- sapply(1:length(sig),
                       function(x) all(sig[x:(x + nsig - 1)]))
  minTakeoff <- min(which(selTakeoff == TRUE), na.rm = TRUE)
  takeoffC <- as.numeric(names(sig[minTakeoff]))
  takeoffF <- fPoints[takeoffC]
  return(list(takeoffC = takeoffC, takeoffF = takeoffF))
}
detectGlobalStartPoint <- function(RFU, takeoffF) {
  FBeforeTakeoffFs <- sort(RFU[RFU <= max(takeoffF)])
  meanFBeforeTakeoffFs <- mean(FBeforeTakeoffFs)
  tail(FBeforeTakeoffFs[FBeforeTakeoffFs <
                          (meanFBeforeTakeoffFs * 4)], 1) * 2
}

detectIndividualEndPoints <- function(cycles, fPoints) {
  der <- c(sapply(1:(length(cycles) - 1),
                  function(i) fPoints[i + 1] - fPoints[i]), NA)
  i <- which(der == max(der, na.rm = TRUE))
  list(derMaxC = cycles[i], derMaxF = fPoints[i])
}
detectGlobalEndPoint <- function(derMaxF) {
  mean(derMaxF)
}

#' Calculates pairwise efficiency
#'
#' Calculates pairwise efficiency
#'
#' @import magrittr data.table
#' @importFrom readxl read_excel
#' @importFrom moments kurtosis
#' @export
pairEff2 <- function(filename) {
  # read table
  inptw <- data.table(read_excel(filename))
  # convert to long format
  wells <- colnames(inptw)[-1]
  plateStr <- data.table(
    well = wells,
    conc = rep(c(100, 50, 25, 12, 6, 3), 16),
    set = as.vector(sapply(1:16, function(i)
      rep(paste0(wells[i * 6 - 5], "-", wells[i * 6]), 6))),
    gene = c(sapply(1:4, function(i) rep(paste0("Gene", i), 12)),
             sapply(1:4, function(i) rep(paste0("Gene", i), 12))),
    stringsAsFactors = FALSE
  )


  inptl <- melt(inptw, id.vars = "Cycle", variable.name = "well", value.name = "RFU")
  inptl <- inptl[plateStr, on = .(well)]
  inptl[, (c("takeoffC", "takeoffF")) := detectIndividualStartPoint(RFU),
        by = well]
  inptl[, minStartPointC := min(takeoffC, na.rm = TRUE),
        by = gene]
  inptl[, geneStartPoint := detectGlobalStartPoint(RFU[Cycle > 5], takeoffF),
        by = gene]
  inptl[, (c("derMaxC", "derMaxF")) := detectIndividualEndPoints(Cycle, RFU),
        by = well]
  inptl[, geneEndPoint := detectGlobalEndPoint(derMaxF),
        by = gene]

  inptlf <- inptl[Cycle > 5 &
                    Cycle >= minStartPointC &
                    RFU <= geneEndPoint & RFU >= geneStartPoint]

  gTbl <- rbindlist(
    lapply(unique(inptlf$set), function(setName) {
      do.call(cbind.data.frame, Map(expand.grid,
                                   i = inptlf[set == setName],
                                   j = inptlf[set == setName]))
    }
    )
  )
  gTbl <- data.table(gTbl)[Cycle.i < Cycle.j][
    , -c("set.j", "gene.j", "geneStartPoint.j", "geneEndPoint.j")
    ]
  colnames(gTbl)[colnames(gTbl) %in% c("set.i", "gene.i", "geneStartPoint.i", "geneEndPoint.i")] <-
    c("set", "gene", "geneStartPoint", "geneEndPoint")

  gTbl[, ':='(dilutionRate = ifelse(conc.i == conc.j, # one dilution
                                    2,
                                    ifelse(conc.i > conc.j,
                                           conc.i / conc.j,
                                           conc.j / conc.i)),
              dilutionDelta = ifelse(conc.i == conc.j,
                                     0,
                                     ifelse(conc.i > conc.j,
                                            1, -1)))]

  gTbl <- gTbl[, pE := pE(dilutionRate,
                          dilutionDelta,
                          RFU.i, RFU.j,
                          Cycle.i, Cycle.j)][
                            , ':='(mean_pE = mean(pE),
                                   totalN = .N), by = set
                            ][
                              #pE > 0 &
                              pE <= 2
                              ]

  gTbl[, pEgroup := round(pE / 0.05)]
  gTbl[, NpEgroup := .N, by = list(set, pEgroup)]
  gTbl <- gTbl[NpEgroup >= 5]
  gTbl[, ':='(mean_pE2 = mean(pE),
              sd_pE2 = sd(pE),
              includedN = .N), by = set]
  gTbl[, F0 := RFU.i / ((1 + pE) ^ Cycle.i)][
    , ':='(mean_F0 = mean(F0),
           sd_F0 = sd(F0),
           kurtosis = kurtosis(pE)),
    by = "set"
    ]
  # browser()
  gTbl[, errors := {
    sorted <- sort(c(mean_pE[1], mean_pE2[1]))
    paste0(
      {if (sorted[2]/sorted[1] > 1.15)
        "Possibly bad data: mean pE diff > 15%;"
        else ""
     }, {if (any(mean_F0 > geneStartPoint[1]))
       "Possibly bad data: F0 > working range start point;"
       else ""
     }, {if (kurtosis[1] < 1)
       "Possibly bad data: excess kurtosis < 1;"
       else ""
     })
  }, by = set]

  res <- unique(gTbl[, .(gene, set,
                         mean_pE, mean_pE2, sd_pE2,
                         mean_F0, sd_F0,
                         totalN, includedN,
                         kurtosis,
                         errors)],
                by = "set")
  list(
    input = as.data.table(inptl),
    pE = gTbl,
    result = res
  )
}

# pairEff2("/Users/kablag/Dropbox/r/pairEff/inst/extdata/Gapdh-Got1-Gtf2-Gusb -  rfu.xls")
# ggplot(gTbl) +
#   geom_point(aes(x = set, y = mean_pE), color = "red") +
#   geom_point(aes(x = set, y = mean_pE2), color = "green")
#
# ggplot(gTbl) +
#   geom_point(aes(x = set, y = mean_F0), color = "red")
#
# ggplot(inptlf) +
#   geom_line(aes(x = Cycle, y = RFU, color = set, group = well)) +
#   geom_point(aes(x = takeoffC, y = takeoffF, color = well)) +
#   geom_hline(aes(yintercept = geneStartPoint)) +
#   geom_hline(aes(yintercept = geneEndPoint)) +
#   facet_grid(. ~ gene) +
#   theme(legend.position = "none")
#
# gTbl[pE == -1 & set.i == "A1-A6"] %>%
#   ggplot() +
#   geom_line(aes(x = Cycle.i, y = RFU.i, color = well.i)) +
#   geom_line(aes(x = Cycle.j, y = RFU.j, color = well.j), linetype = 2)
