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
pairEff <- function(filename, regionStart = "gene", regionEnd = "gene") {
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
  inptl[, plate := "plate1"]
  if (is.numeric(regionStart)) {
    inptl[, (c("takeoffC", "takeoffF")) := list(1, regionStart)]
    inptl[, minStartPointC := 1]
    inptl[, regionStart := regionStart]
  } else {
    inptl[, (c("takeoffC", "takeoffF")) := detectIndividualStartPoint(RFU),
          by = well]
    inptl[, minStartPointC := min(takeoffC, na.rm = TRUE),
          by = regionStart]
    inptl[, regionStart := detectGlobalStartPoint(RFU[Cycle > 5],
                                                  takeoffF[Cycle > 5]),
          by = regionStart]
  }

  if (is.numeric(regionEnd)) {
    inptl[, (c("derMaxC", "derMaxF")) := list(1, regionEnd)]
    inptl[, regionEnd := regionEnd]
  } else {
    inptl[, (c("derMaxC", "derMaxF")) := detectIndividualEndPoints(Cycle, RFU),
          by = well]
    inptl[, regionEnd := detectGlobalEndPoint(derMaxF),
          by = regionEnd]
  }

  inptlf <- inptl[Cycle > 5 &
                    Cycle >= minStartPointC &
                    RFU <= regionEnd & RFU >= regionStart]

  gTbl <- rbindlist(
    lapply(unique(inptlf$set), function(setName) {
      do.call(cbind.data.frame, Map(expand.grid,
                                    i = inptlf[set == setName],
                                    j = inptlf[set == setName]))
    }
    )
  )
  gTbl <- data.table(gTbl)[Cycle.i < Cycle.j][
    , -c("set.j", "gene.j", "regionStart.j", "regionEnd.j")
    ]
  colnames(gTbl)[colnames(gTbl) %in% c("set.i", "gene.i", "regionStart.i", "regionEnd.i")] <-
    c("set", "gene", "regionStart", "regionEnd")

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
  gTbl[, ':='(mean_pEf = mean(pE),
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
    sorted <- sort(c(mean_pE[1], mean_pEf[1]))
    paste0(
      {if (sorted[2]/sorted[1] > 1.15)
        "Possibly bad data: mean pE diff > 15%;"
        else ""
      }, {if (any(mean_F0 > regionStart[1]))
        "Possibly bad data: F0 > working range start point;"
        else ""
      }, {if (kurtosis[1] < 1)
        "Possibly bad data: excess kurtosis < 1;"
        else ""
      })
  }, by = set]

  res <- unique(gTbl[, .(gene, set,
                         mean_pE, mean_pEf, sd_pE2,
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
