pE <- function(dilutionRate,
               dilutionDelta,
               fluor.i, fluor.j,
               cyc.i, cyc.j) {
  (dilutionRate ^ ((log(fluor.j / fluor.i, base = dilutionRate) +
                      dilutionDelta) /
                     (cyc.j - cyc.i))) - 1
}

detectIndividualStartPoint <- function(fPoints, pval = 0.05, nsig = 3)
{
  cycs <- seq_along(fPoints)
  res <- sapply(5:length(cycs),
                function(i) {
                  mod <- lm(fPoints[1:i] ~ cycs[1:i], na.action = na.exclude)
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
detectGlobalStartPoint <- function(fluor, takeoffF) {
  FBeforeTakeoffFs <- sort(fluor[fluor <= max(takeoffF)])
  meanFBeforeTakeoffFs <- mean(FBeforeTakeoffFs)
  tail(FBeforeTakeoffFs[FBeforeTakeoffFs <
                          (meanFBeforeTakeoffFs * 4)], 1) * 2
}

detectIndividualEndPoints <- function(cycs, fPoints) {
  der <- c(sapply(1:(length(cycs) - 1),
                  function(i) fPoints[i + 1] - fPoints[i]), NA)
  i <- which(der == max(der, na.rm = TRUE))
  list(derMaxC = cycs[i], derMaxF = fPoints[i])
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
pairEff <- function(rdml, regionStart = "gene", regionEnd = "gene") {
  inptl <- rdml$GetFData(
    rdml$AsTable(conc = sample[[react$sample$id]]$quantity$value)[sample.type == "std"],
    long.table = TRUE)
  # inptl <- melt(inptw, id.vars = "cyc", variable.name = "position", value.name = "fluor")
  # inptl <- inptl[plateStr, on = .(position)]
  inptl[, plate := exp.id]
  tmp <- strsplit(inptl$target, "@")
  inptl[, set := sapply(tmp, function(el) el[1])]
  inptl[, gene := sapply(tmp, function(el) el[2])]
  if (is.numeric(regionStart)) {
    inptl[, (c("takeoffC", "takeoffF")) := list(1, regionStart)]
    inptl[, minStartPointC := 1]
    inptl[, regionStart := regionStart]
  } else {
    inptl[, (c("takeoffC", "takeoffF")) := detectIndividualStartPoint(fluor),
          by = position]
    inptl[, minStartPointC := min(takeoffC, na.rm = TRUE),
          by = regionStart]
    inptl[, regionStart := detectGlobalStartPoint(fluor[cyc > 5],
                                                  takeoffF[cyc > 5]),
          by = regionStart]
  }

  if (is.numeric(regionEnd)) {
    inptl[, (c("derMaxC", "derMaxF")) := list(1, regionEnd)]
    inptl[, regionEnd := regionEnd]
  } else {
    inptl[, (c("derMaxC", "derMaxF")) := detectIndividualEndPoints(cyc, fluor),
          by = position]
    inptl[, regionEnd := detectGlobalEndPoint(derMaxF),
          by = regionEnd]
  }

  inptlf <- inptl[cyc > 5 &
                    cyc >= minStartPointC &
                    fluor <= regionEnd & fluor >= regionStart]

  gTbl <- rbindlist(
    lapply(unique(inptlf$set), function(setName) {
      do.call(cbind.data.frame, Map(expand.grid,
                                    i = inptlf[set == setName],
                                    j = inptlf[set == setName]))
    }
    )
  )
  gTbl <- data.table(gTbl)[cyc.i < cyc.j][
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
                          fluor.i, fluor.j,
                          cyc.i, cyc.j)][
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
  gTbl[, F0 := fluor.i / ((1 + pE) ^ cyc.i)][
    , ':='(mean_F0 = mean(F0),
           sd_F0 = sd(F0),
           kurtosis = kurtosis(pE)),
    by = "set"
    ]
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
