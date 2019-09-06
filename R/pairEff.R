# pE <- function(dilutionRate,
#                 dilutionDelta,
#                RFU_i, RFU_j,
#                cycle_i, cycle_j) {
#   (dilutionRate ^ (((log(RFU_j, base = dilutionRate)
#                      - log(RFU_i, base = dilutionRate)) +
#                       dilutionDelta) /
#                      (cycle_j - cycle_i))) - 1
# }

pE <- function(dilutionRate,
               dilutionDelta,
               RFU_i, RFU_j,
               cycle_i, cycle_j) {
  (dilutionRate ^ ((log(RFU_j/RFU_i, base = dilutionRate) +
                      dilutionDelta) /
                     (cycle_j - cycle_i))) - 1
}

#' Calculates pairwise efficiency
#'
#' Calculates pairwise efficiency
#'
#' @import tidyverse
#' @importFrom readxl read_excel
#' @export
pairEff <- function(filename, modtype) {
  inptw <- read_excel(filename)

  cat("Fitting models\n")
  fits <- modlist(inptw, 1, 2:ncol(inptw),
                  model = modtype, verbose = FALSE)
  names(fits) <- colnames(inptw)[-1]

  for (cname in names(fits)) {
    inptw[[cname]] <- fits[[cname]]$m$predict()
  }

  cat("Add conc sets\n")
  inptl <- gather(inptw, "Well", "RFU", -Cycle)
  inptl$conc <- rep(rep(c(100, 50, 25, 12, 6, 3), each = max(inptl$Cycle)), 16)
  # maxConc <- max(inptl$conc)
  # inptl$dilution <- log(maxConc / inptl$conc, base = 2)
  inptl$set <- rep(sprintf("Set%02i", c(1:16)), each = max(inptl$Cycle) * 6)

  cat("Filter curves regions\n")
  # browser()
  inptl <- inptl %>%
    group_by(Well) %>%
    mutate(#fTop = RFU[takeoff(fits[[Well[1]]])$top + 1],
      #top = takeoff(fits[[Well[1]]])$top,
      # fmidp = midpoint(fits[[Well[1]]])$f.mp,
      fcpD2 = efficiency(fits[[Well[1]]], type = "cpD2", plot = FALSE)$fluo,
      fcpD1 = efficiency(fits[[Well[1]]], type = "cpD1", plot = FALSE)$fluo
      # ,
      # cpD2 = efficiency(fits[[Well[1]]], type = "cpD2", plot = FALSE)$cpD2,
      # cpD1 = efficiency(fits[[Well[1]]], type = "cpD1", plot = FALSE)$cpD1
      # ,
      # usePoint = RFU >= fcpD2[1] & RFU <= fcpD1[1]) %>%
      # usePoint =
        # Cycle >= (cpD2[1] - 1) & Cycle <= (cpD1[1] + 1)
      ) %>%
    group_by(set) %>%
    mutate(
      # usePoint = RFU >= (min(fmidp) * 1.1) &
        # RFU <= (max(fcpD1) * 0.9)
      # usePoint = RFU >= min(fmidp) & RFU <= max(fcpD1)
      usePoint = RFU >= (min(fcpD2) - (max(fcpD1) - min(fcpD2)) * 0.1) &
        RFU <= (max(fcpD1) + (max(fcpD1) - min(fcpD2)) * 0.1)
      # usePoint = Cycle >= min(cpD2) & Cycle <= max(cpD1)
    )

  finptl <- inptl %>%
    filter(usePoint)

  # finptl <- inptl %>%
    # filter(RFU > 100 & RFU < 300)
    # filter(RFU > 20 & RFU < 180)
  cat("Filter pEs\n")
  toInclude <- function(pes) {
    pEsummary <- summary(pes)
    pes >= pEsummary[2] & pes <= pEsummary[5]
  }

  cat("Calc pEs\n")
  pEtbl <-
    map_dfr(finptl$set %>% unique,
            function(cset) {
              fd <- finptl %>%
                filter(set == cset) %>%
                mutate(i = seq(n()),
                       j = i)
              fdi <- fd %>%
                rename_at(vars(-i, -j), ~paste(., "i", sep = "_"))
              fdj <- fd %>%
                rename_at(vars(-i, -j), ~paste(., "j", sep = "_"))

              expand.grid(i = seq(nrow(finptl)),
                          j = seq(nrow(finptl))) %>%
                filter(j > i) %>%
                left_join(fdi, by = "i") %>%
                rename(j = j.x) %>%
                left_join(fdj, by = "j") %>%
                rename(i = i.x, set = set_i) %>%
                filter(Cycle_i != Cycle_j) %>%
                dplyr::select(-j.y, -i.y, -set_j) %>%
                mutate(
                  dilutionRate = ifelse(conc_i == conc_j, # one dilution
                                        2,
                                        ifelse(conc_i > conc_j,
                                               conc_i / conc_j,
                                               conc_j / conc_i)),
                  dilutionDelta = ifelse(conc_i == conc_j,
                                         0,
                                         ifelse(conc_i > conc_j,
                                                1, -1)),
                  pE = pE(dilutionRate,
                            dilutionDelta,
                            RFU_i, RFU_j,
                            Cycle_i, Cycle_j),
                  F0 = RFU_i / ((1 + pE) ^ Cycle_i),
                  include = toInclude(pE)
                )
              # %>%
              #   group_by(Well_i) %>%
              #   mutate(F0 = RFU_i[pE == max(pE)] /
              #            ((1 + max(pE)) ^ Cycle_i[pE == max(pE)]))
            })


  cat("Calc results\n")
  result <- pEtbl %>%
    # mutate(include = pE > 0.6 & pE < 1.05) %>%
    group_by(set) %>%
    # filter(include) %>%
    summarise(mean_pE = mean(pE[include]),
              sd_pE = sd(pE[include]),
              mean_F0 = mean(F0[include]),
              sd_F0 = sd(F0[include]),
              totalN = n(),
              includeN = sum(include)
    )
  list(pEtbl = pEtbl,
       result = result,
       inptl = inptl,
       mods = fits)
}
