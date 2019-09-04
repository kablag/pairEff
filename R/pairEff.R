# pE <- function(dilutionRate,
#                RFU_i, RFU_j,
#                dilution_i, dilution_j,
#                cycle_i, cycle_j) {
#   (dilutionRate ^ (((log(RFU_j, base = dilutionRate)
#                       - log(RFU_i, base = dilutionRate)) +
#            (dilution_j - dilution_i)) /
#           (cycle_j - cycle_i))) - 1
# }

pE <- function(dilutionRate,
                dilutionDelta,
               RFU_i, RFU_j,
               cycle_i, cycle_j) {
  (dilutionRate ^ (((log(RFU_j, base = dilutionRate)
                     - log(RFU_i, base = dilutionRate)) +
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
pairEff <- function(filename) {
  inptw <- read_excel(filename)
  inptl <- gather(inptw, "Well", "RFU", -Cycle)
  inptl$conc <- rep(rep(c(100, 50, 25, 12, 6, 3), each = max(inptl$Cycle)), 16)
  # maxConc <- max(inptl$conc)
  # inptl$dilution <- log(maxConc / inptl$conc, base = 2)
  inptl$set <- rep(sprintf("Set%02i", c(1:16)), each = max(inptl$Cycle) * 6)
  finptl <- inptl %>%
    filter(RFU > 20 & RFU < 180)
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
                select(-j.y, -i.y, -set_j) %>%
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
                  F0 = RFU_i / ((1 + pE) ^ Cycle_i)
                )
            })
  result <- pEtbl %>%
    mutate(include = pE > 0.6 & pE < 1.05) %>%
    group_by(set) %>%
    summarise(mean_pE = mean(pE[include]),
              sd_pE = sd(pE[include]),
              mean_F0 = mean(F0[include]),
              sd_F0 = sd(F0[include]),
              totalN = n(),
              includeN = sum(include)
    )
  list(pEtbl = pEtbl,
       result = result)
}
