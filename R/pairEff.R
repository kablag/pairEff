pE <- function(dilutionRate,
               RFU_i, RFU_j,
               dilution_i, dilution_j,
               cycle_i, cycle_j) {
  (dilutionRate ^ (((log(RFU_j, base = dilutionRate)
                      - log(RFU_i, base = dilutionRate)) +
           (dilution_j - dilution_i)) /
          (cycle_j - cycle_i))) - 1
}

pE2 <- function(dilutionRate,
                dilutionDelta,
               RFU_i, RFU_j,
               cycle_i, cycle_j) {
  (dilutionRate ^ (((log(RFU_j, base = dilutionRate)
                     - log(RFU_i, base = dilutionRate)) +
                      dilutionDelta) /
                     (cycle_j - cycle_i))) - 1
}
