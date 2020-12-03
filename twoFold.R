twoFold <- function(startCon, 
                    compound = c("Butyrate", "Propionate","Acetate"), 
                    series = 8, 
                    prefix = "CC") {
  require(assertive.types)
  require(assertive.base)
  require(tidyverse)
  # arg match
  compound <- match.arg(compound)
  # checks
  assert_is_numeric(startCon)
  startCon <- use_first(startCon)
  # body
  conc <- NULL
  conc[1] <- cur <- startCon
  i = 2
  while( i <= series){ 
    conc[i] <- cur/2
    cur <- conc[i]
    i <- i + 1
  }
  # combine to tibble and return
  odf <- tibble(compound, conc, curveLab = paste0(prefix,c(series:1)))
  return(odf)
}
