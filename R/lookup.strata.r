
  lookup.strata = function(db) {

    Vn = as.character(c(440:442))   # these are alphanumeric codes
    Vs = as.character(c(443:452))   # these are alphanumeric codes
    W  = as.character(c(453:466))   # these are alphanumeric codes
    X  = as.character(c(470:495))   # these are alphanumeric codes
    Ge = c("5Z1","5Z2","5Z3","5Z4","5Z5","5Z6","5Z7","5Z8")

    out = data.frame(strat=c(Vn, Vs, W, X, Ge), region=NA)
    out$strat = as.character(out$strat)

    out$region[out$strat %in% Vn] = "Vn"
    out$region[out$strat %in% Vs] = "Vs"
    out$region[out$strat %in% W ] = "W"
    out$region[out$strat %in% X ] = "X"
    out$region[out$strat %in% Ge] = "Ge"

    return(out)
  }


