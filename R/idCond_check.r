#' @name idcond_check
#' @description  Check the identifiability conditions C1-C2 of Park and Oh (2015).
#' @title Check the identifiability conditions
#' @usage idCond_check(P)
#' @param P  the source composition matrix in multivariate receptor model


idCond_check = function(P){
  q = base::nrow(P)

  idCond_1 = 1
  idCond_2 = 1
  for (k in 1:q){
    idCol = base::which(P[k,]==0)
    idCond_1 = idCond_1*(base::length(idCol)  >= q-1)  # Condition C1
    idRow = base::c(1:q)[-k]
    P_k = P[idRow, idCol]
    idCond_2 = idCond_2*(base::qr(P_k)$rank  == q-1)   # Condition C2
  }

  idCond = idCond_1 * idCond_2

  return(idCond) # if all the conditions are satisfied, it is TRUE (1).
}
