library(dlm)
library(magrittr)
source("R/models.R")
model <- poly.dm(order = 2, delta = 0.985) %>%
  superposition.dm(
    trig.dm(s = 4, delta = 0.98)
  )
delta <- model$delta[1]
dim_p <- c(2, 3)[1]
sapply(1:length(delta), function(x) matrix(1/delta[x],
                                           nrow=dim_p[x],ncol=dim_p[x]))
C0 <- diag(100,5)
D1 <- matrix(1/delta[1],2,2)
D2 <- matrix(1/delta[2],3,3)
##
DD <- bdiag(D1, D2)
DD <- diag(rep(1/delta, c(2,2)), 5)
DD[which(DD==0)] <- 1
GG <- bdiag(model$GG[[1]], model$GG[[2]])
PP <- tcrossprod(GG %*% C0, GG)
RR <- DD * PP
##
P_1 <- tcrossprod(model$GG[[1]] %*% C0[1:2,1:2], model$GG[[1]])
R_1 <- P_1/delta[1]
P_2 <- tcrossprod(model$GG[[2]] %*% C0[-(1:2),-(1:2)], model$GG[[2]])
R_2 <- P_2/delta[2]
R <- bdiag(R_1, R_2)

RR;R
