set.seed(1)

unif.matrix = function(nrow, ncol, beta){
  M = matrix(runif(nrow*ncol), nrow, ncol)
  M[,1] = M %*% beta_star
  y = (M[,1] > 0.5)
  y = as.numeric(y)
}

W_star.matrix = function(matrix, beta_star){
  support = beta_star
  W_star = diag(matrix %*% bata_star )
}

Cmin = function(X_star.matrix, W_star.matrix, obs.num ){
  ##minimum eigen value of
  #X_star.matrix * W_star.matrix * X_star.matrix / n
}

Dmax = function(X_star.matrix, obs.num){
  #maximum eigen value
  #X_star.matrix * X_star.matrix / n
}

b = function(X_star.matrix, W_star.matrix){
  #inf norm [(X_star.matrix * W_star.matrix * X_star.matrix)^(-1)]*[X_star.matrix*W_star.matrix * X_star]
}

L = function(X.matrix){
  #max(X.matrix)
}

lambda_max = function() {
  (b*Cmin)/(100*(2-b)*Dmax*L)
}

#C_bar greater or equal to (3a/2Cmin)
tunPar = function(lambda_max, beta, C_bar) {
  #indicator function t = 
  lambda_grid = seq(0, lambda_max, 0.1)
  j = N
  while(t != 0 & j > 1){
    j = j - 1
  }
}