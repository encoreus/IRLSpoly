#' Polychoric Correlation
#' @description
#' Estimate the polychoric correlation coefficient.
#'
#' @param X a matrix or dataframe contains two polychoric variable, or a contingency table with both columns and rows names.
#' @param maxn the maximum iterations times.
#' @param e the maximum tolerance of convergence.
#' @param ct \code{TRUE} for contingency table, \code{FALSE} for matrix or dataframe
#'
#' @return \item{rho}{estimated value of polychoric correlation coefficient.}
#' @return \item{std}{standard deviation of rho.}
#' @return \item{iter}{times of iteration convergence.}
#' @return \item{Ex,Ey}{the support points series of regression model}
#'
#' @seealso
#' [esti_polyserial]
#'
#' @export
#'
#' @examples
#' X = gen_polychoric(1000,0.5,0:1,-1:0)
#' result = esti_polychoric(X)
#' print(c(result$rho,result$std,result$iter))
esti_polychoric <- function(X,maxn=100,e=1e-8,ct = FALSE){
  if(ct){
    N = sum(X)
    P = X/N}
  else{
    x = X[1,]
    y = X[2,]
    N = sum(table(y,x))
    P = table(y,x)/N
  }
  if(N<100){
    print('estimated std may smaller than true std.')
  }

  Pa = colSums(P);cPa = cumsum(Pa);s = length(Pa)
  Pb = rowSums(P);cPb = cumsum(Pb);r = length(Pb)
  a = qnorm(cPa);a[s] = 1e8;a = c(-1e8,a)
  b = qnorm(cPb);b[r] = 1e8;b = c(-1e8,b)

  Pl = matrix(P,1,s*r)
  B = t(Pl)%*%Pl
  B = diag(Pl[1,])-B
  B = B/N

  if(ct){
    rho0 = 0
  }else{
    rho0 = cor(x,y)
  }
  ex = (dnorm(a[1:s])-dnorm(a[2:(s+1)]))/Pa
  Pw1 = t(t(P)/Pa)
  Pw2 = P/Pb

  iter = 0
  dif = 1
  n = 100
  e = 1e-8

  Ex = c()
  Ey = c()
  while((iter<n)&(dif>e)){
    midy = matrix(rep(0,s*r),r,s)
    for (i in 1:r) {
      mid1 = (b[i]-rho0*ex)/sqrt(1-rho0^2)
      mid2 = (b[i+1]-rho0*ex)/sqrt(1-rho0^2)
      midy[i,] = rho0*ex + sqrt(1-rho0^2)*(dnorm(mid1)-dnorm(mid2))/(pnorm(mid2)-pnorm(mid1))
    }
    EY = colSums(Pw1*midy)
    Ex = cbind(Ex,ex)
    Ey = cbind(Ey,EY)

    D = matrix(rep(0,s*s*r),s,s*r)
    for (i in 1:s) {
      P_i = diag(rep(Pa[i],r))-t(matrix(rep(P[,i],r),r,r))
      D[i,(1+(i-1)*r):(i*r)] = P_i%*%midy[,i]/(Pa[i]^2)
    }
    Sigma = D%*%B%*%t(D)
    #if(is.singular.matrix(Sigma)){Sigma = diag(s)/length(x)}

    rho = solve(t(ex)%*%solve(Sigma)%*%ex)%*%t(ex)%*%solve(Sigma)%*%EY
    rho = rho[1,1]
    if(is.na(rho)){rho = cor(x,y)}
    if(rho>=1){
      rho=1-1e-4
    }else if(rho<=-1){
      rho = -1+1e-4}
    dif = abs(rho-rho0)

    EYx = rowSums(Pw2*midy)
    midx = matrix(rep(0,r*s),r,s)
    for (i in 1:s) {
      mid1 = (a[i]-rho*EYx)/sqrt(1-rho^2)
      mid2 = (a[i+1]-rho*EYx)/sqrt(1-rho^2)
      midx[,i] = rho*EYx + sqrt(1-rho^2)*(dnorm(mid1)-dnorm(mid2))/(pnorm(mid2)-pnorm(mid1))
    }
    ex = colSums(Pw1*midx)
    rho0 = rho
    iter = iter + 1
  }
  varrho = solve(t(ex)%*%solve(Sigma)%*%ex)

  return(list(rho=rho,
              std = sqrt(varrho[1,1]),
              iter = iter,
              Ex = Ex,
              Ey = Ey))
}
