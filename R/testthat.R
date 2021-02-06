if(FALSE){
n = 100
x = matrix(rnorm(n*3, 2, 1), nrow = n)
expit = function(xi) {
  exp(sum(xi))/(1 + exp(sum(xi)))
}
y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.4*xi), expit(-0.6 + 0.4*xi)))}))
y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - expit(-0.6 + 0.3*xi), expit(-0.6 + 0.3*xi)))}))
datacrc = cbind(y1, y2, exp(x/2))

p = psinhat(List_matrix = datacrc, funcname = c("logit", "gam", "sl"))
plot(psinhat = p)

ss = sample(c('a','b','c','d','e','f'), nrow(datacrc), replace = TRUE, prob = (1:6)/sum(1:6))

datacrc1 = data.frame(datacrc, ss)

p = psinhatcond(List_matrix = datacrc1, condvar = 'ss')
plot(psinhatcond = p)
}
