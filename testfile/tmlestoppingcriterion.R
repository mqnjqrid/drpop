set.seed(9)
nsum = matrix(NA, nrow = 25*3, ncol = 3)
colnames(nsum) = c("epiter", "logit.TMLE", "gam.TMLE")
nsum = as.data.frame(nsum)
nsum$epiter = rep(#c(0.1, 0.01, 0.001)
  c(100, 500, 1000)/20
  , each = 25)
for(s in 1:25){
  print(s)
  listdata = simuldata2(n = 3000, l = 3, K = 3, ep = -4)$List_matrix_xstar
  listdata = reformat(listdata, capturelists = c('y1', 'y2', 'y3'))
  listdata = listdata[,-2]
  qhat = qhateval(List_matrix = listdata, funcname = c("logit", "gam"), nfolds = 2)
  for(epiter in
      #c(0.1, 0.01, 0.001)
      c(100, 500, 1000)/20
      ){
    #listdata = simuldata(n = 6000, l = 2, ep = -3)$List_matrix
    nhat = psinhatgivenq2(List_matrix = listdata, qhateval = qhat, eps_stop = 0.0005, iter = epiter, eps = 0)
    nsum[nsum$epiter == epiter,][s,-1] = nhat$n[,c("logit.TMLE", "gam.TMLE")]
    #print(c(nhat$iterations, nhat$epsilon_error))
  }
}
nsum2 = reshape2::melt(nsum, id = "epiter")
library(ggplot2)
ggplot(nsum2) + geom_density(aes(x = value, fill = factor(epiter)), alpha = 0.3) + geom_vline(xintercept = 3000) + facet_wrap(~variable)
print.psinhat(nhat)
