set.seed(1234)
funcname = c("rangerlogit", "logit", "ranger")
result = numeric(0)
n = 5000
for(s in 1:30){
  print(s)
  ld = simuldata2(n = n, l = 6, ep = -6)
  listdata = ld$List_matrix
  listdata2 = ld$List_matrix_xstar
  p1 = psinhat(listdata, funcname = funcname, nfolds = 2, iter = 250, eps_stop = 0.01)
  p2 = psinhat(listdata2, funcname = funcname, nfolds = 2, iter = 250, eps_stop = 0.01)
  result = rbind(result, cbind(t(p1$n), t(p1$varn), "x"),
                 cbind(t(p1$n), t(p1$varn), "x*"))
}

result2 = cbind(rownames(result), result)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
result2 = as.data.frame(result2)
colnames(result2) = c("Var1", "n", "varn", "covariate")
result2$n = as.numeric(as.character(result2$n))
result2$varn = as.numeric(as.character(result2$varn))
result2$biasn = result2$n - n
ggplot(result2, aes(x = biasn, fill = Var1)) + geom_density(alpha = 0.3)
result3 = ddply(result2, c("Var1", "covariate"), summarise,
                bias = median(abs(biasn)),
                rmse = sqrt(median(biasn^2 + varn)))
result3 = result3 %>% separate(Var1, c("model", "method"))

ggplot(result3 %>% filter(model != "rangerj"), aes(x = model, y = bias, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~covariate) +
  scale_fill_manual("Estimation method", values = c("PI" = "red", "DR" = "#E69F00", "TMLE" = "#56B4E9"))
ggplot(result3 %>% filter(model != "rangerj"), aes(x = model, y = rmse, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~covariate) +
  scale_fill_manual("Estimation method", values = c("PI" = "red", "DR" = "#E69F00", "TMLE" = "#56B4E9"))
