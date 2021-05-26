listdata = simuldata(n = 8000, l = 4, ep = -1.6)$List_matrix
simuldata(n = 2000, l = 3, ep = -1.6)$psi

split = sort(sample(1:nrow(listdata), nrow(listdata)/2))
####### fitting the models to obtain q1, q2 and q12
Sys.time()
qhat1 = qhat_sl(List.train = listdata[split,],
                   List.test = listdata[-split,], K = 2, i = 1, j = 2, num_cores = 1)
Sys.time()
qhat2 = qhat_sl(List.train = listdata[split,],
                 List.test = listdata[-split,], K = 2, i = 1, j = 2, num_cores = 7)
Sys.time()


plot(qhat1$q12, qhat2$q12); abline(0,1)
plot(qhat1$q1, qhat2$q1); abline(0,1)
plot(qhat1$q2, qhat2$q2); abline(0,1)

