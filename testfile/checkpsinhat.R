ld = simuldata(n = 3000, l = 3, ep = -4)
List_matrix = ld$List_matrix
p1 = apply(as.matrix(List_matrix[,-c(1:2)]), 1, ld$pi1)
p2 = apply(as.matrix(List_matrix[,-c(1:2)]), 1, ld$pi2)
gam = 1 - (1 - p1)*(1 - p2)
q10 = p1/gam
q20 = p2/gam
q120 = p1*p2/gam

K = 2
l = ncol(List_matrix) - K
n = nrow(List_matrix)
nfolds = 2
eps = 0.00
func = "rangerlogit"

List_matrix = na.omit(List_matrix)


List_matrix = as.data.frame(List_matrix)
#N = number of observed or captured units
N = nrow(List_matrix)

#renaming the columns of List_matrix for ease of use
colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix) - K), sep = ''))

psiinv_summary = matrix(0, nrow = K*(K - 1)/2, ncol = 3*length(funcname))
rownames(psiinv_summary) = unlist(sapply(1:(K - 1), function(k) {
  sapply((k + 1):K, function(s) {
    return(paste0(k, ",", s))
  })}))
colnames(psiinv_summary) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
var_summary = psiinv_summary

ifvals = matrix(NA, nrow = N*K*(K-1)/2, ncol = length(funcname))
colnames(ifvals) = funcname
rownames(ifvals) = rep(rownames(psiinv_summary), each = N)

nuis = matrix(NA, nrow = N*K*(K-1)/2, ncol = 3*length(funcname))
colnames(nuis) = paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.')
rownames(nuis) = rownames(ifvals)
nuistmle = nuis

permutset = sample(1:N, N, replace = FALSE)
i = 1; j = 2

psiinvmat = matrix(numeric(0), nrow = nfolds, ncol = 3*length(funcname))
colnames(psiinvmat) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
varmat = psiinvmat

ifvalsfold = matrix(numeric(0), nrow = N, ncol = length(funcname))
colnames(ifvalsfold) = funcname

nuisfold = matrix(numeric(0), nrow = N, ncol = 3*length(funcname))
colnames(nuisfold) = paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.')
nuistmlefold = nuisfold

folds = 1
if(nfolds == 1){
  List1 = List_matrix
  List2 = List1
  sbset = 1:N
}else{
  sbset = ((folds - 1)*ceiling(N/nfolds) + 1):(folds*ceiling(N/nfolds))
  sbset = sbset[sbset <= N]
  List1 = List_matrix[permutset[-sbset],]
  List2 = List_matrix[permutset[sbset],]
}

yi = List2[,paste("L", i, sep = '')]
yj = List2[,paste("L", j, sep = '')]

qhat = try(get(paste0("qhat_", func))(List.train = List1, List.test = List2, K, i, j, eps = eps), silent = TRUE)

qhat = qhateval(List_matrix, funcname = c("logit", "ranger", "rangerlogit"), eps = 0)
class(qhat)
summary(qhat$q12mat)

q12 = qhat$q12
q1 = pmin(pmax(q12, qhat$q1), 1)
q2 = pmax(q12/q1, pmin(qhat$q2, 1 + q12 - q1, 1))

plot(q120[permutset[sbset]], q12)
plot(q10[permutset[sbset]], q1)
plot(q20[permutset[sbset]], q2)

gammainvhat = q1*q2/q12
psiinvhat = mean(gammainvhat)
phihat = gammainvhat*(yj/q2 + yi/q1 - yi*yj/q12) - psiinvhat

Qnphihat = mean(phihat + gammainvhat, na.rm = TRUE)
psiinvhat
Qnphihat
summary(phihat)
summary(cbind(gammainvhat*yj/q2, gammainvhat*yi/q1, - gammainvhat*yi*yj/q12))

