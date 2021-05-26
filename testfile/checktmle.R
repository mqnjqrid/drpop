
tmle <- function(datmat, iter = 100, eps = 0.005, eps_stop = 0.01, twolist = FALSE, K = 2,...){

  if(!prod(c("yi", "yj", "yij", "q10", "q02", "q12") %in% colnames(datmat))){
    stop("datmat misses one or more of the following columns: \t (yi, yj, yij, q10, q02, q12).")
    return(list(error = TRUE))
  }

  ifval = (datmat$q10+datmat$q12)*(datmat$q10+datmat$q12)/datmat$q12 *(
    datmat$yi/(datmat$q10+datmat$q12) + datmat$yj/(datmat$q02+datmat$q12) - datmat$yij/(datmat$q12) - 1)
  eps_stop = sqrt(mean((ifval)^2))/max(log(nrow(datmat)), 100)/sqrt(nrow(datmat))
  expit = function(x) {
    exp(x)/(1 + exp(x))
  }
  logit = function(x) {
    log(x/(1 - x))
  }

  epsilon_error = eps_stop + 1
  cnt = 0

  while (abs(epsilon_error) > eps_stop){
    cnt = cnt + 1
    if (cnt > iter){break}

    ########################### model 1 for q12
    dat1 = cbind(datmat$yij, logit(datmat$q12), (datmat$q10 + datmat$q12)/datmat$q12
                 + (datmat$q02 + datmat$q12)/datmat$q12
                 - (datmat$q10 + datmat$q12)*(datmat$q02 + datmat$q12)/datmat$q12^2 )
    colnames(dat1) = c("yij", "logitq12", "ratio")
    dat1 = as.data.frame(dat1)
    mod1 = try(glm(yij ~ -1 + offset(logitq12) + ratio
                   , family = binomial(link = logit), data = dat1, na.action = na.omit), silent = TRUE)
    if (!("try-error" %in% class(mod1))){
      datmat[,"q12"] = predict(mod1, newdata = dat1, type = "response")

    }
    datmat$q12 = pmax(pmin(datmat$q12, 1), eps)

    ########################### model 2 for q1
    dat2 = cbind(datmat$yi*(1 - datmat$yj), logit(datmat$q10), (datmat$q02 + datmat$q12)/datmat$q12)
    colnames(dat2) = c("yi0", "logitq10", "ratio")
    dat2 = as.data.frame(dat2)
    mod2 = try(glm(yi0 ~ -1 + offset(logitq10) + ratio, family = binomial(link = logit), data = dat2, na.action = na.omit), silent = TRUE)
    if (!("try-error" %in% class(mod2))){
      datmat$q10 = predict(mod2, newdata = dat2, type = "response")
      datmat[,"q10"] = pmin(datmat[,"q10"], 1 - datmat$q12)
    }
    datmat$q10 = pmax(pmin(datmat$q10, 1), eps)

    ########################### model 3 for q2
    if (K > 2 | twolist == FALSE){
      dat3 = cbind(datmat$yj*(1 - datmat$yi), logit(datmat$q02), (datmat$q10 + datmat$q12)/datmat$q12)
      colnames(dat3) = c("y0j", "logitq02", "ratio")
      dat3 = as.data.frame(dat3)
      mod3 = try(glm(y0j ~ -1 + offset(logitq02) + ratio, family = binomial(link = logit), data = dat3, na.action = na.omit), silent = TRUE)
      if (!("try-error" %in% class(mod3))){
        datmat$q02 = predict(mod3, newdata = dat3, type = "response")
        datmat[,"q02"] = pmin(datmat[,"q02"], 1 - datmat$q10 - datmat$q12)
      }
    }else{
      mod3 = mod2
      datmat[,"q02"] = pmax(0, 1 - datmat$q10 - datmat$q12)
    }

    datmat$q02 = pmax(pmin(datmat$q02, 1), eps)

     ifval = (datmat$q10+datmat$q12)*(datmat$q10+datmat$q12)/datmat$q12 *(
       datmat$yi/(datmat$q10+datmat$q12) + datmat$yj/(datmat$q02+datmat$q12) - datmat$yij/(datmat$q12) - 1)
     epsilon_error = mean(ifval)
    # print(c(cnt, epsilon_error, eps_stop))
    #epsilon_error = max(abs(c(mod2$coefficients, mod3$coefficients, mod1$coefficients)))
    #epsilon_error = max(abs(expit(dat1$logitq12) - datmat$q12),
    #                    abs(expit(dat2$logitq10) - datmat$q10),
    #                    abs(expit(dat3$logitq02) - datmat$q02))
  }
  return(list(error = abs(epsilon_error) > 1, datmat = datmat, iterations = cnt, epsilon_error = epsilon_error))
}

psinhatgivenq2 <- function(List_matrix, i = 1, j = 2, eps = 0.005, qhateval, q1mat, q2mat, q12mat, idfold, ...){

  K = 2
  n = nrow(List_matrix)

  if(missing(i)){
    i = 1
  }
  if(missing(j)){
    j = 2
  }

  stopifnot(!(missing(qhateval) & missing(q1mat) & missing(q2mat) & missing(q12mat)))

  if(!missing(qhateval)){
    q1mat = qhateval$q1mat
    q2mat = qhateval$q2mat
    q12mat = qhateval$q12mat
    idfold = qhateval$idfold
  }

  stopifnot(!is.null(q1mat) & !is.null(q2mat) & !is.null(q12mat))

  funcname = colnames(q12mat)

  if(missing(idfold) | is.null(idfold)){
    idfold = rep(1, n)
  }
  nfolds = max(idfold)
  stopifnot(!is.null(dim(List_matrix)))

  if(!informat(List_matrix = List_matrix, K = K)){
    List_matrix <- reformat(List_matrix = List_matrix, capturelists = 1:K)
  }

  List_matrix = as.data.frame(List_matrix)
  #N = number of observed or captured units
  N = nrow(List_matrix)

  stopifnot(N > 1)

  #renaming the columns of List_matrix for ease of use
  colnames(List_matrix) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List_matrix) - K), sep = ''))

  psiinv_summary = matrix(0, nrow = K*(K - 1)/2, ncol = 3*length(funcname))
  rownames(psiinv_summary) = paste0(i, ",", j)
  colnames(psiinv_summary) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
  var_summary = psiinv_summary

  ifvals = matrix(NA, nrow = N*K*(K-1)/2, ncol = length(funcname))
  colnames(ifvals) = funcname
  rownames(ifvals) = rep(rownames(psiinv_summary), each = N)

  nuis = matrix(NA, nrow = N*K*(K-1)/2, ncol = 3*length(funcname))
  colnames(nuis) = paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = '.')
  rownames(nuis) = rownames(ifvals)
  nuistmle = nuis

  psiinvmat = matrix(NA, nrow = nfolds, ncol = 3*length(funcname))
  colnames(psiinvmat) = paste(rep(funcname, each = 3), c("PI", "DR", "TMLE"), sep = '.')
  varmat = psiinvmat

  for (folds in 1:nfolds){#print(folds)

    List2 = List_matrix[idfold == folds,]

    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]

    for (func in funcname){

      colsubset = stringr::str_subset(colnames(psiinv_summary), func)

      q12 = q12mat[idfold == folds, func]
      q1 = pmin(pmax(q12, q1mat[idfold == folds, func]), 1)
      q2 = pmax(q12/q1, pmin(q2mat[idfold == folds, func], 1 + q12 - q1, 1))

      nuis[idfold == folds, paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

      gammainvhat = q1*q2/q12
      psiinvhat = mean(gammainvhat, na.rm = TRUE)

      phihat = gammainvhat*(yj/q2 + yi/q1 - yi*yj/q12) - psiinvhat
      ifvals[idfold == folds, func] = phihat

      Qnphihat = mean(phihat, na.rm = TRUE)

      psiinvhat.dr = max(psiinvhat + Qnphihat, 1)

      psiinvmat[folds, colsubset][1:2] = c(psiinvhat, psiinvhat.dr)

      sigmasq = var(phihat, na.rm = TRUE)
      varmat[folds, colsubset][1:2] = sigmasq/N

      datmat = as.data.frame(cbind(yi, yj, yi*yj, q1 - q12, q2 - q12, q12))
      datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))
      colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12")

      tmle = tmle(datmat = datmat, eps = eps, K = 2, ...)
      print(c(tmle$iterations, tmle$epsilon_error))
      if(tmle$error){
        warning("TMLE did not run or converge.")
        psiinvmat[folds,colsubset][3] = NA
        varmat[folds,colsubset][3] = NA
      }else{
        datmat = tmle$datmat
        q12 = pmax(datmat$q12, eps)
        q1 = pmin(datmat$q12 + datmat$q10, 1)
        q2 = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1), q12/q1)

        nuistmle[idfold == folds, paste(func, c("q12", "q1", "q2"), sep = '.')] = cbind(q12, q1, q2)

        gammainvhat = q1*q2/q12
        psiinvhat.tmle = mean(gammainvhat, na.rm = TRUE)

        phihat = gammainvhat*(yi/q1 + yj/q2 - yi*yj/q12) - psiinvhat.tmle

        Qnphihat = mean(phihat, na.rm = TRUE)
        cat(func, Qnphihat, '\n')
        psiinvmat[folds, colsubset][3] = psiinvhat.tmle
        sigmasq = var(phihat, na.rm = TRUE)
        varmat[folds, colsubset][3] = sigmasq/N
      }
    }
  }

  psiinv_summary[paste0(i, ",", j),] = colMeans(psiinvmat, na.rm = TRUE)
  var_summary[paste0(i, ",", j),] = colMeans(varmat, na.rm = TRUE)

  result <- list(psi = 1/psiinv_summary, sigma2 = N*var_summary, n = round(N*psiinv_summary),
                 varn = N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1), N = N,
                 ifvals = ifvals, nuis = nuis, nuistmle = nuistmle,
                 cin.l = round(pmax(N*psiinv_summary - 1.96*sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1)), N)),
                 cin.u = round(N*psiinv_summary + 1.96 *sqrt(N^2*var_summary + N*psiinv_summary*(psiinv_summary - 1))))
  class(result) <- "psinhat"
  return(result)
}
