#' fit crr model
#'
#' Wrapper function to fit fine and gray competing risk model using function crr
#' from package cmprsk
#'
#' @param f formula for the model. Currently the formula only works by using the name
#' of the column in a dataframe. It does not work by using $ or [] notation.
#' @param data dataframe containing data
#' @keywords model
#' @export
crrRx<-function(f,data){
  k<-as.character(f)[3]
  covs<-removedollar(k)
  ff<-modelmatrix(f,data)
  m1<-cmprsk::crr(ff[[1]][,1],ff[[1]][,2],ff[[2]])
  m1$call<-paste("~",covs)
  return(m1)
}

#' fit box cox transformed linear model
#'
#' Wrapper function to fit fine and gray competing risk model using function crr
#' from package cmprsk
#'
#' @param f formula for the model. Currently the formula only works by using the name
#' of the column in a dataframe. It does not work by using $ or [] notation.
#' @param data dataframe containing data
#' @param lambda boolean indicating if you want to output the lamda used in the boxcox transformation. If so the function will return a list of length 2 with the model as the first element and a vector of length 2 as the second.
#' @keywords model
#' @export
boxcoxfitRx<-function(f,data,lambda=F){
  x<-as.character(f)[3]
  y<-as.character(f)[2]
  time<- gsub("\\s","",unlist(strsplit(y,"+",fixed=T))[1])
  covs<-removedollar(x)
  tempindexboxcoxfitRx<-seq_len(nrow(data))
  df1<-data.frame(tempindexboxcoxfitRx,data)
  f2<-as.formula(paste("tempindexboxcoxfitRx+",y,"~",x))
  temp<-modelmatrix(f2,df1)
  ff<-list(temp[[1]][,-1,drop=F],temp[[2,drop=F]])
  temp<-temp[[1]][,1,drop=F]
  lambda1<-unlist(unlist(geoR::boxcoxfit(ff[[1]],ff[[2]],lambda2=T))[1:2])
  ff[[1]]<-((ff[[1]]+lambda1[2])^lambda1[1]-1)/lambda1[1]
  df<-merge(df1,temp,by="tempindexboxcoxfitRx")[,-1,drop=F]
  df[,time]<-ff[[1]]
  out<-lm(f,data=df)
  out$call<-paste("~",covs)
  if(lambda)  return(list(out,lambda1))
  return(out)
}


# bwselect<-function(data,response,covs,strata=1,type,force=NULL,p=0.05,test=F,boxcox=F){
#   while(T){
#
#     if(type=="logistic"){
#       pvalues<-sapply(covs,function(cov){
#         m0<-glm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                 data=subset(data,!is.na(data[,cov])),family="binomial")
#         m1<-update(m0,as.formula(paste(".~.-",cov,sep="")))
#         1-pchisq(summary(m1)$deviance-summary(m0)$deviance,
#                  summary(m1)$df.residual-summary(m0)$df.residual)
#       })
#       if(test) print(pvalues)
#       if(length(covs)==1){
#         if(max(pvalues)>p){
#           return(glm(as.formula(paste(response,"~1",sep="")),
#                      data=data,family="binomial"))
#         }else{
#           return(glm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                      data=data,family="binomial"))
#         }}else{
#           if(max(pvalues)>p){
#             covs<-covs[-which.max(pvalues)]
#           }else{
#             return(glm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                        data=data,family="binomial"))
#           }}
#     }else if(type=="linear"){
#       pvalues<-sapply(covs,function(cov){
#         if(!boxcox){
#           m0<-lm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                  data=subset(data,!is.na(data[,cov])))
#           m1<-update(m0,as.formula(paste(".~.-",cov,sep="")))
#         }else{
#           m0<-boxcoxlm(data[,response],data[,covs])[[1]]
#           m1<-boxcoxlm(data[!is.na(data[,cov]),response],data[!is.na(data[,cov]),setdiff(covs,cov)])[[1]]
#         }
#         1-pchisq(2*(logLik(m0)-logLik(m1)),
#                  summary(m0)$df[1]-summary(m1)$df[1])
#       })
#       if(test) print(pvalues)
#       if(length(covs)==1){
#         if(max(pvalues)>p){
#           return(lm(as.formula(paste(response,"~1",sep="")),
#                     data=data))
#         }else{
#           return(lm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                     data=data))
#         }}else if(length(covs)==2 & boxcox & max(pvalues)>p){
#           return(boxcoxlm(data[,response],data[,covs[which.min(pvalues)]]))
#
#         }
#       else{
#         if(max(pvalues)>p){
#           covs<-covs[-which.max(pvalues)]
#         }else{
#           if(!boxcox){
#             return(lm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                       data=data))
#           }else{
#             return(boxcoxlm(data[,response],data[,covs]))
#           }
#         }}}
#
#   }}

# WE NEED TO BE ATTRITBUTE TO THE ORIGINAL AUTHOR IF THIS IS POSTED ON CRAN!
# Benjamin Schlegel, Marco Steenbergen
# https://cran.r-project.org/web/packages/brant/brant.pdf
# This is code from the brant package, with  this line of code:
# result.matrix = print.testresult(model, X2, df.v, by.var)
# replaced by the contents of brant:::print.testresult) to not print to screen
# Changed: LA 16 December 2020 return the result matrix, and the number of empty cells separately in a list
# the number of empty cells
modified_brant <- function (model, by.var = F)
{
  m_model <- model$call
  if (is.matrix(eval.parent(m_model$data)))
    m_model$data <- as.data.frame(data)
  m_model[-which(names(m_model) %in% c("", "formula", "data"))] <- NULL
  m_model[[1L]] <- quote(stats::model.frame)
  m_model <- eval.parent(m_model)
  Terms <- attr(m_model, "terms")
  x <- model.matrix(Terms, m_model)
  xint <- match("(Intercept)", colnames(x), nomatch = 0L)
  x <- x[, -xint, drop = FALSE]
  y <- as.numeric(model.response(m_model))
  x.variables = names(m_model)[-1]
  temp.data = data.frame(m_model, y)
  if (grepl(":", paste0(colnames(x), collapse = "")) & by.var) {
    by.var = FALSE
    warning("by.var = TRUE currently not supported for interactions, setting by.var to FALSE")
  }
  x.factors = c()
  for (name in x.variables) {
    if (!is.numeric(m_model[, name])) {
      x.factors = c(x.factors, name)
    }
  }
  if (length(x.factors) > 0) {
    tab = table(data.frame(temp.data$y, m_model[, x.factors]))
    count0 = sum(tab==0)
    # count0 = colSums(tab==0)
  }
  else {
    count0 = 0
  }
  J = max(y, na.rm = T)
  K = length(coef(model))
  for (m in 1:(J - 1)) {
    temp.data[[paste0("z", m)]] = ifelse(y > m, 1, 0)
  }
  binary.models = list()
  beta.hat = matrix(NA, nrow = J - 1, ncol = K + 1, byrow = T)
  var.hat = list()
  for (m in 1:(J - 1)) {
    mod = glm(paste0("z", m, " ~ ", as.character(formula(model)[3])),
              data = temp.data, family = "binomial")
    binary.models[[paste0("model", m)]] = mod
    beta.hat[m, ] = coef(mod)
    var.hat[[m]] = vcov(mod)
  }
  X = cbind(1, x)
  tau = matrix(model$zeta, nrow = 1, ncol = J - 1, byrow = T)
  pi.hat = matrix(NA, nrow = length(model$model[, 1]), ncol = J -
                    1, byrow = T)
  for (m in 1:(J - 1)) {
    pi.hat[, m] = binary.models[[m]]$fitted.values
  }
  varBeta = matrix(NA, nrow = (J - 1) * K, ncol = (J - 1) *
                     K)
  for (m in 1:(J - 2)) {
    for (l in (m + 1):(J - 1)) {
      Wml = Matrix::Diagonal(x = pi.hat[, l] - pi.hat[,
                                                      m] * pi.hat[, l])
      Wm = Matrix::Diagonal(x = pi.hat[, m] - pi.hat[,
                                                     m] * pi.hat[, m])
      Wl = Matrix::Diagonal(x = pi.hat[, l] - pi.hat[,
                                                     l] * pi.hat[, l])
      Xt = t(X)
      varBeta[((m - 1) * K + 1):(m * K), ((l - 1) * K +
                                            1):(l * K)] = as.matrix((solve(Xt %*% Wm %*%
                                                                             X) %*% (Xt %*% Wml %*% X) %*% solve(Xt %*% Wl %*%
                                                                                                                   X))[-1, -1])
      varBeta[((l - 1) * K + 1):(l * K), ((m - 1) * K +
                                            1):(m * K)] = varBeta[((m - 1) * K + 1):(m *
                                                                                       K), ((l - 1) * K + 1):(l * K)]
    }
  }
  betaStar = c()
  for (m in 1:(J - 1)) {
    betaStar = c(betaStar, beta.hat[m, -1])
  }
  for (m in 1:(J - 1)) {
    varBeta[((m - 1) * K + 1):(m * K), ((m - 1) * K + 1):(m *
                                                            K)] = var.hat[[m]][-1, -1]
  }
  I = diag(1, K)
  E0 = diag(0, K)
  for (i in 1:(J - 2)) {
    for (j in 1:(J - 1)) {
      if (j == 1) {
        temp = I
      }
      else if (j == i + 1) {
        temp = cbind(temp, -I)
      }
      else {
        temp = cbind(temp, E0)
      }
    }
    if (i == 1) {
      D = temp
    }
    else {
      D = rbind(D, temp)
    }
  }
  X2 = t(D %*% betaStar) %*% solve(D %*% varBeta %*% t(D)) %*%
    (D %*% betaStar)
  df.v = (J - 2) * K
  if (by.var) {
    combinations = getCombiCoefs(model)
    for (v in unique(combinations$var)) {
      k = subset(combinations, var == v)$i
      s = c()
      df.v.temp = 0
      for (e in k) {
        s = c(s, seq(from = e, to = K * (J - 1), by = K))
        df.v.temp = df.v.temp + J - 2
      }
      s = sort(s)
      Ds = D[, s]
      if (!is.null(dim(Ds))) {
        Ds = Ds[which(!apply(Ds == 0, 1, all)), ]
      }
      if (!is.null(dim(Ds)))
        X2 = c(X2, t(Ds %*% betaStar[s]) %*% solve(Ds %*%
                                                     varBeta[s, s] %*% t(Ds)) %*% (Ds %*% betaStar[s]))
      else X2 = c(X2, t(Ds %*% betaStar[s]) %*% solve(Ds %*%
                                                        varBeta[s, s] %*% t(t(Ds))) %*% (Ds %*% betaStar[s]))
      df.v = c(df.v, df.v.temp)
    }
  }
  else {
    for (k in 1:K) {
      s = seq(from = k, to = K * (J - 1), by = K)
      Ds = D[, s]
      if (!is.null(dim(Ds))) {
        Ds = Ds[which(!apply(Ds == 0, 1, all)), ]
      }
      if (!is.null(dim(Ds)))
        X2 = c(X2, t(Ds %*% betaStar[s]) %*% solve(Ds %*%
                                                     varBeta[s, s] %*% t(Ds)) %*% (Ds %*% betaStar[s]))
      else X2 = c(X2, t(Ds %*% betaStar[s]) %*% solve(Ds %*%
                                                        varBeta[s, s] %*% t(t(Ds))) %*% (Ds %*% betaStar[s]))
      df.v = c(df.v, J - 2)
    }
  }
  # result.matrix = print.testresult(model, X2, df.v, by.var)
  p.values = pchisq(X2, df.v, lower.tail = FALSE)
  if (by.var) {
    var.names = unlist(strsplit(as.character(formula(model))[3],
                                split = " \\+ "))
  }
  else {
    var.names = names(coef(model))
  }
  longest.char = max(nchar(var.names))
  n.tabs = ceiling(longest.char/7)
  n.tabs = ifelse(n.tabs < 2, 2, n.tabs)
  result.matrix = matrix(c(X2, df.v, p.values), ncol = 3)
  rownames(result.matrix) = c("Omnibus", var.names)
  colnames(result.matrix) = c("X2", "df", "probability")
  # invisible(result.matrix)
  # Changed: LA 16 December 2020 return the result matrix, and the number of empty cells separately in a list
  list(result = result.matrix,zero_count_cells=count0)
}
