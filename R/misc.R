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


