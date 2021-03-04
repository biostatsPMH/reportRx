#' Plot KM curve
#'
#'This function will plot a KM curve with possible stratification. You can specifyif you want
#'a legend or confidence bands as well as the units of time used.
#'
#' @param data dataframe containing your data
#' @param response character vector with names of columns to use for response
#' @param group string specifiying the column name of stratification variable
#' @param pos what position you want the legend to be. Current option are bottomleft and topright
#' @param units string specifying what the unit of time is use lower case and plural
#' @param CI boolean to specify if you want confidence intervals
#' @param legend boolean to specify if you want a legend
#' @param title title of plot
#' @importFrom graphics axis legend lines mtext par plot
#' @keywords plot
#' @export
#' @examples
#' require(survival)
#' lung$sex<-factor(lung$sex)
#' plotkm(lung,c("time","status"))
#' plotkm(lung,c("time","status"),"sex")
plotkm<-function(data,response,group=1,pos="bottomleft",units="months",CI=F,legend=T,title=""){
  if(class(group)=="numeric"){
    kfit<-survfit(as.formula(paste("Surv(",response[1],",",response[2],")~1",sep="")),data=data)
    sk<-summary(kfit)$table
    levelnames<-paste("N=",sk[1],",Events=",sk[4]," (",round(sk[4]/sk[1],2)*100,"%)",sep="")
    if(title=="")  title<-paste("KM-Curve for ",nicename(response[2]),sep="")

  }else if(length(group)>1){
    return("Currently you can only stratify by 1 variable")
  }else{
    if(class(data[,group])!="factor")
      stop("group must be a vactor variable. (Or leave unspecified for no group)")
    lr<-survdiff(as.formula(paste("Surv(",response[1],",",response[2],")~",paste(group,collapse="+"),sep="")),data=data)
    lrpv<-1-pchisq(lr$chisq,length(lr$n)- 1)
    levelnames<-levels(data[,group])
    kfit<-survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",paste(group,collapse="+"),sep="")),data=data)
    if(title=="") title<-paste("KM-Curve for ",nicename(response[2])," stratified by ",nicename(group),sep="")
    levelnames<-sapply(1:length(levelnames),function(x){paste(levelnames[x]," n=",lr$n[x],sep="")})

  }


  plot(kfit,mark.time=T,lty=1:length(levelnames),xlab=paste("Time (",cap(units),")",sep=""),
       ylab="Suvival Probability ",cex=1.1,conf.int=CI,
       main=title)


  if(legend){
    if(class(group)=="numeric"){legend(pos,levelnames,lty=1:length(levelnames),bty="n")
    }else{ legend(pos,c(levelnames,paste("p-value=",pvalue(lrpv)," (Log Rank)",sep="")),
                  col=c(rep(1,length(levelnames)),"white"),lty=1:(length(levelnames)+1),bty="n")}
  }
}

#'Get event time summary dataframe
#'
#'This function will output a dataframe with usefull summary statistics from a coxph model
#'
#'@param data dataframe containing data
#'@param response character vector with names of columns to use for response
#'@param group string specifiying the column name of stratification variable
#'@param times numeric vector of times you want survival time provbabilities for.
#'@keywords dataframe
#'@export
#'@examples
#'require(survival)
#'lung$sex<-factor(lung$sex)
#'etsum(lung,c("time","status"),"sex")
#'etsum(lung,c("time","status"))
#'etsum(lung,c("time","status"),"sex",c(1,2,3))
etsum<- function(data,response,group=1,times=c(12,24)){
  if(class(group)=="numeric"){
    kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep=""))  ,data=data))
    maxtime=max(kfit$time)
    times[times>maxtime]=maxtime
    kfit2<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep="")) ,data=data),times=times)
    tab<-as.data.frame(cbind(strata=as.character(kfit2$strata),times=kfit2$time,SR=paste(round(kfit2$surv*100,0)," (",round(kfit2$lower*100,0),"-",round(kfit2$upper*100,0),")",sep="")))
    tbl<-kfit2$table
  }else{
    if(class(data[,group])!="factor")
      stop("group variable must be factor or leave unspecified for no group")
    tab<-lapply(levels(data[,group]),function(level){
      subdata<-subset(data,data[,group]==level)
      kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",1,sep=""))  ,data=subdata))
      maxtime=max(kfit$time)
      times[times>maxtime]=maxtime
      kfit2<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",1,sep="")) ,data=subdata),times=times)
      list(cbind(strata=paste0(group,"=",level),times=kfit2$time,SR=paste(round(kfit2$surv*100,0)," (",round(kfit2$lower*100,0),"-",round(kfit2$upper*100,0),")",sep="")),kfit2$table)})
    tbl=t(sapply(tab,"[[",2))
    rownames(tbl)=sapply(levels(data[,group]),function(level)paste0(group,"=",level))
    tab=do.call(rbind.data.frame,lapply(tab,"[[",1))
  }

  if(class(group)!="numeric"){
    kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep=""))  ,data=data))
    med=by(data,data[,group],function(x) median(x[,response[1]],na.rm=T))
    min=by(data,data[,group],function(x) min(x[,response[1]],na.rm=T))
    max=by(data,data[,group],function(x) max(x[,response[1]],na.rm=T))
    survtimes<-data.frame(strata=as.character(kfit$strata),kfit$time)
    minst<-round(as.numeric(by(survtimes,survtimes$strata,function(x) min (x[,2]))),1)
    maxst<-round(as.numeric(by(survtimes,survtimes$strata,function(x) max (x[,2]))),1)
    tab<-reshape::cast(tab,strata ~ times)
    names<-names(tab)
    tab<-data.frame(tab)
    names(tab)<-names
    tab[,1]<-levels(data[,group])
    if(length(times)>1){
      indx<-c(0,sapply(sort(as.numeric(names(tab)[-1])),function(x){which(as.numeric(names(tab)[-1])==x)}))+1
      tab<-tab[,indx]
      tab<-tab[c(2:length(tab),1)]
    }else{
      tab<-tab[c(2:length(tab),1)]
    }
    noeventsindx<-ifelse(length(which(tbl[,4]==0))!=0,
                         which(tbl[,4]==0),NA)
    if(!is.na(noeventsindx)){
      for(i in noeventsindx){
        if(i==1){
          minst<-c(0,minst)
          maxst<-c(0,maxst)
        }else if(i>length(minst)){
          minst<-c(minst,0)
          maxst<-c(maxst,0)
        }else{
          minst<-c(minst[1:i-1],0,minst[i:length(minst)])
          maxst<-c(maxst[1:i-1],0,maxst[i:length(maxst)])
        }}}


    tab<-cbind("n"=tbl[,1],"Events"=tbl[,4],"MedKM"=round(tbl[,5],1),
               "LCI"=round(tbl[,6],1),"UCI"=round(tbl[,7],1),
               "MedFU"=round(as.numeric(med),1),
               "MinFU"=round(as.numeric(min),1),"MaxFU"=round(as.numeric(max),1),
               "MinET"=minst,"MaxET"=maxst,tab)
    rownames(tab)<-NULL
  }else{
    med=median(data[,response[1]],na.rm=T)
    min=min(data[,response[1]],na.rm=T)
    max=max(data[,response[1]],na.rm=T)
    if(length(times)>1){
      tab<-data.frame(t(tab))
      rownames(tab)<-NULL
      names(tab)<-as.numeric(as.matrix(tab[1,]))
      tab<-tab[-1,]
    }else{
      rownames(tab)<-NULL
      names(tab)[2]<-times
      tab<-tab[-1]
    }
    tab<-cbind("n"=tbl[1],"Events"=tbl[4],"MedKM"=round(tbl[5],1),"LCI"=round(tbl[6],1),"UCI"=round(tbl[7],1),
               "MedFU"=round(as.numeric(med),1),"MinFU"=round(as.numeric(min),1),"MaxFU"=round(as.numeric(max),1),
               "MinET"=round(min(kfit$time),1),"MaxET"=round(max(kfit$time),1),tab)
    rownames(tab)<-NULL
  }
  return(tab)
}

#'Print LaTeX event time summary
#'
#'Wrapper for the etsum function that prints paragraphs of text in LaTeX
#'
#'@param data dataframe containing data
#'@param response character vector with names of columns to use for response
#'@param group string specifiying the column name of stratification variable
#'@param times numeric vector of times you want survival time provbabilities for.
#'@param units string indicating the unit of time. Use lower case and plural.
#'@keywords print
#'@export
#'@examples
#'require(survival)
#'lung$sex<-factor(lung$sex)
#'petsum(lung,c("time","status"),"sex")
#'petsum(lung,c("time","status"))
#'petsum(lung,c("time","status"),"sex",c(1,2,3),"months")
petsum<-function(data,response,group=1,times=c(12,14),units="months"){
  t<-etsum(data,response,group,times)

  #plotkm(nona,response,group)

  names<-names(t)
  if("strata"%in% names){
    strta<-sapply(t[,"strata"],function(x) paste(x,": ",sep=""))
    offset<-2
    ofst<-1
  }else{
    strta=matrix(c("",""))
    offset<-1
    ofst<-0
  }


  out<-sapply(seq_len(nrow(t)),function(i){

    if(is.na(t[i,3])) {km<-paste("The KM median event time has not been achieved due to lack of events.",sep="")
    }else if (!is.na(t[i,5])){km<-paste("The KM median event time is ",t[i,3]," with 95",sanitizestr("%")," confidence Interval (",t[i,4],",",t[i,5],").",sep="")
    }else{km<-paste("The KM median event time is ",t[i,3]," ",units," with 95",sanitizestr("%")," confidence Interval (",t[i,4],",",t[i,10],").",sep="")}

    # if at least one event
    if(t[i,2]!=0){
      flet<-paste(" The first and last event times occurred at ",t[i,9],
                  " and ",t[i,10]," ",units," respectively. ",sep="")

      psindex=11:(ncol(t)-ofst)
      psindex=psindex[which(!is.na(t[i,psindex]))]
      if(length(psindex)>1){
        lastindex=psindex[length(psindex)]
        firstindex=psindex[-length(psindex)]
        ps<-paste("The ",paste(names[firstindex],collapse=",")," and ",names[lastindex]," " ,substring(units,1,nchar(units)-1),
                  " probabilities of 'survival' and their 95",sanitizestr("%")," confidence intervals are ",
                  paste(sapply(t[i,firstindex],function(x) paste(x)),collapse=",")," and ",t[i,lastindex]," percent.",sep="")

      }else{
        ps<-paste("The ",names[psindex]," ",substring(units,1,nchar(units)-1),
                  " probability of 'survival' and 95",sanitizestr("%")," confidence interval is ",
                  t[i,psindex]," percent.",sep="")
      }
      #if no events
    }else{
      km=""
      ps=""
      flet=""
    }


    out<-paste(lbld(sanitizestr(nicename(strta[i])))," There are ",t[i,1]," patients. There were ",t[i,2],
               " (",round(100*t[i,2]/t[i,1],0),sanitizestr("%"),") events. The median and range of the follow-up times is ",
               t[i,6]," (",t[i,7],"-",t[i,8],") ",units,". ",km,flet,ps,sep="")
    cat("\n",out,"\n")
  })
}

#'Get covariate summary dataframe
#'
#'Returns a dataframe corresponding to a descriptive table
#'
#'@param data dataframe containing data
#'@param covs character vector with the names of columns to include in table
#'@param maincov covariate to stratify table by
#'@param numobs named list overriding the number of people you expect to have the covariate
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param IQR boolean indicating if you want to display the inter quantile range (Q1,Q3) as opposed to (min,max) in the summary for continuous variables
#'@param digits.cat number of digits for the proportions when summarizing categorical data (default: 0)
#'@param testcont test of choice for continuous variables,one of \emph{rank-sum} (default) or \emph{ANOVA}
#'@param testcat test of choice for categorical variables,one of \emph{Chi-squared} (default) or \emph{Fisher}
#'@keywords dataframe
#'@export
#'@seealso \code{\link{fisher.test}},\code{\link{chisq.test}},\code{\link{wilcox.test}},\code{\link{kruskal.test}},and \code{\link{anova}}
covsum<-function(data,covs,maincov=NULL,numobs=NULL,markup=TRUE,sanitize=TRUE,nicenames=TRUE,IQR = FALSE,digits.cat = 0,
                 testcont = c('rank-sum test','ANOVA'),testcat = c('Chi-squared','Fisher')){
  # New LA 18 Feb, test for presence of variables in data and convert character to factor
  missing_vars = setdiff(covs,names(data))
  if (length(missing_vars)>0){  stop(paste('These covariates are not in the data:',missing_vars))  }
  for (v in c(maincov,covs)) if (class(data[[v]])[1]=='character') data[[v]] <- factor(data[[v]])

  testcont <- match.arg(testcont)
  testcat <- match.arg(testcat)
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity
  }
  digits.cat <- as.integer(digits.cat)
  if( digits.cat<0 ) stop("parameter 'digits.cat' cannot be negative!")
  if(!sanitize) sanitizestr<-identity
  if(!nicenames) nicename<-identity
  if(!is.null(maincov)){
    levels<-names(table(data[[maincov]]))
    levels<-c(list(levels),as.list(levels))
  }else{
    levels<-"NOMAINCOVNULLNA"
  }
  N=nrow(data)
  if(!is.null(maincov)){
    nmaincov<-c(sum(table(data[[maincov]])),table(data[[maincov]]))
  }else{
    nmaincov<-N
    p<-NULL
  }
  out<-lapply(covs,function(cov){
    ismiss=F
    n<-sum(table(data[[cov]]))

    #Set up the first coulmn
    factornames<-NULL
    if(is.null(numobs[[cov]]))  numobs[[cov]]<-nmaincov
    if(numobs[[cov]][1]-n>0) {ismiss=T
    factornames<-c(factornames,"Missing")
    }
    #if the covariate is a factor
    if(is.factor(data[[cov]])){
      factornames<-c(levels(data[[cov]]),factornames)
      if (!is.null(maincov)) {
        p <- if( testcat=='Fisher'){ try(fisher.test(data[[maincov]],data[[cov]])$p.value,silent=T)
        } else try(chisq.test(data[[maincov]],data[[cov]])$p.value,silent=T)
        if (class(p) == "try-error")
          p <- NA
        p <- lpvalue(p)
      }


      #set up the main columns
      onetbl<-mapply(function(sublevel,N){
        missing<-NULL
        if(sublevel[1]!="NOMAINCOVNULLNA"){
          subdata<-subset(data,subset=data[[maincov]]%in%sublevel)
        }else{
          subdata<-data
        }
        table<-table(subdata[[cov]])
        tbl<-table(subdata[[cov]])
        n<-sum(tbl)
        prop <- round(tbl/n,2+digits.cat)*100
        prop <- sapply(prop,function(x){if(!is.nan(x)){x} else{0}})
        prop.fmt <- sprintf(paste0("%.",digits.cat,"f"),prop)
        tbl<-mapply(function(num,prop){paste(num," (",prop,")",sep="")},tbl,prop.fmt)
        if(ismiss) missing<-N-n
        tbl<-c(tbl,lbld(missing))
        return(tbl)
      },levels,numobs[[cov]])

      #if the covariate is not a factor
    }else{
      #setup the first column
      factornames <- c("Mean (sd)",ifelse(IQR,"Median (Q1,Q3)","Median (Min,Max)"),factornames)
      if (!is.null(maincov)) {
        p <- if( testcont=='rank-sum test'){
          if( length(unique(data[[maincov]]))==2 ){
            try( wilcox.test(data[[cov]] ~ data[[maincov]])$p.value )
          } else try( kruskal.test(data[[cov]] ~ data[[maincov]])$p.value )
        } else try(anova(lm(data[[cov]] ~ data[[maincov]]))[5][[1]][1])
        if (class(p) == "try-error")
          p <- NA
        p <- lpvalue(p)
      }
      #set up the main columns
      onetbl <- mapply(function(sublevel,N){
        missing <- NULL
        if(sublevel[1]!="NOMAINCOVNULLNA"){
          subdata<-subset(data,subset=data[[maincov]]%in%sublevel)
        }else{subdata<-data}
        summary<-round(summary(subdata[[cov]]),1)
        meansd<-paste(summary[4]," (",round(sd(subdata[[cov]],na.rm=T),1),")",sep="")
        mmm <- if( IQR ){ paste(summary[3]," (",summary[2],",",summary[5],")",sep = "")
        }else paste(summary[3]," (",summary[1],",",summary[6],")",sep = "")

        #if there is a missing in the whole data
        if(ismiss){
          n<-sum(table(subdata[[cov]]))
          missing<-N-n
        }
        tbl<-c(meansd,mmm,lbld(missing))

        return(tbl)}
        ,levels,numobs[[cov]])}

    #Add the first column to the main columns and get the matrix ready for later
    factornames<-addspace(sanitizestr(nicename(factornames)))
    # LA Added 20 Jan 2021 to deal with one-level factors
    if (is.null(nrow(onetbl))){onetbl <- matrix(data=onetbl,ncol=length(onetbl),nrow=1) }

    onetbl<-cbind(factornames,onetbl)

    if(!is.null(maincov)){
      onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),rep("",length(levels[[1]])+1)),onetbl)
      onetbl<-cbind(onetbl,c(p,rep("",nrow(onetbl)-1)))
    }else{
      onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),""),onetbl)
    }
    rownames(onetbl)<-NULL
    colnames(onetbl)<-NULL
    return(onetbl)})
  table <- do.call("rbind",lapply(out,data.frame,stringsAsFactors = FALSE))
  ### unlist each column of the table
  table <- data.frame(apply(table,2,unlist),stringsAsFactors = FALSE)
  rownames(table)<-NULL
  if(!is.null(maincov)){
    colnames(table)<-c("Covariate",paste("Full Sample (n=",N,")",sep=""),
                       mapply(function(x,y){paste(x," (n=",y,")",sep="")},
                              names(table(data[[maincov]])),table(data[[maincov]])),"p-value")
  }else{
    colnames(table)<-c("Covariate",paste("n=",N,sep=""))

  }
  colnames(table)<-sanitizestr(colnames(table))
  return(table)
}

#'Print covariate summary Latex
#'
#'Returns a dataframe corresponding to a descriptive table
#'
#'@param data dataframe containing data
#'@param covs character vector with the names of columns to include in table
#'@param maincov covariate to stratify table by
#'@param TeX boolean indicating if you want to be able to view extra long tables in the LaTeX pdf. If TeX is T then the table will not convert properly to docx
#'@param ... additional options passed to function  \code{\link{covsum}}
#'@keywords print
#'@export
pcovsum<-function(data,covs,maincov=NULL,TeX=FALSE,...){
  if(!TeX){
    print.xtable(xtable(covsum(data,covs,maincov,...)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
  }else{

    print.xtable(xtable(covsum(data,covs,maincov,...)),include.rownames=F,sanitize.text.function=identity,table.placement="H",floating=FALSE,tabular.environment="longtable")
  }}

#'Get univariate summary dataframe
#'
#'Returns a dataframe corresponding to a univariate table
#'
#'@param response string vector with name of response
#'@param covs character vector with the names of columns to fit univariate models to
#'@param data dataframe containing data
#'@param type string indicating he type of univariate model to fit. The function will try and guess what type you want based on your response. If you want to override this you can manually specify the type.
#'Options in clude "linear","logistic","coxph","crr","boxcox","logistic"
#'@param strata character vector of covariates to stratify by. Only used for coxph and crr
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param CIwidth width of confidence interval, default is 0.95
#'@param testing boolean to indicate if you want to print out the covariates before the model fits.
#'This will allow you to see which model is not fitting if the function throws an error
#'@keywords dataframe
#'@export
uvsum<-function(response,covs,data,type=NULL,strata=1,markup=T,sanitize=T,nicenames=T,CIwidth=0.95,testing=F){
  # New LA 24 Feb, test for presence of variables in data and convert character to factor
  missing_vars = setdiff(c(response,covs),names(data))
  if (length(missing_vars)>0){
    stop(paste('These covariates are not in the data:',missing_vars))
  }
  for (v in c(response,covs)) if (class(data[[v]])[1]=='character') data[[v]] <- factor(data[[v]])

  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity
  }
  if(!sanitize) sanitizestr<-identity
  if(!nicenames) nicename<-identity

  if(class(strata)!="numeric") {strata<-sapply(strata,function(stra){paste("strata(",stra,")",sep="")})
  }else{strata<-""}

  if(!is.null(type)){
    if(type=="logistic"){beta<-"OR"
    }else if(type=="linear"|type=="boxcox"){beta<-"Estimate"
    }else if(type=="coxph"|type=="crr"){beta<-"HR"
    }else{stop("type must be either coxph,logisitc,linear,coxbox,crr (or NULL)")
    }}else
    {if(length(response)==2) {
      if(length(unique(data[,response[2]]))<3){type<-"coxph"
      }else{type<-"crr"}
      beta<-"HR"
    }else if (length(unique(data[,response] ))==2) {type<-"logistic"
    beta<-"OR"
    }else {type<-"linear"
    beta<-"Estimate"
    }
    }
  beta = betaWithCI(beta,CIwidth)

  if(strata!="" & type!="coxph") stop("strata can only be used with coxph")



  out<-lapply(covs,function(cov){
    cov2<-cov
    if(testing) print(cov)
    if(is.factor(data[[cov]])){
      levelnames<-sapply(sapply(sapply(levels(factor(data[[cov]])),nicename),sanitizestr),addspace)
      cov<-lbld(sanitizestr(nicename(cov)))
      title<-NULL
      body<-NULL
      if(type=="coxph"){
        m2<-coxph(as.formula(paste(paste("Surv(",response[1],",",response[2],")",sep=""),"~",cov2,ifelse(strata=="","","+"),paste(strata,collapse="+"),sep="")),data=data)
        hazardratio<-c("Reference",apply(matrix(summary(m2,conf.int=CIwidth)$conf.int[,c(1,3,4)],ncol=3),1,psthr))
        pvalue<-c("",sapply(summary(m2)$coef[,5],lpvalue))
        title<-c(cov,"","",lpvalue(summary(m2,conf.int=CIwidth)$waldtest[3]))
      }else if(type=="crr"){
        m2<-crrRx(as.formula(paste(paste(response,collapse="+"),"~",cov2,sep="")),data=data)
        hazardratio<-c("Reference",apply(matrix(summary(m2,conf.int=CIwidth)$conf.int[,c(1,3,4)],ncol=3),1,psthr))
        pvalue<-c("",sapply(summary(m2)$coef[,5],lpvalue))
        globalpvalue<-try(aod::wald.test(b=m2$coef,Sigma=m2$var,Terms=seq_len(length(m2$coef)))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        title<-c(cov,"","",lpvalue(globalpvalue))

      }else if(type=="logistic"){
        m2<-glm(as.formula(paste(response,"~",cov2,sep="")),family="binomial",data=data)
        #globalpvalue<-1-pchisq(2*(summary(m2)$null.deviance-summary(m2)$deviance),summary(m2)$df.null-summary(m2)$df.residual)
        globalpvalue<-try(aod::wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"

        m<-summary(m2)$coefficients
        Z_mult = qnorm(1-(1-CIwidth)/2)
        hazardratio <- c("Reference", apply(cbind(exp(m[-1,1]), exp(m[-1, 1] - Z_mult * m[-1, 2]), exp(m[-1,1] + Z_mult * m[-1, 2])), 1, psthr))
        pvalue <- c("", sapply(m[-1, 4], lpvalue))
      }else if(type=="linear"|type=="boxcox"){

        if(type=="linear"){m2<-lm(as.formula(paste(response,"~",cov2,sep="")),data=data)
        }else{m2<-boxcoxfitRx(as.formula(paste(response,"~",cov2,sep="")),data=data)}
        m<-summary(m2)$coefficients
        #globalpvalue<-anova(m2)[5][[1]][1])
        globalpvalue<-try(aod::wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"

        T_mult = stats::qt(1-(1-CIwidth)/2,m2$df.residual)
        hazardratio <- c("Reference", apply(cbind(m[-1,1], m[-1, 1] - T_mult * m[-1, 2], m[-1, 1] +T_mult * m[-1, 2]), 1, psthr))
        pvalue<-c("",sapply(m[-1,4],lpvalue))
        title<-c(cov,"","",lpvalue(globalpvalue))
      }
      if(length(levelnames)==2){
        body<-cbind(levelnames,hazardratio,c("",""),c("",""))
      }else{
        body<-cbind(levelnames,hazardratio,pvalue,rep("",length(levelnames)))
      }
      out<-rbind(title,body)
      rownames(out)<-NULL
      colnames(out)<-NULL
      return(list(out,nrow(out)))
    }else
    {
      cov<-lbld(sanitizestr(nicename(cov)))
      if(type=="coxph"){
        m2<-coxph(as.formula(paste(paste("Surv(",response[1],",",response[2],")",sep=""),"~",cov2,ifelse(strata=="","","+"),paste(strata,collapse="+"),sep="")),data=data)
        out<-matrix(c(cov,psthr(summary(m2,conf.int=CIwidth)$conf.int[,c(1,3,4)]),"",lpvalue(summary(m2)$waldtest[3])),ncol=4)

      }else if(type=="crr"){

        m2<-crrRx(as.formula(paste(paste(response,collapse="+"),"~",cov2,sep="")),data=data)
        globalpvalue<-try(aod::wald.test(b=m2$coef,Sigma=m2$var,Terms=seq_len(length(m2$coef)))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        out<-matrix(c(cov,psthr(summary(m2,conf.int=CIwidth)$conf.int[,c(1,3,4)]),"",lpvalue(globalpvalue)),ncol=4)

      }else if(type=="logistic"){
        m2<-glm(as.formula(paste(response,"~",cov2,sep="")),family="binomial",data=data)

        m<-summary(m2)$coefficients

        #globalpvalue<-1-pchisq(2*(summary(m2)$null.deviance-summary(m2)$deviance),summary(m2)$df.null-summary(m2)$df.residual)
        globalpvalue<-try(aod::wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"

        Z_mult = qnorm(1-(1-CIwidth)/2)
        out <- matrix(c(cov, psthr(c(exp(m[-1, 1]), exp(m[-1,1] - Z_mult * m[-1, 2]), exp(m[-1, 1] + Z_mult *m[-1, 2]))), "", lpvalue(globalpvalue)), ncol = 4)


      }else if(type=="linear"|type=="boxcox"){
        if(type=="linear"){m2<-lm(as.formula(paste(response,"~",cov2,sep="")),data=data)
        }else{m2<-boxcoxfitRx(as.formula(paste(response,"~",cov2,sep="")),data=data)}
        #globalpvalue<-anova(m2)[5][[1]][1])
        globalpvalue<-try(aod::wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"

        m<-summary(m2)$coefficients

        T_mult = stats::qt(1-(1-CIwidth)/2,m2$df.residual)
        out <- matrix(c(cov, psthr(c(m[-1, 1], m[-1,1] - T_mult * m[-1, 2], m[-1, 1] + T_mult * m[-1,2])), "", lpvalue(globalpvalue)), ncol = 4)

      }
      return(list(out,nrow(out)))}})
  table<-lapply(out,function(x){return(x[[1]])})
  table<-do.call("rbind",lapply(table,data.frame,stringsAsFactors = FALSE))

  colnames(table)<-sapply(c("Covariate",sanitizestr(beta),"p-value","Global p-value"),lbld)
  return(table)
}

#'Print univariate summary LaTeX table
#'
#'Returns a LaTeX table of the univariate summary
#'
#'@param response string vector with name of response
#'@param covs character vector with the names of columns to fit univariate models to
#'@param data dataframe containing data
#'@param type string indicating he type of univariate model to fit. The function will try and guess what type you want based on your response. If you want to override this you can manually specify the type. Options in clude "linear","logistic","coxph","crr","boxcox","logistic"
#'@param strata character vector of covariates to stratify by. Only used for coxph and crr
#'@param TeX boolean indicating if you want to be able to view extra long tables in the LaTeX pdf. If TeX is T then the table will not convert properly to docx
#'@importFrom stats anova as.formula chisq.test coef fisher.test glm kruskal.test lm median model.matrix pchisq qnorm sd time vcov wilcox.test
#'@importFrom xtable xtable print.xtable
#'@keywords dataframe
#'@export
puvsum<-function(response,covs,data,type=NULL,strata=1,TeX=F){
  if(!TeX){
    print.xtable(xtable(uvsum(response,covs,data,type,strata)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
  }else{
    print.xtable(xtable(uvsum(response,covs,data,type,strata)),include.rownames=F,sanitize.text.function=identity,table.placement="H",floating=FALSE,tabular.environment="longtable")
  }

}

# TODO: Add support for svyglm and svycoxph functions. May need to check the feasibility of the global p-value here
#'Get multivariate summary dataframe
#'
#'Returns a dataframe corresponding to a univariate table
#'
#'@param model fitted model object
#'@param data dataframe containing data
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@keywords dataframe
#'@export
mvsum<-function(model,data,markup=T,sanitize=T,nicenames=T){
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity
  }
  if(!sanitize) sanitizestr<-identity
  if(!nicenames) nicename<-identity

  call<-paste(deparse(summary(model)$call),collapse="")
  call<-unlist(strsplit(call,"~",fixed=T))[2]
  call<-unlist(strsplit(call,",",fixed=T))[1]
  if(substr(call,nchar(call),nchar(call))=="\"") call<-substr(call,1,nchar(call)-1)
  call<-unlist(strsplit(call,"\"",fixed=T))[1]
  call<-unlist(strsplit(call,"+",fixed=T))
  call<-unlist(strsplit(call,"*",fixed=T))
  call<-unique(call)
  call<-call[which(is.na(sapply(call,function(cov){charmatch("strata(",cov)}))==T)]
  call<-gsub("\\s","",call)
  type<-class(model)[1]

  if(type!="lme"){
    betanames<-attributes(summary(model)$coef)$dimnames[[1]]
  }else{
    betanames<-names(model$coef$fixed)
  }



  if(type=="glm"){ beta<-"OR(95%CI)"
  betanames<-betanames[-1]
  }else if(type=="lm"|type=="lm"|type=="lme"){ beta<-"Estimate(95%CI)"
  betanames<-betanames[-1]
  }else if(type=="coxph"|type=="crr"){beta<-"HR(95%CI)"
  }else{ stop("type must be either coxph,logistic,lm,crr,lme (or NULL)")}

  ucall=unique(call)

  #   indx<-as.vector(sapply(betanames,function(string){
  #
  #     indx<-which(sapply(call,function(cov){charmatch(cov,string)})==1)
  #     if(length(indx)==1) return(indx)
  #     #If one  facorname is a subset of another
  #     indx2<-which.max(sapply(call[indx],nchar))
  #     if(length(indx2)==1) return(indx[indx2])
  #     indx3<-which(sapply(call[indx2],function(c){substr(betaname,1,nchar(c))==c}))
  #     if(length(indx3)==1)  return(call[indx[indx2[indx3]]])
  #     return(-1)
  #   }))

  indx=matchcovariate(betanames,ucall)
  if(min(indx)==-1) stop("Factor name + level name is the same as another factor name. Please change. Will fix this issue later")





  y<-betaindx(indx)
  if(type%in%c("lm","glm","lm","lme")){
    y<-lapply(y,function(x){ x+1})
    betanames<-c("intercept",betanames)
  }

  out<-lapply(y,function(covariateindex){

    #Get attribute names and split by ineractions
    betaname<-betanames[covariateindex]
    betaname<-strsplit(betaname,":",fixed=T)
    #get the covariate names
    oldcovname<-covnm(betaname[[1]],call)

    #get the levelnames
    levelnames<-unlist(lapply(betaname,function(level){
      paste(mapply(function(lvl,cn){result<-unlist(strsplit(lvl,cn,fixed=T))[2]
      out<-ifelse(is.na(result),cn,result)},level,oldcovname),collapse=":")}))
    levelnames<-addspace(sanitizestr(nicename(levelnames)))
    covariatename<-lbld(sanitizestr(nicename(paste(oldcovname,collapse=":"))))
    reference=NULL
    title=NULL
    body=NULL
    if(type=="lme"){
      globalpvalue<-try(aod::wald.test(b=model$coef$fixed[covariateindex],Sigma=vcov(model)[covariateindex,covariateindex],Terms=seq_along(covariateindex))$result$chi2[3]);
    }
    else if(type!="crr"){
      globalpvalue<-try(aod::wald.test(b=coef(model)[covariateindex],Sigma=vcov(model)[covariateindex,covariateindex],Terms=seq_along(covariateindex))$result$chi2[3]);
    }else{
      globalpvalue<-try(aod::wald.test(b=model$coef[covariateindex],Sigma=model$var[covariateindex,covariateindex],Terms=seq_along(covariateindex))$result$chi2[3]);

    }
    if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
    globalpvalue<-lpvalue(globalpvalue)

    if(type=="coxph"|type=="crr"){
      hazardratio<-c(apply(matrix(summary(model)$conf.int[covariateindex,c(1,3,4)],ncol=3),1,psthr))
      pvalues<-c(sapply(summary(model)$coef[covariateindex,5],lpvalue))

    }else if (type=="glm"){
      m<-summary(model)$coefficients
      hazardratio<-apply(cbind(exp(m[covariateindex,1]),exp(m[covariateindex,1]-1.96*m[covariateindex,2]),
                               exp(m[covariateindex,1]+1.96*m[covariateindex,2])),1,psthr)
      pvalues<-c(sapply(m[covariateindex,4],lpvalue))
    }else if (type=="lm"|type=="lm"){
      m<-summary(model)$coefficients

      hazardratio<-apply(cbind(m[covariateindex,1],m[covariateindex,1]-1.96*m[covariateindex,2],
                               m[covariateindex,1]+1.96*m[covariateindex,2]),1,psthr)
      pvalues<-sapply(m[covariateindex,4],lpvalue)
    }else if(type=="lme"){
      m<-summary(model)$tTable
      hazardratio<-apply(cbind(m[covariateindex,1],m[covariateindex,1]-1.96*m[covariateindex,2],
                               m[covariateindex,1]+1.96*m[covariateindex,2]),1,psthr)

      pvalues<-c(sapply(m[covariateindex,5],lpvalue))
    }

    #if not interaction

    if(length(betaname[[1]])==1){
      #if cts
      if(!is.factor(data[,oldcovname])){

        title<-c(nicename(covariatename),hazardratio,"",globalpvalue)
      }else if(length(levelnames)==1){
        title<-c(covariatename,"","",globalpvalue)
        if(!is.null(data)) reference<-c(addspace(sanitizestr(names(table(data[,which(names(data)==oldcovname)]))[1])),"reference","","")
        body<-c(levelnames,hazardratio,"","")


      }else{
        if(!is.null(data)) reference<-c(addspace(sanitizestr(names(table(data[,which(names(data)==oldcovname)]))[1])),"reference","","")
        title<-c(covariatename,"","",globalpvalue)
        body<-cbind(levelnames,hazardratio,pvalues,rep("",length(levelnames)))

        #if interaction
      }}else{
        if(length(levelnames)!=1){
          title<-c(covariatename,"","",globalpvalue)
          body<-cbind(levelnames,hazardratio,pvalues,rep("",length(levelnames)))
        }else{
          title<-c(covariatename,hazardratio,"",globalpvalue)
        }
      }

    out<-rbind(title,reference,body)
    rownames(out)<-NULL
    colnames(out)<-NULL
    return(list(out,nrow(out)))
  })

  table<-lapply(out,function(x){return(x[[1]])})
  index<-unlist(lapply(out,function(x){return(x[[2]])}))
  table<-do.call("rbind",lapply(table,data.frame,stringsAsFactors = FALSE))

  colnames(table)<-sapply(c("Covariate",sanitizestr(beta),"p-value","Global p-value"),lbld)
  return(table)
}

#'Print multivariate summary LaTeX table
#'
#'Returns a LaTeX table of the multivariate summary.
#'
#'@param model fitted model object
#'@param data dataframe containing data
#'@keywords print
#'@export

pmvsum<-function(model,data){
  print.xtable(xtable(mvsum(model,data)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
}

#' Convert .TeX to .docx
#'
#' Converts the knitr-compiled .TeX file to a .docx file. The function calls ImageMagick to convert .pdf images into .png and then calls pandoc (included with RStudio) to convert the .Tex file into .docx.
#'
#' @param dir full path of .TeX file directory
#' @param fname .TeX file filename. Do not include extension.
#' @param pdwd full path to pandoc. The default value is the usual installation path for Windows systems if RStudio is installed.
#' @param imwd full path to image magick. Only include if there is at least one graphic.
#' @keywords print
#' @export

makedocx<-function(dir,fname,pdwd="C:/PROGRA~1/RStudio/bin/pandoc",imwd=""){
  oldwd<-getwd()
  ### convert path to forward slashes if any
  dir <- gsub("\\\\","/",dir)
  setwd(dir)
  if( grepl(" ",fname) ){## check whether fname has spaces,if it does rename it with a warning.
    file.rename(from = fname,to = gsub(" ","-",fname))
    fname <- gsub(" ","-",fname)
    warning("Your fname has spaces,the new filename is: ",fname)
  }

  if(imwd!="" & dir.exists(paste0(dir,"/figure"))){
    setwd(imwd)
    if( grepl("*ImageMagick-6[.]*",imwd) ){
      exec <- "mogrify"
    }else if ( grepl("*ImageMagick-7[.]*",imwd) ){
      exec <- "magick mogrify"
    }else warning("The installation path for ImageMagick has an unrecognized format.")

    command <- paste0(exec,' -path "',dir,'/figure/" ','-format png "',dir,'/figure/*.pdf"')
    shell(command)
  }
  setwd(dir)
  ### The lines below ammend the issue with newer versions of pandoc with the kframe environment
  ## read the original .tex file
  tx  <- readLines(paste0(fname,'.tex'),warn=FALSE)
  ## rename the environment to something simpler as suggested by John MacFarlane in the pandoc-discuss thread
  tx2 <- gsub(pattern = "\\begin{document}",
              replace = "\\renewenvironment{kframe}{}{}\\begin{document}",
              x = tx,fixed = TRUE)
  ## create a file with the workaround for the kframe environment and use it in the pandoc call below
  zz <- file(paste0(fname,'_cp.tex'),"wb")
  writeLines(tx2,con=zz)
  close(zz)

  command <- paste0('"',pdwd,'/pandoc" -o ',fname,'.docx ',fname,'_cp.tex ',
                    "--default-image-extension=png ")
  shell(command)
  # remove the files created
  file.remove(paste0(fname,'_cp.tex'))
  pngfiles = list.files(paste0(dir,'/figure/'),pattern = "*.png",full.names = TRUE)
  file.remove(pngfiles)
  setwd(oldwd)
}


#'Plot CI curve
#'
#'Plots a CI curve. Currently not very powerful. Only plots a single curve
#'
#'@param data dataframe containing data
#'@param response character vector or list of character vector. If a list it plot the '1' event for all outcomes on
#'the same plot
#'@param group string of the group want to stratify by
#'@param units units of time
#'@param main String corresponding to title
#'@param CI Bool If True will plot CI and only the '1' event. if F will plot all events except for the final one
#'@param legpos string indicating which position to put legend choies are "topright" etc
#'@param xlim numeric vector corresponding to xlimits. Default is NULL
#'@param outcomes character vector of the names of the different competing outcomes
#'@importFrom cmprsk cuminc
#'@keywords print
#'@export

plotci<-function (data,response,group=NULL,units = "months",main="Viral Infections",CI=F,legpos="topleft",xlim=NULL,outcomes=NULL){
  if(!is.null(group)){
    groups=levels(data[,group])
  }
  #If response is a list plot the '1' event for all outcomes on same plot
  if(class(response)!="list"){
    if(!is.null(group)){
      groups=levels(data[,group])
      fita <- cuminc(data[,response[1]],data[,response[2]],data[,group])
    }else{
      fita <- cuminc(data[,response[1]],data[,response[2]])
    }
    if(CI){
      plot(fita[[1]]$time,sapply(fita[[1]]$est + 1.96 * sqrt(fita[[1]]$var),
                                  function(x) min(x,1)),type = "l",lty = 2,main = paste("CI plot for ",
                                                                                            sanitizestr(nicename(response[2])),sep = ""),xlab = paste("Time (",
                                                                                                                                                        cap(units),")",sep = ""),ylim = c(0,1),ylab = paste("Incidence of ",
                                                                                                                                                                                                                 sanitizestr(nicename(response[2])),sep = ""),xlim=xlim)

      lines(fita[[1]]$time,fita[[1]]$est)
      lines(fita[[1]]$time,sapply(fita[[1]]$est - 1.96 * sqrt(fita[[1]]$var),
                                   function(x) max(x,0)),lty = 2)
    }else{
      plot(fita[[1]]$time,fita[[1]]$est,
           type = "l", main = paste("CI plot for ",
                                     sanitizestr(nicename(response[2])),sep = ""),xlab = paste("Time (",
                                                                                                 cap(units),")",sep = ""),ylim = c(0,1),ylab = paste("Incidence of ",
                                                                                                                                                          sanitizestr(nicename(response[2])),sep = ""),xlim=xlim)
      numoutcomes<-length(fita)-1
      if(numoutcomes>1){
        for (i in 2:numoutcomes){
          lines(fita[[i]]$time,fita[[i]]$est,lty=i,lwd=2)
        }
        legend(legpos,outcomes,lty = 1:numoutcomes,bty = "n",lwd=2)
      }
    }

  }else{
    d<-lapply(response,function(respons){
      fita <- cuminc(data[,respons[1]],data[,respons[2]])
      list(fita[[1]]$time,fita[[1]]$est)})
    if(is.null(xlim)) xlim=c(0,ceiling(max(sapply(d,function(x) max(x[[1]])))))
    plot(1,type="n",xlim=xlim,ylim=c(0,1),
         ylab="Cumulative Incidence",xlab = paste("Time (",cap(units),")",sep = ""),main=paste("Cumulative Incidence plot for",main))
    for(i in 1:length(d)){
      lines(d[[i]][[1]],d[[i]][[2]],lty=i,lwd=2)
    }
    legend(legpos,sapply(response,function(x) x[2]) ,col =rep(1,length(response)),lty = 1:length(response),bty = "n",lwd=2)

  }
}




#' Get CI cinfidence interval
#'
#' Returns the confidence interval of a CI at a specified time. Currently not very powerful. Only works on single strata.
#'
#' @param data dataframe containing data
#' @param response character vector of response
#' @param times numeric vector specifying single time to get CI for
#' @param units string specifying the unit of times
#' @param outcomes character vector specifying names of competing outcomes.
#' Leave NULL if there is only one outcome
#' @param decimals positive integer corresponding to the number of decimals
#' @keywords print
#' @export
citime<-function (data,response,times,units="Years",outcomes=NULL,decimals=2)
{
  out<-sapply(times,function(time){
    fita <- cuminc(data[,response[1]],data[,response[2]])
    numoutcomes<-length(fita)-1
    sapply(1:numoutcomes,function(i){
      index <- max(which(fita[[i]]$time <= time))
      est <- fita[[i]]$est[index]
      pm <- 1.96 * sqrt(fita[[i]]$var[index])
      psthr(c(est,max(est - pm,0),min(est + pm,1)),decimals)
    })
  })
  if (class(out)!="matrix")
    out<-t(out)
  out<-data.frame(out,stringsAsFactors=F)
  rownames(out)<-NULL
  if(!is.null(outcomes)){
    out<-cbind(outcomes,out)
    colnames(out)<-c("Outcome",paste(times,units))
  }else{
    colnames(out)<-paste(times,units)
  }
  return(out)
}


#' Create a forest plot
#'
#' Create a forest plot. All entires with cutoff=T will be plotted with an NA
#' rather than their original value.
#'
#' @param data dataframe containing data
#' @param xlab String corresponding to xlabel. By default is set to names(data)[2]
#' @param ylab String corresponding to ylabel. By default is set to names(data)[1]
#' @param main String corresponding to main title. By default is set to "Forest plot for subgroup analysis"
#' @param space numeric corresponding to offset of y label. Should be positive if y label is on top of the names of the y axis
#' @param bool A boolean vector. All entries with T will be invisible in the plot
#' @param xlim vector of length 2 corresponding to limits of x-axis. Default to NULL.
#' @importFrom graphics abline axis legend lines mtext par plot
#' @keywords print
#' @export
forestplot<-function (data,xlab = NULL,ylab = NULL,main = NULL,space = 0,bool=F,xlim=NULL)
{
  if (is.null(xlab))
    xlab <- names(data)[2]
  if (is.null(ylab))
    ylab <- names(data)[1]
  if (is.null(main))
    main <- "Forest plot for subgroup analysis"
  par(oma = c(0,space,0,0))
  l1 <- nrow(data)
  colors<- ifelse(bool,"white","black")
  if(is.null(xlim)) xlim<-c(0,max(data[!bool,4]))
  plot(data[,2],c(1:l1),col = colors,pch = "|",bg = colors,
       yaxt = "n",xlim = xlim,ylab = "",xlab = "",
       main = main)
  abline(v = 1,col = "red",lty = 2)
  segments(data[,3],c(1:l1),data[,4],c(1:l1),col=colors)
  axis(2,at = c(1:l1),labels = data[,1],las = 1,cex.axis = 0.8)
  mtext(side = 1,xlab,line = 2)
  mtext(side = 2,ylab,line = space + 2.5)
}


#'Plot KM and CI curves using ggplot2
#'
#'Generates Kaplan-Meier and Cumulative Incidence plots including a numbers at risk table via ggplot2 with high level of customization.
#'
#'@param response covariate to stratify table by
#'@param cov character vector with the names of columns to include in table
#'@param data dataframe containing data
#'@param type ...
#'@param times ...
#'@param table boolean indicating if you want latex markup
#'@param returns boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param xlabs booling indicating if you want to replace . and _ in strings with a space
#'@param ylabs boolean indicating if you want to display the inter quantile range (Q1,Q3) as opposed to (min,max) in the summary for continuous variables
#'@param main number of digits for the proportions when summarizing categorical data (default: 0)
#'@param ystratalabs test of choice for continuous variables, one of \emph{rank-sum} (default) or \emph{ANOVA}
#'@param ystrataname test of choice for categorical variables, one of \emph{Chi-squared} (default) or \emph{Fisher}
#'@param censor.marks ...
#'@param HR ...
#'@param HR.pval ...
#'@param conf.type ...
#'@param fsize ...
#'@param nsize ...
#'@param lsize ...
#'@param psize ...
#'@param ylims ...
#'@param xlims ......
#'@param col ...
#'@param linetype ...
#'@param legend.pos ...
#'@param pval ...
#'@param pval.pos ...
#'@param plot.event ...
#'@param event ...
#'@param flip.CIF ...
#'@param cut ...
#'@param ticklabs ...
#'@param logse ...
#'@param eventlabs ...
#'@param event.name ...
#'@param ... additional options passed to function
#'@import ggplot2
#'@import survival
#'@importFrom gridExtra grid.arrange arrangeGrob
#'@keywords dataframe
#'@export
#'@seealso \code{\link{fisher.test}}
#' @examples
#' require(survival)
#' lung$sex<-factor(lung$sex)
#' ggsurv(c("time","status"), data=lung)
#' ggsurv(c("time","status"),"sex", data=lung)
ggsurv <- function(response, cov=NULL, data, type=NULL, times = NULL, table = TRUE, returns = FALSE, xlabs = "Time", ylabs=NULL, main = NULL, ystratalabs = NULL, ystrataname = nicename(cov), censor.marks = FALSE, HR = FALSE, HR.pval= FALSE, conf.type = "none", fsize = 15, nsize = 4.5, lsize = 1.5, psize = 4.5, ylims=c(0,1), xlims=NULL, col=NULL, linetype=NULL, legend.pos = NULL, pval = TRUE, pval.pos=NULL, plot.event=1, event=c("col","linetype"), flip.CIF = FALSE, cut=NULL, ticklabs=NULL, logse=FALSE, eventlabs=NULL, event.name=NULL,...){

  event <- match.arg(event)
  if (!is.factor(data[,cov])  & !is.numeric(data[,cov]) & !is.null(cov)){ message("Coercing the cov variable to factor"); data[,cov] <- factor(data[,cov])}


  if (is.numeric(data[,cov])){
    numeric = T
    if(is.null(cut)) cut <- median(data[,cov],na.rm = T)
    data[,cov] <- factor(ifelse(data[,cov]<=cut,paste0("<=",round(cut,2)),paste0(">",round(cut,2))),levels = c(paste0("<=",round(cut,2)),paste0(">",round(cut,2))))
    if(is.null(ystratalabs)) ystratalabs <- levels(data[,cov])
  }

  # Specifing the type of plot ----------------------------------------------


  #Specifying KM or CIF & is.null(type)
  if (length(unique(data[, response[2]]))< 3 & is.null(type)) type = "KM"
  if (length(unique(data[, response[2]]))>= 3 & is.null(type)) type = "CIF"
  if(type=="KM") {
    if(is.null(main)) main <- "Kaplan-Meier Plot"
    if(is.null(ylabs)) ylabs = "Survival probability"
  }else if(type=="CIF"){
    #For the number at risk calcuation
    #data_sfit[,response[2]][data_sfit[,response[2]]!=0] <- 1
    if(is.null(main)) main <- "Cumulative Incidence Plot"
    if(is.null(ylabs)) ylabs = "Probability of an event"
  }else(stop("Type must be either KM or CIF"))


  # Labels ------------------------------------------------------------------

  multiple_lines <- !is.null(cov)

  if(is.null(ystratalabs) & multiple_lines) ystratalabs <- nicename(levels(factor(data[, cov])))
  if(is.null(ticklabs)) ticklabs <- ystratalabs

  # HR and p-val cox----------------------------------------------------------------------
  if(type=="KM" &multiple_lines & (HR|HR.pval)){
    coxfit <- coxph(as.formula(paste(paste("Surv(", response[1],
                                           ",", response[2], ")", sep = ""), "~", cov,
                                     sep = "")), data = data)


    HR_vals <- paste0("HR=",sapply(seq(length(ystratalabs)-1),function(i){
      return(psthr0(summary(coxfit)$conf.int[i,c(1, 3, 4)]))
    }))

    if(HR)ystratalabs[-1] <- paste(ystratalabs[-1],HR_vals)
    if(HR.pval) ystratalabs[-1] <- paste(ystratalabs[-1],"p=",sapply(summary(coxfit)$coef[,5],lpvalue2))
    ystratalabs[1] <- paste(ystratalabs[1],"REF")
  }


  # HR and p-val crr --------------------------------------------------------
  if(type=="CIF" & multiple_lines & (HR|HR.pval)&length(plot.event==1) & plot.event[1]==1){
    crrfit <- crrRx(as.formula(paste(paste(response,
                                           collapse = "+"), "~", cov, sep = "")),
                    data = data)

    HR_vals <- paste0("HR=",sapply(seq(length(ystratalabs)-1),function(i){
      return(psthr0(summary(crrfit)$conf.int[i,c(1, 3, 4)]))
    }))

    if(HR)ystratalabs[-1] <- paste(ystratalabs[-1],HR_vals)
    if(HR.pval) ystratalabs[-1] <- paste(ystratalabs[-1],"p=",sapply(summary(crrfit)$coef[,5],lpvalue2))
    ystratalabs[1] <- paste(ystratalabs[1],"REF")

  }
  # Model fitting KM and creating a dataframe of times--------------------------------------------------------

  if(type=="KM"){
    if(!multiple_lines){

      sfit <- survfit(as.formula(paste(paste("Surv(", response[1],
                                             ",", response[2], ")", sep = ""), "~", 1,
                                       sep = "")), data = data,conf.type=conf.type)
      ystratalabs <- "All"
    }else{
      sfit <- survfit(as.formula(paste(paste("Surv(", response[1],
                                             ",", response[2], ")", sep = ""), "~", cov,
                                       sep = "")), data = data,conf.type=conf.type)
    }



    df <- NULL
    df <- data.frame(time = sfit$time,
                     n.risk = sfit$n.risk,  n.censor = sfit$n.censor,
                     n.event = sfit$n.event, surv = sfit$surv,
                     strata = if(multiple_lines){
                       summary(sfit, censored = T)$strata
                     }else factor("All"),
                     upper = if(conf.type != "none"){
                       sfit$upper
                     }else factor(NA),
                     lower = if(conf.type != "none"){
                       sfit$lower
                     }else factor(NA))
    levels(df$strata) <- ystratalabs
    zeros <- data.frame(time = 0, surv = 1,
                        strata = if(multiple_lines){
                          levels(df$strata)
                        }else factor("All"),
                        upper = 1, lower = 1)
    df <- plyr::rbind.fill(zeros, df) # Forcing the curves to start at 1

    df$strata <- factor(df$strata,levels=ystratalabs)
  }


  # Model fitting CIF -------------------------------------------------------
  if(type=="CIF"){



    if(!multiple_lines){

      invisible(capture.output(fit <-  cuminc( data[,response[1]],data[,response[2]] )))
      ystratalabs <- " "
      gsep = " "

      if(table){ #Sfit is for the numbers at risk so both events are counted the same way

        temp <- data
        temp[,response[2]][temp[,response[2]] > 0] <- 1
        sfit <- survfit(as.formula(paste(paste("Surv(", response[1],
                                               ",", response[2], ")", sep = ""), "~", 1,
                                         sep = "")), data = temp)
      }

    }else{
      newgpvar <- paste0(data[,cov],":")
      newgpvar <- factor(newgpvar, levels = paste0(levels(data[,cov]),":") )
      invisible(capture.output(fit <- cuminc(data[,response[1]],data[,response[2]], newgpvar)))
      gsep = ": "

      if(table){ #Sfit is for the numbers at risk so both events are counted the same way

        temp <- data
        temp[,response[2]][temp[,response[2]] > 0] <- 1
        sfit <- survfit(as.formula(paste(paste("Surv(", response[1],
                                               ",", response[2], ")", sep = ""), "~", cov,
                                         sep = "")), data = temp)
      }

    }

    ##Setting up the data frame
    if (!is.null(fit$Tests)){
      test <- fit$Test
      fit <- fit[names(fit) != "Tests"]
    }
    fit2 <- lapply(fit, `[`, 1:3)
    gnames <- names(fit2)

    fit2_list <- lapply(seq_along(gnames), function(ind) {
      df <- as.data.frame(fit2[[ind]])
      df$name <- gnames[ind]
      df
    })

    df <- do.call(rbind, fit2_list)
    df$event <- sapply(strsplit(df$name, split = gsep), `[`, 2)
    df$strata <- sapply(strsplit(df$name, split = gsep), `[`, 1)
    df$strata <- factor(df$strata, levels = levels(data[,cov]) )
    levels(df$strata) <- ystratalabs

    if(multiple_lines){df$strata <- factor(df$strata,levels=ystratalabs)
    }else df$strata <- "ALL"

    df$std <- std <- sqrt(df$var)

    names(df)[names(df)=="est"] <- 'surv'

    df <- df[df$event %in% plot.event,]
    df$upper <- NA
    df$lower <- NA



    if(conf.type!="none") {
      CIF_confint <- survfit_confint(p=df$surv, se=df$std, conf.type=conf.type, conf.int=0.95,ulimit=TRUE,logse=logse)
      df$upper <- CIF_confint$upper
      df$lower <- CIF_confint$lower
    }


    if(flip.CIF) {
      df$surv <- 1-df$surv
      df$upper <- 1-df$upper
      df$lower <- 1-df$lower
    }

    if(!is.null(eventlabs)) {
      df$event <- factor(df$event)
      levels(df$event) <- eventlabs
    }
  }


  # axis and legend positions --------------------------------------------------------------------
  m <- max(nchar(ystratalabs))
  maxxval = max(df$time,times[length(times)])
  if( is.null(xlims) ){
    maxxlim =  maxxval
  }else {
    if( length(xlims)!=2 | !is.numeric(xlims) ) stop("xlims must be a numeric vector of length 2.")
    maxxlim = xlims[2]
  }


  linetype_name <- col_name <- ystrataname

  leg.pos <- legend.pos
  d <- length(levels(df$strata))


  if(!multiple_lines & (type=="KM" | (type=="CIF" & length(plot.event)==1))){
    leg.pos <- 'none'
  }else if(is.null(legend.pos)) {
    if(type=="CIF" & flip.CIF==F){ leg.pos <- c(min(0.05+m/200,0.5), 0.95-d*0.05)
    }else leg.pos <- c(min(0.05+m/200,0.5), 0.05+d*0.05)
  }


  #leg.pos <- if(){ "none" }else leg.pos ## don't plot if no strata and only one

  if(is.null(times)) times <- break_function(maxxlim)

  # plotting ----------------------------------------------------------------

  if(type=="CIF" & length(plot.event)>1){
    if(is.null(event.name)) event.name <- 'event'
    if(event=='linetype') {p <- ggplot(df)+ geom_step( aes(time, surv, color = strata,linetype=event),size = lsize); linetype_name = event.name}
    if(event=='col') {p <- ggplot(df)+ geom_step( aes(time, surv, color = event,linetype=strata),size = lsize); col_name=event.name}
  }else{
    p <- ggplot(df) + geom_step(aes(time, surv, group = strata,linetype = strata, col=strata), size = lsize)
  }


  # Confidence intervals to the plot ----------------------------------------


  if(conf.type!="none" ){
    if(type=="KM")p <- p +  geom_ribbon(data=df[!is.na(df$upper) & !is.na(df$lower),], aes(x=time, fill=strata, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, show.legend = F)

    if(type=="CIF" & event == "linetype"){
      for( evnt in unique(df$event) ){
        p <- p + geom_ribbon(data=df[!is.na(df$upper) & !is.na(df$lower) & df$event == evnt,], aes(x=time, fill=strata, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, show.legend = F)
      }
    }

    if(type=="CIF" & event == "col" & length(plot.event)>1){
      for( stra in unique(df$strata) ){
        p <- p + geom_ribbon(data=df[!is.na(df$upper) & !is.na(df$lower) & df$strata == stra,], aes(x=time, fill=event, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, show.legend = F)
      }
    }

  }

  # Modyfing axis, titles, fonts, and theming -------------------------------


  p <- p+
    theme_classic() +
    theme(axis.title.x = element_text(vjust = 0.5)) +
    theme(text = element_text(size = fsize),
          axis.text=element_text(size=fsize)) +  ##increase font size
    scale_x_continuous(paste0("\n",xlabs), breaks = times,
                       limits = c(0, maxxval)) +
    coord_cartesian(xlim=c(0,maxxlim)) + ### changes the actual plotted limits if needed
    scale_y_continuous(paste0(ylabs,"\n"), limits = ylims) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.key = element_rect(colour = "transparent", fill = "transparent")) +
    theme(legend.background=element_blank()) +
    labs(linetype = linetype_name, col=col_name) +
    ggtitle(main) + theme(plot.title = element_text(face="bold",hjust = 0.5))+
    theme(legend.position = leg.pos)


  # censoring ---------------------------------------------------------------


  if( censor.marks & type=="KM") p <- p + geom_point(data = subset(df, n.censor>0),
                                                     aes(x=time, y=surv, group = strata, col=strata),
                                                     shape=3, size=3, stroke = 1.5,show.legend =F)


  # Log rank p-val ----------------------------------------------------------
  if(pval & type=="KM" & multiple_lines) {
    sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    pval <- pchisq(sdiff$chisq, length(sdiff$n)-1, lower.tail = FALSE)
    pvaltxt <- ifelse(pval < 0.0001, "p < 0.0001", paste("p =", signif(pval, 3)))
    pvaltxt <- paste(pvaltxt,"(Log Rank)")
    #pvaltxt <- paste("p =", signif(pval, 3), "(Log Rank)")
    if(is.null(pval.pos)){
      p <- p + annotate("text", x = 0.9 * max(times), y = ylims[1], label = pvaltxt, size = psize)
    }else p <- p + annotate("text", x = pval.pos[1], y = pval.pos[2], label = pvaltxt, size = psize)
  }


  # Gray's pvalue -----------------------------------------------------------

  if( !is.null(cov) & pval & type=="CIF") {
    if(length(plot.event)==1 ){
      test <- test[rownames(test)==plot.event,]
      pval <- test[2]
      pvaltxt <- ifelse(pval < 0.0001, "p < 0.0001", paste("p =", signif(pval, 3)))
      pvaltxt <- paste(pvaltxt,"(Gray's test)")
      if(is.null(pval.pos)){
        p <- p + annotate("text", x = 0.9 * max(times), y = ylims[1], label = pvaltxt, size = psize)
      }else p <- p + annotate("text", x = pval.pos[1], y = c(pval.pos[2]), label = pvaltxt)
    }else{

      pval <- test[,2]
      pvaltxt <- ifelse(pval < 0.0001, "p < 0.0001", paste("p =", signif(pval, 3)))
      pvaltxt <- c("Gray's test",paste("Event", rownames(test), pvaltxt))
      if(is.null(pval.pos)){

        p <- p + annotate("text", x = 0.9 * max(df$time), y = c(0.12,0.08,0.04), label = pvaltxt)
      }else p <- p + annotate("text", x = pval.pos[1], y = c(pval.pos[2],pval.pos[2]-0.04, pval.pos[2]-0.08), label = pvaltxt)
    }
  }




  # Colour, linetype and fill -----------------------------------------------



  if(!is.null(col)) {
    p <- p + scale_colour_manual(values = col)
    p <- p + scale_fill_manual(values = col)
  }
  if(!is.null(linetype)) {
    p <- p + scale_linetype_manual(values = linetype)

  }




  # Numbers at risk ---------------------------------------------------------

  if(table) {
    ## Create table graphic to include at-risk numbers

    blank.pic <- ggplot(df, aes(time, surv)) +
      geom_blank() +
      theme_bw() +
      scale_x_continuous(breaks = times, limits = c(0, maxxval)) +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.ticks = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), panel.background = element_rect(fill = "transparent"))

    sfit.summary <- summary(sfit, times = times, extend = TRUE)
    risk.data <- data.frame(strata = if(multiple_lines){
      sfit.summary$strata
    }else factor("All"),
    time = sfit.summary$time,
    n.risk = sfit.summary$n.risk)
    # if risk and event do paste0(sfit.summary$n.risk, "(",sfit.summary$n.events,")")
    risk.data$strata <- factor(risk.data$strata, levels=rev(levels(risk.data$strata)))

    cols1 <- .extract_ggplot_colors(p, grp.levels = ystratalabs)
    if(multiple_lines==F) cols1 = "black"
    ### TODO: check length of m, if too long use dashes as thickmarks

    n_strata <- length(ystratalabs)
    #yticklabs <- rep("-", n_strata)
    yticklabs <- unname(rev(ticklabs))
    strataylab = ifelse(n_strata==1, "\n", ystrataname)

    data.table <- ggplot(risk.data, aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
      geom_text(hjust="middle", vjust="center", size = nsize) +
      theme_bw() +
      scale_x_continuous("Numbers at risk", breaks = times, limits = c(0,maxxval)) +
      coord_cartesian(xlim=c(0,maxxlim)) + ### changes the actual plotted limits if needed
      theme(legend.position = "none") +
      theme(text = element_text(size = fsize)) + ##increase font size
      theme(plot.margin = unit(c(-1.5, 1, 0.1, 0.2), "lines"))+
      scale_y_discrete(strataylab, breaks = as.character(levels(risk.data$strata)), labels = yticklabs)

    data.table <- data.table +suppressWarnings(theme(axis.title.x = element_text(size = fsize, vjust = 1), panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(), panel.border = element_blank(),
                                                     axis.text.x = element_blank(), axis.ticks = element_blank(),
                                                     axis.text.y = element_text(face = "bold", hjust = 1,colour = rev(cols1))))




    if(all(yticklabs=="-")) data.table <- .set_large_dash_as_ytext(data.table)



    gA <- ggplotGrob(p)
    gB <- ggplotGrob(blank.pic)
    gC <- ggplotGrob(data.table)
    maxWidth = grid::unit.pmax(gA$widths[2:5], gC$widths[2:5])
    gA$widths[2:5] <- as.list(maxWidth)
    gB$widths[2:5] <- as.list(maxWidth)
    gC$widths[2:5] <- as.list(maxWidth)

    grid.arrange(gA, gB, gC,
                 clip = FALSE, nrow = 3, ncol = 1,
                 heights = unit(c(2, .1, .25), c("null", "null", "null")))

    if(returns) {
      a <- arrangeGrob(p, blank.pic, data.table,
                       clip = FALSE, nrow = 3, ncol = 1,
                       heights = unit(c(2, .1, .25), c("null", "null", "null")))
      return(a)
    }
  }


}

#' Plot univariate relationships all on one plot
#' Need a hierarchy of how to display plots sensibly
#' If response is continuous
#'   For a numeric predictor -> scatterplot
#'   For a categorical predictor -> boxplot
#' If response is a factor
#'   For a numeric predictor -> boxplot
#'   For a categorical predictor -> barplot
#' @param response character vector with names of columns to use for response
#' @param covs character vector with names of columns to use for covariates
#' @param data dataframe containing your data
#' @param showN boolean indicating whether sample sizes should be shown on the plots
#' @param na.rm boolean indicating whether na values should be shown or removed
#' @param response_title character value with title of the plot
#' @keywords plot
#' @importFrom ggplot2 ggplot aes_string geom_boxplot  geom_point geom_text stat_summary scale_x_discrete stat theme labs .data
#' @importFrom ggpubr ggarrange
#' @export
plot_univariate <- function(response,covs,data,showN=FALSE,na.rm=TRUE,response_title=NULL){
  for (v in c(response,covs)){
    if (class(data[[v]])=='character') data[[v]] <- factor(data[[v]])
  }

  if (is.null(response_title)) response_title = response
  response_title = niceStr(response_title)
  plist <- NULL
  if (class(data[[response]])[1] %in% c('factor','ordered')){
    levels(data[[response]]) = niceStr(levels(data[[response]]))
    for (x_var in covs){
      # remove missing data, if requested
      if (na.rm) pdata = stats::na.omit(data[,c(response,x_var)]) else pdata = data[,c(response,x_var)]

      if (class(pdata[[x_var]])[1] =='numeric' ){
        p <- ggplot(data=pdata, aes(y=.data[[response]],x=.data[[x_var]],fill=.data[[response]])) +
          geom_boxplot()
        if (showN){
          p=  p+
            stat_summary(geom='text',fun.data = lbl_count,vjust=-0.5,hjust=1)
        }
      } else {
        p <- ggplot(data=pdata, aes(x=.data[[x_var]],fill=.data[[response]])) +
          geom_bar(position='fill') +
          scale_x_discrete(labels= function(x) wrp_lbl(x))
        if (showN){
          p <- p +
            geom_text(aes(label=stat(count)),stat='count',position='fill',vjust=1)
        }
        if (length(unique(pdata[[x_var]]))>8){
          p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
      }
      plist[[x_var]] <- p  +
        theme(axis.text.y=element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = 'bottom',
              plot.title = element_text(size=10),
              plot.margin = unit(c(0,1,0,1), "lines")) +
        labs(title=niceStr(x_var),x='',y='',fill=response_title)
    }
  } else{
    for (x_var in covs){
      # remove missing data, if requested
      if (na.rm) pdata = stats::na.omit(data[,c(response,x_var)]) else pdata = data[,c(response,x_var)]

      if (class(pdata[[x_var]])[1] =='numeric' ){
        p <- ggplot(data=pdata, aes(y=.data[[response]],x=.data[[x_var]])) +
          geom_point()
      } else
        p <- ggplot(data=pdata, aes(y=.data[[response]],x=.data[[x_var]],fill=.data[[response]])) +
          geom_boxplot() +
          scale_x_discrete(labels= function(x) wrp_lbl(x))
      if (showN){
        p=  p+
          stat_summary(geom='text',fun.data = lbl_count,vjust=-0.5,hjust=1)
      }
      if (length(unique(pdata[[x_var]]))>8){
        p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
      plist[[x_var]] <- p  +
        theme(
          legend.position = 'bottom',
          plot.title = element_text(size=10),
          plot.margin = unit(c(0,1,0,1), "lines")) +
        labs(title=niceStr(x_var),x='',y='',fill=response_title)
    }

  }
  ggarrange(plotlist=plist,
            common.legend = T,
            ncol=2,
            nrow=ceiling(length(plist)/2))
}

# Rmarkdown Reporting --------------------------------------------------------------
#' The output function for the print methods
#' Table output defaults to kable, but the kableExtra package doesn't work well with Word.
#' To export nice tables to Word use options('doc_type'='doc')
#' @param tab a table to format
#' @param to_indent numeric vector the length of nrow(tab) indicating which rows to indent
#' @param to_bold numeric vector the length of nrow(tab) indicating which rows to bold
#' @param caption table caption
#' @param chunk_label only used if options("doc_type"='doc') to allow cross-referencing
#' @param ... other variables passed to covsum and the table output function
#' @export
outTable <- function(tab,to_indent=numeric(0),to_bold=numeric(0),caption=NULL,chunk_label,...){
  tab = as.data.frame(tab) #strips any tibble aspects 
  out_fmt = ifelse(is.null(getOption('doc_type')),'pdf',getOption('doc_type'))
  out_fmt =ifelse(out_fmt%in%c('doc','docx'),'doc','pdf')
  chunk_label = ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label)

  # set NA to empty in kable
  options(knitr.kable.NA = '')

  if (is.null(to_indent)) to_indent = numeric(0)
  to_indent = as.vector(to_indent)
  rownames(tab) <- NULL


  if (out_fmt=='doc'){
    caption = if (!is.null(caption)) {ifelse(chunk_label=='NOLABELTOADD',caption,paste0('(\\#tab:',chunk_label,')',caption))}
    tab[is.na(tab)] <-'&nbsp;' # This is necessary to assign the 'Compact' style to empty cells
    tab[tab==''] <-'&nbsp;'

    tab[[1]][to_indent] <- sapply(tab[[1]][to_indent],function(x) paste('&nbsp;&nbsp;',x))
    if (length(to_bold)>0) {
      pander::pander(tab,
                     caption=caption,
                     emphasize.strong.rows=to_bold,
                     split.table=Inf, split.cells=15,
                     justify = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))

    } else {
      pander::pander(tab,
                     caption=caption,
                     split.table=Inf, split.cells=15,
                     justify = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))
    }
  } else {
    if (nrow(tab)>30){
      kout <- knitr::kable(tab, booktabs=TRUE,
                           longtable=TRUE,
                           linesep='',
                           caption=caption,
                           align = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))
      if (ncol(tab)>4) {
        kout <- kableExtra::kable_styling(kout,full_width = T,latex_options = c('repeat_header'))
      } else {
        kout <- kableExtra::kable_styling(kout,latex_options = c('repeat_header'))
      }
    } else {
      kout <- knitr::kable(tab, booktabs=TRUE,
                           longtable=FALSE,
                           linesep='',
                           caption=caption,
                           align = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))
      if (ncol(tab)>4) kout <- kableExtra::kable_styling(kout,full_width = T)
    }
    kout <- kableExtra::add_indent(kout,positions = to_indent)
    if (length(to_bold)>0){
      kout<- kableExtra::row_spec(kout,to_bold,bold=TRUE)
    }
    kout
  }

}

#'
#'Returns a dataframe corresponding to a descriptive table for printing in Rmarkdown
#' The default output is a kable table for use in pdfs or html, but pander tables can be produced
#' for Word documents by specifying options('doc_type'='doc') in the setup chunk of the markdown document.
#'
#'@param data dataframe containing data
#'@param covs character vector with the names of columns to include in table
#'@param maincov covariate to stratify table by
#'@param caption character containing table caption
#'@param excludeLevels a named list of levels to exclude in the form list(varname =c('level1','level2')) these levels will be excluded from association tests
#'@param showP boolean indicating if p values should be displayed (may only want corrected p-values)
#'@param IQR boolean indicating if you want to display the inter quantile range (Q1,Q3) as opposed to (min,max) in the summary for continuous variables
#'@param tableOnly should a dataframe or a formatted print object be returned
#'@param covTitle character with the names of the covariate column
#'@param chunk_label only used if options("doc_type"='doc') to allow cross-referencing
#'@param ... other variables passed to covsum and the table output function
#'@keywords dataframe
#'@return A formatted table displaying a summary of the covariates stratified by maincov
#'@export
#'@seealso \code{\link{fisher.test}}, \code{\link{chisq.test}}, \code{\link{wilcox.test}}, \code{\link{kruskal.test}}, and \code{\link{anova}}
rm_covsum <- function(data,covs,maincov=NULL,caption=NULL,excludeLevels=NULL,showP=TRUE,IQR=FALSE,tableOnly=FALSE,covTitle='Covariate',
                        chunk_label,...){

  tab <- covsum(data,covs,maincov,markup=FALSE,IQR=IQR,sanitize=FALSE,...)
  if (is.null(maincov))
    showP=FALSE
  to_bold = numeric(0)
  if (showP) {
    # format p-values nicely
    p_vals <- tab[,ncol(tab)]
    new_p <- sapply(p_vals,formatp)
    tab[,ncol(tab)] <- new_p
    to_bold <- which(as.numeric(new_p)<0.05)
  } else {
    tab <- tab[,!names(tab) %in%'p-value']
  }
  nice_var_names = gsub('[_.]',' ',covs)
  to_indent <- which(!tab$Covariate %in% nice_var_names )
  if (IQR) {
    tab$Covariate <- gsub('[(]Q1,Q3[)]','(IQR)',tab$Covariate)
  }
  if (covTitle !='Covariate') names(tab[1]) <-covTitle
  if (tableOnly){
    return(tab)
  }
  if (is.null(caption)) {
    if (!is.null(maincov)){
      caption = paste0('Summary sample statistics by ',nicename(maincov),'.')
    } else
      caption = 'Summary sample statistics.'
  }

  outTable(tab=tab,to_indent=to_indent,to_bold=to_bold,
           caption=caption,
           chunk_label=ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label))

}

#' Output several univariate models nicely in a single table in Rmarkdown
#' The default output is a kable table for use in pdfs or html, but pander tables can be produced
#' for Word documents by specifying options('doc_type'='doc') in the setup chunk of the markdown document.
#'
#' @param response string vector with name of response
#' @param covs character vector with the names of columns to fit univariate models to
#' @param data dataframe containing data
#' @param caption table caption
#' @param CIwidth width of confidence interval, default is 0.95
#' @param showP boolean indicating if p-values should be shown (may only want to show holm-corrected values)
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param removeInf boolean indicating if infinite estimates should be removed from the table
#' @param HolmGlobalp boolean indicting if a Holm-corrected p-value should be presented
#' @param chunk_label only used if options("doc_type"='doc') to allow cross-referencing
#' @param ... other variables passed to uvsum and the table output functions
#' @export
#'
rm_uvsum <- function(response, covs , data ,caption=NULL,CIwidth=0.95,showP=T,tableOnly=FALSE,removeInf=T,HolmGlobalp=FALSE,chunk_label,...){

  # get the table
  tab <- uvsum(response,covs,data,CIwidth=CIwidth,markup = FALSE,sanitize=FALSE,...)

  cap_warn <- character(0)
  if (removeInf){
    # Do not display unstable estimates
    inf_values =  grep('Inf',tab[,2])
    if (length(inf_values)>0){
      tab[inf_values,2:4] <-NA
      cap_warn <- paste0(cap_warn,'Covariates with unstable estimates:',paste(tab$Covariate[inf_values],collapse=','),'.')

    }
  }

  # If HolmGlobalp = T then report an extra column with the adjusted p and only bold these values
  if (HolmGlobalp){
    p_sig <- stats::p.adjust(tab$`Global p-value`,method='holm')
    tab$"Holm Adj p" = p_sig
  } else {
    p_sig <- tab$`Global p-value`
  }

  to_bold <- which(as.numeric(p_sig)<0.05)
  nice_var_names = gsub('[_.]',' ',covs)
  to_indent <- which(!tab$Covariate %in% nice_var_names )

  tab[["p-value"]] <- formatp(tab[["p-value"]])
  tab[["Global p-value"]] <- formatp(tab[["Global p-value"]])
  if (HolmGlobalp){
    tab[["Holm Adj p"]] <- formatp(tab[["Holm Adj p"]])
  }

  # If all outcomes are continunous (and so all p-values are NA), remove this column & rename Global p-value to p-value
  if (sum(is.na(tab[["p-value"]]))==nrow(tab)) {
    tab <- tab[,-which(names(tab)=="p-value")]
    names(tab) <- gsub('Global p-value','p-value',names(tab))
  }

  if (tableOnly){
    if (nchar(cap_warn)>0) warning(cap_warn)
    return(tab)
  }
  if(is.null(caption)){
    if (length(response)==2) {outcome = 'survival'} else{ outcome = niceStr(response)}
    caption = paste0('Univariate analysis of predictors of ',outcome,'.',cap_warn)
  } else if (caption=='none' & identical(cap_warn,character(0))) {
    caption=NULL
  } else caption = paste0(caption,cap_warn)

  outTable(tab=tab,to_indent=to_indent,to_bold=to_bold,
           caption=caption,
           chunk_label=ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label))
}

#' Output a multivariable model nicely in Rmarkdown
#' The default output is a kable table for use in pdfs or html, but pander tables can be produced
#' for Word documents by specifying options('doc_type'='doc') in the setup chunk of the markdown document.
#'
#' @param model model fit
#' @param data data that model was fit on
#' @param caption table caption
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param HolmGlobalp boolean indicting if a Holm-corrected p-value should be presented
#' @param chunk_label only used if options("doc_type"='doc') to allow cross-referencing
#' @export
rm_mvsum <- function(model , data ,caption=NULL,tableOnly=FALSE,HolmGlobalp=FALSE,chunk_label){

  # get the table
  tab <- mvsum(model=model,data=data,markup = FALSE, sanitize = FALSE, nicenames = T)


  # Reduce the number of significant digits in p-values
  p_val <-  formatp(tab$`p-value`)
  g_p_val = formatp(tab$`Global p-value`)
  # If HolmGlobalp = T then report an extra column with the adjusted p and only bold these values
  if (HolmGlobalp){
    gp <- stats::p.adjust(tab$`Global p-value`,method='holm')
  } else {
    gp <- tab$`Global p-value`
  }
  to_bold <- which(as.numeric(gp)<0.05)
  to_indent <- which(is.na(g_p_val))

  tab$`p-value` <- p_val
  tab$`Global p-value` <- g_p_val

  # Reorder
  if (HolmGlobalp){
    tab$`Holm Adj p` <- formatp(gp)
  }
  # TO DO: possibly automate this... need to extract response from mvsum
  # if(is.null(caption)){
  #   caption = paste0('Multivariable analysis of predictors of ',niceStr(response),'.')
  # } else if (caption=='none') {
  #   caption=NULL
  # }


  # If all outcomes are contunous (and so all p-values are NA), remove this column
  if (sum(is.na(tab[["p-value"]]))==nrow(tab)) tab <- tab[,-which(names(tab)=="p-value")]

  if (tableOnly){
    return(tab)
  }
  outTable(tab=tab,to_indent=to_indent,to_bold=to_bold,
           caption=caption,
           chunk_label=ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label))

}
