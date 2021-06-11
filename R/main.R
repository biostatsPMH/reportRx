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
#'Returns a dataframe corresponding to a descriptive table. 
#'
#'Comparisons for categorical variables default to chi-square tests, but if there are counts of <5 then the Fisher Exact 
#'test will be used and if this is unsuccessful then a second attempt will be made computing p-values using MC simulation. 
#'If testcont='ANOVA' then the t-test with unequal variance will be used for two groups and an ANOVA will be used for three or more.
#'The statistical test used can be displayed by specifying show.tests=TRUE.
#'
#'@param data dataframe containing data
#'@param covs character vector with the names of columns to include in table
#'@param maincov covariate to stratify table by
#'@param digits number of digits for summarizing mean data
#'@param numobs named list overriding the number of people you expect to have the covariate
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param IQR boolean indicating if you want to display the inter quantile range (Q1,Q3) as opposed to (min,max) in the summary for continuous variables
#'@param all.stats boolean indicating if all summary statistics (Q1,Q3 + min,max on a separate line) should be displayed. Overrides IQR.
#'@param pvalue boolean indicating if you want p-values included in the table
#'@param show.tests boolean indicating if the type of statistical used should be shown in a column beside the pvalues. Ignored if pvalue=FALSE.
#'@param excludeLevels a named list of covariate levels to exclude from statistical tests in the form list(varname =c('level1','level2')). These levels will be excluded from association tests, but not the table. This can be useful for levels where there is a logical skip (ie not missing, but not presented). Ignored if pvalue=FALSE.
#'@param full boolean indicating if you want the full sample included in the table, ignored if maincov is NULL
#'@param digits.cat number of digits for the proportions when summarizing categorical data (default: 0)
#'@param testcont test of choice for continuous variables,one of \emph{rank-sum} (default) or \emph{ANOVA}
#'@param testcat test of choice for categorical variables,one of \emph{Chi-squared} (default) or \emph{Fisher}
#'@param include_missing Option to include NA values of maincov. NAs will not be included in statistical tests
#'@param percentage choice of how percentages are presented ,one of \emph{column} (default) or \emph{row}
#'@keywords dataframe
#'@export
#'@seealso \code{\link{fisher.test}},\code{\link{chisq.test}},\code{\link{wilcox.test}},\code{\link{kruskal.test}},and \code{\link{anova}}
covsum <- function(data,covs,maincov=NULL,digits=1,numobs=NULL,markup=TRUE,sanitize=TRUE,nicenames=TRUE,IQR = FALSE,all.stats=FALSE,pvalue=TRUE,show.tests=FALSE,excludeLevels=NULL,full=TRUE,
                   digits.cat = 0,testcont = c('rank-sum test','ANOVA'),testcat = c('Chi-squared','Fisher'),include_missing=FALSE,percentage=c('column','row'))
{
  # New LA 18 Feb, test for presence of variables in data and convert character to factor
  missing_vars = setdiff(covs,names(data))
  if (length(missing_vars)>0){  stop(paste('These covariates are not in the data:',missing_vars))  }
  for (v in c(maincov,covs)) if (class(data[[v]])[1]=='character') data[[v]] <- factor(data[[v]])
  
  testcont <- match.arg(testcont)
  testcat <- match.arg(testcat)
  percentage <- match.arg(percentage)
  
  if (!pvalue) {
    show.tests<- FALSE
    excludeLevels<- NULL
    }
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity
  }
  digits <- as.integer(digits)
  digits.cat <- as.integer(digits.cat)
  if( digits<0 ) stop("parameter 'digits' cannot be negative!")
  if( digits.cat<0 ) stop("parameter 'digits.cat' cannot be negative!")
  if(!sanitize) sanitizestr<-identity
  if(!nicenames) nicename<-identity
  if(!is.null(maincov)){
    
    ##JW Removes missing of maincov
    if(include_missing==FALSE)  data <- data[!is.na(data[[maincov]]),]
    
    #JW May 2021 keeps the NAs (all instances of useNA = 'ifany' have been added) 
    levels <- names(table(data[[maincov]],useNA = 'ifany'))
    levels<-c(list(levels),as.list(levels))
  }else{
    full = TRUE # don't allow users to specify full = FALSE when there is no main covariate
    levels<-"NOMAINCOVNULLNA"
  }
  N=nrow(data)
  if(!is.null(maincov)){
    nmaincov<-c(sum(table(data[[maincov]],useNA = 'ifany')),table(data[[maincov]],useNA = 'ifany'))
  }else{
    nmaincov<-N
    p<-NULL
  }
  out<-lapply(covs,function(cov){
    ismiss=F
    n<-sum(table(data[[cov]]))
    
    # Exclude specified levels
    if (!is.null(excludeLevels[[cov]])){
      excludeLevel = excludeLevels[[cov]]
    } else excludeLevel = ''
    
    # Set up the first column
    factornames<-NULL
    if(is.null(numobs[[cov]]))  numobs[[cov]]<-nmaincov
    if(numobs[[cov]][1]-n>0) {ismiss=T
    factornames<-c(factornames,"Missing")
    }
    #if the covariate is a factor
    if(is.factor(data[[cov]])){
      factornames<-c(levels(data[[cov]]),factornames)
      if (!is.null(maincov)) {
        if(pvalue){
          pdata = data[!(data[[cov]] %in% excludeLevel),]
          # Check for low counts and if found perform Fisher test
          if( testcat[1]=='Fisher' | sum(table(pdata[[maincov]], pdata[[cov]],exclude = excludeLevel)<5)>0){
            p_type <- 'Fisher Exact'
            p <- try(fisher.test(pdata[[maincov]], pdata[[cov]])$p.value,silent=T)
            if (class(p)[1]=='try-error'){
              p <- try(fisher.test(pdata[[maincov]], pdata[[cov]],simulate.p.value =  T)$p.value,silent=T)
              p_type <- 'MC sim'
            }
          } else {
            p_type = 'Chi Sq'
            p = try(chisq.test(pdata[[maincov]], pdata[[cov]])$p.value,silent=T)
          }
          # p <- if( testcat=='Fisher' | sum(table(pdata[[maincov]], pdata[[cov]]))<5){ try(fisher.test(pdata[[maincov]], pdata[[cov]])$p.value,silent=T)
          # } else try(chisq.test(pdata[[maincov]], pdata[[cov]])$p.value,silent=T)
          if (class(p) == "try-error") p <- NA
          p <- lpvalue(p)
        }
      }
      #set up the main columns
      if (percentage == "column")   {
        onetbl<-mapply(function(sublevel,N){
          missing<-NULL
          if(is.na(sublevel[1])| sublevel[1]!="NOMAINCOVNULLNA"){
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
      }
      
      if(percentage=='row') {
        onetbl<-mapply(function(sublevel,N){
          missing<-NULL
          if(is.na(sublevel[1])| sublevel[1]!="NOMAINCOVNULLNA"){
            subdata<-subset(data,subset=data[[maincov]]%in%sublevel)
          }else{
            subdata<-data
          }
          table<-table(subdata[[cov]])
          tbl<-table(subdata[[cov]])
          n<-sum(tbl)
          if(ismiss) missing<-N-n
          tbl<-c(tbl,lbld(missing))
          return(tbl)
        },levels,numobs[[cov]])
        
        if(ismiss){
          #More than two rows with missing 
          if(dim(onetbl)[1]>2){
            onetbl[-nrow(onetbl),-1] <- t(apply(onetbl[-nrow(onetbl),-1],1,function(x){
              x <- as.numeric(x)
              prop <- round(x/sum(x),2+digits.cat)*100
              
              prop.fmt <- sprintf(paste0("%.",digits.cat,"f"),prop)
              return(paste(x," (",prop.fmt,")",sep=""))
            }))
          }else{
            #Only one row with missing 
            onetbl[-nrow(onetbl),-1] <- (function(x){
              x <- as.numeric(x)
              prop <- round(x/sum(x),2+digits.cat)*100
              
              prop.fmt <- sprintf(paste0("%.",digits.cat,"f"),prop)
              return(paste(x," (",prop.fmt,")",sep=""))
            })( onetbl[-nrow(onetbl),-1])
          }
          
          
        }else{
          #multiple rows and no missing
          if(!is.null(dim(onetbl))){
            onetbl[,-1] <- t(apply( onetbl[,-1],1,function(x){
              x <- as.numeric(x)
              prop <- round(x/sum(x),2+digits.cat)*100
              
              prop.fmt <- sprintf(paste0("%.",digits.cat,"f"),prop)
              return(paste(x," (",prop.fmt,")",sep=""))
            }))
          }else{
            #one row and no missing
            onetbl[-1] <- (function(x){
              x <- as.numeric(x)
              prop <- round(x/sum(x),2+digits.cat)*100
              
              prop.fmt <- sprintf(paste0("%.",digits.cat,"f"),prop)
              return(paste(x," (",prop.fmt,")",sep=""))
            })(onetbl[-1])
          }
        }
      }

      #if the covariate is not a factor
    }else{
      #setup the first column
      #factornames <- c("Mean (sd)",ifelse(IQR,"Median (Q1,Q3)","Median (Min,Max)"),factornames)
      if (all.stats) {factornames <- c("Mean (sd)", "Median (Q1,Q3)", "Range (min, max)", factornames)
      } else factornames <- c("Mean (sd)",ifelse(IQR,"Median (Q1,Q3)","Median (Min,Max)"),factornames)
      if (!is.null(maincov)) {
        if(pvalue){
          # p <- if( testcont=='rank-sum test'){
          #   if( length(unique(data[[maincov]]))==2 ){
          #     try( wilcox.test(data[[cov]] ~ data[[maincov]])$p.value )
          #   } else try( kruskal.test(data[[cov]] ~ data[[maincov]])$p.value )
          # } else try(anova(lm(data[[cov]] ~ data[[maincov]]))[5][[1]][1])
          #LA enable output of statistical tests
          if( testcont=='rank-sum test'){
            if( length(unique(data[[maincov]]))==2 ){
              p_type = 'Wilcoxon Rank Sum'
              p <- try( wilcox.test(data[[cov]] ~ data[[maincov]])$p.value )
            } else {
              p_type='Kruskal Wallis'
              p <-try( kruskal.test(data[[cov]] ~ data[[maincov]])$p.value )
            }
          } else {
            if( length(unique(data[[maincov]]))==2 ){
              p_type = 't-test'
              p <-try( t.test(data[[cov]] ~ data[[maincov]])$p.value )
            } else {
              p_type = 'ANOVA'
              p <-try(anova(lm(data[[cov]] ~ data[[maincov]]))[5][[1]][1])
            }}
          
          if (class(p) == "try-error")
            p <- NA
          p <- lpvalue(p)
        }
      }
      #set up the main columns
      onetbl <- mapply(function(sublevel,N){
        missing <- NULL
        if(is.na(sublevel[1])| sublevel[1]!="NOMAINCOVNULLNA"){
          subdata<-subset(data,subset=data[[maincov]]%in%sublevel)
        }else{subdata<-data}
        #if there is a missing in the whole data
        if(ismiss){
          n<-sum(table(subdata[[cov]]))
          missing<-N-n
        }
        # Updated LA to remove NaN from tables
        sumCov <-round(summary(subdata[[cov]]), digits)
        if (sumCov[4]=="NaN"){
          meansd <-''
          mmm <-''
          if (all.stats) mmm=c('','')
        } else {
          meansd <- paste(niceNum(sumCov["Mean"],digits), " (", niceNum(sd(subdata[[cov]], na.rm = T),digits), ")", sep = "")
          mmm <- if (IQR |all.stats) {
            if(all(c(sumCov['Median'],sumCov["1st Qu."],sumCov["3rd Qu."]) ==floor(c(sumCov['Median'],sumCov["1st Qu."],sumCov["3rd Qu."])))){
              paste(sumCov['Median'], " (", sumCov["1st Qu."], ",", sumCov["3rd Qu."],")", sep = "")
            } else {paste(niceNum(sumCov['Median'],digits), " (", niceNum(sumCov["1st Qu."],digits), ",", niceNum(sumCov["3rd Qu."],digits),")", sep = "")}
          } else {
            if(all(c(sumCov['Median'],sumCov["Min."],sumCov["Max."]) ==floor(c(sumCov['Median'],sumCov["Min."],sumCov["Max."])))){
              paste(sumCov["Median"], " (", sumCov["Min."], ",",sumCov["Max."], ")", sep = "")
            } else {paste(niceNum(sumCov['Median'],digits), " (", niceNum(sumCov["Min."],digits), ",", niceNum(sumCov["Max."],digits),")", sep = "")}
          }}
        if (all.stats) {
          mmm <- c(mmm,if(all(c(sumCov["Min."],sumCov["Max."]) ==floor(c(sumCov["Min."],sumCov["Max."])))){
            paste("(", sumCov["Min."], ",",sumCov["Max."], ")", sep = "")
          } else {paste("(", niceNum(sumCov["Min."],digits), ",", niceNum(sumCov["Max."],digits),")", sep = "")})
        }
        tbl <- c(meansd, mmm, lbld(missing))
        return(tbl)}
        ,levels,numobs[[cov]])}
    
    #Add the first column to the main columns and get the matrix ready for later
    factornames<-addspace(sanitizestr(nicename(factornames)))
    # LA Added 20 Jan 2021 to deal with one-level factors
    if (is.null(nrow(onetbl))){onetbl <- matrix(data=onetbl,ncol=length(onetbl),nrow=1) }
    
    onetbl<-cbind(factornames,onetbl)
    
    if(!is.null(maincov)){
      onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),rep("",length(levels[[1]])+1)),onetbl)
#      if(pvalue) onetbl<-cbind(onetbl,c(p,rep("",nrow(onetbl)-1)))
      if (pvalue){
        p_NA = rep("", nrow(onetbl) - 1)
        p_NA[levels(data[[cov]]) %in% excludeLevel] <-'excl'
        onetbl <- cbind(onetbl, c(p,p_NA))
        if (show.tests) onetbl <- cbind(onetbl, c(p_type,rep("", nrow(onetbl) - 1)))
      }
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
    colnm_table <- c("Covariate",paste("Full Sample (n=",N,")",sep=""),
                       mapply(function(x,y){paste(x," (n=",y,")",sep="")},
                              names(table(data[[maincov]],useNA='ifany')),table(data[[maincov]],useNA='ifany')))
    if(pvalue) colnm_table <- c(colnm_table,"p-value")
    if (show.tests) colnm_table <- c(colnm_table,'StatTest')
    colnames(table) <- colnm_table
    
  }else{
    colnames(table)<-c("Covariate",paste("n=",N,sep=""))
    
  }
  colnames(table)<-sanitizestr(colnames(table))
  
  #JW May 20 2021 adding option to remove full sample
  if(!full) table <- table[,-2] 
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
#'Options include "linear", "logistic", "coxph", "crr", "boxcox","logistic"
#'@param strata character vector of covariates to stratify by. Only used for coxph and crr
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param testing boolean to indicate if you want to print out the covariates before the model fits.
#'@param showN boolean indicating if you want to show sample sizes
#'@param CIwidth width of confidence interval, default is 0.95
#'This will allow you to see which model is not fitting if the function throws an error
#'@keywords dataframe
#' @importFrom survival coxph Surv
#'@export
uvsum <- function (response, covs, data, type = NULL, strata = 1, markup = T,
                   sanitize = T, nicenames = T, testing = F,showN=T,CIwidth=0.95)
{
  # New LA 24 Feb, test for presence of variables in data and convert character to factor
  missing_vars = setdiff(c(response,covs),names(data))
  if (length(missing_vars)>0){
    stop(paste('These covariates are not in the data:',missing_vars))
  }
  for (v in c(response,covs)) if (class(data[[v]])[1]=='character') data[[v]] <- factor(data[[v]])
  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  if (class(strata) != "numeric") {
    strataVar=strata
    strata <- sapply(strata, function(stra) {
      paste("strata(", stra, ")", sep = "")
    })
  } else {
    strataVar <-""
    strata <- ""
  }
  if (!is.null(type)) {
    if (type == "logistic") {
      beta <- "OR"
    } else if (type == "linear" | type == "boxcox") {
      beta <- "Estimate"
    } else if (type == "coxph" | type == "crr") {
      beta <- "HR"
    } else {
      stop("type must be either coxph, logisitc, linear, coxbox, crr (or NULL)")
    }
  } else {
    if (length(response) == 2) {
      if (length(unique(data[[response[2]]])) < 3) {
        type <- "coxph"
      } else {
        type <- "crr"
      }
      beta <- "HR"
    } else if (length(unique(data[[response]])) == 2) {
      type <- "logistic"
      beta <- "OR"
    } else {
      type <- "linear"
      beta <- "Estimate"
    }
  }
  beta = betaWithCI(beta,CIwidth)
  
  if (strata != "" & type != "coxph") {
    stop("strata can only be used with coxph")
  }
  out <- lapply(covs, function(x_var) {
    data <- dplyr::select(data,dplyr::any_of(c(response,x_var,strataVar)))
    data <-na.omit(data)
    x_var_str <- x_var
    if (testing)
      print(x_var)
    if (is.factor(data[[x_var]])) {
      
      data[[x_var]] <- factor(data[[x_var]],ordered = FALSE)
      levelnames = sapply(sapply(sapply(levels(data[[x_var]]),nicename),sanitizestr),addspace)
      #levelnames <- sapply(sapply(sapply(levels(factor(data[[x_var]])), nicename), sanitizestr), addspace)
      x_var_str <- lbld(sanitizestr(nicename(x_var)))
      title <- NULL
      body <- NULL
      if (type == "coxph") {
        m2 <- survival::coxph(as.formula(paste(paste("survival::Surv(", response[1],
                                                     ",", response[2], ")", sep = ""), "~", x_var,
                                               ifelse(strata == "", "", "+"), paste(strata,
                                                                                    collapse = "+"), sep = "")), data = data)
        hazardratio <- c("Reference", apply(matrix(summary(m2,conf.int=CIwidth)$conf.int[,
                                                                                         c(1, 3, 4)], ncol = 3), 1, psthr))
        pvalue <- c("", sapply(summary(m2,conf.int=CIwidth)$coef[, 5],
                               lpvalue))
        title <- c(x_var_str, "", "", lpvalue(summary(m2,conf.int=CIwidth)$waldtest[3]))
      } else if (type == "crr") {
        m2 <- crrRx(as.formula(paste(paste(response,
                                           collapse = "+"), "~", x_var, sep = "")), data = data)
        hazardratio <- c("Reference", apply(matrix(summary(m2,conf.int=CIwidth)$conf.int[,
                                                                                         c(1, 3, 4)], ncol = 3), 1, psthr))
        pvalue <- c("", sapply(summary(m2)$coef[, 5],
                               lpvalue))
        globalpvalue <- try(aod::wald.test(b = m2$coef, Sigma = m2$var,
                                           Terms = seq_len(length(m2$coef)))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        title <- c(x_var_str, "", "", lpvalue(globalpvalue))
      } else if (type == "logistic") {
        m2 <- glm(as.formula(paste(response, "~", x_var, 
                                   sep = "")), family = "binomial", data = data)
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                           Sigma = vcov(m2)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        m <- summary(m2)$coefficients
        Z_mult = qnorm(1-(1-CIwidth)/2)
        hazardratio <- c("Reference", apply(cbind(exp(m[-1,
                                                        1]), exp(m[-1, 1] - Z_mult * m[-1, 2]), exp(m[-1,
                                                                                                      1] + Z_mult * m[-1, 2])), 1, psthr))
        pvalue <- c("", sapply(m[-1, 4], lpvalue))
        title <- c(x_var_str, "", "", lpvalue(globalpvalue))
      } else if (type == "linear" | type == "boxcox") {
        if (type == "linear") {
          m2 <- lm(as.formula(paste(response, "~", x_var,
                                    sep = "")), data = data)
        } else {
          m2 <- boxcoxfitRx(as.formula(paste(response,
                                             "~", x_var, 
                                             sep = "")), data = data)
        }
        m <- summary(m2)$coefficients
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                           Sigma = vcov(m2)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        T_mult = qt(1-(1-CIwidth)/2,m2$df.residual)
        hazardratio <- c("Reference", apply(cbind(m[-1,
                                                    1], m[-1, 1] - T_mult * m[-1, 2], m[-1, 1] +
                                                    T_mult * m[-1, 2]), 1, psthr))
        pvalue <- c("", sapply(m[-1, 4], lpvalue))
        title <- c(x_var_str, "", "", lpvalue(globalpvalue))
      }
      if (length(levelnames) == 2) {
        body <- cbind(levelnames, hazardratio, c("",
                                                 ""), c("", ""))
      } else {
        body <- cbind(levelnames, hazardratio, pvalue,
                      rep("", length(levelnames)))
      }
      out <- rbind(title, body)
      if (showN){
        n_by_level = c(nrow(data),
                       as.vector(table(data[[x_var]])))
        out <- cbind(out,n_by_level)
      }
      rownames(out) <- NULL
      colnames(out) <- NULL
      return(list(out, nrow(out)))
    } else {
      x_var_str <- lbld(sanitizestr(nicename(x_var)))
      if (type == "coxph") {
        m2 <- survival::coxph(as.formula(paste(paste("survival::Surv(", response[1],
                                                     ",", response[2], ")", sep = ""), "~", x_var,
                                               ifelse(strata == "", "", "+"), paste(strata,
                                                                                    collapse = "+"), sep = "")), data = data)
        out <- matrix(c(x_var_str, psthr(summary(m2,conf.int=CIwidth)$conf.int[,
                                                                               c(1, 3, 4)]), "", lpvalue(summary(m2)$waldtest[3])),
                      ncol = 4)
        
      } else if (type == "crr") {
        m2 <- crrRx(as.formula(paste(paste(response,
                                           collapse = "+"), "~", x_var, sep = "")), data = data)
        globalpvalue <- try(aod::wald.test(b = m2$coef, Sigma = m2$var,
                                           Terms = seq_len(length(m2$coef)))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        out <- matrix(c(x_var_str, psthr(summary(m2,conf.int=CIwidth)$conf.int[,
                                                                               c(1, 3, 4)]), "", lpvalue(globalpvalue)), ncol = 4)
      } else if (type == "logistic") {
        m2 <- glm(as.formula(paste(response, "~", x_var,
                                   sep = "")), family = "binomial", data = data)
        m <- summary(m2)$coefficients
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                           Sigma = vcov(m2)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        Z_mult = qnorm(1-(1-CIwidth)/2)
        out <- matrix(c(x_var_str, psthr(c(exp(m[-1, 1]), exp(m[-1,
                                                                1] - Z_mult * m[-1, 2]), exp(m[-1, 1] + Z_mult *
                                                                                               m[-1, 2]))), "", lpvalue(globalpvalue)), ncol = 4)
      } else if (type == "linear" | type == "boxcox") {
        if (type == "linear") {
          m2 <- lm(as.formula(paste(response, "~", x_var,
                                    sep = "")), data = data)
        } else {
          m2 <- boxcoxfitRx(as.formula(paste(response,
                                             "~", x_var, sep = "")), data = data)
        }
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                           Sigma = vcov(m2)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        m <- summary(m2)$coefficients
        T_mult = qt(1-(1-CIwidth)/2,m2$df.residual)
        out <- matrix(c(x_var_str, psthr(c(m[-1, 1], m[-1,
                                                       1] - T_mult * m[-1, 2], m[-1, 1] + T_mult * m[-1,
                                                                                                     2])), "", lpvalue(globalpvalue)), ncol = 4)
      }
      if (showN){
        out<- cbind(out,nrow(data))
      }
      
      return(list(out, nrow(out)))
    }
  })
  table <- lapply(out, function(x) {
    return(x[[1]])
  })
  table <- do.call("rbind", lapply(table, data.frame, stringsAsFactors = FALSE))
  if (showN){
    colnames(table) <- sapply(c("Covariate", sanitizestr(beta),
                                "p-value", "Global p-value","N"), lbld)
    
  } else{
    colnames(table) <- sapply(c("Covariate", sanitizestr(beta),
                                "p-value", "Global p-value"), lbld)
    
  }
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
#'@param showN boolean indicating if you want to show sample sizes
#'@param CIwidth width of confidence interval, default is 0.95
#'@importFrom stats anova as.formula chisq.test coef fisher.test glm kruskal.test lm median model.matrix pchisq qnorm sd time vcov wilcox.test
#'@importFrom xtable xtable print.xtable
#'@keywords dataframe
#'@export
puvsum<-function(response,covs,data,type=NULL,strata=1,TeX=FALSE,showN=FALSE,CIwidth=0.95){
  if(!TeX){
    print.xtable(xtable(uvsum(response,covs,data,type,strata,showN = showN,CIwidth = CIwidth)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
  }else{
    print.xtable(xtable(uvsum(response,covs,data,type,strata,showN = showN,CIwidth = CIwidth)),include.rownames=F,sanitize.text.function=identity,table.placement="H",floating=FALSE,tabular.environment="longtable")
  }

}

# TODO: Add support for svyglm and svycoxph functions. May need to check the feasibility of the global p-value here
#'Get multivariate summary dataframe
#'
#'Returns a dataframe corresponding to a univariate table
#'
#'@param model fitted model object
#'@param data dataframe containing data
#'@param showN boolean indicating sample sizes should be shown for each comparison, can be useful for interactions
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param CIwidth width for confidence intervals, defaults to 0.95
#'@keywords dataframe
#'@export
mvsum <-function(model, data, showN = F, markup = T, sanitize = T, nicenames = T,CIwidth=0.95)
{
  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  call <- paste(deparse(summary(model)$call), collapse = "")
  call <- unlist(strsplit(call, "~", fixed = T))[2]
  call <- unlist(strsplit(call, ",", fixed = T))[1]
  if (substr(call, nchar(call), nchar(call)) == "\"")    call <- substr(call, 1, nchar(call) - 1)
  call <- unlist(strsplit(call, "\"", fixed = T))[1]
  call <- unlist(strsplit(call, "+", fixed = T))
  call <- unlist(strsplit(call, "*", fixed = T))
  call <- unique(call)
  call <- call[which(is.na(sapply(call, function(cov) {charmatch("strata(", cov)})) == T)]
  call <- gsub("\\s", "", call)
  type <- class(model)[1]
  
  # THis needs to be changed for multinom
  # Still to fix: doesn't work when interactions are specified as x1:x2 instead of x1*x2
  if (type=='lm'){
    betanames <- attributes(summary(model)$coef)$dimnames[[1]][-1]
    beta <- 'Estimate'
    expnt = FALSE
    ss_data <- model$model
  } else if (type=='polr'){
    betanames <- names(model$coefficients)
    beta <- "OR"
    ss_data <- model$model
  } else if (type=='lme'){
    betanames <- names(model$coef$fixed)[-1]
    beta <- 'Estimate'
    ss_data <- model$data
  } else if (type =='glm'){
    if (model$family$link=='logit'){
      beta <- "OR"
      expnt = TRUE
    } else if (model$family$link=='log') {
      beta <- "RR"
      expnt = TRUE
    } else {
      beta <- "Estimate"
      expnt = FALSE
    }
    betanames <- attributes(summary(model)$coef)$dimnames[[1]][-1]
    ss_data <- model$model
  } else if (type == "coxph" | type == "crr") {
    beta <- "HR"
    betanames <- attributes(summary(model)$coef)$dimnames[[1]]
    ss_data <- model.frame(model$call$formula,eval(parse(text=paste('data=',deparse(model$call$data)))))
  } else {
    stop("type must be either polr, coxph, logistic, lm, crr, lme (or NULL)")
  }
  
  if (missing(data)) if(class(ss_data)=='data.frame') {data=ss_data} else{ stop('data can not be derived from model')}
  
  beta = betaWithCI(beta,CIwidth)
  
  ucall = unique(call)
  indx = matchcovariate(betanames, ucall)
  data = as.data.frame(data) # LA to enable tibbles to work
  for (v in ucall[-1]){ if(class(data[[v]])[1]=='character') data[[v]]<-factor(data[[v]])} # LA change char to factor
  if (min(indx) == -1)  stop("Factor name + level name is the same as another factor name. Please change. Will fix this issue later")
  
  y <- betaindx(indx)
  if (type %in% c("lm", "glm", "lme")) {
    y <- lapply(y, function(x) {x + 1 })
    betanames <- c("intercept", betanames)
  }
  
  out <- lapply(y, function(covariateindex) {
    #Get attribute names and split by interactions
    betaname <- betanames[covariateindex]
    betaname <- strsplit(betaname, ":", fixed = T)
    # betaname <- gsub(' - ','-',betaname)  # Added LA, 14 Dec 2020
    # betaname <- gsub(' + ','-',betaname)  # Added LA, 14 Dec 2020
    # Get covariate names
    oldcovname <- covnm(betaname[[1]], call)
    oldcovname <- getvarname(oldcovname)
    # # Changed Mar 2021 to enable sample size calculations
    # levelnames <- unlist(lapply(betaname, function(level) {
    #   paste(mapply(function(lvl, cn) {
    #     # # Changed LA , Dec 15 for level names that also include the varname
    #     #   result <- unlist(strsplit(lvl, cn, fixed = T))[2]
    #     #   out <- ifelse(is.na(result), cn, result)
    #     result <- ifelse(length(grep(paste0(cn,cn),lvl))>0, unlist(sub(paste0(cn,cn),cn,lvl)),unlist(sub(cn,'',lvl)))
    #     out <- ifelse(result=='',cn,result)
    #   }, level, oldcovname), collapse = ":")
    # }))
    levelnameslist <- lapply(betaname, function(level) {
      mapply(function(lvl, cn) {
        result <- ifelse(length(grep(paste0(cn,cn),lvl))>0, unlist(sub(paste0(cn,cn),cn,lvl)),unlist(sub(cn,'',lvl)))
        out <- ifelse(result=='',cn,result)
      }, level, oldcovname)
    })
    levelnames <-unlist(lapply(levelnameslist,function(x) paste(x,collapse=':'))    )
    levelnames <- addspace(sanitizestr(nicename(levelnames)))
    covariatename <- lbld(sanitizestr(nicename(paste(oldcovname,collapse = ":"))))
    reference = NULL
    title = NULL
    body = NULL
    if (type == "lme") {
      globalpvalue <- try(aod::wald.test(b = model$coef$fixed[covariateindex],
                                         Sigma = vcov(model)[covariateindex, covariateindex],
                                         Terms = seq_along(covariateindex))$result$chi2[3])
    } else if (type == "polr") {
      globalpvalue <- try(aod::wald.test(b = model$coefficients[covariateindex],
                                         Sigma = vcov(model)[covariateindex, covariateindex],
                                         Terms = seq_along(covariateindex))$result$chi2[3])
    } else if (type != "crr") {
      globalpvalue <- try(aod::wald.test(b = coef(model)[covariateindex],
                                         Sigma = vcov(model)[covariateindex, covariateindex],
                                         Terms = seq_along(covariateindex))$result$chi2[3])
    } else {
      globalpvalue <- try(aod::wald.test(b = model$coef[covariateindex],
                                         Sigma = model$var[covariateindex, covariateindex],
                                         Terms = seq_along(covariateindex))$result$chi2[3])
    }
    if (class(globalpvalue) == "try-error")
      globalpvalue <- "NA"
    globalpvalue <- lpvalue(globalpvalue)
    if (type == "coxph" | type == "crr") {
      hazardratio <- c(apply(matrix(summary(model,conf.int=CIwidth)$conf.int[covariateindex,c(1, 3, 4)], ncol = 3), 1, psthr))
      pvalues <- c(sapply(summary(model)$coef[covariateindex,5], lpvalue))
    } else if (type == "glm" & expnt) {
      m <- summary(model,conf.int=CIwidth)$coefficients
      Z_mult = qnorm(1-(1-CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex,1]), exp(m[covariateindex, 1] - Z_mult * m[covariateindex,2]), exp(m[covariateindex, 1] + Z_mult * m[covariateindex,2])), 1, psthr)
      pvalues <- c(sapply(m[covariateindex, 4], lpvalue))
    } else if (type == "polr" ) {
      m <- summary(model)$coefficients
      Z_mult = qnorm(1-(1-CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex,1]), exp(m[covariateindex, 1] - Z_mult * m[covariateindex,2]), exp(m[covariateindex, 1] + Z_mult * m[covariateindex,2])), 1, psthr)
      pvalues =  pnorm(abs(m[covariateindex, "Value"]/m[covariateindex, "Std. Error"]),lower.tail = FALSE) * 2
      pvalues <- c(sapply(pvalues, lpvalue))
    } else if (type == "lm" | type == "glm" & !expnt) {
      T_mult = abs(qt((1-CIwidth)/2,model$df.residual))
      m <- summary(model,conf.int=CIwidth)$coefficients
      hazardratio <- apply(cbind(m[covariateindex, 'Estimate'],m[covariateindex, 'Estimate'] - T_mult * m[covariateindex,'Std. Error'], m[covariateindex, 'Estimate'] + T_mult * m[covariateindex,'Std. Error']), 1, psthr)
      pvalues <- sapply(m[covariateindex, 4], lpvalue)
    } else if (type == "lme") {
      T_mult = abs(qt((1-CIwidth)/2,summary(model)$fixDF$X))[covariateindex]
      m <- summary(model,conf.int=CIwidth)$tTable
      hazardratio <- apply(cbind(m[covariateindex, 1],m[covariateindex, 1] - T_mult * m[covariateindex,2], m[covariateindex, 1] + T_mult * m[covariateindex,2]), 1, psthr)
      pvalues <- c(sapply(m[covariateindex, 5], lpvalue))
    }
    if (length(betaname[[1]]) == 1) {
      if (!is.factor(data[, oldcovname])) {
        title <- c(nicename(covariatename), hazardratio,"", globalpvalue)
      } else if (length(levelnames) == 1) {
        title <- c(covariatename, "", "", globalpvalue)
        if (!is.null(data))
          reference <- c(addspace(sanitizestr(names(table(data[,which(names(data) == oldcovname)]))[1])),"reference", "", "")
        body <- c(levelnames, hazardratio, "", "")
      } else {
        if (!is.null(data)){
          reference <- c(addspace(sanitizestr(names(table(data[,which(names(data) == oldcovname)]))[1])),"reference", "", "")
        }
        title <- c(covariatename, "", "", globalpvalue)
        body <- cbind(levelnames, hazardratio, pvalues,rep("", length(levelnames)))
      }
    } else {
      if (length(levelnames) != 1) {
        title <- c(covariatename, "", "", globalpvalue)
        body <- cbind(levelnames, hazardratio, pvalues,
                      rep("", length(levelnames)))
      } else {
        title <- c(covariatename, hazardratio, "", globalpvalue)
      }
    }
    out <- rbind(title, reference, body)
    # New sample size work
    if (out[1,2]=="") {
      if (length(grep(':',title[1]))>0){
        ss_N = c('',unlist(lapply(levelnameslist, function(level){
          N<-mapply(function(cn,lvl){
            if (cn==lvl) {nrow(ss_data)} else {sum(ss_data[[cn]]==lvl)}
          },oldcovname,level)
          return(min(N))
        })))
      } else{ss_N = c('',table(data[[oldcovname]]))}
    } else {ss_N = nrow(ss_data)}
    out <- cbind(out,ss_N)
    #    print(out)
    rownames(out) <- NULL
    colnames(out) <- NULL
    return(list(out, nrow(out)))
  })
  table <- lapply(out, function(x) {
    return(x[[1]])
  })
  index <- unlist(lapply(out, function(x) {
    return(x[[2]])
  }))
  table<-lapply(out,function(x){return(x[[1]])})
  index<-unlist(lapply(out,function(x){return(x[[2]])}))
  table<-do.call("rbind",lapply(table,data.frame,stringsAsFactors = FALSE))
  colName = c("Covariate",sanitizestr(beta),"p-value","Global p-value")
  if (!showN) {table <-table[,-5]} else {colName = c(colName,'N')}
  colnames(table)<-sapply(colName,lbld)
  return(table)
}

#'Print multivariate summary LaTeX table
#'
#'Returns a LaTeX table of the multivariate summary.
#'
#'@param model fitted model object
#'@param data dataframe containing data
#'@param showN boolean indicating sample sizes should be shown for each comparison, can be useful for interactions
#'@param CIwidth width for confidence intervals, defaults to 0.95
#'@keywords print
#'@export

pmvsum<-function(model,data,showN=FALSE,CIwidth=0.95){
  print.xtable(xtable(mvsum(model=model,data=data,showN=showN,CIwidth=CIwidth)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
}

#' Fit and format an ordinal logistic regression using polr from the {MASS} package. The parallel regression assumption can
#' be tested using the Brant test (modfied from the the Brant package).
#'@param data dataframe containing data [REQUIRED]
#'@param covs character vector with the names of columns to include in table [REQUIRED]
#'@param response ordinal outcome variable [REQUIRED]
#'@param reflevel manual specification of the reference level, must match level exactly
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param excludeLevels a named list of levels to exclude from the response variable
#'@param testPO logical, should the proportional odds (parallel regression) assumption be tested with the Brant test, defaults to TRUE, values greater than alpha are desirable.
#'@param showN logical, should sample sizes be shown for each lvel, defaults to TRUE
#'@param digits number of digits to display, defaults to
#'@param CIwidth level of significance for computing the confidence intervals, default is 0.95
#'@return A formatted table displaying the odds ratio associated with each covariate
#'@keywords ordinal regression, Brant test
#'@importFrom MASS polr
#'@export
#'
ordsum  <- function(data, covs, response,reflevel,markup=FALSE,sanitize=TRUE,nicenames=TRUE,
                    excludeLevels,testPO=TRUE,showN=TRUE,digits=1,CIwidth=0.95){
  
  if (!markup) {
    lbld <- identity # not yet used
    addspace <- identity  # not yet used
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  
  missing_covs = setdiff(covs,names(data))
  if (length(missing_covs)>0) {
    stop(paste('Check the covarariates, the following variables are not in the data:',missing_covs))
  }
  if (!class(data[[response]])[1] %in% c('factor','ordered')) {
    warning('Response variable is not a factor, will be converted to an ordered factor')
    data[[response]] <- factor(data[[response]],ordered=T)
  }
  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  if (!missing(excludeLevels)){
    if (response %in% names(excludeLevels)){
      to_remove = sapply(data[[response]],function (x) {x %in% excludeLevels[[response]]})
      data = data[!to_remove,]
    }
  }
  if (!missing(reflevel)){
    data[[response]] <- relevel(data[[response]], ref=reflevel)
  }
  for (v in covs) { if (class(data[[v]])[1] %in% c('character','ordered')) data[[v]] <- factor(data[[v]],ordered=F)}
  out <- lapply(covs, function(x_var) {
    polr.form = as.formula(paste(response,'~',x_var))
    fit = MASS::polr(data=data,polr.form,method='logistic',Hess = TRUE)
    brant_test = try(modified_brant(fit,by.var=T),silent = T)
    if (class(brant_test)[1] == "try-error") {
      po_test_omni = data.frame(Covariate = x_var,"PO Test" = 'Not Tested')
      zero_counts <- NA
    } else {
      po_test_omni = data.frame(Covariate=rownames(brant_test$result),brant_test$result)
      po_test_omni$"PO Test" = formatp(po_test_omni$probability)
      zero_counts <-brant_test$zero_count_cells
    }
    coef_tbl <- data.frame(summary(fit)$coef)
    coef_tbl <- coef_tbl[grep(x_var,rownames(summary(fit)$coef)),]
    coef_tbl$p_value <- pnorm(abs(coef_tbl$"t.value"), lower.tail = FALSE) * 2
    coef_tbl$"p-value" = sapply(coef_tbl$p_value,lpvalue)
    coef_tbl$OR <- exp(coef_tbl$"Value")
    qstdN <- qnorm(p=1-(1-CIwidth)/2)
    coef_tbl$LB <- exp(coef_tbl$"Value"-qstdN*coef_tbl$"Std..Error")
    coef_tbl$UB <- exp(coef_tbl$"Value"+qstdN*coef_tbl$"Std..Error")
    coef_tbl$"OR_CI" = paste0(niceNum(coef_tbl$OR,digits = digits),
                              ' (',niceNum(coef_tbl$LB,digits = digits),
                              ',',niceNum(coef_tbl$UB,digits = digits),')')
    tbl <- cbind(Covariate=rownames(coef_tbl),data.frame(coef_tbl[,c('OR_CI','p-value')]))
    
    if (class(data[[x_var]])[1] %in% c("ordered", "factor" )){
      brant_test_level = try(modified_brant(model=fit,by.var=F),silent = T)
      if (class(brant_test_level)[1] == "try-error") {
        po_test_level = data.frame(Covariate = tbl,"PO Test" = 'NT')
      } else {
        po_test_level = data.frame(cbind(brant_test_level$result,Covariate=rownames(brant_test_level$result)))
        po_test_level$"PO Test" <- formatp(po_test_level$probability)
        #          po_test_level$"PO Test"[po_test_level$"PO Test"=='1'] <- 'NT' # WHY DID I ADD THIS?
      }
      
      tbl$"PO Test" = sapply(tbl$Covariate, function(x) {po_test_level$"PO Test"[po_test_level$Covariate==x]})
      tbl$Covariate = sub(x_var,'',tbl$Covariate)
      reflevel=setdiff(levels(data[[x_var]]),c(excluded_xvar,tbl$Covariate))
      tbl <- rbind(c(reflevel,"Reference","",""),tbl)
      nterms=length(fit$coefficients)
      globalp <- try(aod::wald.test(Sigma=vcov(fit)[1:nterms,1:nterms],
                                    b=fit$coefficients,Terms=1:nterms)$result$chi2[3],silent = T)
      if (class(globalp)[1]=='try-error') globalp <- NA
      tbl$"globalPval" = ''
      tbl <- rbind(c(x_var,"","",
                     po_test_omni$"PO Test"[po_test_omni$Covariate==x_var],
                     lpvalue(globalp)),tbl)
      
      n_by_level <- as.vector(sapply(tbl$Covariate[-1], function(x){ sum(fit$model[[x_var]]==x)}))
      n_by_level <- c(sum(n_by_level),n_by_level)
      tbl <- cbind(tbl,N=n_by_level)
      
      
    } else {
      tbl$"PO Test" = po_test_omni$"PO Test"[po_test_omni$Covariate==x_var]
      tbl$"globalPval" = tbl$p.value
      tbl$p.value <- ''
      
      tbl <- cbind(tbl,N=nrow(fit$fitted.values))
      
    }
    tbl <- cbind(tbl,ZeroCount=c(zero_counts,rep(NA,nrow(tbl)-1)))
    tbl <- tbl[,c("Covariate","N","OR_CI","globalPval","p.value","PO Test","ZeroCount")]
    
    names(tbl) <- c("Covariate","N","OR_CI","Global p-value","p-value","PO Test","ZeroCount")
    return(tbl)
  })
  onetbl = do.call("rbind",out)
  onetbl$Covariate <- nicename(onetbl$Covariate)
  
  if (sum(onetbl$ZeroCount>0,na.rm=T)!=0){    warning("Caution, proportional odds test performed with empty cells.")  }
  
  if (showN){
    onetbl <- onetbl[,c("Covariate","N","OR_CI","Global p-value","p-value","PO Test")]
  } else {     onetbl <- onetbl[,c("Covariate","OR_CI","Global p-value","p-value","PO Test")]}
  if (sum(onetbl$`p-value`=='')==nrow(onetbl)) {
    onetbl <- onetbl[,which(names(onetbl)=='p-value')]
  }
  
  beta = sanitizestr(betaWithCI('OR',CIwidth))
  names(onetbl) <- gsub('OR_CI',beta,names(onetbl))
  rownames(onetbl) <- NULL
  return(onetbl)
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


#' Create a forest plot using ggplot2
#' 
#'This function will accept a log or logistic regression fit from glm, and display the
#'OR or RR for each variable on the appropriate log scale.
#'
#' @param glm_fit an object output from the glm function, must be from a logistic regression
#' @param conf.level controls the width of the confidence interval
#' @param orderByRisk logical, should the plot be ordered by risk
#' @param colours can specify colours for risks less than, 1 and greater tham 1.0. Default is red, black, green
#' @param showEst logical, should the risks be displayed on the plot in text
#' @param rmRef logical, should the reference levels be removed for the plot?
#' @param logScale logical, should OR/RR be shown on log scale, defaults to TRUE
#' @param nxTicks Number of tick marks supplied to the log_breaks function to produce
#' @import ggplot2
#' @importFrom scales log_breaks
#' @keywords plot
#' @export
#'
forestplot2 = function(glm_fit,conf.level=0.95,orderByRisk=T,colours='default',showEst=TRUE,rmRef=FALSE,logScale=TRUE,nxTicks=5){
  
  if (class(glm_fit)[1]=='glm'){
    if(glm_fit$family$link=='log'){
      x_lab = 'Relative Risk'
    } else if (glm_fit$family$link=='logit'){
      x_lab='Odds Ratio'
    } else stop('glm_fit must be a logit or log link fit')
  } else {
    x_lab='Odds Ratio'
  }
  
  tab = format_glm(glm_fit,conf.level = conf.level,orderByRisk=orderByRisk)
  if (rmRef) tab = tab[setdiff(1:nrow(tab),which(tab$estimate.label=='1.0 (Reference)')),]
  
  
  yvals=1
  for (i in 2:nrow(tab)) {
    yvals = c(yvals, ifelse(tab$var.order[i]==tab$var.order[i-1],
                            yvals[i-1]+.5,
                            yvals[i-1]+1))
  }
  tab$estimate.label = ifelse(is.na(tab$estimate.label),'',tab$estimate.label)
  tab$estimate.label = ifelse(tab$estimate.label == '1.0 (Reference)','(Reference)',tab$estimate.label)
  
  if (showEst){
    yLabels = data.frame(y.pos=yvals,
                         labels=ifelse(is.na(tab$level.name),
                                       paste(tab$variable,tab$estimate.label),
                                       paste(tab$level.name,tab$estimate.label)))
  } else {
    yLabels = data.frame(y.pos=yvals,
                         labels=ifelse(is.na(tab$level.name),
                                       tab$variable,
                                       ifelse(tab$estimate.label == '(Reference)',
                                              paste(tab$level.name,tab$estimate.label),
                                              tab$level.name)
                         ))
  }
  yLabels$labels <- gsub('_',' ',yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos),]
  # TODO: add indents using '  ' to category labels, omit RR estimates?, make hjust=0
  
  tab$x.val = ifelse(tab$estimate.label == '(Reference)',1,tab$estimate)
  tab$y.val = yLabels$y.pos
  
  # set colours
  tab$colour <- ifelse(tab$x.val<1,'a',ifelse(tab$x.val==1,'b','c'))
  
  if (colours=='default'){
    colours = c(a='red',b='black',c='darkgreen')
  }  else {
    names(colours) = c('a','b','c')
  }
  
  # ensure that colours are always red, black, green
  colours <- colours[sort(unique(tab$colour))]
  
  p = ggplot(tab, aes_(x=~x.val,y=~y.val,colour=~colour))+
    geom_point(na.rm=TRUE,size=2) +
    geom_errorbarh(aes_(xmin = ~conf.low, xmax = ~conf.high),
                   height  = 0,
                   size   = 0.9,
                   na.rm=TRUE) +
    geom_vline(xintercept = 1.0) +
    labs(y='',x=x_lab) +
    guides(colour='none')+
    scale_y_continuous(breaks = yLabels$y.pos,labels=yLabels$labels) +
    scale_colour_manual(values=colours)+
    theme_bw() +
    theme(axis.text.y = element_text(hjust=0),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
  
  
  if (logScale) p + scale_x_log10(breaks=scales::log_breaks(n=nxTicks)) else p
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

      invisible(utils::capture.output(fit <-  cuminc( data[,response[1]],data[,response[2]] )))
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
      invisible(utils::capture.output(fit <- cuminc(data[,response[1]],data[,response[2]], newgpvar)))
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

#' Plot bivariate relationships in a single plot
#' 
#' This function is designed to accompany \code{\link{uvsum}} as a means
#' of visualising the results. It is not designed for publication, but rather
#' to facilitate explanations of the uvsum output.
#' 
#' The following hierarchy of how to display plots is used:
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
plotuv <- function(response,covs,data,showN=FALSE,na.rm=TRUE,response_title=NULL){
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
        theme_bw() +
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
  ggpubr::ggarrange(plotlist=plist,
            common.legend = T,
            ncol=2,
            nrow=ceiling(length(plist)/2))
}

# Rmarkdown Reporting --------------------------------------------------------------

#' Generate a bibfile for an Rmd document from a master *.bib file
#'
#' This function will search through the current Rmd document for citations and
#' extract these from a central bib file to save as a file-specific bibfile.
#' It will also search through for citations to R packages and include those,
#' which can be useful to keep track of which version of a package was used.
#' @param bibfile a master bib file contianing all the references in the document (Zotero or Mendeley).
#' @param outfile a filename to write the bibfile to. If missing this will be the rmd file name with a bib extension
#' @importFrom bib2df df2bib
#' @importFrom knitr write_bib
#' @importFrom rstudioapi getSourceEditorContext
#' @keywords citations
#' @export
rmdBibfile <- function(bibfile,outfile){
  # First sanitise the specified bibfile
  if (!file.exists(bibfile)) stop(paste(bibfile,'does not exist.'))
  
  if (file.access(bibfile)==0){
    bib <- bib_ReadGatherTidy(bibfile)
  } else stop(paste('Can not access',bibfile))
  
  # Replace some troublesome characters
  bib$URL <- gsub("\\{\\\\_\\}","_", bib$URL)
  
  # remove abstract, keywords
  rm <- c(which(names(bib)=='ABSTRACT'),which(names(bib)=='KEYWORDS'))
  if (length(rm) >0) bib <- bib[,-rm]
  
  thisFile = rstudioapi::getSourceEditorContext()$path
  if (thisFile=="") stop('Please save the current file before running writeBibFile.')
  thisDir = dirname(thisFile)
  
  # open current file and search for [@xx] entries
  fileWords = scan(thisFile,what='character')
  refs = fileWords[grep("\\[@.*\\]",fileWords)]
  stripped = gsub(".*\\[|\\].*","",refs)
  allrefs = unique(unlist(strsplit(stripped,';')))
  rRefs = allrefs[grep("R-",allrefs)]
  otherRefs = setdiff(allrefs,rRefs)
  Rpckgs = gsub("@R-","",rRefs)
  sepRefs = unlist(strsplit(otherRefs,";"))
  
  tidyRefs <-gsub(".*@","",sepRefs)
  tidyRefs <-gsub("[.]| .*|,","",tidyRefs)
  libRefs = unique(tidyRefs)
  
  # Extract the references from the master library
  paperBib <- NULL
  for ( p in libRefs){
    ind = which(bib$BIBTEXKEY==p)
    if (length(ind)>0) {
      paperBib <- rbind(paperBib,bib[ind,])
    } else warning(paste(p,'not in',bibfile,'\n'))
  }
  
  if (missing(outfile)){
    bib_out = paste0(thisDir,gsub(".Rmd","",gsub(thisDir,"",thisFile)),".bib")
  } else bib_out = paste0(thisDir,"/",gsub(".bib","",outfile),".bib")
  
  # write R packages to bibfile and append remaining references
  if (length(Rpckgs)>0){
    add = TRUE
    knitr::write_bib(Rpckgs,file =bib_out)
  } else {add=FALSE}
  
  if (!is.null(paperBib)){
    bib2df::df2bib(paperBib,file = bib_out,append = add)
  } else { if (add==FALSE) warning('No citations in current file to write')}
  
}


#' Print tables to PDF/Latex HTML or Word
#'
#' Output the table nicely to whatever format is appropriate. 
#' This is the output function used by the rm_* printing functions.
#' @param tab a table to format
#' @param to_indent numeric vector the length of nrow(tab) indicating which rows to indent
#' @param to_bold numeric vector the length of nrow(tab) indicating which rows to bold
#' @param caption table caption
#' @param chunk_label only used knitting to Word docs to allow cross-referencing
#' @export
outTable <- function(tab,to_indent=numeric(0),to_bold=numeric(0),caption=NULL,chunk_label){
  
  # strip tibble aspects
  tab=as.data.frame(tab)
  rownames(tab) <- NULL
  
#  out_fmt = ifelse(is.null(getOption('doc_type')),'pdf',getOption('doc_type'))
  out_fmt = ifelse(is.null(knitr::pandoc_to()),'html',
                           ifelse(knitr::pandoc_to(c('doc','docx')),'doc',
                           ifelse(knitr::is_latex_output(),'latex','html')))
  # if (out_fmt=='tblOnly'){
  #   tab
  # } else{
  #   
  #    out_fmt =ifelse(out_fmt%in%c('doc','docx'),'doc','pdf')
    chunk_label = ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label)
    
    if (is.null(to_indent)) to_indent = numeric(0)
    to_indent = as.vector(to_indent)
    
    
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
      # set NA to empty in kable
      options(knitr.kable.NA = '')
      if (nrow(tab)>30){
        kout <- knitr::kable(tab, format = out_fmt,
                             booktabs=TRUE,
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
        kout <- knitr::kable(tab, format = out_fmt,
                             booktabs=TRUE,
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
    # }
  
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
#'@param tableOnly should a dataframe or a formatted print object be returned
#'@param covTitle character with the names of the covariate column
#'@param chunk_label only used if output is to Word to allow cross-referencing
#'@param ... additional options passed to function  \code{\link{covsum}}
#'@keywords dataframe
#'@return A formatted table displaying a summary of the covariates stratified by maincov
#'@export
#'@seealso \code{\link{fisher.test}}, \code{\link{chisq.test}}, \code{\link{wilcox.test}}, \code{\link{kruskal.test}}, and \code{\link{anova}}
rm_covsum <- function(data,covs,maincov=NULL,caption=NULL,tableOnly=FALSE,covTitle='Covariate',chunk_label,...){

  tab <- covsum(data,covs,maincov,markup=FALSE,sanitize=FALSE,...)
  to_bold = numeric(0)
  if ('p-value' %in% names(tab)) {
    # format p-values nicely
    p_vals <- tab[['p-value']]
    new_p <- sapply(p_vals,formatp)
    tab[['p-value']] <- new_p
    to_bold <- which(suppressWarnings(as.numeric(p_vals))<0.05)
  } 
  nice_var_names = gsub('[_.]',' ',covs)
  to_indent <- which(!tab$Covariate %in% nice_var_names )
  # if (IQR) {
  #   tab$Covariate <- gsub('[(]Q1,Q3[)]','(IQR)',tab$Covariate)
  # }
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

#' Output several univariate models nicely in a single table
#' 
#' Wrapper for the uvsum function for use with Rmarkdown.
#'
#' @param response string vector with name of response
#' @param covs character vector with the names of columns to fit univariate models to
#' @param data dataframe containing data
#' @param caption table caption
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param removeInf boolean indicating if infinite estimates should be removed from the table
#' @param HolmGlobalp boolean indicting if a Holm-corrected p-value should be presented
#' @param chunk_label only used if output is to Word to allow cross-referencing
#'@param ... additional options passed to function  \code{\link{uvsum}}
#' @export
#'
rm_uvsum <- function(response, covs , data ,caption=NULL,tableOnly=FALSE,removeInf=T,HolmGlobalp=FALSE,chunk_label,...){

  # get the table
  tab <- uvsum(response,covs,data,markup = FALSE,sanitize=FALSE,...)

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

  to_bold <- which(suppressWarnings(as.numeric(p_sig))<0.05)
  nice_var_names = gsub('[_.]',' ',covs)
  to_indent <- which(!tab$Covariate %in% nice_var_names )

  tab[["p-value"]] <- formatp(tab[["p-value"]])
  tab[["Global p-value"]] <- formatp(tab[["Global p-value"]])
  if (HolmGlobalp){
    tab[["Holm Adj p"]] <- formatp(tab[["Holm Adj p"]])
  }

  # If all outcomes are continuous (and so all p-values are NA), remove this column & rename Global p-value to p-value
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


#' Fit and format an ordinal logistic regression using polr from the {MASS} package. The parallel regression assumption can
#' be tested using the Brant test in the Brant package and visually. Only logistic ordinal regression is supported currently.
#'@param data dataframe containing data [REQUIRED]
#'@param covs character vector with the names of columns to include in table [REQUIRED]
#'@param response ordinal outcome variable [REQUIRED]
#'@param reflevel manual specification of the reference level, must match level exactly
#'@param caption Table caption
#'@param showN logical, should sample sizes be shown for each lvel, defaults to TRUE
#'@param excludeLevels a named list of levels to exclude from factor variables. Currently, this has only been implemented for the response variable.
#'@param testPO logical, should the proportional odds (parallel regression) assumption be tested with the Brant test, defaults to TRUE
#'@param digits number of digits to display, defaults to
#'@param CIwidth level of significance for computing the confidence intervals, default is 0.95
#'@param excludeLevels a named list of levels to exclude from the response variable
#'@param chunk_label only used if output is to Word to allow cross-referencing
#'@return A formatted table displaying the odds ratio associated with each covariate
#'@keywords ordinal regression, Brant test
#'@export
#'
rm_ordsum <- function(data, covs, response, reflevel,caption = NULL, showN=T,
                      testPO=TRUE,digits=2,CIwidth=0.95,excludeLevels=NULL,chunk_label){
  
  
  tab <- ordsum(data=data,
                covs=covs,
                response=response,
                reflevel=reflevel,
                markup=FALSE,
                sanitize=FALSE,
                nicenames=T,
                testPO=testPO,
                showN=showN,
                digits = digits,
                CIwidth=CIwidth,
                excludeLevels=excludeLevels)
  
  if (is.null(caption))
    caption = paste('Univariate ordinal logistic regression analysis of predictors of',nicename(response),'.')
  
  # format p-values nicely
  tab$`Global p-value` <- sapply(tab$`Global p-value`,formatp)
  if (length(which(names(tab)=='p-value'))>0)
    tab[,which(names(tab)=='p-value')] <- sapply(tab[[which(names(tab)=='p-value')]],formatp)
  
  nice_var_names = gsub('_',' ',covs)
  to_indent <- which(!tab$Covariate %in% nice_var_names )
  to_bold <- which(as.numeric(tab[["Global p-value"]])<(1-CIwidth))
  
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
#' @param showN boolean indicating sample sizes should be shown for each comparison, can be useful for interactions
#' @param CIwidth width for confidence intervals, defaults to 0.95
#' @param caption table caption
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param HolmGlobalp boolean indicting if a Holm-corrected p-value should be presented
#' @param chunk_label only used if output is to Word to allow cross-referencing
#' @export
rm_mvsum <- function(model , data ,showN=FALSE,CIwidth=0.95,caption=NULL,tableOnly=FALSE,HolmGlobalp=FALSE,chunk_label){

  # get the table
  tab <- mvsum(model=model,data=data,markup = FALSE, sanitize = FALSE, nicenames = T,showN=showN,CIwidth = CIwidth)


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

#'Print Event time summary
#'
#'Wrapper for the etsum function that prints paragraphs of text in R Markdown
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
rm_etsum<-function(data,response,group=1,times=c(12,14),units="months"){
  t<-etsum(data,response,group,times)
  
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
    
    
    out<-paste0(ifelse(strta[i]=='','',paste('**',sanitizestr(nicename(strta[i])),'**'))," There are ",t[i,1]," patients. There were ",t[i,2],
                " (",round(100*t[i,2]/t[i,1],0),sanitizestr("%"),") events. The median and range of the follow-up times is ",
                t[i,6]," (",t[i,7],"-",t[i,8],") ",units,". ",km,flet,ps,sep="")
    cat("\n",out,"\n")
  })
}


#' Plot KM and CIF curves with ggplot
#'
#'This function will plot a KM or CIF curve with option to add the number at risk. You can specify if you want
#'confidence bands, the hazard ratio, and pvalues, as well as the units of time used.
#'
#'
#'
#' @param response character vector with names of columns to use for response
#' @param cov String specifying the column name of stratification variable
#' @param data dataframe containing your data
#' @param type string indicating he type of univariate model to fit. The function will try and guess what type you want based on your response. If you want to override this you can manually specify the type.
#'Options include "KM", and ,"CIF"

#' @param pval boolean to specify if you want p-values in the plot (Log Rank test for KM and Gray's test for CIF)
#' @param HR boolean to specify if you want hazard ratios included in the plot
#' @param HR_pval boolean to specify if you want HR p-values in the plot
#' @param conf.type One of "none"(the default), "plain", "log" , "log-log" or "logit". Only enough of the string to uniquely identify it is necessary. The first option causes confidence intervals not to be generated. The second causes the standard intervals curve +- k *se(curve), where k is determined from conf.int. The log option calculates intervals based on the cumulative hazard or log(survival). The log-log option bases the intervals on the log hazard or log(-log(survival)), and the logit option on log(survival/(1-survival)).
#' @param table Logical value. If TRUE, includes the number at risk table
#' @param times Numeric vector of times for the x-axis
#' @param xlab String corresponding to xlabel. By default is "Time"
#' @param ylab String corresponding to ylabel. When NULL uses "Survival probability" for KM cuves, and "Probability of an event" for CIF
#' @param main String corresponding to main title. When NULL uses Kaplan-Meier Plot s, and "Cumulative Incidence Plot for CIF"

#' @param stratalabs string corresponding to the labels of the covariate, when NULL will use the levels of the covariate
#' @param strataname String of the covariate name default is  nicename(cov)
#' @param stratalabs.table String corresponding to the levels of the covariate for the number at risk table, when NULL will use the levels of the covariate. Can use a string of "-" when the labels are long
#' @param strataname.table String of the covariate name for the number at risk table default is  nicename(cov)

#' @param median.text boolean to specify if you want the median values added to the legend (or as added text if there are no covariates)
#' @param median.lines boolean to specify if you want the median values added as lines to the plot
#' @param set.time.text string for the text to add survival at a specified time (eg. 5year OS)
#' @param set.time.line boolean to specify if you want the survival added as lines to the plot at a specified point
#' @param set.time Numeric values of the specific time of intrest, default is 5



#' @param censor.marks logical value. If TRUE, includes censor marks (only for KM curves)
#' @param censor.size size of censor marks, default is 3
#' @param censor.stroke stroke of censor marks, default is 1.5
#' @param fsize font size
#' @param nsize font size for numbers in the numbers at risk table
#' @param lsize line size
#' @param psize size of the pvalue
#' @param median.size size of the median text (Only when there are no covariates)
#' @param median.pos vector of length 2 corresponding to the median position (Only when there are no covariates)
#' @param median.lsize line size of the median lines
#' @param set.size size of the survival at a set time text (Only when there are no covariates)
#' @param set.pos  vector of length 2 corresponding to the survival at a set point position (Only when there are no covariates)
#' @param set.lsize line size of the survival at set points
#' @param ylim vector of length 2 corresponding to limits of y-axis. Default to NULL
#' @param col vector of colours
#' @param linetype vector of line types
#' @param xlim  vector of length 2 corresponding to limits of x-axis. Default to NULL
#' @param legend.pos Can be either a string corresponding to the legend position ("left","top", "right", "bottom", "none") or a vector of length 2 corresponding to the legend position (uses normalized units (ie the c(0.5,0.5) is the middle of the plot))
#' @param pval.pos  vector of length 2 corresponding to the p-value position
#' @param plot.event  Which event(s) to plot (1,2, or c(1,2))
#' @param event String specifying if the event should be mapped to the colour, or linetype when plotting both events to colour = "col", line type
#' @param flip.CIF boolean to flip the CIF curve to start at 1
#' @param cut numeric value indicating where to divide a continuous covariate (default is the median)
#' @param eventlabs String corresponding to the event type names
#' @param event.name String corresponding to the label of the event types
#' @param Numbers_at_risk_text String for the label of the number at risk
#' @param HR.digits Number of digits printed of the  hazard ratio
#' @param HR.pval.digits Number of digits printed of the hazard ratio pvalue
#' @param pval.digits Number of digits printed of the Gray's/log rank pvalue

#' @param median.digits Number of digits printed of the median pvalue
#' @param set.time.digits Number of digits printed of the probability at a specified time
#' @param print.n.missing Logical, should the number of missing be shown !Needs to be checked
#' @param returns Logical value returns a list with all ggplot objects in a list
#'@export
ggkmcif <- function(response,cov=NULL,data,type=NULL,
                    pval = TRUE,HR=F,HR_pval=F, conf.type = "none",table = TRUE,times = NULL,xlab = "Time",ylab=NULL ,
                    main = NULL,stratalabs = NULL,strataname = nicename(cov),
                    stratalabs.table=NULL,strataname.table=strataname,
                    median.text=F,median.lines=F,set.time.text=NULL,set.time.line=F,set.time=5,
                    censor.marks = TRUE,censor.size = 3,censor.stroke = 1.5,
                    fsize = 15, nsize = 4.5, lsize = 1.5, psize = 4.5,
                    median.size=4.5,median.pos=NULL,median.lsize=1.1,
                    set.size=4.5,set.pos=NULL,set.lsize=1.1,
                    ylim=c(0,1), col=NULL,linetype=NULL, xlim=NULL,
                    legend.pos = NULL,  pval.pos=NULL,plot.event=1,event=c("col","linetype"),flip.CIF =F,
                    cut=NULL,eventlabs=NULL,event.name=NULL,Numbers_at_risk_text="Numbers at risk", HR.digits = 2,HR.pval.digits=3, pval.digits=3,
                    median.digits=3,set.time.digits=3,returns = FALSE,print.n.missing=T){
  event <- match.arg(event)
  
  
  
  ##Removing all missing data
  #rowSums(is.na(final[ , 5:6])) == 0,
  if(!is.null(cov)) remove <- rowSums(is.na(data[ ,c(response,cov)])) > 0
  if(is.null(cov)) remove <- rowSums(is.na(data[ ,c(response)])) > 0
  
  data <- data[!remove,]
  
  if(print.n.missing==T & sum(remove) > 0) print(paste(sum(remove),"observations have been removed due to missing data"))
  
  
  if (!is.factor(data[,cov])  & !is.numeric(data[,cov]) & !is.null(cov)){ message("Coercing the cov variable to factor"); data[,cov] <- factor(data[,cov])}
  
  
  if (is.numeric(data[,cov])){
    numeric = T
    if(is.null(cut)) cut <- median(data[,cov],na.rm = T)
    data[,cov] <- factor(ifelse(data[,cov]<=cut,paste0("<=",round(cut,2)),paste0(">",round(cut,2))),levels = c(paste0("<=",round(cut,2)),paste0(">",round(cut,2))))
    if(is.null(stratalabs)) stratalabs = c(paste0("\u2264",round(cut,2)),paste0(">",round(cut,2)))
  }
  
  data[,cov] <- droplevels(data[,cov])
  if(is.null(stratalabs)) stratalabs <- levels(data[,cov])
  if(is.null(stratalabs.table)) stratalabs.table <- stratalabs
  
  if(!is.null(col) & !is.null(cov) & (event != "col" | length(plot.event)==1)) col <- rep(col,length.out=length(levels(data[,cov])))
  if(!is.null(linetype) & !is.null(cov) & (event != "linetype" | length(plot.event)==1)) linetype <- rep(linetype,length.out=length(levels(data[,cov])))
  
  
  if(is.null(col)){
    if(length(plot.event)==1|event != "col"){
      col_length <- ifelse(is.null(cov),1,length(levels(data[,cov])))
    }else col_length = 2
    
    col <- color_palette_surv_ggplot(col_length)
  }
  
  # Specifing the type of plot ----------------------------------------------
  
  
  #Specifying KM or CIF & is.null(type)
  if (length(unique(data[, response[2]]))< 3 & is.null(type)) type = "KM"
  if (length(unique(data[, response[2]]))>= 3 & is.null(type)) type = "CIF"
  if(type=="KM") {
    if(is.null(main)) main <- "Kaplan-Meier Plot"
    if(is.null(ylab)) ylab = "Survival Probability"
  }else if(type=="CIF"){
    #For the number at risk calcuation
    #data_sfit[,response[2]][data_sfit[,response[2]]!=0] <- 1
    if(is.null(main)) main <- "Cumulative Incidence Plot"
    if(is.null(ylab)) ylab = "Probability of an Event"
  }else(stop("Type must be either KM or CIF"))
  
  
  if(flip.CIF==T & type=="CIF") {
    median.text=F
    median.lines=F
    set.time.text=NULL
    set.time.line=F
    set.time=NULL
  }
  
  # Labels ------------------------------------------------------------------
  
  multiple_lines <- !is.null(cov)
  
  
  
  # HR and p-val cox----------------------------------------------------------------------
  if(type=="KM" & multiple_lines & (HR|HR_pval)){
    coxfit <- survival::coxph(as.formula(paste(paste("survival::(", response[1],
                                           ",", response[2], ")", sep = ""), "~", cov,
                                     sep = "")), data = data)
    
    
    HR_vals <- paste0("HR=",sapply(seq(length(stratalabs)-1),function(i){
      return(psthr0(summary(coxfit)$conf.int[i,c(1, 3, 4)],digits=HR.digits))
    }))
    
    if(HR)stratalabs[-1] <- paste(stratalabs[-1],HR_vals)
    if(HR_pval) stratalabs[-1] <- paste(stratalabs[-1],sapply(summary(coxfit)$coef[,5],lpvalue2,digits=HR.pval.digits))
    stratalabs[1] <- paste(stratalabs[1],"REF")
  }
  
  
  # HR and p-val crr --------------------------------------------------------
  if(type=="CIF" & multiple_lines & (HR|HR_pval)&length(plot.event==1) & plot.event[1]==1){
    crrfit <- crrRx(as.formula(paste(paste(response,
                                           collapse = "+"), "~", cov, sep = "")),
                    data = data)
    
    HR_vals <- paste0("HR=",sapply(seq(length(stratalabs)-1),function(i){
      return(psthr0(summary(crrfit)$conf.int[i,c(1, 3, 4)],digits=HR.digits))
    }))
    
    if(HR)stratalabs[-1] <- paste(stratalabs[-1],HR_vals)
    if(HR_pval) stratalabs[-1] <- paste(stratalabs[-1],sapply(summary(crrfit)$coef[,5],lpvalue2,digits=HR.pval.digits))
    stratalabs[1] <- paste(stratalabs[1],"REF")
    
  }
  # Model fitting KM and creating a dataframe of times--------------------------------------------------------
  
  if(type=="KM"){
    if(!multiple_lines){
      
      sfit <- survfit(as.formula(paste(paste("survival::Surv(", response[1],
                                             ",", response[2], ")", sep = ""), "~", 1,
                                       sep = "")), data = data,conf.type=conf.type)
      
      if(median.lines==T|median.text==T) median_vals <- summary(sfit)$table['median']
      if(!is.null(set.time.text)|set.time.line==T){
        
        sum <- summary(sfit,time=c(0,set.time))
        set.surv <- sum$surv[sum$time==set.time]
        set.surv <- ifelse(length(set.surv)==0,'NA',set.surv)
        
        
      }
      stratalabs <- "All"
      
    }else{
      sfit <- survfit(as.formula(paste(paste("survival::Surv(", response[1],
                                             ",", response[2], ")", sep = ""), "~", cov,
                                       sep = "")), data = data,conf.type=conf.type)
      
      if(median.lines==T|median.text==T) {
        median_vals <- summary(sfit)$table[,'median']
        if(median.text==T) stratalabs <- paste(stratalabs,", Median=",round_sprintf(median_vals,digits=median.digits))
      }
      
      if(!is.null(set.time.text)|set.time.line==T) {
        sum <- summary(sfit,time=c(0,set.time)) ##The zero is to make sure all levels are included
        df_text <- data.frame(strat=sum$strata,surv=sum$surv,time=sum$time)
        df_text$surv[df_text$time==0] <- -1
        
        
        set.surv <- stats::aggregate(df_text$surv~df_text$strat,FUN=max)[,2]
        set.surv[set.surv<0] <- NA 
        
        if(!is.null(set.time.text)) stratalabs <- paste0(stratalabs,", ",set.time,"",set.time.text,"=",round_sprintf(set.surv,digits = set.time.digits))
      }
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
    levels(df$strata) <- stratalabs
    zeros <- data.frame(time = 0, surv = 1,
                        strata = if(multiple_lines){
                          levels(df$strata)
                        }else factor("All"),
                        upper = 1, lower = 1)
    df <- plyr::rbind.fill(zeros, df) # Forcing the curves to start at 1
    
    df$strata <- factor(df$strata,levels=stratalabs)
  }
  
  
  # Model fitting CIF -------------------------------------------------------
  if(type=="CIF"){
    
    ### Functions for getting the median values
    median_time_to_event <- function(name){
      estall=get(name,fit)$est
      timall=get(name,fit)$time
      return(timall[estall>=0.5][1])
      
    }
    
    if(!multiple_lines){
      
      invisible(utils::capture.output(fit <-  cuminc( data[,response[1]],data[,response[2]] )))
      stratalabs <- " "
      gsep = " "
      
      last_character <- substr(names(fit), nchar(names(fit)), nchar(names(fit)))
      get_values <- names(fit)[last_character %in% plot.event]
      if(median.lines==T|median.text==T) median_vals <- sapply(get_values, median_time_to_event)
      if(!is.null(set.time.text)|set.time.line==T) set.surv <- cmprsk::timepoints(fit,times=set.time)$est[rownames(cmprsk::timepoints(fit,times=set.time)$est) %in% get_values,]
      
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
      invisible(utils::capture.output(fit <- cuminc(data[,response[1]],data[,response[2]], newgpvar)))
      gsep = ": "
      
      
      last_character <- substr(names(fit), nchar(names(fit)), nchar(names(fit)))
      get_values <- names(fit)[last_character %in% plot.event]
      
      if(median.lines==T|median.text==T) median_vals <- sapply(get_values, median_time_to_event)
      if(median.text==T & length(plot.event)==1) stratalabs <- paste(stratalabs,", Median=",round_sprintf(median_vals,digits=median.digits))
      
      if(!is.null(set.time.text)|set.time.line==T) set.surv <- cmprsk::timepoints(fit,times=set.time)$est[rownames(cmprsk::timepoints(fit,times=set.time)$est) %in% get_values,]
      if(!is.null(set.time.text)& length(plot.event)==1) stratalabs <- paste0(stratalabs,", ",set.time,"",set.time.text,"=",round_sprintf(set.surv,digits=set.time.digits))
      
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
    levels(df$strata) <- stratalabs
    
    if(multiple_lines){df$strata <- factor(df$strata,levels=stratalabs)
    }else df$strata <- "ALL"
    
    df$std <- std <- sqrt(df$var)
    
    names(df)[names(df)=="est"] <- 'surv'
    
    df <- df[df$event %in% plot.event,]
    df$upper <- NA
    df$lower <- NA
    
    
    
    if(conf.type!="none") {
      CIF_confint <- survfit_confint(p=df$surv, se=df$std, conf.type=conf.type, conf.int=0.95,ulimit=TRUE,logse=F)
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
    }else eventlabs <- c(1,2)
    
    
    
  }
  
  
  # axis and legend positions --------------------------------------------------------------------
  m <- max(nchar(stratalabs))
  maxxval = max(df$time,times[length(times)])
  if( is.null(xlim) ){
    maxxlim =  maxxval
  }else {
    if( length(xlim)!=2 | !is.numeric(xlim) ) stop("xlim must be a numeric vector of length 2.")
    maxxlim = xlim[2]
  }
  
  
  linetype_name <- col_name <- strataname
  
  leg.pos <- legend.pos
  d <- length(levels(df$strata))
  
  
  if(!multiple_lines & (type=="KM" | (type=="CIF" & length(plot.event)==1))){
    leg.pos <- 'none'
  }else if(is.null(legend.pos)) {
    if(type=="CIF" & flip.CIF==F){ leg.pos <- c(min(0.1+m/200,0.5), 0.9-d*0.05)
    }else leg.pos <- c(min(0.1+m/200,0.5), 0.05+d*0.05)
  }
  
  
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
    
    if(type=="CIF" & (event == "linetype"|length(plot.event)==1)){
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
    scale_x_continuous(paste0("\n",xlab), breaks = times,
                       limits = c(0, maxxval)) +
    coord_cartesian(xlim=c(0,maxxlim)) + ### changes the actual plotted limits if needed
    scale_y_continuous(paste0(ylab,"\n"), limits = ylim) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.key = element_rect(colour = "transparent", fill = "transparent")) +
    theme(legend.background=element_blank()) +
    labs(linetype = linetype_name, col=col_name) +
    ggtitle(main) + theme(plot.title = element_text(face="bold",hjust = 0.5))+
    theme(legend.position = leg.pos)
  
  
  # censoring ---------------------------------------------------------------
  
  
  if( censor.marks & type=="KM") p <- p + geom_point(data = subset(df, n.censor>0),
                                                     aes(x=time, y=surv, group = strata, col=strata),
                                                     shape=3, size=censor.size, stroke = censor.stroke,show.legend =F)
  
  
  # Log rank p-val ----------------------------------------------------------
  if(pval & type=="KM" & multiple_lines) {
    sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    pval <- pchisq(sdiff$chisq, length(sdiff$n)-1, lower.tail = FALSE)
    pvaltxt <- lpvalue2(pval,pval.digits)
    pvaltxt <- paste(pvaltxt,"(Log Rank)")
    #pvaltxt <- paste("p =", signif(pval, 3), "(Log Rank)")
    if(is.null(pval.pos)){
      p <- p + annotate("text", x = 0.85 * max(times), y = ylim[1], label = pvaltxt, size = psize)
    }else p <- p + annotate("text", x = pval.pos[1], y = pval.pos[2], label = pvaltxt, size = psize)
  }
  
  
  # Gray's pvalue -----------------------------------------------------------
  
  if( !is.null(cov) & pval & type=="CIF") {
    if(length(plot.event)==1 ){
      test <- test[rownames(test)==plot.event,]
      pval <- test[2]
      pvaltxt <- lpvalue2(pval,pval.digits)
      pvaltxt <- paste(pvaltxt,"(Gray's test)")
      if(is.null(pval.pos)){
        p <- p + annotate("text", x = 0.85 * max(times), y = ylim[1], label = pvaltxt, size = psize)
      }else p <- p + annotate("text", x = pval.pos[1], y = c(pval.pos[2]), label = pvaltxt, size = psize)
    }else{
      
      
      
      pval <- test[,2]
      pvaltxt <- sapply(pval,lpvalue2,pval.digits)
      pvaltxt <- c("Gray's test",paste(eventlabs, pvaltxt))
      if(is.null(pval.pos)){
        
        p <- p + annotate("text", x = 0.85 * max(df$time), y = c(0.12,0.08,0.04), label = pvaltxt)
      }else p <- p + annotate("text", x = pval.pos[1], y = c(pval.pos[2],pval.pos[2]-0.04, pval.pos[2]-0.08), label = pvaltxt)
    }
  }
  
  
  # median with one level text ---------------------------------------------------
  
  if(median.text &!multiple_lines) {
    median_txt <- paste0("Median=",round_sprintf(median_vals,digits=median.digits))
    if(length(plot.event)==2) median_txt <- paste(paste0(eventlabs,":",median_txt),collapse="\n")
    
    
    if(is.null(median.pos) & type=="KM")  median.pos <- c( 0.1 * max(times),ylim[1])
    if(is.null(median.pos) & type=="CIF")  median.pos <- c( 0.1 * max(times),ylim[2]*.95)
    
    p <- p + annotate("text", x =median.pos[1], y = median.pos[2], label = median_txt, size = median.size)
    
  }
  
  if(!is.null(set.time.text) &!multiple_lines) {
    set.surv.text <- paste0(set.time,"",set.time.text,"=",round_sprintf(set.surv,digits=set.time.digits))
    if(length(plot.event)==2) set.surv.text <- paste(paste0(eventlabs,":",set.surv.text),collapse="\n")
    
    if(is.null(set.pos) & type=="KM")  set.pos <- c( 0.1 * max(times),ylim[1]+0.1)
    if(is.null(set.pos) & type=="CIF")  set.pos <- c( 0.1 * max(times),ylim[2]*.85)
    
    
    p <- p + annotate("text", x =set.pos[1], y = set.pos[2], label = set.surv.text, size = set.size)
    
  }
  
  
  # median lines ------------------------------------------------------------
  
  if(median.lines){
    temp <- data.frame(x=median_vals,y=0.5)
    p <- p + geom_segment(data=temp,aes(x=x,xend=x,y=0,yend=0.5),lty=2,lwd=median.lsize)
    temp2 <- data.frame(x=max(median_vals,na.rm=T),y=0.5)
    p <- p + geom_segment(data=temp2,aes(x=0,xend=x,y=0.5,yend=0.5),lty=2,lwd=median.lsize)
  }
  
  if(set.time.line){
    temp <- data.frame(x=set.time,y=set.surv)
    p <- p + geom_segment(data=temp,aes(x=0,xend=x,y=set.surv,yend=set.surv),lty=2,lwd=set.lsize)
    temp2 <- data.frame(x=set.time,y=max(set.surv,na.rm=T))
    p <- p + geom_segment(data=temp2,aes(x=x,xend=x,y=0,yend=y),lty=2,lwd=set.lsize)
  }
  
  
  # Colour, linetype and fill -----------------------------------------------
  
  
  ##Getting the colour labels
  if(multiple_lines==T & (length(plot.event)==1|event=="linetype")){
    col_labs <- stratalabs
  }else col_labs <- eventlabs
  
  p <- p + scale_colour_manual(values = col,labels=col_labs)
  p <- p + scale_fill_manual(values = col,labels=col_labs)
  
  ##Getting the linetype labels
  if(multiple_lines==T & (length(plot.event)==1|event=="col")){
    linetype_labs <- stratalabs
  }else linetype_labs <- eventlabs
  
  
  if(!is.null(linetype)) p <- p + scale_linetype_manual(labels=linetype_labs,values = linetype)
  if(is.null(linetype)) p <- p + scale_linetype_discrete(labels=linetype_labs)
  
  
  if(event == "linetype" & length(plot.event)==2 & multiple_lines==F){
    
    p <- p + guides(col = F)
  }
  
  if(event == "col" & length(plot.event)==2 & multiple_lines==F){
    
    p <- p + guides(linetype = F)
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
    
    if(!is.null(col)){cols1 <- col
    }else{cols1 <- .extract_ggplot_colors(p, grp.levels = stratalabs)}
    
    
    
    if(multiple_lines==F|(event == "col" & length(plot.event)>1)) cols1 = "black"
    
    
    n_strata <- length(stratalabs)
    #yticklabs <- rep("-", n_strata)
    
    
    yticklabs <- unname(rev(stratalabs.table))
    #strataname.table = ifelse(n_strata==1, "\n", strataname)
    
    data.table <- ggplot(risk.data, aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
      geom_text(hjust="middle", vjust="center", size = nsize) +
      theme_bw() +
      scale_x_continuous(Numbers_at_risk_text, breaks = times, limits = c(0,maxxval)) +
      coord_cartesian(xlim=c(0,maxxlim)) + ### changes the actual plotted limits if needed
      theme(legend.position = "none") +
      theme(text = element_text(size = fsize)) + ##increase font size
      theme(plot.margin = unit(c(-1.5, 1, 0.1, 0.2), "lines"))+
      scale_y_discrete(strataname.table, breaks = as.character(levels(risk.data$strata)), labels = yticklabs)
    
    data.table <- data.table +suppressWarnings(theme(axis.title.x = element_text(size = fsize, vjust = 1), panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(), panel.border = element_blank(),
                                                     axis.text.x = element_blank(), axis.ticks = element_blank(),
                                                     axis.text.y = element_text(face = "bold", hjust = 1,colour = rev(cols1))))
    
    
    
    
    if(all(as.character(yticklabs)=="-")) data.table <- .set_large_dash_as_ytext(data.table)
    
    
    
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
      # a <- arrangeGrob(p, blank.pic, data.table,
      #                  clip = FALSE, nrow = 3, ncol = 1,
      #                  heights = unit(c(2, .1, .25), c("null", "null", "null")))
      
      if(table==T){
        a <- list(p,blank.pic,data.table)
      } else a <- p
      
      return(a)
    }
  } else(return(p))
  
  
}

#' Plot KM and CIF curves with ggplot
#'
#'This function put together a survival curve, and a number at risk table
#'
#'
#'
#' @param list_gg list containing the results of ggkmcif
#'@export
modify_ggkmcif <- function(list_gg){
  gA <- ggplotGrob(list_gg[[1]])
  gB <- ggplotGrob(list_gg[[2]])
  gC <- ggplotGrob(list_gg[[3]])
  maxWidth = grid::unit.pmax(gA$widths[2:5], gC$widths[2:5])
  gA$widths[2:5] <- as.list(maxWidth)
  gB$widths[2:5] <- as.list(maxWidth)
  gC$widths[2:5] <- as.list(maxWidth)
  
  gridExtra::grid.arrange(gA, gB, gC,
                          clip = FALSE, nrow = 3, ncol = 1,
                          heights = unit(c(2, .1, .25), c("null", "null", "null")))
}