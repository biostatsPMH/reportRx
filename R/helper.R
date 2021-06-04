#' Retrieve columns number from Excel columns specified as unquoted letters
#' @param ... unquoted excel column headers (i.e. excelCol(A,CG,AA)) separated by commas
#' @importFrom rlang as_string
#' @export
excelCol<- function(...){
  args <- as.list(match.call())[-1]
  args <-unname(unlist(lapply(args,function(x) {rlang::as_string(x)})))
  rtn<-sapply(args, function(x){
    colHead <- toupper(x)
    if (nchar(colHead)>1){
      l1 = substr(colHead,1,1)
      l2 = substr(colHead,2,2)
      rtn <- 26*which(LETTERS==l1)+which(LETTERS==l2)
    } else {
      rtn <- which(LETTERS==colHead)
    }
  })
  names(rtn) <- toupper(names(rtn))
  return(rtn)
}

niceNum <- function(x,digits=2){
  rndx = sapply(x, function(x) {format(round(as.numeric(x),digits),nsmall=digits)})
  return(gsub(" ","",rndx))
}

#' Return a data frame of codes
#' 
#' Accepts a string input in the form "code1=label1,code2=label2,.." and 
#' returns a data frame with a column of codes and a column of labels
#' 
#' @param labelStr in the format code1=label1,code2=label2
#' @param delim delimeter separating codes in labelStr, defaults to ','
#' @param codeName column name for codes, defaults to code
#' @param lblName column name for labels, defaults to label
#' @export
importCodes<-function(labelStr,delim=',',codeName='code',lblName='label'){
  x=strsplit(labelStr,split=delim)[[1]]
  codeLst=strsplit(x,"=")
  tbl <- NULL
  for (i in seq_along(codeLst)) tbl<-rbind(tbl,cbind(code=codeLst[[i]][1],label=codeLst[[i]][2]))
  tbl <- as.data.frame(tbl)
  if (isTRUE(all.equal(as.character(suppressWarnings(as.numeric(tbl[[1]]))),tbl[[1]]))) tbl[[1]] <- as.numeric(tbl[[1]])
  names(tbl)=c(codeName,lblName)
  return(tbl)
}

#' Paste with parentheses
#'
#' Paste with parentheses
#'
#' @param x a vector
#'@keywords helper
#'@export
#'@examples
#'pstprn(c(1,2,3,4,5))
#'pstprn(c("Hello","Hi",2))
pstprn<-function(x){paste(x[1]," (",paste(x[-1],collapse=","),")",sep="")}

#' Round and paste with parentheses
#'
#' Round and paste with parentheses
#'
#' @param x a numeic vector
#' @param y integer corresponding to the number of digits to round by
#'@keywords helper
#'@export
#'@examples
#'psthr(c(1.111111,2.2222222,3.333333))
# LA updated to always return a formatted string
psthr<- function (x, y = 2)
{
  x <- sapply(x, function(x) {
    ifelse(abs(x) < 0.01 | abs(x) > 1000, format(x, scientific = TRUE,digits = y), format(round(x, y),nsmall = y))
  })
  pstprn(x)
}
covnm<-function(betanames,call){
  sapply(betanames,function(betaname){
    
    # indx<-which(sapply(call,function(cov){charmatch(cov,betaname)})==1)
    indx=which(sapply(call,function(cov)grepl(cov,betaname,fixed=TRUE))) ## changed on Feb 21, 2019 for checkings
    if(length(indx)==1) return(call[indx])
    #If one  facorname is a subset of another
    indx2<-which.max(sapply(call[indx],nchar))
    if(length(indx2)==1) return(call[indx[indx2]])
    indx3<-which(sapply(call[indx2],function(c){substr(betaname,1,nchar(c))==c}))
    if(length(indx3)==1)  return(call[indx[indx2[indx3]]])
  })
}

alleql<-function(x,y){
  !any((x==y)==F)
}


betaindx<-function(x){
  i=1
  out<-1
  result<-NULL
  while(TRUE){
    if(i+1>length(x)){
      result<-c(result,list(out))
      return(result)
    }
    else if(alleql(x[[i+1]],x[[i]])){
      out<-c(out,i+1)
    }
    else{
      result<-c(result,list(out))
      out<-i+1
    }
    i=i+1
  }
}


#' Capitalize a string
#'
#' Calitalize a string
#'
#' @param x string
#' @keywords helper
#' @export
cap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

#'Clean strings for printing
#'
#' Returns strings with . and _ replaced by a space. This is nice when printing column names of your dataframe in a report
#' @param strings vector of strings to give a nice name
#' @keywords helper
#' @export
nicename<-function(strings){
  out<-sapply(strings,function(x){
    x<-chartr(".", " ",x)
    x<-chartr("_", " ",x)
    return(x)})
  return(out)
}

#' Formats p-values
#'
#' Returns <0.001 if pvalue is <0.001. Else rounds the pvalue to 2 significant digits
#'
#' @param x an integer
#' @export
pvalue<-function(x){
  if(is.na(x)|class(x)=="character") return(x)
  else if (x<=0.001) return("<0.001")
  else return(signif(x,2))
}

sanitize <- function(str) {
  result <- str
  result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
  result <- gsub("$", "\\$", result, fixed = TRUE)
  result <- gsub(">", "$>$", result, fixed = TRUE)
  result <- gsub("<", "$<$", result, fixed = TRUE)
  result <- gsub("|", "$|$", result, fixed = TRUE)
  result <- gsub("{", "\\{", result, fixed = TRUE)
  result <- gsub("}", "\\}", result, fixed = TRUE)
  result <- gsub("%", "\\%", result, fixed = TRUE)
  result <- gsub("&", "\\&", result, fixed = TRUE)
  result <- gsub("_", "\\_", result, fixed = TRUE)
  result <- gsub("#", "\\#", result, fixed = TRUE)
  result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
  result <- gsub("~", "\\~{}", result, fixed = TRUE)
  result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$",
                 result, fixed = TRUE)
  return(result)
}

#' Sanitizes strings to not break LaTeX
#'
#' Strings with special charaters will break LaTeX if returned 'asis' by knitr. This happens every time we use one of the main reportRx functions. We first sanitize our strings with this function to stop LaTeX from breaking.
#'
#'@param str a vector of strings to sanitize
#'@export
sanitizestr<-function(str){
  as.vector(sapply(str,function(char){sanitize(char)}))
}

#'Bold strings in LaTeX
#'
#'Bold strings in LaTeX.
#'
#'@param strings A vector of strings to bold.
#'@export
lbld<-function(strings){sapply(strings,function(x){
  if(is.null(x)) return(x)
  if(is.na(x)) return(x)
  return(paste("\\textbf{",x,"}",sep=""))})}

#'Add spaces to strings in LaTeX
#'
#'Add spaces to strings in LaTeX. Returns appends ~~~ before the string
#'
#'@param x string
#'@export
addspace<-function(x){
  paste("~~~",x,sep="")
}
#' Formats p-values for LaTeX
#'
#' Returns <0.001 if pvalue is <0.001. Else rounds the pvalue to 2 significant digits. Will bold the p-value if it is <= 0.05
#' @param x an integer
#' @export
lpvalue<-function(x){
  if(is.na(x)|class(x)=="character") return(x)
  else if (x<=0.001) return("\\textbf{$<$0.001}")
  else x=signif(x,2)
  if(x<=0.05) return(paste("\\textbf{",x,"}",sep=""))
  else return(x)
}



removedollar<-function(x){
  colnms<-strsplit(x,":")
  indx<-unlist(lapply(colnms,function(colnm) sapply(colnm, function(coln) regexpr("$",coln,fixed=T)[1]+1)))
  if(length(unique(indx))==1){
    if(unique(indx)!=0) x<-unlist(lapply(colnms,function(colnm) paste(substring(colnm,indx[1]),collapse=":")))
  }
  return(x)
}

modelmatrix<-function(f,data=NULL){
  k<-as.character(f)
  y<-NULL
  if(!length(k)%in%c(2,3)) stop("formula not properly formed")
  if(length(k)==3) {
    f<-as.formula(paste("~",k[2],"+",k[3],sep=""))
    y<-model.matrix(as.formula(paste("~",k[2],sep="")),data)[,-1,drop=F]}
  x<-model.matrix(f,data)[,-1,drop=F]
  colnames(x)<-removedollar(colnames(x))
  if(!is.null(y)){
    return(list(x[,1:ncol(y),drop=F],x[,(ncol(y)+1):ncol(x),drop=F]))
  }else{
    return(x)
  }}

matchcovariate=function(betanames,ucall){
  out=as.vector(sapply(betanames,function(betaname){
    splitbetaname=unlist(strsplit(betaname,":",fixed=T))
    out=sapply(splitbetaname,function(bname){
      bname=gsub(" ","",bname) # added 14 Dec 2020 to allow matching with centred variables
      #indx=which(sapply(ucall,function(cov)charmatch(cov,bname))==1)
      indx=which(sapply(ucall,function(cov)grepl(cov,bname,fixed=TRUE))) ## changed on Feb 21, 2019 for checkings
      if(length(indx)==1)return(indx)
      #If one  facorname is a subset of another
      indx2<-which.max(sapply(ucall[indx],nchar))
      if(length(indx2)==1) return(indx[indx2])
      indx3<-which(sapply(ucall[indx2],function(c){substr(betaname,1,nchar(c))==c}))
      if(length(indx3)==1)  return(ucall[indx[indx2[indx3]]])
      return(-1)
    })
    if(-1 %in% out) return(-1)
    result=0
    n=length(out)
    for(i in 1:length(out)){
      result=result+out[i]*100^(n-1)
      n=n-1
    }
    return(result)}))
  if(-1 %in% out) return(-1)
  return (out)
}

# (ggsurv) ---------------------------------------------------------

round_sprintf <- function(value,digits){
  sprintf( paste0("%.",digits,"f"), round(value,digits))
}

pstprn0 <- function (x)
{
  paste0(x[1], "(", paste0(x[-1], collapse = ","),
         ")", sep = "")
}

psthr0 <- function (x, digits = 2)
{
  x <- sapply(x, function(x) {
    ifelse(abs(x) < 0.01 | abs(x) > 1000, format(x, scientific = TRUE,
                                                 digits = digits), round_sprintf(x, digits))
  })
  pstprn0(x)
}

break_function <- function(xmax){
  
  xmax_length <- ifelse(xmax>1,nchar(round(xmax)),round(abs(log10(xmax))))
  
  byx <- if(xmax>1) {round(xmax/10,digits = 2-xmax_length)
  }else round(xmax/10,digits = xmax_length+1)
  
  breaks <- seq(0,xmax,by=byx)
  if(max(breaks)<byx) breaks <- c(breaks,max(breaks)+byx)
  return(breaks)
}

lpvalue2 <- function (x,digits)
{
  if (is.na(x) | class(x) == "character")
    return(x)
  else if (x < 10^-(digits))
    return(paste0("p < ",10^-(digits)))
  else return(paste0("p = ",round_sprintf(x, digits)))
  
}

.extract_ggplot_colors <- function(p, grp.levels){
  g <- ggplot_build(p)
  .cols <- unlist(unique(g$data[[1]]["colour"]))
  if(!is.null(grp.levels)){
    if(length(.cols)==1) .cols <- rep(.cols, length(grp.levels))
    names(.cols) <- grp.levels
  }
  .cols
}

.set_large_dash_as_ytext <- function(ggp){
  ggp + theme(axis.text.y = element_text(size = 50, vjust = 0.35),
              axis.ticks.y = element_blank())
}

##This function is used by the survfit package
survfit_confint <- function(p, se, logse=TRUE, conf.type, conf.int=0.95,
                            selow, ulimit=TRUE) {
  zval <- qnorm(1- (1-conf.int)/2, 0,1)
  if (missing(selow)) scale <- 1.0
  else scale <- ifelse(selow==0, 1.0, selow/se)  # avoid 0/0 at the origin
  if (!logse) se <- ifelse(se==0, 0, se/p)   # se of log(survival) = log(p)
  
  if (conf.type=='plain') {
    se2 <- se* p * zval  # matches equation 4.3.1 in Klein & Moeschberger
    if (ulimit) list(lower= pmax(p -se2*scale, 0), upper = pmin(p + se2, 1))
    else  list(lower= pmax(p -se2*scale, 0), upper = p + se2)
  }
  else if (conf.type=='log') {
    #avoid some "log(0)" messages
    xx <- ifelse(p==0, NA, p)
    se2 <- zval* se
    temp1 <- exp(log(xx) - se2*scale)
    temp2 <- exp(log(xx) + se2)
    if (ulimit) list(lower= temp1, upper= pmin(temp2, 1))
    else  list(lower= temp1, upper= temp2)
  }
  else if (conf.type=='log-log') {
    xx <- ifelse(p==0 | p==1, NA, p)
    se2 <- zval * se/log(xx)
    temp1 <- exp(-exp(log(-log(xx)) - se2*scale))
    temp2 <- exp(-exp(log(-log(xx)) + se2))
    list(lower = temp1 , upper = temp2)
  }
  else if (conf.type=='logit') {
    xx <- ifelse(p==0, NA, p)  # avoid log(0) messages
    se2 <- zval * se *(1 + xx/(1-xx))
    
    temp1 <- 1- 1/(1+exp(log(p/(1-p)) - se2*scale))
    temp2 <- 1- 1/(1+exp(log(p/(1-p)) + se2))
    list(lower = temp1, upper=temp2)
  }
  else if (conf.type=="arcsin") {
    xx <- ifelse(p==0, NA, p)
    se2 <- .5 *zval*se * sqrt(xx/(1-xx))
    list(lower= (sin(pmax(0, asin(sqrt(xx)) - se2*scale)))^2,
         upper= (sin(pmin(pi/2, asin(sqrt(xx)) + se2)))^2)
  }
  else stop("invalid conf.int type")
}


color_palette_surv_ggplot <- function(length){
  if(length==1) return("black")
  if(length==2) return(c("#D53E4F","#3288BD"))
  if(length==3) return(c("#D53E4F","#ABDDA4","#3288BD"))
  if(length==4) return(c("#D53E4F","#FDAE61","#ABDDA4","#3288BD"))
  if(length==5) return(c("#D53E4F","#FDAE61","#FEE08B","#ABDDA4","#3288BD"))
  if(length==6) return(c("#D53E4F","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD"))
  if(length==7) return(c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD"))
  if(length==8) return(c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))
  if(length==9) return(c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))
  if(length==10) return(c("black","#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))
  if(length>10) {message("10 colours maximum in default")}
  return(rep(c("black","#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"),length.out=length))
}


# (forestplot2) ---------------------------------------------------------
format_glm = function(glm_fit,conf.level = 0.95,digits=c(2,3),orderByRisk=TRUE){
  
  if (! class(glm_fit)[1] %in% c('glm','polr')) stop('Only objects of class glm and polr are accepted.')
  
  #extracting ORs and p values
  Z = qnorm(1-(1-conf.level)/2)
  tab <- as.data.frame(summary(glm_fit)$coefficients)
  names(tab) =  c("estimate",  "std.error" ,"statistic", "p.value")
  tab <- cbind(term= rownames(tab),tab)
  rownames(tab) <- NULL
  tab$conf.low=exp(tab$estimate-Z*tab$std.error)
  tab$conf.high=exp(tab$estimate+Z*tab$std.error)
  tab$estimate = exp(tab$estimate)
  tab$estimate.label = paste0(niceNum(tab$estimate), ' (',niceNum(tab$conf.low),', ',niceNum(tab$conf.high),')')
  
  if (class(glm_fit)[1]=='glm'){
    tab = tab[-which(tab$term=='(Intercept)'),]
  }  else {
    tab$coef.type = ifelse(str_detect(rownames(tab),"[|]"),"scale","coefficient")
    tab <- tab[tab$coef.type=='coefficient',]
    tab$p.value = pnorm(abs(tab$`t value`),lower.tail = FALSE) * 2
  }
  

  tab$p.label = ifelse(tab$p.value<0.001, '<0.001', niceNum(tab$p.value,digits[2]))
  names(tab)[1] = 'variable'
  
  tab = tab[,c('variable', 'estimate', 'p.label', 'p.value', 'conf.low', 'conf.high')]
  
  
  if (orderByRisk){
    tab$var.order = rank(tab$estimate)
  } else{
    tab$var.order = 1:nrow(tab)
  }
  
  # Extract the reference levels if needed
  if (length(glm_fit$xlevels)!=0){
    ref_levels <- NULL
    for (i in seq_along(glm_fit$xlevels)){
      ref_levels <- rbind(ref_levels,
                          data.frame(var.name=rep(names(glm_fit$xlevels)[i],length(glm_fit$xlevels[[i]])+1),
                                     level.name = c(names(glm_fit$xlevels)[i],glm_fit$xlevels[[i]]),
                                     level.order=1:(length(glm_fit$xlevels[[i]])+1),
                                     variable=paste0(names(glm_fit$xlevels)[i],c('',glm_fit$xlevels[[i]]))))
    }
    
    
    tab = merge(ref_levels, tab, by='variable',all = T)
    
    tab$estimate.label = ifelse(is.na(tab$estimate), '1.0 (Reference)',
                                paste0(niceNum(tab$estimate), ' (',niceNum(tab$conf.low),', ',niceNum(tab$conf.high),')'))
    
    varOrders <- tapply(X = tab$var.order,
                        INDEX=tab$var.name,
                        FUN = function(x) min(x,na.rm=T))
    varOrderLookup <- data.frame(var.name=names(varOrders),var.order=varOrders)
    
    
    varOrderLookup <- stats::na.omit(tab[,c("var.name","var.order")])
    
    for (i in 1:nrow(varOrderLookup)){
      tab$var.order[tab$var.name==varOrderLookup$var.name[i]] <- varOrderLookup$var.order[i]
    }
    
    tab$estimate.label = ifelse(tab$level.name %in% names(glm_fit$xlevels),NA_character_,tab$estimate.label)
    tab[order(tab$var.order,tab$level.order,decreasing=c(F,T)),]
  } else {
    tab$estimate.label = paste0(niceNum(tab$estimate), ' (',niceNum(tab$conf.low),', ',niceNum(tab$conf.high),')')
    tab$level.order=1
    tab$var.name=tab$variable
    tab$level.name=tab$variable
    tab[order(tab$var.order),]
  }
  
}

# LA new 2021 ---------------------------------------------------------
# New function to strip centering from a covariate
getvarname = function(betaname){
  sapply(betaname,function(x){
    x = gsub('I[(]','',x)
    x = gsub('[-+].*','',x)
    x = trimws(x)
    return(x)
  })
}

lbl_count <- function(y){
  q75 <- summary(y)[5]
  return(data.frame(y=max(y),  label=paste('n =',length(y))))
}

betaWithCI <-function(betaname,CIwidth=0.95){
  paste0(betaname,"(",100*CIwidth,"%CI)")
}

niceStr <- function (strings)
{
  out <- sapply(strings, function(x) {
    x <- chartr('/',' ',x)
    x <- chartr(".", " ", x)
    x <- chartr("_", " ", x)
    return(x)
  })
  return(out)
}

wrp_lbl <- function(x,width = 10){
  x <- niceStr(x)
  stringr::str_wrap(x,width = width)
}


label_wrap_reportRx <- function (width = 25, multi_line = TRUE) {
  fun <- function(labels) {
    labels <- label_value(labels, multi_line = multi_line)
    lapply(labels, function(x) {
      x <- niceStr(x)
      x <- strwrap(x, width = width, simplify = FALSE)
      vapply(x, paste, character(1), collapse = "\n")
    })
  }
  structure(fun, class = "labeller")
}

formatp<- function(pvalues,sigdigits=2){
  p_out <- sapply(pvalues, function(x){
    xsig <- suppressWarnings(signif(as.numeric(x),sigdigits))
    x <- ifelse(x=='excl','excl',ifelse(is.na(xsig),NA_character_,ifelse(xsig<0.001,"<0.001",format(xsig))))})
  p_out = unname(p_out)
  return(p_out)
}

bib_ReadGatherTidy <- function(file){
  #
  # Note that this is an amalgamation of the bib2df read, gather and tidy functions which are included here because they are not exported from the bib2df namespace
  # changes have been made to enable the code to work outside the bib2df package.
  bib <- readLines(file)
  bib <- stringr::str_replace_all(bib, "[^[:graph:]]", " ")
  
  from <- which( stringr::str_extract(bib, "[:graph:]") == "@")
  to <- c(from[-1] - 1, length(bib))
  if (!length(from)) {
    return(empty)
  }
  itemslist <- mapply(function(x, y) return(bib[x:y]), x = from,
                      y = to - 1, SIMPLIFY = FALSE)
  keys <- lapply(itemslist, function(x) {
    stringr::str_extract(x[1], "(?<=\\{)[^,]+")
  })
  fields <- lapply(itemslist, function(x) {
    stringr::str_extract(x[1], "(?<=@)[^\\{]+")
  })
  fields <- lapply(fields, toupper)
  categories <- lapply(itemslist, function(x) {
    stringr::str_extract(x, "[:graph:]+")
  })
  dupl <- sum(unlist(lapply(categories, function(x) sum(duplicated(x[!is.na(x)])))))
  if (dupl > 0) {
    message("Some BibTeX entries may have been dropped.\n            The result could be malformed.\n            Review the .bib file and make sure every single entry starts\n            with a '@'.")
  }
  values <- lapply(itemslist, function(x) {
    stringr::str_extract(x, "(?<==).*")
  })
  values <- lapply(values, function(x) {
    stringr::str_extract(x, "(?![\"\\{\\s]).*")
  })
  values <- lapply(values, function(x) {
    gsub("?(^[\\{\"])", "", x)
  })
  values <- lapply(values, function(x) {
    gsub("?([\\}\"]\\,$)", "", x)
  })
  values <- lapply(values, function(x) {
    gsub("?([\\}\"]$)", "", x)
  })
  values <- lapply(values, function(x) {
    gsub("?(\\,$)", "", x)
  })
  values <- lapply(values, trimws)
  items <- mapply(cbind, categories, values, SIMPLIFY = FALSE)
  items <- lapply(items, function(x) {
    x <- cbind(toupper(x[, 1]), x[, 2])
  })
  items <- lapply(items, function(x) {
    x[stats::complete.cases(x), ]
  })
  items <- mapply(function(x, y) {
    rbind(x, c("CATEGORY", y))
  }, x = items, y = fields, SIMPLIFY = FALSE)
  items <- lapply(items, t)
  items <- lapply(items, function(x) {
    colnames(x) <- x[1, ]
    x <- x[-1, ]
    return(x)
  })
  items <- lapply(items, function(x) {
    x <- t(x)
    x <- data.frame(x, stringsAsFactors = FALSE)
    return(x)
  })
  empty = data.frame(CATEGORY=character(), BIBTEXKEY=character(), ADDRESS=character(), ANNOTE=character(), AUTHOR=character(), BOOKTITLE=character(), CHAPTER=character(), CROSSREF=character(), EDITION=character(), EDITOR=character(), HOWPUBLISHED=character(), INSTITUTION=character(), JOURNAL=character(), KEY=character(), MONTH=character(), NOTE=character(), NUMBER=character(), ORGANIZATION=character(), PAGES=character(), PUBLISHER=character(), SCHOOL=character(), SERIES=character(), TITLE=character(), TYPE=character(), VOLUME=character(), YEAR=character())
  dat <- dplyr::bind_rows(c(list(empty), items))
  dat <- as.data.frame(dat)
  dat$BIBTEXKEY <- unlist(keys)
  bib <- dat
  if (dim(bib)[1] == 0) {
    return(bib)
  }
  AUTHOR <- EDITOR <- YEAR <- CATEGORY <- NULL
  if ("AUTHOR" %in% colnames(bib)) {
    bib$AUTHOR = strsplit(bib$AUTHOR, " and ",fixed = TRUE)
  }
  if ("EDITOR" %in% colnames(bib)) {
    bib$EDITOR = strsplit(bib$EDITOR, " and ",fixed = TRUE)
  }
  if ("YEAR" %in% colnames(bib)) {
    if (sum(is.na(as.numeric(bib$YEAR))) == 0) {
      bib$YEAR = as.numeric(bib$YEAR)
    } else { warning('Check YEAR in bibfile, may be some missing or character strings.')}
  }
  return(bib)
}



checkOutput <- function(fmt){
  if (missing(fmt)){
    message(paste("Is Latex:", knitr::is_latex_output()))
    message(paste("Is HTML:",    knitr::is_html_output()))
    message(paste("Is HTML excl markdoen epub:",  knitr::is_html_output(excludes = c("markdown", "epub"))))
  }
  # Get current formats
  message(paste('Pandoc To:',knitr::pandoc_to()))
}