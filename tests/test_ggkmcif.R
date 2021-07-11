library(survival)
data('lung')
lung$status <- lung$status-1 
lung$sex <- ifelse(lung$sex==1,"F","M")


lung$status2 <- lung$status
lung$status2[1:160] <- 0

lung$status3 <- lung$status
lung$status2[1:160] <- NA
# KM ----------------------------------------------------------------------

###########################NO COVARIATE #########################
survfit(Surv(time,status)~1,data=lung)
summary(survfit(Surv(time,status)~1,data=lung),times = c(100,500))

#Confidence intervals
ggkmcif(response = c('time','status'),data=lung,median.text = T,set.time.text = 'Y OS',
        set.time = c(500,100,1100),conf.type = 'log-log',conf.curves=TRUE)
plot(survfit(Surv(time,status)~1,data=lung),mark.time=T,conf.int = T,conf.type = 'log-log')



#MEDIAN and set time
ggkmcif(response = c('time','status'),data=lung,median.text = T,set.time.text = 'Y OS',
        set.time = c(500,100,1100),
        xlab='xlab',ylab='ylab',main='main',
        median.size=6,median.CI = T,median.pos=c(500,0.75),
        median.lines = T,median.lsize = 0.75,
        set.size = 3,set.pos = c(500,0.9),set.time.line = T,
        set.lsize=4)

##aethetics1 
ggkmcif(response = c('time','status'),data=lung,set.time = c(500,100,1100),
        times = c(100,150,200,700,1000),
        xlab='xlab',ylab='ylab',main='main',censor.marks = T,censor.size = 5,
        censor.stroke = 2)


##aethetics2 
ggkmcif(response = c('time','status'),data=lung,
        xlab='xlab',ylab='ylab',main='main',censor.marks = F,
        fsize=5,nsize=10,lsize=4,ylim=c(0,1.5))

##aethetics3
ggkmcif(response = c('time','status'),data=lung,
        xlab='xlab',ylab='ylab',main='main',censor.marks = F,col='#2ca25f',linetype = 2,
        xlim=c(0,1200),Numbers_at_risk_text = "Number")

##Testing missing info
lung$status2 <- lung$status
lung$status2[1:50] <- NA

lung$time2 <- lung$time
lung$time2[160:190] <- NA

summary(survfit(Surv(time2,status)~1,data=lung),times = c(100,500))
ggkmcif(response = c('time2','status'),data=lung,median.text = T,set.time.text = 'Y OS',
        set.time = c(500),conf.type = 'log',conf.curves=TRUE,median.CI = TRUE,set.time.CI = T)

summary(survfit(Surv(time,status2)~1,data=lung),times = c(100,500))
ggkmcif(response = c('time','status2'),data=lung,median.text = T,set.time.text = 'Y OS',
        set.time = c(500),conf.type = 'log',conf.curves=TRUE,median.CI = TRUE,set.time.CI = T,
        times = c(0,80,100,160,240,500))

summary(survfit(Surv(time2,status2)~1,data=lung),times = c(100,500))
ggkmcif(response = c('time2','status2'),data=lung,median.text = T,set.time.text = 'Y OS',
        set.time = c(500),conf.type = 'log',conf.curves=TRUE,median.CI = TRUE,set.time.CI = T)

###########################COVARIATE #########################
lung$sex2 <- lung$sex
lung$sex2[lung$status==0 & lung$age > 50] <- "Missing"
lung$sex2[1:3] <- "Missing"
lung$sex2 <- factor(lung$sex2,levels = c('M',"F","Missing"))

survfit(Surv(time,status)~sex2,data=lung)
summary(survfit(Surv(time,status)~sex2,data=lung),times = c(100,500))

#Confidence intervals
ggkmcif(response = c('time','status'),cov='sex2',data=lung,median.text = T,set.time.text = 'Y OS',
        set.time = c(500,100,1100),conf.type = 'plain',conf.curves=TRUE)
plot(survfit(Surv(time,status)~sex2,data=lung),mark.time=T,conf.int = T,conf.type='plain')


#MEDIAN and set time

ggkmcif(response = c('time','status'),cov='sex2',data=lung,median.text = T,median.digits=1,
        median.CI=TRUE,median.lines = T,set.time.text = 'year OS',set.time=c(500,200,900))

ggkmcif(response = c('time','status'),cov='sex2',data=lung,median.text = T,median.digits=1,
        median.CI=TRUE,median.lines = T,set.time.text = 'year OS',set.time=c(500,200,900),
        set.time.line=T,set.lsize=3,median.size=.5)


summary(survfit(Surv(time,status)~sex2,data=lung),times = c(200,500,900))
ggkmcif(response = c('time','status'),cov='sex2',data=lung,median.text = F,median.digits=1,
        median.CI=TRUE,median.lines = F,set.time.text = 'year OS',set.time=c(500,200,900),
        set.time.CI = T,set.time.digits = 2)

#log rank p-values
ggkmcif(response = c('time','status'),cov='sex',data=lung,median.text = T,median.digits=1,
        median.CI=TRUE,median.lines = T,psize = 3,pval.pos = c(100,0.2))

#continuous
ggkmcif(response = c('time','status'),cov='ph.karno',data=lung,median.text = T,median.digits=1,
        median.CI=TRUE,median.lines = T,cut=60,pval.digits = 5)

lung$temp <- lung$ph.karno > 60
table(lung$temp)

survdiff(Surv(time, status) ~ temp, data = lung)


ggkmcif(response = c('time','status'),cov='sex2',data=lung,median.text = T,median.digits=1,
        median.CI=TRUE,median.lines = T,cut=60,pval.digits = 3,pval.pos = c(600,0.5),psize=10)

lung$temp <- lung$ph.karno > 60
table(lung$temp)

survdiff(Surv(time, status) ~ sex2, data = lung)


#HR
ggkmcif(response = c('time','status'),cov='sex2',HR=T,data=lung,HR_pval = TRUE,HR.pval.digits = 4,
        HR.digits = 4)
summary(coxph(Surv(time, status) ~ sex2, data = lung))

#aestetics
summary(survfit(Surv(time,status)~sex2,data=lung),times = c(250))
ggkmcif(response = c('time','status'),cov='sex2',HR=T,data=lung,table=TRUE,
        times=c(0,250,500,750,1000),xlab='X',ylab="Y",main="MAIN",
        stratalabs=c("Male","Female","Prefer not to answer"),
        strataname = "Sex",stratalabs.table=c("M","F","POA"))

ggkmcif(response = c('time','status'),cov='sex2',HR=T,data=lung,table=TRUE,
        times=c(0,250,500,750,1000),xlab='X',ylab="Y",main="MAIN",
        stratalabs=c("Male","Female","Prefer not to answer"),
        strataname = "Sex",stratalabs.table=c("-","-","-"),strataname.table="strataname.table",
        censor.marks=FALSE)


ggkmcif(response = c('time','status'),cov='sex2',data=lung,censor.size = 5,
        censor.stroke = 2,censor.marks=TRUE,fsize=3,nsize=5,lsize=0.5,linetype=c(1))


ggkmcif(response = c('time','status'),cov='sex2',data=lung,linetype=c(1),col=c(1,2,3))

ggkmcif(response = c('time','status'),cov='sex2',data=lung,linetype=c(1,2,3),col=c(1),
        ylim=c(0,2),xlim=c(0,2000),legend.pos="top",Numbers_at_risk_text="TEST")


##MISSING
lung$sex3 <- lung$sex
lung$sex3[1:50] <- NA

summary(survfit(Surv(time,status)~sex3,data=lung),times = c(80,160))

#Confidence intervals
ggkmcif(response = c('time','status'),cov='sex2',data=lung,median.text = T,set.time.text = 'Y OS',
        set.time = c(80),conf.type = 'log',conf.curves=TRUE,)





# CIF ---------------------------------------------------------------------

lung$crr_status <- lung$status*1
lung$crr_status[seq(1,200,by=2)] <- 2


#Confidence intervals
ggkmcif(response = c('time','crr_status'),data=lung,conf.type = 'log',conf.curves = T)
sum(lung$time>500)

#median and set time

#Median was removed 
ggkmcif(response = c('time','crr_status'),data=lung,conf.curves = T,
        median.text = T,set.time=c(700,300,1200),
        set.time.line = T,set.time.text = "day rate",set.pos = c(500,0.75),
        plot.event=c(1,2),set.lsize=10)




fit <- cuminc(lung$time,lung$crr_status)
cmprsk::timepoints(fit,times=c(300,700,1200))

#Aestetics

ggkmcif(response = c('time','crr_status'),data=lung,conf.curves = T,
        plot.event=c(1),table=FALSE)

ggkmcif(response = c('time','crr_status'),data=lung,conf.curves = T,
        plot.event=c(1),table=TRUE,times=c(0,400,1000),xlab='x',ylab='y',
        main='main',Numbers_at_risk_text='text',fsize=20,nsize = 6,lsize=2,ylim=c(0,2),flip.CIF = T)


ggkmcif(response = c('time','crr_status'),data=lung,col=c(1),plot.event=c(1,2),event='linetype',xlim=c(0,500),
        eventlabs=c("rec","death"),event.name='type',legend.pos = 'bottom')

###########################COVARIATE #########################
lung$sex2 <- lung$sex
lung$sex2[lung$status==0 & lung$age > 50] <- "Missing"
lung$sex2[1:3] <- "Missing"
lung$sex2 <- factor(lung$sex2,levels = c('M',"F","Missing"))

#Confidence intervals
ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,conf.curves = T)

##Set times 
ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,set.time =c(1000,1100),
        set.time.text='day rate',set.time.line = T)

fit <- cuminc(lung$time,lung$crr_status,lung$sex2)
cmprsk::timepoints(fit,times=c(1100,1000))


ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,set.time =c(700,400),
        set.time.text='day rate',set.time.line = T,set.lsize = .25)

fit <- cuminc(lung$time,lung$crr_status,lung$sex2)

cmprsk::timepoints(fit,times=c(600))
ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,set.time =c(600),
        set.time.text='day rate',set.time.line = T,plot.event = c(1),set.time.digits=5,set.time.CI = T)

ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,set.time =c(600),
        set.time.text='day rate',set.time.line = T,plot.event = c(2),set.time.digits=5,set.time.CI = T)


#pvalue

cuminc(lung$time,lung$crr_status,lung$sex)$Tests
ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,set.time =c(600),
        set.time.text='day rate',set.time.line = T)

ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,set.time =c(600),
        set.time.text='day rate',set.time.line = T,plot.event = 2)

ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,
        set.time.text='day rate',set.time.line = T,plot.event = c(1,2),
        eventlabs = c('death by A','death by other'),legend.pos = 'bottom',
        event.name = 'Type of death',psize=4,
        pval.pos=c(200,0.75))

ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,
        set.time.text='day rate',set.time.line = T,plot.event = c(1,2),
        eventlabs = c('death by A','death by other'),legend.pos = 'bottom',
        event.name = 'Type of death',psize=3,
        pval.pos=c(200,0.75))

##HR
ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,HR_pval=T,plot.event = c(1))
ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,HR=T,HR.pval=T,plot.event = c(2))
ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,HR=T,HR.pval=T,plot.event = c(1,2))


ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,HR=T,HR.pval=T,plot.event = c(2))
ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,HR=T,HR.pval=T,plot.event = c(1,2))

uvsum(response=c('time','crr_status'),cov='sex',data=lung)

ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,HR_pval=T,HR=T,plot.event = c(1),HR.digits = 1,
        HR.pval.digits = 3,table=F)
uvsum(response=c('time','crr_status'),cov='sex2',data=lung)

ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,HR_pval=T,HR=T,plot.event = c(1),HR.digits = 1,
        HR.pval.digits = 3)


ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,HR_pval=T,HR=T,plot.event = c(1),HR.digits = 1,
        HR.pval.digits = 3,table=F,flip.CIF = T)

ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,HR_pval=T,HR=T,plot.event = c(1),HR.digits = 1,
        HR.pval.digits = 3,table=F,flip.CIF = T,conf.curves = T)

ggkmcif(response = c('time','crr_status'),cov='sex2',data=lung,HR_pval=T,HR=T,plot.event = c(1),HR.digits = 1,
        HR.pval.digits = 3,times=c(100,1000),xlab="x",ylab="y",main = "main",fsize=10,nsize=7,lsize = 0.4,ylim=c(0,1.23))

ggkmcif(response = c('time','crr_status'),cov='sex',data=lung,plot.event=c(1,2),
        stratalabs = c('Female',"Male"),strataname='SEX',stratalabs.table=c("F","M"),
        strataname.table="Gender",event='col',event.name = c("Cause of death"),eventlabs =c("Cancer","Other"))


# -------------------------------------------------------------------------



