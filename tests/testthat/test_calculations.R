#-------------------------------------------------------------------------------------
# Data Definition
library(survival)
data(lung)
lung$Status=factor(lung$status-1)
lung$Sex = factor(lung$sex,labels = c('Male','Female'))
lung$AgeGroup = factor(cut(lung$age, breaks=seq(0,100,10)), labels = c('30s','40s','50s','60s','70s','80s'))
lung$OneLevelFactor = factor(x='one level')
lung = lung[order(lung$Status),]

lung$x_null = rnorm(nrow(lung))
lung$x_pred = c(rnorm(sum(lung$Status==0),0,1),
                rnorm(sum(lung$Status==1),1,1))

test_data = data.frame(
  y= rnorm(100),
  x0= geoR::rboxcox(100, lambda=0.5, mean=10, sd=2))
test_data$x1=  test_data$x0*test_data$y
#-------------------------------------------------------------------------------------


test_that("covsum calculates correctly with no maincov", {
  output = covsum(data=lung,
                  covs=c('Status','Sex','wt.loss','OneLevelFactor'),
                  markup=F)
  expect_equal(names(output) , c("Covariate",'n=228'))
  expect_equal(output$Covariate , c("Status","0","1","Sex","Male","Female","wt loss","Mean (sd)","Median (Min,Max)","Missing","OneLevelFactor","one level"))
  expect_equal(output$`n=228`,c("","63 (28)","165 (72)","","138 (61)","90 (39)","","9.8 (13.1)","7 (-24,68)","14","","228 (100)"))

})

test_that("covsum calculates correctly with maincov", {
  output = covsum(data=lung,
                  maincov='Status',
                  covs=c('Sex','wt.loss','OneLevelFactor'),
                  markup=F)
  expect_equal(names(output) ,c("Covariate","Full Sample (n=228)","0 (n=63)","1 (n=165)","p-value") )
  expect_equal(output$Covariate , c("Sex","Male","Female","wt loss","Mean (sd)","Median (Min,Max)","Missing","OneLevelFactor","one level"))
  expect_equal(output[,3],c("","26 (41)","37 (59)","","9.1 (12.9)","4 (-10,49)","1","","63 (100)"))

})


test_that("uvsum logistic regression CIS are correct",{
  digits = 2 # TODO: add function flexibility
  for (ci_width in c(0.9,.95,.99,.995)){
    m1 = glm(Status~x_null,data=lung,family='binomial')
    x1=summary(m1)$coefficients
    m2 = glm(Status~x_pred,data=lung,family='binomial')
    x2=summary(m2)$coefficients
    # TODO: update the wald CIs to likelihood profile CIs?
    # expected = c(  paste(format(round(exp(x1[2,1]),digits),nsmall=digits),
    #                      paste0("(",paste0(format(round(exp(confint(m1,level=ci_width)[2,]),digits),nsmall=digits),collapse=","),")")),
    #                paste(format(round(exp(x2[2,1]),digits),nsmall=digits),
    #                       paste0("(",paste0(format(round(exp(confint(m2,level=ci_width)[2,]),digits),nsmall=digits),collapse=","),")")))
    expected = c(  paste(format(round(exp(x1[2,1]),digits),nsmall=digits),
                         paste0("(",paste0(format(round(c(exp(x1[2,1]-qnorm(1-(1-ci_width)/2)*x1[2,2]),exp(x1[2,1]+qnorm(1-(1-ci_width)/2)*x1[2,2])),digits),nsmall=digits),collapse=","),")")),
                   paste(format(round(exp(x2[2,1]),digits),nsmall=digits),
                         paste0("(",paste0(format(round(c(exp(x2[2,1]-qnorm(1-(1-ci_width)/2)*x2[2,2]),exp(x2[2,1]+qnorm(1-(1-ci_width)/2)*x2[2,2])),digits),nsmall=digits),collapse=","),")")))
    output = uvsum(response = 'Status',
                   covs=c('x_null','x_pred'),
                   data=lung,
                   type='logistic',
                   CIwidth = ci_width,
                   markup=F)

    expect_equal(output[,2],expected)
    expect_equal(names(output)[2],paste0("OR(",ci_width*100,"\\%CI)"))

  }

})

test_that("uvsum linear regression CIS are correct",{
  digits = 2 # TODO: add function flexibility
  for (ci_width in c(0.9,.95,.99,.995)){
    m1 = lm(age~wt.loss,data=lung)
    x1=summary(m1)$coefficients
    m2 = lm(age~Sex,data=lung)
    x2=summary(m2)$coefficients
    expected = c(  paste(format(round(x1[2,1],digits),nsmall=digits),
                         paste0("(",paste0(format(round(c(x1[2,1]-qt(1-(1-ci_width)/2,m1$df.residual)*x1[2,2],x1[2,1]+qt(1-(1-ci_width)/2,m1$df.residual)*x1[2,2]),digits),nsmall=digits),collapse=","),")")),
                   paste(format(round(x2[2,1],digits),nsmall=digits),
                         paste0("(",paste0(format(round(c(x2[2,1]-qt(1-(1-ci_width)/2,m2$df.residual)*x2[2,2],x2[2,1]+qt(1-(1-ci_width)/2,m2$df.residual)*x2[2,2]),digits),nsmall=digits),collapse=","),")")))
    expected = gsub(', ',',',expected)
    output = uvsum(response = 'age',
                   covs=c('wt.loss','Sex'),
                   data=lung,
                   type='linear',
                   CIwidth = ci_width,
                   markup=F)

    expect_equal(output[c(1,4),2],expected)
    expect_equal(names(output)[2],paste0("Estimate(",ci_width*100,"\\%CI)"))
  }

})

# # Uncomment this if you need to ensure the tests are being run
# test_that("this script is being executed",{
#   expect_equal("This script was run","YES!")
# })
# TODO: Finishing writing CI test scripts for uvsum (boxcos, coxph, crr)
# Write script to check sample size calculations

# TODO: test_that('CI works correctly ')

