################################################################################################
####          Maternal Preconceptional and Prenatal Exposure to El Nino                     ####
####     Southern Oscillation Levels and Child Mortality: A Multi-country Study             ####
####                                                                                        ####
####       Hongbing XU, Castiel Chen ZHUANG*, Vanessa M. ODDO, Espoir Bwenge MALEMBAKA,     ####
####                     Xinghou HE, Qinghong ZHANG, Wei HUANG                              ####
####                                                                                        ####
####                              First Version: October 2022                               ####
####                              First Revision: August 2023                               ####
####                              Current Version: March 2024                               ####
################################################################################################



#### Data Preparation ####

## The dataset used in this study, "data_IPUMS", included the following main variables:
##
## country, dhsid, neonatal_death, infant_death, under_five_death, survival_time_neonatal,
## survival_time_infant, survival_time_under_five, kid_sex, maternal_age, place_of_residence, 
## education, household_wealth, toilet, water_source, birth_order, delivery_location, 
## marital_status, MEI_lag0, MEI_lag1, MEI_lag2, MEI_lag3, MEI_lag4, MEI_lag5, MEI_lag6,
## MEI_lag7, MEI_lag8, MEI_lag9, MEI_lag10, MEI_lag11, MEI_lag12, ESPI_lag0, ESPI_lag1,
## ESPI_lag2, ESPI_lag3, ESPI_lag4, ESPI_lag5, ESPI_lag6, ESPI_lag7, ESPI_lag8, ESPI_lag9,
## ESPI_lag10, ESPI_lag11, ESPI_lag12, ONI_lag0, ONI_lag1, ONI_lag2, ONI_lag3, ONI_lag4,
## ONI_lag5, ONI_lag6, ONI_lag7, ONI_lag8, ONI_lag9, ONI_lag10, ONI_lag11, ONI_lag12,
## Nino12_lag0, Nino12_lag1, Nino12_lag2, Nino12_lag3, Nino12_lag4, Nino12_lag5,
## Nino12_lag6, Nino12_lag7, Nino12_lag8, Nino12_lag9, Nino12_lag10, Nino12_lag11,
## Nino12_lag12, Nino34_lag0, Nino34_lag1, Nino34_lag2, Nino34_lag3, Nino34_lag4,
## Nino34_lag5, Nino34_lag6, Nino34_lag7, Nino34_lag8, Nino34_lag9, Nino34_lag10, 
## Nino34_lag11, Nino34_lag12


# Load Package 
library("survival");library(splines);library(dlnm);library(ggplot2);library(glmnet);library(meta)

# Create data frame with outcomes and covariates
data <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                      "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                      "place_of_residence","education","household_wealth","toilet", 
                      "water_source","birth_order","delivery_location","marital_status")]
colnames(data)



#### Step 1. Selecting covariates for each mortality outcome based on LASSO regression ####

# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# before the LASSO regression using package "neuralnet"

# LASSO for neonatal mortality 
# This sub-section produces the following: 
# (1) Extended Data Figure S8 (1)(A)
# (2) Extended Data Figure S8 (1)(B)

X<-data.frame(data$maternal_age,data$kid_sex,data$place_of_residence,data$education_dum,data$household_wealth_dum,
              data$toilet,data$water_source,data$birth_order_dum,data$delivery_location,data$marital_status,
              stringsAsFactors = FALSE)
y<-data$neonatal_death
                      
# LASSO procedure 
la.eq <- glmnet(X, y, 
                family='binomial', 
                intercept = F, alpha=1)

# Lambda selection
print(la.eq ) 

# Plot Extended Data Figure S8 (1)(A)
la.eq <- glmnet(X, y, family="binomial", 
                intercept = F, alpha=1) 

tiff(file = "D:/DHS/neonatal_death_plot.tiff",width = 3500, height = 2000, res = 400)
plot(la.eq,xvar = "lambda", label = F)
dev.off()

la.eq <- glmnet(X, y,lambda="0.001443", family="binomial", 
                intercept = F, alpha=1) # Value for 0.001443 was derived from procedure of lambda selection
cof<-coef(la.eq)
cof
X <- as.matrix(X)
Y <- y

# Run cross-validation & select lambda
mod_cv <- cv.glmnet(x=X, y=Y, family="binomial", nfolds = 10,
                    intercept = F, alpha=1)
plot(mod_cv) 

# Plot Extended Data Figure S8 (1)(B)
tiff(file = "D:/DHS/neonatal_death_plotcvfit.tiff",width = 3500, height = 2000, res = 400)
plot(mod_cv)
dev.off()

# Lambda.min: the λ at which the minimal MSE is achieved.
# Lambda.1se: the largest λ at which the MSE is within one standard error of the minimal MSE.
print(paste(mod_cv$lambda.min,
            log(mod_cv$lambda.min)))
print(paste(mod_cv$lambda.1se,
            log(mod_cv$lambda.1se)))

best_lambda <- mod_cv$lambda.min
best_lambda

# Find coefficients of best model
best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
coef(best_model)


# LASSO for infant mortality 
# This sub-section produces the following: 
# (1) Extended Data Figure S8 (2)(C)
# (2) Extended Data Figure S8 (2)(D)

X<-data.frame(data$maternal_age,data$kid_sex,data$place_of_residence,data$education_dum,data$household_wealth_dum,
              data$toilet,data$water_source,data$birth_order_dum,data$delivery_location,data$marital_status,
              stringsAsFactors = FALSE)
y<-data$infant_death

# LASSO procedure 
la.eq <- glmnet(X, y, 
                family='binomial', 
                intercept = F, alpha=1)

# Lambda selection
print(la.eq ) 

# Plot Extended Data Figure S8 (2)(C)
la.eq <- glmnet(X, y, family="binomial", 
                intercept = F, alpha=1) 

tiff(file = "D:/DHS/infant_death_plot.tiff",width = 3500, height = 2000, res = 400)
plot(la.eq,xvar = "lambda", label = F)
dev.off()

la.eq <- glmnet(X, y,lambda="0.001255", family="binomial", 
                intercept = F, alpha=1) # Value for 0.001255 was derived from procedure of lambda selection
cof<-coef(la.eq)
cof
X <- as.matrix(X)
Y <- y

# Run cross-validation & select lambda
mod_cv <- cv.glmnet(x=X, y=Y, family="binomial", nfolds = 10,
                    intercept = F, alpha=1)
plot(mod_cv) 

# Plot Extended Data Figure S8 (2)(D)
tiff(file = "D:/DHS/infant_death_plotcvfit.tiff",width = 3500, height = 2000, res = 400)
plot(mod_cv)
dev.off()

# Lambda.min: the λ at which the minimal MSE is achieved.
# Lambda.1se: the largest λ at which the MSE is within one standard error of the minimal MSE.
print(paste(mod_cv$lambda.min,
            log(mod_cv$lambda.min)))
print(paste(mod_cv$lambda.1se,
            log(mod_cv$lambda.1se)))

best_lambda <- mod_cv$lambda.min
best_lambda

# Find coefficients of best model
best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
coef(best_model)


# LASSO for under-five mortality 
# This sub-section produces the following: 
# (1) Extended Data Figure S8 (3)(E)
# (2) Extended Data Figure S8 (3)(F)

X<-data.frame(data$maternal_age,data$kid_sex,data$place_of_residence,data$education_dum,data$household_wealth_dum,
              data$toilet,data$water_source,data$birth_order_dum,data$delivery_location,data$marital_status,
              stringsAsFactors = FALSE)
y<-data$under_five_death

# LASSO procedure 
la.eq <- glmnet(X, y, 
                family='binomial', 
                intercept = F, alpha=1)

# Lambda selection
print(la.eq ) 

# Plot Extended Data Figure S8 (3)(E)
la.eq <- glmnet(X, y, family="binomial", 
                intercept = F, alpha=1) 

tiff(file = "D:/DHS/under_five_death_plot.tiff",width = 3500, height = 2000, res = 400)
plot(la.eq,xvar = "lambda", label = F)
dev.off()

la.eq <- glmnet(X, y,lambda="0.001135", family="binomial", 
                intercept = F, alpha=1) # Value for 0.001135 was derived from procedure of lambda selection
cof<-coef(la.eq)
cof
X <- as.matrix(X)
Y <- y

# Run cross-validation & select lambda
mod_cv <- cv.glmnet(x=X, y=Y, family="binomial", nfolds = 10,
                    intercept = F, alpha=1)
plot(mod_cv) 

# Plot Extended Data Figure S8 (3)(F)
tiff(file = "D:/DHS/under_five_death_plotcvfit.tiff",width = 3500, height = 2000, res = 400)
plot(mod_cv)
dev.off()

# Lambda.min : the λ at which the minimal MSE is achieved.
# Lambda.1se : the largest λ at which the MSE is within one standard error of the minimal MSE.
print(paste(mod_cv$lambda.min,
            log(mod_cv$lambda.min)))
print(paste(mod_cv$lambda.1se,
            log(mod_cv$lambda.1se)))

best_lambda <- mod_cv$lambda.min
best_lambda

# Find coefficients of best model
best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
coef(best_model)



#### Step 2.1. Main models: Deriving results for the associations between ENSO levels and child survival (1) ####

##  This section produces the following: 
##  (1) Figure 2 (overall cumulative association), MEI
##  (2) Figure 3 (contour plot), MEI
##  NOTE: These sample codes do not include style (e.g., font size, color) modifications

## Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"

outcome <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                         "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                         "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                         "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with MEI
MEI<-data_IPUMS[,c("MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6",
                   "MEI_lag7","MEI_lag8","MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12")]
colnames(MEI)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbMEI <- crossbasis(MEI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between MEI and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbMEI+factor(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)


# Plot the contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="MEI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="MEI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

## Association between MEI and infant mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbMEI+factor(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text 
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="MEI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="MEI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

## Association between MEI and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbMEI+factor(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 2 in the main text
plot(pred,"contour", xlab="MEI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="MEI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))



#### Step 2.2. Main models: Deriving results for the associations between ENSO levels and child survival (2) ####

##  This section produces the following: 
##  (1) Figure 2 (overall cumulative association), ESPI
##  (2) Figure 3 (contour plot), ESPI
##  NOTE: These sample codes do not include style (e.g., font size, color) modifications

# Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
outcome <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                         "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                         "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                         "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with ESPI
ESPI<-data_IPUMS[,c("ESPI_lag0","ESPI_lag1","ESPI_lag2","ESPI_lag3","ESPI_lag4","ESPI_lag5","ESPI_lag6",
                    "ESPI_lag7","ESPI_lag8","ESPI_lag9","ESPI_lag10","ESPI_lag11","ESPI_lag12")]
colnames(ESPI)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbESPI <- crossbasis(ESPI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between ESPI and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbESPI+factor(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbESPI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbESPI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="ESPI",ylab="",cex.lab=1.5,cex.axis=1.5)


# Plot the contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="ESPI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="ESPI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))


## Association between ESPI and infant mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbESPI+factor(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbESPI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text 
redcum <- crossreduce(cbESPI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="ESPI",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="ESPI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="ESPI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

## Association between ESPI and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbESPI+factor(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbESPI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbESPI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="ESPI",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 2 in the main text
plot(pred,"contour", xlab="ESPI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="ESPI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))



#### Step 2.3. Main models: Deriving results for the associations between ENSO levels and child survival (3) ####

##  This section produces the following: 
##  (1) Figure 2 (overall cumulative association), ONI
##  (2) Figure 3 (contour plot), ONI
##  NOTE: These sample codes do not include style (e.g., font size, color) modifications

# Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
outcome <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                         "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                         "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                         "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with ONI measure
ONI<-data_IPUMS[,c("ONI_lag0","ONI_lag1","ONI_lag2","ONI_lag3","ONI_lag4","ONI_lag5","ONI_lag6",
                   "ONI_lag7","ONI_lag8","ONI_lag9","ONI_lag10","ONI_lag11","ONI_lag12")]
colnames(ONI)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbONI <- crossbasis(ONI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between ONI and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbONI+factor(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbONI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbONI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="ONI",ylab="",cex.lab=1.5,cex.axis=1.5)


# Plot the contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="ONI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="ONI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))


## Association between ONI and infant mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbONI+factor(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbONI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text 
redcum <- crossreduce(cbONI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="ONI",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="ONI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="ONI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

## Association between ONI and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbONI+factor(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbONI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbONI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="ONI",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 2 in the main text
plot(pred,"contour", xlab="ONI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="ONI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))



#### Step 2.4. Main models: Deriving results for the associations between ENSO levels and child survival (4) ####

##  This section produces the following: 
##  (1) Figure 2 (overall cumulative association), Nino 1+2
##  (2) Figure 3 (contour plot), Nino 1+2
##  NOTE: These sample codes do not include style (e.g., font size, color) modifications

# Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
outcome <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                         "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                         "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                         "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with Nino12 measure
Nino12<-data_IPUMS[,c("Nino12_lag0","Nino12_lag1","Nino12_lag2","Nino12_lag3","Nino12_lag4","Nino12_lag5","Nino12_lag6",
                      "Nino12_lag7","Nino12_lag8","Nino12_lag9","Nino12_lag10","Nino12_lag11","Nino12_lag12")]
colnames(Nino12)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbNino12 <- crossbasis(Nino12,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between Nino12 and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbNino12+factor(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbNino12,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbNino12, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="Nino12",ylab="",cex.lab=1.5,cex.axis=1.5)


# Plot contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="Nino12", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="Nino12",ylab="Lag (months)",cex.main=2,cex.lab=1.5))


## Association between Nino12 and infant mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbNino12+factor(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbNino12,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text 
redcum <- crossreduce(cbNino12, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="Nino12",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="Nino12", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="Nino12",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

## Association between Nino12 and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbNino12+factor(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbNino12,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbNino12, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="Nino12",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 2 in the main text
plot(pred,"contour", xlab="Nino12", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="Nino12",ylab="Lag (months)",cex.main=2,cex.lab=1.5))



#### Step 2.5. Main models: Deriving results for the associations between ENSO levels and child survival (5) ####

##  This section produces the following: 
##  (1) Figure 2 (overall cumulative association), Nino 3.4
##  (2) Figure 3 (contour plot), Nino 3.4
##  NOTE: These sample codes do not include style (e.g., font size, color) modifications

# Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
outcome <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                         "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                         "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                         "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with Nino34 measure
Nino34<-data_IPUMS[,c("Nino34_lag0","Nino34_lag1","Nino34_lag2","Nino34_lag3","Nino34_lag4","Nino34_lag5","Nino34_lag6",
                      "Nino34_lag7","Nino34_lag8","Nino34_lag9","Nino34_lag10","Nino34_lag11","Nino34_lag12")]
colnames(Nino34)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbNino34 <- crossbasis(Nino34,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between Nino34 and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbNino34+factor(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbNino34,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbNino34, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="Nino34",ylab="",cex.lab=1.5,cex.axis=1.5)


# Plot the contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="Nino34", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="Nino34",ylab="Lag (months)",cex.main=2,cex.lab=1.5))


## Association between Nino34 and infant mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbNino34+factor(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbNino34,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text 
redcum <- crossreduce(cbNino34, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="Nino34",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 3 in the main text
plot(pred,"contour", xlab="Nino34", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="Nino34",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

## Association between Nino34 and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbNino34+factor(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbNino34,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Figure 2 in the main text
redcum <- crossreduce(cbNino34, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="Nino34",ylab="",cex.lab=1.5,cex.axis=1.5)

# Plot the contour presented in Figure 2 in the main text
plot(pred,"contour", xlab="Nino34", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="Nino34",ylab="Lag (months)",cex.main=2,cex.lab=1.5))



#### Step 2.6. Attributable fractions: Estimating the proportion of child deaths that can be attributed to ENSO exposure level ####

##  This section produces the following: 
##  Extended Data Table S6
##  NOTE: Here we take MEI as an example; codes for other ENSO indices are similar

## Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
data <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                      "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                      "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                      "birth_order","delivery_location","marital_status","MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6",
                      "MEI_lag7","MEI_lag8","MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12")]
colnames(data)

# Create categorical variables for ENSO measures at 0th-12th lagged months
data$MEI_mean0_12<-(data$MEIlag0+data$MEIlag1+data$MEIlag2+data$MEIlag3+data$MEIlag4+data$MEIlag5+
                      data$MEIlag6+data$MEIlag7+data$MEIlag8+data$MEIlag9+data$MEIlag10+data$MEIlag11+data$MEIlag12)/13

data$MEI_group<-ifelse(data$MEI_mean0_12 <= -0.5,"1",  # La Niña condition
                       ifelse(data$MEI_mean0_12>-0.5 & data$MEI_mean0_12<0.5,"0", # Neutral condition
                              ifelse(data$MEI_mean0_12>=0.5 & data$MEI_mean0_12<1,"2", # Weak El Niño condition
                                     ifelse(data$MEI_mean0_12>=1,"3",NA)))) # Mderate El Niño condition

data$MEI_group<-as.factor(data$MEI_group)

## Association between MEI_group and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~MEI_group+factor(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=data)    
summary(res.cox)

## Association between MEI_group and infant mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~MEI_group+factor(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=data)     
summary(res.cox)

## Association between MEI_group and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~MEI_group+factor(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=data)                                                        
summary(res.cox)



#### Step 3. Stratified analyses: Results by subgroups defined by characteristics of study participants ####

##  This section produces the following: 
##  (1) Figure 4; Extended Data Figure S2
##  (2) Extended Data Tables S4 and S5
##  NOTE: Here we take subgroups by sex as an illustrative example to avoid repetition.

# NOTE: The illustration below only takes MEI as a example; codes for other ENSO indices are similar

### Associations with MEI among boys
subgroup<-subset(data_IPUMS,data_IPUMS$kid_sex==1) # 1: boy; 0: girl 

outcome <- subgroup[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                       "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                       "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                       "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with MEI
MEI<-subgroup[,c("MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6",
                 "MEI_lag7","MEI_lag8","MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12")]
colnames(MEI)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbMEI <- crossbasis(MEI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between MEI and neonatal mortality presented in Extended Data Figure S2
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbMEI+factor(country)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

write.csv(pred$cumfit,file="D:/DHS/neonatal_death_cumest.csv")
write.csv(pred$cumse,file="D:/DHS/neonatal_death_cumse.csv")

## Association between MEI and infant mortality presented in Figure 4
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbMEI+factor(country)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

write.csv(pred$cumfit,file="D:/DHS/infant_death_cumest.csv")
write.csv(pred$cumse,file="D:/DHS/infant_death_cumse.csv")

## Association between MEI and under-five mortality presented in Figure 4
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbMEI+factor(country)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

write.csv(pred$cumfit,file="D:/DHS/iunder_five_death_cumest.csv")
write.csv(pred$cumse,file="D:/DHS/under_five_death_cumse.csv")

### Associations with MEI among girls
subgroup<-subset(data_IPUMS,data_IPUMS$kid_sex==0) # 1: boy; 0: girl 

outcome <- subgroup[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                       "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                       "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                       "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with MEI
MEI<-subgroup[,c("MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6",
                 "MEI_lag7","MEI_lag8","MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12")]
colnames(MEI)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbMEI <- crossbasis(MEI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between MEI and neonatal mortality presented in Extended Data Figure S2
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbMEI+factor(country)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

write.csv(pred$cumfit,file="D:/DHS/neonatal_death_cumest.csv")
write.csv(pred$cumse,file="D:/DHS/neonatal_death_cumse.csv")

## Association between MEI and infant mortality presented in Figure 4
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbMEI+factor(country)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

write.csv(pred$cumfit,file="D:/DHS/infant_death_cumest.csv")
write.csv(pred$cumse,file="D:/DHS/infant_death_cumse.csv")

## Association between MEI and under-five mortality presented in Figure 4
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbMEI+factor(country)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

write.csv(pred$cumfit,file="D:/DHS/iunder_five_death_cumest.csv")
write.csv(pred$cumse,file="D:/DHS/under_five_death_cumse.csv")



#### Step 4. Meta-analyses: Country-specific analyses for associations between ENSO levels and child survival ####

##  This section produces the following: 
##  (1) Figure 5 
##  (2) Extended Data Figures S3 and S4
##  NOTE: Here we take one country as an illustrative example to avoid repetition.

# NOTE: The illustration below only takes Angola and MEI as a example; codes for other countries and indices are similar

### Country-specific analyses
Angola<-subset(data_IPUMS,data_IPUMS$country==24) # 24 is identifying code for Angola 
# NOTE: For more details of each country code, please visit the DHS program website: https://dhsprogram.com/)

# Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
outcome <- Angola[,c("neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                     "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                     "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                     "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with MEI
MEI<-Angola[,c("MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6",
               "MEI_lag7","MEI_lag8","MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12")]
colnames(MEI)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbMEI <- crossbasis(MEI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between MEI and neonatal mortality presented in Extended Data Figures S3 and S4
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbMEI+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

write.csv(pred$cumfit,file="D:/DHS/iunder_five_death_Angola_cumest.csv")
write.csv(pred$cumse,file="D:/DHS/under_five_death_Angola_cumse.csv")

## Association between MEI and infant mortality presented in Extended Data Figures S3 and S4
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbMEI+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

write.csv(pred$cumfit,file="D:/DHS/iunder_five_death_Angola_cumest.csv")
write.csv(pred$cumse,file="D:/DHS/under_five_death_Angola_cumse.csv")

## Association between MEI and under-five mortality presented in Figure 5 and Extended Data Figure S4
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbMEI+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

write.csv(pred$cumfit,file="D:/DHS/iunder_five_death_Angola_cumest.csv")
write.csv(pred$cumse,file="D:/DHS/under_five_death_Angola_cumse.csv")


# Meta-analyses are conducted by pooling 
# the estimates of the above 
# country-specific associations using 
# random effects models 
#
# This sub-section produces the following: 
# (1) Figure 5 
# (2) Extended Data Figures S3 and S4

## Results of meta-analysis for association between MEI and neonatal mortality 
#  Based on 90th percentile of MEI
data_meta<-read.csv("D:/DHS/Meta/neonatal_MEI_P90.csv")
lnHR<-log(data_meta[,"HR"]) #estimate
lnHRL<-log(data_meta[,"HR.L"]) # 95% confidence interval (lower)
lnHRH<-log(data_meta[,"HR.H"]) # 95% confidence interval (upper)
selnHR<-(lnHRH-lnHRL)/(2*1.96)
tiff(file="D:/DHS/Meta/neonatal_MEI_P90_plot.tiff",width = 4000,height = 3800,res=400)
m<-metagen(lnHR,selnHR,sm="HR",data=data_meta,fixed=F,random=T,studlab=paste(data_meta$Region,data_meta$Country,sep="  |  "))
forest(m,col.square ="DeepSkyBlue3", col.diamond = "LightCoral",
       digits = 2)
dev.off()

#  Based on 10th percentile of MEI
data_meta<-read.csv("D:/DHS/Meta/neonatal_MEI_P10.csv")
lnHR<-log(data_meta[,"HR"]) #estimate
lnHRL<-log(data_meta[,"HR.L"]) # 95% confidence interval (lower)
lnHRH<-log(data_meta[,"HR.H"]) # 95% confidence interval (upper)
selnHR<-(lnHRH-lnHRL)/(2*1.96)
tiff(file="D:/DHS/Meta/neonatal_MEI_P10_plot.tiff",width = 4000,height = 3800,res=400)
m<-metagen(lnHR,selnHR,sm="HR",data=data_meta,fixed=F,random=T,studlab=paste(data_meta$Region,data_meta$Country,sep="  |  "))
forest(m,col.square ="DeepSkyBlue3", col.diamond = "LightCoral",
       digits = 2)
dev.off()

## Results of meta-analysis for association between MEI and infant mortality 
#  Based on 90th percentile of MEI
data_meta<-read.csv("D:/DHS/Meta/infant_MEI_P90.csv")
lnHR<-log(data_meta[,"HR"]) #estimate
lnHRL<-log(data_meta[,"HR.L"]) # 95% confidence interval (lower)
lnHRH<-log(data_meta[,"HR.H"]) # 95% confidence interval (upper)
selnHR<-(lnHRH-lnHRL)/(2*1.96)
tiff(file="D:/DHS/Meta/infant_MEI_P90_plot.tiff",width = 4000,height = 3800,res=400)
m<-metagen(lnHR,selnHR,sm="HR",data=data_meta,fixed=F,random=T,studlab=paste(data_meta$Region,data_meta$Country,sep="  |  "))
forest(m,col.square ="DeepSkyBlue3", col.diamond = "LightCoral",
       digits = 2)
dev.off()

#  Based on 10th percentile of MEI
data_meta<-read.csv("D:/DHS/Meta/infant_MEI_P10.csv")
lnHR<-log(data_meta[,"HR"]) #estimate
lnHRL<-log(data_meta[,"HR.L"]) # 95% confidence interval (lower)
lnHRH<-log(data_meta[,"HR.H"]) # 95% confidence interval (upper)
selnHR<-(lnHRH-lnHRL)/(2*1.96)
tiff(file="D:/DHS/Meta/infant_MEI_P10_plot.tiff",width = 4000,height = 3800,res=400)
m<-metagen(lnHR,selnHR,sm="HR",data=data_meta,fixed=F,random=T,studlab=paste(data_meta$Region,data_meta$Country,sep="  |  "))
forest(m,col.square ="DeepSkyBlue3", col.diamond = "LightCoral",
       digits = 2)
dev.off()

## Results of meta-analysis for association between MEI and Under-five mortality 
#  Based on 90th percentile of MEI
data_meta<-read.csv("D:/DHS/Meta/Under_five_MEI_P90.csv")
colnames(a1)
lnHR<-log(data_meta[,"HR"]) #estimate
lnHRL<-log(data_meta[,"HR.L"]) # 95% confidence interval (lower)
lnHRH<-log(data_meta[,"HR.H"]) # 95% confidence interval (upper)
selnHR<-(lnHRH-lnHRL)/(2*1.96)
tiff(file="D:/DHS/Meta/Under_five_MEI_P90_plot.tiff",width = 4000,height = 3800,res=400)
m<-metagen(lnHR,selnHR,sm="HR",data=data_meta,fixed=F,random=T,studlab=paste(data_meta$Region,data_meta$Country,sep="  |  "))
forest(m,col.square ="DeepSkyBlue3", col.diamond = "LightCoral",
       digits = 2)
dev.off()

#  Based on 10th percentile of MEI
data_meta<-read.csv("D:/DHS/Meta/Under_five_MEI_P10.csv")
colnames(a1)
lnHR<-log(data_meta[,"HR"]) #estimate
lnHRL<-log(data_meta[,"HR.L"]) # 95% confidence interval (lower)
lnHRH<-log(data_meta[,"HR.H"]) # 95% confidence interval (upper)
selnHR<-(lnHRH-lnHRL)/(2*1.96)
tiff(file="D:/DHS/Meta/Under_five_MEI_P10_plot.tiff",width = 4000,height = 3800,res=400)
m<-metagen(lnHR,selnHR,sm="HR",data=data_meta,fixed=F,random=T,studlab=paste(data_meta$Region,data_meta$Country,sep="  |  "))
forest(m,col.square ="DeepSkyBlue3", col.diamond = "LightCoral",
       digits = 2)
dev.off()



#### Step 5.1. Sensitivity analysis (1) ####
##  (a) Assessing the impact of ENSO exposure during the mothers’ preconceptional periods
##  (b) Assessing the impact of ENSO exposure during the mothers’ prenatal periods

##  This section produces the following: 
##  (1) Extended Data Figure S5
##  (2) Extended Data Figure S6
##  NOTE: Here we take the analysis for association between MEI during the mothers’ 
##  preconceptional and prenatal periods of exposure and child survival as a example, 
##  and the codes for other ENSO measures are similar.

### Preconceptional periods
# Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
outcome <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                         "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                         "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                         "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with MEI (lag 9-12 months)
MEI_lag9_12<-data_IPUMS[,c("MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12")]
MEI<-rename(MEI_lag9_12,c("MEIlag0"="MEIlag9","MEIlag1"="MEIlag10","MEIlag2"="MEIlag11","MEIlag3"="MEIlag12"))
colnames(MEI)
colnames(MEI)

# Generate crossbasis matrix 
lk = logknots(c(0,3),2)
cbMEI <- crossbasis(MEI,lag=c(0,3),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between MEI and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbMEI+factor(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S5 
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=3,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)


## Association between MEI and infant mortality 
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbMEI+factor(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S5 
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=3,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

## Association between MEI and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbMEI+factor(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S5 
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=3,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)


### Prenatal periods
# Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
outcome <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                         "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                         "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                         "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with MEI (lag 0-8 months)
MEI<-data_IPUMS[,c("MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6",
                   "MEI_lag7","MEI_lag8")]
colnames(MEI)

# Generate crossbasis matrix 
lk = logknots(c(0,8),2)
cbMEI <- crossbasis(MEI,lag=c(0,8),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between MEI and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbMEI+factor(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S6
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=8,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

## Association between MEI and infant mortality 
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbMEI+factor(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S6 
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=8,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

## Association between MEI and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbMEI+factor(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S6 
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=8,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)




#### Step 5.2. Sensitivity analysis (2) ####
##  Repeating the main analysis by including a random effect for each country

##  This section produces the following: 
##  Extended Data Figure S7A
##  NOTE: We ake the analysis for the association between MEI and child survival as a example, 
##  and the codes for other ENSO measures are similar.

# Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
outcome <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                         "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                         "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                         "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with MEI
MEI<-data_IPUMS[,c("MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6",
                   "MEI_lag7","MEI_lag8","MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12")]
colnames(MEI)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbMEI <- crossbasis(MEI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between MEI and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbMEI+frailty(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S7A
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

## Association between MEI and infant mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbMEI+frailty(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S7A 
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

## Association between MEI and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbMEI+frailty(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S7A
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)



#### Step 5.3. Sensitivity analysis (3) #### 
##  Repeating the main analysis by including a random effect for each DHS survey cluster

##  This section produces the following: 
##  Extended Data Figure S7B
##  NOTE: We ake the analysis for the association between MEI and child survival as a example, 
##  and the codes for other ENSO measures are similar.

# Create data frame with outcomes and covariates
# NOTE: Covariates with three or more categories were firstly converted to categorical dummy variables
# using package "neuralnet"
outcome <- data_IPUMS[,c("dhsid","neonatal_death","infant_death","under_five_death","survival_time_neonatal",
                         "survival_time_infant","survival_time_under_five","maternal_age","kid_sex",
                         "place_of_residence", "education", "household_wealth", "toilet", "water_source",
                         "birth_order","delivery_location","marital_status")]
colnames(outcome)

# Create data frame with MEI
MEI<-data_IPUMS[,c("MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6",
                   "MEI_lag7","MEI_lag8","MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12")]
colnames(MEI)

# Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbMEI <- crossbasis(MEI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  

## Association between MEI and neonatal mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbMEI+frailty(dhsid)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=outcome)    
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S7B
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

## Association between MEI and infant mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbMEI+frailty(dhsid)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S7B 
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

## Association between MEI and under-five mortality
# Fit the Cox model by including selected covariates based on LASSO 
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbMEI+frailty(dhsid)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=outcome)                                                        
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

# Plot the overall cumulative association presented in Extended Data Figure S7B
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

