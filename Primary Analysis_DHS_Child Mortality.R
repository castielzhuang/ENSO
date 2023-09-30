################################################################################################
####            Preconceptional and Prenatal Exposure to Climate Variability                ####
####   driven by  El Ni??o Southern Oscillation and Child Mortality: A Multi-country Study   ####
####                                                                                        ####
####       Hongbing XU, ECastiel Chen ZHUANG, Vanessa M. ODDO, Espoir Bwenge MALEMBAKA,     ####
####                     YXinghou HE, Qinghong ZHANG, Wei HUANG                             ####
####                                                                                        ####
####                              October,2022-August,2023                                  ####
################################################################################################


# Load Package 
library("survival");library(splines);library(dlnm);library(ggplot2)

################# Primary analysis ###############################
## data_IPUMS included the following main variables:
## "country","neonatal_death","infant_death","under_five_death","survival_time_neonatal","survival_time_infant","survival_time_under_five",
## "kid_sex",maternal_age","place_of_residence", "education", "household_wealth", "toilet", "water_source","birth_order","delivery_location","marital_status",
## "MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6","MEI_lag7","MEI_lag8","MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12",
## "ESPI_lag0","ESPI_lag1","ESPI_lag2","ESPI_lag3","ESPI_lag4","ESPI_lag5","ESPI_lag6","ESPI_lag7","ESPI_lag8","ESPI_lag9","ESPI_lag10","ESPI_lag11","ESPI_lag12",
## "ONI_lag0","ONI_lag1","ONI_lag2","ONI_lag3","ONI_lag4","ONI_lag5","ONI_lag6","ONI_lag7","ONI_lag8","ONI_lag9","ONI_lag10","ONI_lag11","ONI_lag12",
## "Ni??o12_lag0","Ni??o12_lag1","Ni??o12_lag2","Ni??o12_lag3","Ni??o12_lag4","Ni??o12_lag5","Ni??o12_lag6","Ni??o12_lag7","Ni??o12_lag8","Ni??o12_lag9","Ni??o12_lag10","Ni??o12_lag11","Ni??o12_lag12",
## "Ni??o34_lag0","Ni??o34_lag1","Ni??o34_lag2","Ni??o34_lag3","Ni??o34_lag4","Ni??o34_lag5","Ni??o34_lag6","Ni??o34_lag7","Ni??o34_lag8","Ni??o34_lag9","Ni??o34_lag10","Ni??o34_lag11","Ni??o34_lag12",


#Create dataframe with outcomes and covariates
data <- data_IPUMS[,c("country","neonatal_death","infant_death","under_five_death","survival_time_neonatal","survival_time_infant","survival_time_under_five",
                     "maternal_age","kid_sex","place_of_residence", "education", "household_wealth", "toilet", "water_source","birth_order","delivery_location","marital_status")]
colnames(data)

#Create dataframe with ENSO indices
MEI<-data_IPUMS[,c("MEI_lag0","MEI_lag1","MEI_lag2","MEI_lag3","MEI_lag4","MEI_lag5","MEI_lag6","MEI_lag7","MEI_lag8","MEI_lag9","MEI_lag10","MEI_lag11","MEI_lag12")]
colnames(MEI)

ESPI<-data_IPUMS[,c("ESPI_lag0","ESPI_lag1","ESPI_lag2","ESPI_lag3","ESPI_lag4","ESPI_lag5","ESPI_lag6","ESPI_lag7","ESPI_lag8","ESPI_lag9","ESPI_lag10","ESPI_lag11","ESPI_lag12")]
colnames(ESPI)

ONI<-data_IPUMS[,c("ONI_lag0","ONI_lag1","ONI_lag2","ONI_lag3","ONI_lag4","ONI_lag5","ONI_lag6","ONI_lag7","ONI_lag8","ONI_lag9","ONI_lag10","ONI_lag11","ONI_lag12")]
colnames(ONI)

Ni??o12<-data_IPUMS[,c("Ni??o12_lag0","Ni??o12_lag1","Ni??o12_lag2","Ni??o12_lag3","Ni??o12_lag4","Ni??o12_lag5","Ni??o12_lag6","Ni??o12_lag7","Ni??o12_lag8","Ni??o12_lag9","Ni??o12_lag10","Ni??o12_lag11","Ni??o12_lag12")]
colnames(Ni??o12)

Ni??o34<-data_IPUMS[,c("Ni??o34_lag0","Ni??o34_lag1","Ni??o34_lag2","Ni??o34_lag3","Ni??o34_lag4","Ni??o34_lag5","Ni??o34_lag6","Ni??o34_lag7","Ni??o34_lag8","Ni??o34_lag9","Ni??o34_lag10","Ni??o34_lag11","Ni??o34_lag12")]
colnames(Ni??o34)


#Generate crossbasis matrix 
lk = logknots(c(0,12),2)
cbMEI <- crossbasis(MEI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk))  
cbESPI <- crossbasis(ESPI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk)) 
cbONI <- crossbasis(ONI,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk)) 
cbNi??o12 <- crossbasis(Ni??o12,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk)) 
cbNi??o34 <- crossbasis(Ni??o34,lag=c(0,12),argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=lk)) 



####################################################################################
#
#  Characterize the association between ENSO exposure (e.g., MEI) and child survival 
#
####################################################################################

### Association between MEI and neonatal mortality
#Fit the cox model and predict
res.cox <-coxph(Surv(survival_time_neonatal,neonatal_death)~cbMEI+factor(country)+factor(kid_sex)+
                  maternal_age+I(maternal_age^2)+factor(delivery_location)+factor(toilet)+
                  factor(water_source)+
                  factor(household_wealth3)+factor(birth_order2)+factor(birth_order3),data=data)    
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

#### Plot overall cumulative association
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)


###Plot contour
plot(pred,"contour", xlab="MEI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="MEI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))




### Association between MEI and infant mortality
#Fit the cox model and predict
res.cox <-coxph(Surv(survival_time_infant,infant_death)~cbMEI+factor(country)+factor(kid_sex)+maternal_age+
                  I(maternal_age^2)+factor(delivery_location)+
                  factor(toilet)+factor(marital_status)+factor(water_source)+factor(education3)+
                  factor(birth_order2)+factor(birth_order3)+factor(birth_order4),data=outcome)     
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

#### Plot overall cumulative association
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)

###Plot contour
plot(pred,"contour", xlab="MEI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="MEI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))




### Association between MEI and under-five mortality
#Fit the cox model and predict
res.cox <-coxph(Surv(survival_time_under_five,under_five_death)~cbMEI+factor(country)+factor(kid_sex)+ maternal_age+
                  I(maternal_age^2)+ factor(delivery_location)+factor(toilet)+factor(marital_status)+factor(water_source)+
                  factor(education2)+factor(education3)+ factor(birth_order2)+factor(birth_order3)+
                  factor(birth_order4),data=data)                                                        
pred <- crosspred(cbMEI,res.cox,cumul=TRUE,cen=0)

#### Plot overall cumulative association
redcum <- crossreduce(cbMEI, res.cox, type="overall", lag=12,cen=0)
plot<-plot(redcum,col=1,lty=1,ci.arg=list(col="Gainsboro"),lwd=3,xlab="MEI",ylab="",cex.lab=1.5,cex.axis=1.5)


###Plot contour
plot(pred,"contour", xlab="MEI", key.title=title("HR"),cex=5,cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title("Contour plot",xlab="MEI",ylab="Lag (months)",cex.main=2,cex.lab=1.5))


