################################################################################
### R code for "Time lags in the regulation of symbiotic nitrogen fixation"
### Authors: Thomas A. Bytnerowicz, Kevin L. Griffin, Duncan N.L. Menge
################################################################################

################################################################################
### R code first performs analyses and then generates figures
### Analyses must be run first for figures to be generated correctly
################################################################################

### Load necessary packages
library(bbmle)
library(MASS)

### Define necessary functions

# Equation 5
sigmoid<-function(K,r,tL,t){
  y<-K/(1+exp(-r*(t-tL)))
  y
}

# Up-regulation curve from 5% up-regulation
up.rate<-function(K,r,t){
  y<-K/(1+exp(-r*(t-log(19)/r)))
  y
}

### Import data
snf.data<-read.csv("snf_time_lags.csv")
co2.data<-read.csv("co2_time_lags.csv")
rinse<-read.csv("rinse_test.csv") #data for Figure S1

### Subset N fixation data
A21DA<-snf.data[snf.data$Individual=="ALRU21DA",]
A21DB<-snf.data[snf.data$Individual=="ALRU21DB",]
A21DC<-snf.data[snf.data$Individual=="ALRU21DC",]
A21DD<-snf.data[snf.data$Individual=="ALRU21DD",]

G21DA<-snf.data[snf.data$Individual=="GLSE21DA",]
G21DB<-snf.data[snf.data$Individual=="GLSE21DB",]
G21DC<-snf.data[snf.data$Individual=="GLSE21DC",]
G21DD<-snf.data[snf.data$Individual=="GLSE21DD",]

M21DA<-snf.data[snf.data$Individual=="MOCE21DA",]
M21DB<-snf.data[snf.data$Individual=="MOCE21DB",]
M21DC<-snf.data[snf.data$Individual=="MOCE21DC",]
M21DD<-snf.data[snf.data$Individual=="MOCE21DD",]

R21DA<-snf.data[snf.data$Individual=="ROPS21DA",]
R21DB<-snf.data[snf.data$Individual=="ROPS21DB",]
R21DC<-snf.data[snf.data$Individual=="ROPS21DC",]
R21DD<-snf.data[snf.data$Individual=="ROPS21DD",]

G21UA<-snf.data[snf.data$Individual=="GLSE21UA",]
G21UB<-snf.data[snf.data$Individual=="GLSE21UB",]
G21UC<-snf.data[snf.data$Individual=="GLSE21UC",]

M21UA<-snf.data[snf.data$Individual=="MOCE21UA",]
M21UB<-snf.data[snf.data$Individual=="MOCE21UB",]
M21UC<-snf.data[snf.data$Individual=="MOCE21UC",]
M21UD<-snf.data[snf.data$Individual=="MOCE21UD",]

R21UA<-snf.data[snf.data$Individual=="ROPS21UA",]
R21UB<-snf.data[snf.data$Individual=="ROPS21UB",]
R21UC<-snf.data[snf.data$Individual=="ROPS21UC",]

A31DA<-snf.data[snf.data$Individual=="ALRU31DA",]
A31DB<-snf.data[snf.data$Individual=="ALRU31DB",]
A31DC<-snf.data[snf.data$Individual=="ALRU31DC",]
A31DD<-snf.data[snf.data$Individual=="ALRU31DD",]

G31DA<-snf.data[snf.data$Individual=="GLSE31DA",]
G31DB<-snf.data[snf.data$Individual=="GLSE31DB",]
G31DC<-snf.data[snf.data$Individual=="GLSE31DC",]
G31DD<-snf.data[snf.data$Individual=="GLSE31DD",]

M31DA<-snf.data[snf.data$Individual=="MOCE31DA",]
M31DB<-snf.data[snf.data$Individual=="MOCE31DB",]
M31DC<-snf.data[snf.data$Individual=="MOCE31DC",]
M31DD<-snf.data[snf.data$Individual=="MOCE31DD",]

R31DA<-snf.data[snf.data$Individual=="ROPS31DA",]
R31DB<-snf.data[snf.data$Individual=="ROPS31DB",]
R31DC<-snf.data[snf.data$Individual=="ROPS31DC",]
R31DD<-snf.data[snf.data$Individual=="ROPS31DD",]

A31UA<-snf.data[snf.data$Individual=="ALRU31UA",]
A31UB<-snf.data[snf.data$Individual=="ALRU31UB",]
A31UC<-snf.data[snf.data$Individual=="ALRU31UC",]

G31UA<-snf.data[snf.data$Individual=="GLSE31UA",]
G31UB<-snf.data[snf.data$Individual=="GLSE31UB",]
G31UC<-snf.data[snf.data$Individual=="GLSE31UC",]
G31UD<-snf.data[snf.data$Individual=="GLSE31UD",]
G31UE<-snf.data[snf.data$Individual=="GLSE31UE",]
G31UF<-snf.data[snf.data$Individual=="GLSE31UF",]

M31UA<-snf.data[snf.data$Individual=="MOCE31UA",]
M31UB<-snf.data[snf.data$Individual=="MOCE31UB",]
M31UC<-snf.data[snf.data$Individual=="MOCE31UC",]
M31UD<-snf.data[snf.data$Individual=="MOCE31UD",]
M31UE<-snf.data[snf.data$Individual=="MOCE31UE",]
M31UF<-snf.data[snf.data$Individual=="MOCE31UF",]
M31UG<-snf.data[snf.data$Individual=="MOCE31UG",]
M31UH<-snf.data[snf.data$Individual=="MOCE31UH",]

R31UA<-snf.data[snf.data$Individual=="ROPS31UA",]
R31UB<-snf.data[snf.data$Individual=="ROPS31UB",]
R31UC<-snf.data[snf.data$Individual=="ROPS31UC",]
R31UD<-snf.data[snf.data$Individual=="ROPS31UD",]
R31UE<-snf.data[snf.data$Individual=="ROPS31UE",]
R31UF<-snf.data[snf.data$Individual=="ROPS31UF",]

### Subset CO2 data
A21DA.co2<-co2.data[co2.data$Individual=="ALRU21DA",]
A21DB.co2<-co2.data[co2.data$Individual=="ALRU21DB",]
A21DC.co2<-co2.data[co2.data$Individual=="ALRU21DC",]
A21DD.co2<-co2.data[co2.data$Individual=="ALRU21DD",]

G21DA.co2<-co2.data[co2.data$Individual=="GLSE21DA",]
G21DB.co2<-co2.data[co2.data$Individual=="GLSE21DB",]
G21DC.co2<-co2.data[co2.data$Individual=="GLSE21DC",]
G21DD.co2<-co2.data[co2.data$Individual=="GLSE21DD",]

M21DA.co2<-co2.data[co2.data$Individual=="MOCE21DA",]
M21DB.co2<-co2.data[co2.data$Individual=="MOCE21DB",]
M21DC.co2<-co2.data[co2.data$Individual=="MOCE21DC",]
M21DD.co2<-co2.data[co2.data$Individual=="MOCE21DD",]

R21DA.co2<-co2.data[co2.data$Individual=="ROPS21DA",]
R21DB.co2<-co2.data[co2.data$Individual=="ROPS21DB",]
R21DC.co2<-co2.data[co2.data$Individual=="ROPS21DC",]
R21DD.co2<-co2.data[co2.data$Individual=="ROPS21DD",]

G21UA.co2<-co2.data[co2.data$Individual=="GLSE21UA",]
G21UB.co2<-co2.data[co2.data$Individual=="GLSE21UB",]
G21UC.co2<-co2.data[co2.data$Individual=="GLSE21UC",]

M21UA.co2<-co2.data[co2.data$Individual=="MOCE21UA",]
M21UB.co2<-co2.data[co2.data$Individual=="MOCE21UB",]
M21UC.co2<-co2.data[co2.data$Individual=="MOCE21UC",]
M21UD.co2<-co2.data[co2.data$Individual=="MOCE21UD",]
M21UE.co2<-co2.data[co2.data$Individual=="MOCE21UE",]

R21UA.co2<-co2.data[co2.data$Individual=="ROPS21UA",]
R21UB.co2<-co2.data[co2.data$Individual=="ROPS21UB",]
R21UC.co2<-co2.data[co2.data$Individual=="ROPS21UC",]
R21UD.co2<-co2.data[co2.data$Individual=="ROPS21UD",]
R21UE.co2<-co2.data[co2.data$Individual=="ROPS21UE",]

A31DA.co2<-co2.data[co2.data$Individual=="ALRU31DA",]
A31DB.co2<-co2.data[co2.data$Individual=="ALRU31DB",]
A31DC.co2<-co2.data[co2.data$Individual=="ALRU31DC",]
A31DD.co2<-co2.data[co2.data$Individual=="ALRU31DD",]

G31DA.co2<-co2.data[co2.data$Individual=="GLSE31DA",]
G31DB.co2<-co2.data[co2.data$Individual=="GLSE31DB",]
G31DC.co2<-co2.data[co2.data$Individual=="GLSE31DC",]
G31DD.co2<-co2.data[co2.data$Individual=="GLSE31DD",]

M31DA.co2<-co2.data[co2.data$Individual=="MOCE31DA",]
M31DB.co2<-co2.data[co2.data$Individual=="MOCE31DB",]
M31DC.co2<-co2.data[co2.data$Individual=="MOCE31DC",]
M31DD.co2<-co2.data[co2.data$Individual=="MOCE31DD",]

R31DA.co2<-co2.data[co2.data$Individual=="ROPS31DA",]
R31DB.co2<-co2.data[co2.data$Individual=="ROPS31DB",]
R31DC.co2<-co2.data[co2.data$Individual=="ROPS31DC",]
R31DD.co2<-co2.data[co2.data$Individual=="ROPS31DD",]

A31UA.co2<-co2.data[co2.data$Individual=="ALRU31UA",]
A31UB.co2<-co2.data[co2.data$Individual=="ALRU31UB",]
A31UC.co2<-co2.data[co2.data$Individual=="ALRU31UC",]

G31UA.co2<-co2.data[co2.data$Individual=="GLSE31UA",]
G31UB.co2<-co2.data[co2.data$Individual=="GLSE31UB",]
G31UC.co2<-co2.data[co2.data$Individual=="GLSE31UC",]
G31UD.co2<-co2.data[co2.data$Individual=="GLSE31UD",]
G31UE.co2<-co2.data[co2.data$Individual=="GLSE31UE",]
G31UF.co2<-co2.data[co2.data$Individual=="GLSE31UF",]

M31UA.co2<-co2.data[co2.data$Individual=="MOCE31UA",]
M31UB.co2<-co2.data[co2.data$Individual=="MOCE31UB",]
M31UC.co2<-co2.data[co2.data$Individual=="MOCE31UC",]
M31UD.co2<-co2.data[co2.data$Individual=="MOCE31UD",]
M31UE.co2<-co2.data[co2.data$Individual=="MOCE31UE",]
M31UF.co2<-co2.data[co2.data$Individual=="MOCE31UF",]
M31UG.co2<-co2.data[co2.data$Individual=="MOCE31UG",]
M31UH.co2<-co2.data[co2.data$Individual=="MOCE31UH",]

R31UA.co2<-co2.data[co2.data$Individual=="ROPS31UA",]
R31UB.co2<-co2.data[co2.data$Individual=="ROPS31UB",]
R31UC.co2<-co2.data[co2.data$Individual=="ROPS31UC",]
R31UD.co2<-co2.data[co2.data$Individual=="ROPS31UD",]
R31UE.co2<-co2.data[co2.data$Individual=="ROPS31UE",]
R31UF.co2<-co2.data[co2.data$Individual=="ROPS31UF",]

### Fit down-regulation curves

## Species x temperature interaction with no tL parameter (r at treatment level)
sigmoid_down_sptemp_no_tL_normNLL <- function(sdNase,
                                         rdownA21,
                                         rdownG21,
                                         rdownM21,
                                         rdownR21,
                                         rdownA31,
                                         rdownG31,
                                         rdownM31,
                                         rdownR31,
                                         K1,K2,K3,K4,K5,K6,K7,K8,
                                         K9,K10,K11,K12,K13,K14,K15,K16,
                                         K17,K18,K19,K20,K21,K22,K23,K24,
                                         K25,K26,K27,K28,K29,K30,K31,K32,
                                         t1,t2,t3,t4,t5,t6,t7,t8,
                                         t9,t10,t11,t12,t13,t14,t15,t16,
                                         t17,t18,t19,t20,t21,t22,t23,t24,
                                         t25,t26,t27,t28,t29,t30,t31,t32,
                                         Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                         Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                         Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                         Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32){
  Nasemean1 <- sigmoid(exp(K1),rdownA21,0,t1)
  Nasemean2 <- sigmoid(exp(K2),rdownA21,0,t2)
  Nasemean3 <- sigmoid(exp(K3),rdownA21,0,t3)
  Nasemean4 <- sigmoid(exp(K4),rdownA21,0,t4)
  Nasemean5 <- sigmoid(exp(K5),rdownG21,0,t5)
  Nasemean6 <- sigmoid(exp(K6),rdownG21,0,t6)
  Nasemean7 <- sigmoid(exp(K7),rdownG21,0,t7)
  Nasemean8 <- sigmoid(exp(K8),rdownG21,0,t8)
  Nasemean9 <- sigmoid(exp(K9),rdownM21,0,t9)
  Nasemean10 <- sigmoid(exp(K10),rdownM21,0,t10)
  Nasemean11 <- sigmoid(exp(K11),rdownM21,0,t11)
  Nasemean12 <- sigmoid(exp(K12),rdownM21,0,t12)
  Nasemean13 <- sigmoid(exp(K13),rdownR21,0,t13)
  Nasemean14 <- sigmoid(exp(K14),rdownR21,0,t14)
  Nasemean15 <- sigmoid(exp(K15),rdownR21,0,t15)
  Nasemean16 <- sigmoid(exp(K16),rdownR21,0,t16)
  Nasemean17 <- sigmoid(exp(K17),rdownA31,0,t17)
  Nasemean18 <- sigmoid(exp(K18),rdownA31,0,t18)
  Nasemean19 <- sigmoid(exp(K19),rdownA31,0,t19)
  Nasemean20 <- sigmoid(exp(K20),rdownA31,0,t20)
  Nasemean21 <- sigmoid(exp(K21),rdownG31,0,t21)
  Nasemean22 <- sigmoid(exp(K22),rdownG31,0,t22)
  Nasemean23 <- sigmoid(exp(K23),rdownG31,0,t23)
  Nasemean24 <- sigmoid(exp(K24),rdownG31,0,t24)
  Nasemean25 <- sigmoid(exp(K25),rdownM31,0,t25)
  Nasemean26 <- sigmoid(exp(K26),rdownM31,0,t26)
  Nasemean27 <- sigmoid(exp(K27),rdownM31,0,t27)
  Nasemean28 <- sigmoid(exp(K28),rdownM31,0,t28)
  Nasemean29 <- sigmoid(exp(K29),rdownR31,0,t29)
  Nasemean30 <- sigmoid(exp(K30),rdownR31,0,t30)
  Nasemean31 <- sigmoid(exp(K31),rdownR31,0,t31)
  Nasemean32 <- sigmoid(exp(K32),rdownR31,0,t32)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_down_sptemp_no_tL <- mle2(sigmoid_down_sptemp_no_tL_normNLL,start=list(sdNase=-1,
                                                                         K1=log(max(A21DA$SNF)),K2=log(max(A21DB$SNF)),
                                                                         K3=log(max(A21DC$SNF)),K4=log(max(A21DD$SNF)),
                                                                         K5=log(max(G21DA$SNF)),K6=log(max(G21DB$SNF)),
                                                                         K7=log(max(G21DC$SNF)),K8=log(max(G21DD$SNF)),
                                                                         K9=log(max(M21DA$SNF)),K10=log(max(M21DB$SNF)),
                                                                         K11=log(max(M21DC$SNF)),K12=log(max(M21DD$SNF)),
                                                                         K13=log(max(R21DA$SNF)),K14=log(max(R21DB$SNF)),
                                                                         K15=log(max(R21DC$SNF)),K16=log(max(R21DD$SNF)),
                                                                         K17=log(max(A31DA$SNF)),K18=log(max(A31DB$SNF)),
                                                                         K19=log(max(A31DC$SNF)),K20=log(max(A31DD$SNF)),
                                                                         K21=log(max(G31DA$SNF)),K22=log(max(G31DB$SNF)),
                                                                         K23=log(max(G31DC$SNF)),K24=log(max(G31DD$SNF)),
                                                                         K25=log(max(M31DA$SNF)),K26=log(max(M31DB$SNF)),
                                                                         K27=log(max(M31DC$SNF)),K28=log(max(M31DD$SNF)),
                                                                         K29=log(max(R31DA$SNF)),K30=log(max(R31DB$SNF)),
                                                                         K31=log(max(R31DC$SNF)),K32=log(max(R31DD$SNF)),
                                                                         rdownA21=-0.1,rdownG21=-0.1,rdownM21=-0.1,rdownR21=-0.1,
                                                                         rdownA31=-0.1,rdownG31=-0.1,rdownM31=-0.1,rdownR31=-0.1),
                                 data=list(t1=A21DA$Day,t2=A21DB$Day,
                                           t3=A21DC$Day,t4=A21DD$Day,
                                           t5=G21DA$Day,t6=G21DB$Day,
                                           t7=G21DC$Day,t8=G21DD$Day,
                                           t9=M21DA$Day,t10=M21DB$Day,
                                           t11=M21DC$Day,t12=M21DD$Day,
                                           t13=R21DA$Day,t14=R21DB$Day,
                                           t15=R21DC$Day,t16=R21DD$Day,
                                           t17=A31DA$Day,t18=A31DB$Day,
                                           t19=A31DC$Day,t20=A31DD$Day,
                                           t21=G31DA$Day,t22=G31DB$Day,
                                           t23=G31DC$Day,t24=G31DD$Day,
                                           t25=M31DA$Day,t26=M31DB$Day,
                                           t27=M31DC$Day,t28=M31DD$Day,
                                           t29=R31DA$Day,t30=R31DB$Day,
                                           t31=R31DC$Day,t32=R31DD$Day,
                                           Nasedat1=A21DA$SNF,Nasedat2=A21DB$SNF,
                                           Nasedat3=A21DC$SNF,Nasedat4=A21DD$SNF,
                                           Nasedat5=G21DA$SNF,Nasedat6=G21DB$SNF,
                                           Nasedat7=G21DC$SNF,Nasedat8=G21DD$SNF,
                                           Nasedat9=M21DA$SNF,Nasedat10=M21DB$SNF,
                                           Nasedat11=M21DC$SNF,Nasedat12=M21DD$SNF,
                                           Nasedat13=R21DA$SNF,Nasedat14=R21DB$SNF,
                                           Nasedat15=R21DC$SNF,Nasedat16=R21DD$SNF,
                                           Nasedat17=A31DA$SNF,Nasedat18=A31DB$SNF,
                                           Nasedat19=A31DC$SNF,Nasedat20=A31DD$SNF,
                                           Nasedat21=G31DA$SNF,Nasedat22=G31DB$SNF,
                                           Nasedat23=G31DC$SNF,Nasedat24=G31DD$SNF,
                                           Nasedat25=M31DA$SNF,Nasedat26=M31DB$SNF,
                                           Nasedat27=M31DC$SNF,Nasedat28=M31DD$SNF,
                                           Nasedat29=R31DA$SNF,Nasedat30=R31DB$SNF,
                                           Nasedat31=R31DC$SNF,Nasedat32=R31DD$SNF),
                                 control=list(maxit=20000))
summary(fit_sigmoid_down_sptemp_no_tL)

## Species x temperature interaction (with tL; r fit at treatment level)
sigmoid_down_sptemp_normNLL <- function(sdNase,
                                        rdownA21,
                                        rdownG21,
                                        rdownM21,
                                        rdownR21,
                                        rdownA31,
                                        rdownG31,
                                        rdownM31,
                                        rdownR31,
                                        tLdownA21A,tLdownA21B,tLdownA21C,tLdownA21D,
                                        tLdownG21A,tLdownG21B,tLdownG21C,tLdownG21D,
                                        tLdownM21A,tLdownM21B,tLdownM21C,tLdownM21D,
                                        tLdownR21A,tLdownR21B,tLdownR21C,tLdownR21D,
                                        tLdownA31A,tLdownA31B,tLdownA31C,tLdownA31D,
                                        tLdownG31A,tLdownG31B,tLdownG31C,tLdownG31D,
                                        tLdownM31A,tLdownM31B,tLdownM31C,tLdownM31D,
                                        tLdownR31A,tLdownR31B,tLdownR31C,tLdownR31D,
                                        K1,K2,K3,K4,K5,K6,K7,K8,
                                        K9,K10,K11,K12,K13,K14,K15,K16,
                                        K17,K18,K19,K20,K21,K22,K23,K24,
                                        K25,K26,K27,K28,K29,K30,K31,K32,
                                        t1,t2,t3,t4,t5,t6,t7,t8,
                                        t9,t10,t11,t12,t13,t14,t15,t16,
                                        t17,t18,t19,t20,t21,t22,t23,t24,
                                        t25,t26,t27,t28,t29,t30,t31,t32,
                                        Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                        Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                        Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                        Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32){
  Nasemean1 <- sigmoid(exp(K1),rdownA21,exp(tLdownA21A),t1)
  Nasemean2 <- sigmoid(exp(K2),rdownA21,exp(tLdownA21B),t2)
  Nasemean3 <- sigmoid(exp(K3),rdownA21,exp(tLdownA21C),t3)
  Nasemean4 <- sigmoid(exp(K4),rdownA21,exp(tLdownA21D),t4)
  Nasemean5 <- sigmoid(exp(K5),rdownG21,exp(tLdownG21A),t5)
  Nasemean6 <- sigmoid(exp(K6),rdownG21,exp(tLdownG21B),t6)
  Nasemean7 <- sigmoid(exp(K7),rdownG21,exp(tLdownG21C),t7)
  Nasemean8 <- sigmoid(exp(K8),rdownG21,exp(tLdownG21D),t8)
  Nasemean9 <- sigmoid(exp(K9),rdownM21,exp(tLdownM21A),t9)
  Nasemean10 <- sigmoid(exp(K10),rdownM21,exp(tLdownM21B),t10)
  Nasemean11 <- sigmoid(exp(K11),rdownM21,exp(tLdownM21C),t11)
  Nasemean12 <- sigmoid(exp(K12),rdownM21,exp(tLdownM21D),t12)
  Nasemean13 <- sigmoid(exp(K13),rdownR21,exp(tLdownR21A),t13)
  Nasemean14 <- sigmoid(exp(K14),rdownR21,exp(tLdownR21B),t14)
  Nasemean15 <- sigmoid(exp(K15),rdownR21,exp(tLdownR21C),t15)
  Nasemean16 <- sigmoid(exp(K16),rdownR21,exp(tLdownR21D),t16)
  Nasemean17 <- sigmoid(exp(K17),rdownA31,exp(tLdownA31A),t17)
  Nasemean18 <- sigmoid(exp(K18),rdownA31,exp(tLdownA31B),t18)
  Nasemean19 <- sigmoid(exp(K19),rdownA31,exp(tLdownA31C),t19)
  Nasemean20 <- sigmoid(exp(K20),rdownA31,exp(tLdownA31D),t20)
  Nasemean21 <- sigmoid(exp(K21),rdownG31,exp(tLdownG31A),t21)
  Nasemean22 <- sigmoid(exp(K22),rdownG31,exp(tLdownG31B),t22)
  Nasemean23 <- sigmoid(exp(K23),rdownG31,exp(tLdownG31C),t23)
  Nasemean24 <- sigmoid(exp(K24),rdownG31,exp(tLdownG31D),t24)
  Nasemean25 <- sigmoid(exp(K25),rdownM31,exp(tLdownM31A),t25)
  Nasemean26 <- sigmoid(exp(K26),rdownM31,exp(tLdownM31B),t26)
  Nasemean27 <- sigmoid(exp(K27),rdownM31,exp(tLdownM31C),t27)
  Nasemean28 <- sigmoid(exp(K28),rdownM31,exp(tLdownM31D),t28)
  Nasemean29 <- sigmoid(exp(K29),rdownR31,exp(tLdownR31A),t29)
  Nasemean30 <- sigmoid(exp(K30),rdownR31,exp(tLdownR31B),t30)
  Nasemean31 <- sigmoid(exp(K31),rdownR31,exp(tLdownR31C),t31)
  Nasemean32 <- sigmoid(exp(K32),rdownR31,exp(tLdownR31D),t32)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_down_sptemp <- mle2(sigmoid_down_sptemp_normNLL,start=list(sdNase=-1,
                                                                       K1=log(max(A21DA$SNF)),K2=log(max(A21DB$SNF)),
                                                                       K3=log(max(A21DC$SNF)),K4=log(max(A21DD$SNF)),
                                                                       K5=log(max(G21DA$SNF)),K6=log(max(G21DB$SNF)),
                                                                       K7=log(max(G21DC$SNF)),K8=log(max(G21DD$SNF)),
                                                                       K9=log(max(M21DA$SNF)),K10=log(max(M21DB$SNF)),
                                                                       K11=log(max(M21DC$SNF)),K12=log(max(M21DD$SNF)),
                                                                       K13=log(max(R21DA$SNF)),K14=log(max(R21DB$SNF)),
                                                                       K15=log(max(R21DC$SNF)),K16=log(max(R21DD$SNF)),
                                                                       K17=log(max(A31DA$SNF)),K18=log(max(A31DB$SNF)),
                                                                       K19=log(max(A31DC$SNF)),K20=log(max(A31DD$SNF)),
                                                                       K21=log(max(G31DA$SNF)),K22=log(max(G31DB$SNF)),
                                                                       K23=log(max(G31DC$SNF)),K24=log(max(G31DD$SNF)),
                                                                       K25=log(max(M31DA$SNF)),K26=log(max(M31DB$SNF)),
                                                                       K27=log(max(M31DC$SNF)),K28=log(max(M31DD$SNF)),
                                                                       K29=log(max(R31DA$SNF)),K30=log(max(R31DB$SNF)),
                                                                       K31=log(max(R31DC$SNF)),K32=log(max(R31DD$SNF)),
                                                                       rdownA21=-0.1,rdownG21=-0.1,rdownM21=-0.1,rdownR21=-0.1,
                                                                       rdownA31=-0.1,rdownG31=-0.1,rdownM31=-0.1,rdownR31=-0.1,
                                                                       tLdownA21A=log(7),tLdownA21B=log(7),tLdownA21C=log(7),tLdownA21D=log(7),
                                                                       tLdownG21A=log(7),tLdownG21B=log(7),tLdownG21C=log(7),tLdownG21D=log(7),
                                                                       tLdownM21A=log(7),tLdownM21B=log(7),tLdownM21C=log(7),tLdownM21D=log(7),
                                                                       tLdownR21A=log(7),tLdownR21B=log(7),tLdownR21C=log(7),tLdownR21D=log(7),
                                                                       tLdownA31A=log(7),tLdownA31B=log(7),tLdownA31C=log(7),tLdownA31D=log(7),
                                                                       tLdownG31A=log(7),tLdownG31B=log(7),tLdownG31C=log(7),tLdownG31D=log(7),
                                                                       tLdownM31A=log(0.2),tLdownM31B=log(7),tLdownM31C=log(7),tLdownM31D=log(7),
                                                                       tLdownR31A=log(7),tLdownR31B=log(7),tLdownR31C=log(7),tLdownR31D=log(7)),
                                data=list(t1=A21DA$Day,t2=A21DB$Day,
                                          t3=A21DC$Day,t4=A21DD$Day,
                                          t5=G21DA$Day,t6=G21DB$Day,
                                          t7=G21DC$Day,t8=G21DD$Day,
                                          t9=M21DA$Day,t10=M21DB$Day,
                                          t11=M21DC$Day,t12=M21DD$Day,
                                          t13=R21DA$Day,t14=R21DB$Day,
                                          t15=R21DC$Day,t16=R21DD$Day,
                                          t17=A31DA$Day,t18=A31DB$Day,
                                          t19=A31DC$Day,t20=A31DD$Day,
                                          t21=G31DA$Day,t22=G31DB$Day,
                                          t23=G31DC$Day,t24=G31DD$Day,
                                          t25=M31DA$Day,t26=M31DB$Day,
                                          t27=M31DC$Day,t28=M31DD$Day,
                                          t29=R31DA$Day,t30=R31DB$Day,
                                          t31=R31DC$Day,t32=R31DD$Day,
                                          Nasedat1=A21DA$SNF,Nasedat2=A21DB$SNF,
                                          Nasedat3=A21DC$SNF,Nasedat4=A21DD$SNF,
                                          Nasedat5=G21DA$SNF,Nasedat6=G21DB$SNF,
                                          Nasedat7=G21DC$SNF,Nasedat8=G21DD$SNF,
                                          Nasedat9=M21DA$SNF,Nasedat10=M21DB$SNF,
                                          Nasedat11=M21DC$SNF,Nasedat12=M21DD$SNF,
                                          Nasedat13=R21DA$SNF,Nasedat14=R21DB$SNF,
                                          Nasedat15=R21DC$SNF,Nasedat16=R21DD$SNF,
                                          Nasedat17=A31DA$SNF,Nasedat18=A31DB$SNF,
                                          Nasedat19=A31DC$SNF,Nasedat20=A31DD$SNF,
                                          Nasedat21=G31DA$SNF,Nasedat22=G31DB$SNF,
                                          Nasedat23=G31DC$SNF,Nasedat24=G31DD$SNF,
                                          Nasedat25=M31DA$SNF,Nasedat26=M31DB$SNF,
                                          Nasedat27=M31DC$SNF,Nasedat28=M31DD$SNF,
                                          Nasedat29=R31DA$SNF,Nasedat30=R31DB$SNF,
                                          Nasedat31=R31DC$SNF,Nasedat32=R31DD$SNF),
                                control=list(maxit=20000))
summary(fit_sigmoid_down_sptemp)

## Calculate AICc to see if tL is a necessary parameter in the model for down-regulation
AICctab(fit_sigmoid_down_sptemp_no_tL,fit_sigmoid_down_sptemp,nobs=272)
#tL is an important parameter (delta AICc 557.4)

## One r for all plants
sigmoid_down_same_normNLL <- function(sdNase,rdown,
                                      tLdownA21A,tLdownA21B,tLdownA21C,tLdownA21D,
                                      tLdownG21A,tLdownG21B,tLdownG21C,tLdownG21D,
                                      tLdownM21A,tLdownM21B,tLdownM21C,tLdownM21D,
                                      tLdownR21A,tLdownR21B,tLdownR21C,tLdownR21D,
                                      tLdownA31A,tLdownA31B,tLdownA31C,tLdownA31D,
                                      tLdownG31A,tLdownG31B,tLdownG31C,tLdownG31D,
                                      tLdownM31A,tLdownM31B,tLdownM31C,tLdownM31D,
                                      tLdownR31A,tLdownR31B,tLdownR31C,tLdownR31D,
                                      K1,K2,K3,K4,K5,K6,K7,K8,
                                      K9,K10,K11,K12,K13,K14,K15,K16,
                                      K17,K18,K19,K20,K21,K22,K23,K24,
                                      K25,K26,K27,K28,K29,K30,K31,K32,
                                      t1,t2,t3,t4,t5,t6,t7,t8,
                                      t9,t10,t11,t12,t13,t14,t15,t16,
                                      t17,t18,t19,t20,t21,t22,t23,t24,
                                      t25,t26,t27,t28,t29,t30,t31,t32,
                                      Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                      Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                      Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                      Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32){
  Nasemean1 <- sigmoid(exp(K1),rdown,exp(tLdownA21A),t1)
  Nasemean2 <- sigmoid(exp(K2),rdown,exp(tLdownA21B),t2)
  Nasemean3 <- sigmoid(exp(K3),rdown,exp(tLdownA21C),t3)
  Nasemean4 <- sigmoid(exp(K4),rdown,exp(tLdownA21D),t4)
  Nasemean5 <- sigmoid(exp(K5),rdown,exp(tLdownG21A),t5)
  Nasemean6 <- sigmoid(exp(K6),rdown,exp(tLdownG21B),t6)
  Nasemean7 <- sigmoid(exp(K7),rdown,exp(tLdownG21C),t7)
  Nasemean8 <- sigmoid(exp(K8),rdown,exp(tLdownG21D),t8)
  Nasemean9 <- sigmoid(exp(K9),rdown,exp(tLdownM21A),t9)
  Nasemean10 <- sigmoid(exp(K10),rdown,exp(tLdownM21B),t10)
  Nasemean11 <- sigmoid(exp(K11),rdown,exp(tLdownM21C),t11)
  Nasemean12 <- sigmoid(exp(K12),rdown,exp(tLdownM21D),t12)
  Nasemean13 <- sigmoid(exp(K13),rdown,exp(tLdownR21A),t13)
  Nasemean14 <- sigmoid(exp(K14),rdown,exp(tLdownR21B),t14)
  Nasemean15 <- sigmoid(exp(K15),rdown,exp(tLdownR21C),t15)
  Nasemean16 <- sigmoid(exp(K16),rdown,exp(tLdownR21D),t16)
  Nasemean17 <- sigmoid(exp(K17),rdown,exp(tLdownA31A),t17)
  Nasemean18 <- sigmoid(exp(K18),rdown,exp(tLdownA31B),t18)
  Nasemean19 <- sigmoid(exp(K19),rdown,exp(tLdownA31C),t19)
  Nasemean20 <- sigmoid(exp(K20),rdown,exp(tLdownA31D),t20)
  Nasemean21 <- sigmoid(exp(K21),rdown,exp(tLdownG31A),t21)
  Nasemean22 <- sigmoid(exp(K22),rdown,exp(tLdownG31B),t22)
  Nasemean23 <- sigmoid(exp(K23),rdown,exp(tLdownG31C),t23)
  Nasemean24 <- sigmoid(exp(K24),rdown,exp(tLdownG31D),t24)
  Nasemean25 <- sigmoid(exp(K25),rdown,exp(tLdownM31A),t25)
  Nasemean26 <- sigmoid(exp(K26),rdown,exp(tLdownM31B),t26)
  Nasemean27 <- sigmoid(exp(K27),rdown,exp(tLdownM31C),t27)
  Nasemean28 <- sigmoid(exp(K28),rdown,exp(tLdownM31D),t28)
  Nasemean29 <- sigmoid(exp(K29),rdown,exp(tLdownR31A),t29)
  Nasemean30 <- sigmoid(exp(K30),rdown,exp(tLdownR31B),t30)
  Nasemean31 <- sigmoid(exp(K31),rdown,exp(tLdownR31C),t31)
  Nasemean32 <- sigmoid(exp(K32),rdown,exp(tLdownR31D),t32)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_down_same <- mle2(sigmoid_down_same_normNLL,start=list(sdNase=-1,
                                                                   K1=log(0.014),K2=log(0.01),
                                                                   K3=log(0.014),K4=log(0.017),
                                                                   K5=log(0.007),K6=log(0.001),
                                                                   K7=log(0.007),K8=log(0.008),
                                                                   K9=log(0.053),K10=log(0.029),
                                                                   K11=log(0.05),K12=log(0.078),
                                                                   K13=log(0.055),K14=log(0.038),
                                                                   K15=log(0.029),K16=log(0.032),
                                                                   K17=log(0.029),K18=log(0.024),
                                                                   K19=log(0.021),K20=log(0.012),
                                                                   K21=log(0.057),K22=log(0.033),
                                                                   K23=log(0.061),K24=log(0.016),
                                                                   K25=log(0.033),K26=log(0.061),
                                                                   K27=log(0.022),K28=log(0.036),
                                                                   K29=log(0.043),K30=log(0.017),
                                                                   K31=log(0.027),K32=log(0.021),
                                                                   rdown=-0.12,
                                                                   tLdownA21A=log(12),tLdownA21B=log(0.4),tLdownA21C=log(0.001),tLdownA21D=log(0.001),
                                                                   tLdownG21A=log(32),tLdownG21B=log(6),tLdownG21C=log(6),tLdownG21D=log(7),
                                                                   tLdownM21A=log(0.001),tLdownM21B=log(16),tLdownM21C=log(14),tLdownM21D=log(0.6),
                                                                   tLdownR21A=log(6),tLdownR21B=log(43),tLdownR21C=log(40),tLdownR21D=log(12),
                                                                   tLdownA31A=log(0.01),tLdownA31B=log(4),tLdownA31C=log(6),tLdownA31D=log(10),
                                                                   tLdownG31A=log(0.1),tLdownG31B=log(0.001),tLdownG31C=log(0.001),tLdownG31D=log(20),
                                                                   tLdownM31A=log(9),tLdownM31B=log(0.002),tLdownM31C=log(0.001),tLdownM31D=log(6),
                                                                   tLdownR31A=log(0.002),tLdownR31B=log(20),tLdownR31C=log(0.001),tLdownR31D=log(13)),
                              data=list(t1=A21DA$Day,t2=A21DB$Day,
                                        t3=A21DC$Day,t4=A21DD$Day,
                                        t5=G21DA$Day,t6=G21DB$Day,
                                        t7=G21DC$Day,t8=G21DD$Day,
                                        t9=M21DA$Day,t10=M21DB$Day,
                                        t11=M21DC$Day,t12=M21DD$Day,
                                        t13=R21DA$Day,t14=R21DB$Day,
                                        t15=R21DC$Day,t16=R21DD$Day,
                                        t17=A31DA$Day,t18=A31DB$Day,
                                        t19=A31DC$Day,t20=A31DD$Day,
                                        t21=G31DA$Day,t22=G31DB$Day,
                                        t23=G31DC$Day,t24=G31DD$Day,
                                        t25=M31DA$Day,t26=M31DB$Day,
                                        t27=M31DC$Day,t28=M31DD$Day,
                                        t29=R31DA$Day,t30=R31DB$Day,
                                        t31=R31DC$Day,t32=R31DD$Day,
                                        Nasedat1=A21DA$SNF,Nasedat2=A21DB$SNF,
                                        Nasedat3=A21DC$SNF,Nasedat4=A21DD$SNF,
                                        Nasedat5=G21DA$SNF,Nasedat6=G21DB$SNF,
                                        Nasedat7=G21DC$SNF,Nasedat8=G21DD$SNF,
                                        Nasedat9=M21DA$SNF,Nasedat10=M21DB$SNF,
                                        Nasedat11=M21DC$SNF,Nasedat12=M21DD$SNF,
                                        Nasedat13=R21DA$SNF,Nasedat14=R21DB$SNF,
                                        Nasedat15=R21DC$SNF,Nasedat16=R21DD$SNF,
                                        Nasedat17=A31DA$SNF,Nasedat18=A31DB$SNF,
                                        Nasedat19=A31DC$SNF,Nasedat20=A31DD$SNF,
                                        Nasedat21=G31DA$SNF,Nasedat22=G31DB$SNF,
                                        Nasedat23=G31DC$SNF,Nasedat24=G31DD$SNF,
                                        Nasedat25=M31DA$SNF,Nasedat26=M31DB$SNF,
                                        Nasedat27=M31DC$SNF,Nasedat28=M31DD$SNF,
                                        Nasedat29=R31DA$SNF,Nasedat30=R31DB$SNF,
                                        Nasedat31=R31DC$SNF,Nasedat32=R31DD$SNF),
                              control=list(maxit=20000))
summary(fit_sigmoid_down_same)

## Symbiosis
sigmoid_down_sym_normNLL <- function(sdNase,rdownRhiz,rdownActin,
                                     tLdownA21A,tLdownA21B,tLdownA21C,tLdownA21D,
                                     tLdownG21A,tLdownG21B,tLdownG21C,tLdownG21D,
                                     tLdownM21A,tLdownM21B,tLdownM21C,tLdownM21D,
                                     tLdownR21A,tLdownR21B,tLdownR21C,tLdownR21D,
                                     tLdownA31A,tLdownA31B,tLdownA31C,tLdownA31D,
                                     tLdownG31A,tLdownG31B,tLdownG31C,tLdownG31D,
                                     tLdownM31A,tLdownM31B,tLdownM31C,tLdownM31D,
                                     tLdownR31A,tLdownR31B,tLdownR31C,tLdownR31D,
                                     K1,K2,K3,K4,K5,K6,K7,K8,
                                     K9,K10,K11,K12,K13,K14,K15,K16,
                                     K17,K18,K19,K20,K21,K22,K23,K24,
                                     K25,K26,K27,K28,K29,K30,K31,K32,
                                     t1,t2,t3,t4,t5,t6,t7,t8,
                                     t9,t10,t11,t12,t13,t14,t15,t16,
                                     t17,t18,t19,t20,t21,t22,t23,t24,
                                     t25,t26,t27,t28,t29,t30,t31,t32,
                                     Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                     Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                     Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                     Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32){
  Nasemean1 <- sigmoid(exp(K1),rdownActin,exp(tLdownA21A),t1)
  Nasemean2 <- sigmoid(exp(K2),rdownActin,exp(tLdownA21B),t2)
  Nasemean3 <- sigmoid(exp(K3),rdownActin,exp(tLdownA21C),t3)
  Nasemean4 <- sigmoid(exp(K4),rdownActin,exp(tLdownA21D),t4)
  Nasemean5 <- sigmoid(exp(K5),rdownRhiz,exp(tLdownG21A),t5)
  Nasemean6 <- sigmoid(exp(K6),rdownRhiz,exp(tLdownG21B),t6)
  Nasemean7 <- sigmoid(exp(K7),rdownRhiz,exp(tLdownG21C),t7)
  Nasemean8 <- sigmoid(exp(K8),rdownRhiz,exp(tLdownG21D),t8)
  Nasemean9 <- sigmoid(exp(K9),rdownActin,exp(tLdownM21A),t9)
  Nasemean10 <- sigmoid(exp(K10),rdownActin,exp(tLdownM21B),t10)
  Nasemean11 <- sigmoid(exp(K11),rdownActin,exp(tLdownM21C),t11)
  Nasemean12 <- sigmoid(exp(K12),rdownActin,exp(tLdownM21D),t12)
  Nasemean13 <- sigmoid(exp(K13),rdownRhiz,exp(tLdownR21A),t13)
  Nasemean14 <- sigmoid(exp(K14),rdownRhiz,exp(tLdownR21B),t14)
  Nasemean15 <- sigmoid(exp(K15),rdownRhiz,exp(tLdownR21C),t15)
  Nasemean16 <- sigmoid(exp(K16),rdownRhiz,exp(tLdownR21D),t16)
  Nasemean17 <- sigmoid(exp(K17),rdownActin,exp(tLdownA31A),t17)
  Nasemean18 <- sigmoid(exp(K18),rdownActin,exp(tLdownA31B),t18)
  Nasemean19 <- sigmoid(exp(K19),rdownActin,exp(tLdownA31C),t19)
  Nasemean20 <- sigmoid(exp(K20),rdownActin,exp(tLdownA31D),t20)
  Nasemean21 <- sigmoid(exp(K21),rdownRhiz,exp(tLdownG31A),t21)
  Nasemean22 <- sigmoid(exp(K22),rdownRhiz,exp(tLdownG31B),t22)
  Nasemean23 <- sigmoid(exp(K23),rdownRhiz,exp(tLdownG31C),t23)
  Nasemean24 <- sigmoid(exp(K24),rdownRhiz,exp(tLdownG31D),t24)
  Nasemean25 <- sigmoid(exp(K25),rdownActin,exp(tLdownM31A),t25)
  Nasemean26 <- sigmoid(exp(K26),rdownActin,exp(tLdownM31B),t26)
  Nasemean27 <- sigmoid(exp(K27),rdownActin,exp(tLdownM31C),t27)
  Nasemean28 <- sigmoid(exp(K28),rdownActin,exp(tLdownM31D),t28)
  Nasemean29 <- sigmoid(exp(K29),rdownRhiz,exp(tLdownR31A),t29)
  Nasemean30 <- sigmoid(exp(K30),rdownRhiz,exp(tLdownR31B),t30)
  Nasemean31 <- sigmoid(exp(K31),rdownRhiz,exp(tLdownR31C),t31)
  Nasemean32 <- sigmoid(exp(K32),rdownRhiz,exp(tLdownR31D),t32)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_down_sym <- mle2(sigmoid_down_sym_normNLL,start=list(sdNase=-1,
                                                                 K1=log(max(A21DA$SNF)),K2=log(max(A21DB$SNF)),
                                                                 K3=log(max(A21DC$SNF)),K4=log(max(A21DD$SNF)),
                                                                 K5=log(max(G21DA$SNF)),K6=log(max(G21DB$SNF)),
                                                                 K7=log(max(G21DC$SNF)),K8=log(max(G21DD$SNF)),
                                                                 K9=log(max(M21DA$SNF)),K10=log(max(M21DB$SNF)),
                                                                 K11=log(max(M21DC$SNF)),K12=log(max(M21DD$SNF)),
                                                                 K13=log(max(R21DA$SNF)),K14=log(max(R21DB$SNF)),
                                                                 K15=log(max(R21DC$SNF)),K16=log(max(R21DD$SNF)),
                                                                 K17=log(max(A31DA$SNF)),K18=log(max(A31DB$SNF)),
                                                                 K19=log(max(A31DC$SNF)),K20=log(max(A31DD$SNF)),
                                                                 K21=log(max(G31DA$SNF)),K22=log(max(G31DB$SNF)),
                                                                 K23=log(max(G31DC$SNF)),K24=log(max(G31DD$SNF)),
                                                                 K25=log(max(M31DA$SNF)),K26=log(max(M31DB$SNF)),
                                                                 K27=log(max(M31DC$SNF)),K28=log(max(M31DD$SNF)),
                                                                 K29=log(max(R31DA$SNF)),K30=log(max(R31DB$SNF)),
                                                                 K31=log(max(R31DC$SNF)),K32=log(max(R31DD$SNF)),
                                                                 rdownRhiz=-0.13,rdownActin=-0.1,
                                                                 tLdownA21A=log(12),tLdownA21B=log(0.4),tLdownA21C=log(0.001),tLdownA21D=log(0.001),
                                                                 tLdownG21A=log(32),tLdownG21B=log(6),tLdownG21C=log(6),tLdownG21D=log(7),
                                                                 tLdownM21A=log(0.001),tLdownM21B=log(16),tLdownM21C=log(14),tLdownM21D=log(0.6),
                                                                 tLdownR21A=log(6),tLdownR21B=log(43),tLdownR21C=log(40),tLdownR21D=log(12),
                                                                 tLdownA31A=log(0.01),tLdownA31B=log(4),tLdownA31C=log(6),tLdownA31D=log(10),
                                                                 tLdownG31A=log(0.1),tLdownG31B=log(0.001),tLdownG31C=log(0.001),tLdownG31D=log(20),
                                                                 tLdownM31A=log(9),tLdownM31B=log(0.002),tLdownM31C=log(0.001),tLdownM31D=log(6),
                                                                 tLdownR31A=log(0.002),tLdownR31B=log(20),tLdownR31C=log(0.001),tLdownR31D=log(13)),
                             data=list(t1=A21DA$Day,t2=A21DB$Day,
                                       t3=A21DC$Day,t4=A21DD$Day,
                                       t5=G21DA$Day,t6=G21DB$Day,
                                       t7=G21DC$Day,t8=G21DD$Day,
                                       t9=M21DA$Day,t10=M21DB$Day,
                                       t11=M21DC$Day,t12=M21DD$Day,
                                       t13=R21DA$Day,t14=R21DB$Day,
                                       t15=R21DC$Day,t16=R21DD$Day,
                                       t17=A31DA$Day,t18=A31DB$Day,
                                       t19=A31DC$Day,t20=A31DD$Day,
                                       t21=G31DA$Day,t22=G31DB$Day,
                                       t23=G31DC$Day,t24=G31DD$Day,
                                       t25=M31DA$Day,t26=M31DB$Day,
                                       t27=M31DC$Day,t28=M31DD$Day,
                                       t29=R31DA$Day,t30=R31DB$Day,
                                       t31=R31DC$Day,t32=R31DD$Day,
                                       Nasedat1=A21DA$SNF,Nasedat2=A21DB$SNF,
                                       Nasedat3=A21DC$SNF,Nasedat4=A21DD$SNF,
                                       Nasedat5=G21DA$SNF,Nasedat6=G21DB$SNF,
                                       Nasedat7=G21DC$SNF,Nasedat8=G21DD$SNF,
                                       Nasedat9=M21DA$SNF,Nasedat10=M21DB$SNF,
                                       Nasedat11=M21DC$SNF,Nasedat12=M21DD$SNF,
                                       Nasedat13=R21DA$SNF,Nasedat14=R21DB$SNF,
                                       Nasedat15=R21DC$SNF,Nasedat16=R21DD$SNF,
                                       Nasedat17=A31DA$SNF,Nasedat18=A31DB$SNF,
                                       Nasedat19=A31DC$SNF,Nasedat20=A31DD$SNF,
                                       Nasedat21=G31DA$SNF,Nasedat22=G31DB$SNF,
                                       Nasedat23=G31DC$SNF,Nasedat24=G31DD$SNF,
                                       Nasedat25=M31DA$SNF,Nasedat26=M31DB$SNF,
                                       Nasedat27=M31DC$SNF,Nasedat28=M31DD$SNF,
                                       Nasedat29=R31DA$SNF,Nasedat30=R31DB$SNF,
                                       Nasedat31=R31DC$SNF,Nasedat32=R31DD$SNF),
                             control=list(maxit=20000))
summary(fit_sigmoid_down_sym)

## Symbiosis x temperature interaction
sigmoid_down_symtemp_normNLL <- function(sdNase,rdownRhiz21,rdownActin21,rdownRhiz31,rdownActin31,
                                         tLdownA21A,tLdownA21B,tLdownA21C,tLdownA21D,
                                         tLdownG21A,tLdownG21B,tLdownG21C,tLdownG21D,
                                         tLdownM21A,tLdownM21B,tLdownM21C,tLdownM21D,
                                         tLdownR21A,tLdownR21B,tLdownR21C,tLdownR21D,
                                         tLdownA31A,tLdownA31B,tLdownA31C,tLdownA31D,
                                         tLdownG31A,tLdownG31B,tLdownG31C,tLdownG31D,
                                         tLdownM31A,tLdownM31B,tLdownM31C,tLdownM31D,
                                         tLdownR31A,tLdownR31B,tLdownR31C,tLdownR31D,
                                         K1,K2,K3,K4,K5,K6,K7,K8,
                                         K9,K10,K11,K12,K13,K14,K15,K16,
                                         K17,K18,K19,K20,K21,K22,K23,K24,
                                         K25,K26,K27,K28,K29,K30,K31,K32,
                                         t1,t2,t3,t4,t5,t6,t7,t8,
                                         t9,t10,t11,t12,t13,t14,t15,t16,
                                         t17,t18,t19,t20,t21,t22,t23,t24,
                                         t25,t26,t27,t28,t29,t30,t31,t32,
                                         Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                         Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                         Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                         Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32){
  Nasemean1 <- sigmoid(exp(K1),rdownActin21,exp(tLdownA21A),t1)
  Nasemean2 <- sigmoid(exp(K2),rdownActin21,exp(tLdownA21B),t2)
  Nasemean3 <- sigmoid(exp(K3),rdownActin21,exp(tLdownA21C),t3)
  Nasemean4 <- sigmoid(exp(K4),rdownActin21,exp(tLdownA21D),t4)
  Nasemean5 <- sigmoid(exp(K5),rdownRhiz21,exp(tLdownG21A),t5)
  Nasemean6 <- sigmoid(exp(K6),rdownRhiz21,exp(tLdownG21B),t6)
  Nasemean7 <- sigmoid(exp(K7),rdownRhiz21,exp(tLdownG21C),t7)
  Nasemean8 <- sigmoid(exp(K8),rdownRhiz21,exp(tLdownG21D),t8)
  Nasemean9 <- sigmoid(exp(K9),rdownActin21,exp(tLdownM21A),t9)
  Nasemean10 <- sigmoid(exp(K10),rdownActin21,exp(tLdownM21B),t10)
  Nasemean11 <- sigmoid(exp(K11),rdownActin21,exp(tLdownM21C),t11)
  Nasemean12 <- sigmoid(exp(K12),rdownActin21,exp(tLdownM21D),t12)
  Nasemean13 <- sigmoid(exp(K13),rdownRhiz21,exp(tLdownR21A),t13)
  Nasemean14 <- sigmoid(exp(K14),rdownRhiz21,exp(tLdownR21B),t14)
  Nasemean15 <- sigmoid(exp(K15),rdownRhiz21,exp(tLdownR21C),t15)
  Nasemean16 <- sigmoid(exp(K16),rdownRhiz21,exp(tLdownR21D),t16)
  Nasemean17 <- sigmoid(exp(K17),rdownActin31,exp(tLdownA31A),t17)
  Nasemean18 <- sigmoid(exp(K18),rdownActin31,exp(tLdownA31B),t18)
  Nasemean19 <- sigmoid(exp(K19),rdownActin31,exp(tLdownA31C),t19)
  Nasemean20 <- sigmoid(exp(K20),rdownActin31,exp(tLdownA31D),t20)
  Nasemean21 <- sigmoid(exp(K21),rdownRhiz31,exp(tLdownG31A),t21)
  Nasemean22 <- sigmoid(exp(K22),rdownRhiz31,exp(tLdownG31B),t22)
  Nasemean23 <- sigmoid(exp(K23),rdownRhiz31,exp(tLdownG31C),t23)
  Nasemean24 <- sigmoid(exp(K24),rdownRhiz31,exp(tLdownG31D),t24)
  Nasemean25 <- sigmoid(exp(K25),rdownActin31,exp(tLdownM31A),t25)
  Nasemean26 <- sigmoid(exp(K26),rdownActin31,exp(tLdownM31B),t26)
  Nasemean27 <- sigmoid(exp(K27),rdownActin31,exp(tLdownM31C),t27)
  Nasemean28 <- sigmoid(exp(K28),rdownActin31,exp(tLdownM31D),t28)
  Nasemean29 <- sigmoid(exp(K29),rdownRhiz31,exp(tLdownR31A),t29)
  Nasemean30 <- sigmoid(exp(K30),rdownRhiz31,exp(tLdownR31B),t30)
  Nasemean31 <- sigmoid(exp(K31),rdownRhiz31,exp(tLdownR31C),t31)
  Nasemean32 <- sigmoid(exp(K32),rdownRhiz31,exp(tLdownR31D),t32)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_down_symtemp <- mle2(sigmoid_down_symtemp_normNLL,start=list(sdNase=-1,
                                                                         K1=log(0.014),K2=log(0.01),
                                                                         K3=log(0.014),K4=log(0.017),
                                                                         K5=log(0.007),K6=log(0.001),
                                                                         K7=log(0.007),K8=log(0.008),
                                                                         K9=log(0.053),K10=log(0.029),
                                                                         K11=log(0.05),K12=log(0.078),
                                                                         K13=log(0.055),K14=log(0.038),
                                                                         K15=log(0.029),K16=log(0.032),
                                                                         K17=log(0.029),K18=log(0.024),
                                                                         K19=log(0.021),K20=log(0.012),
                                                                         K21=log(0.057),K22=log(0.033),
                                                                         K23=log(0.061),K24=log(0.016),
                                                                         K25=log(0.033),K26=log(0.061),
                                                                         K27=log(0.022),K28=log(0.036),
                                                                         K29=log(0.043),K30=log(0.017),
                                                                         K31=log(0.027),K32=log(0.021),
                                                                         rdownActin21=-0.1,
                                                                         rdownRhiz21=-0.1,
                                                                         rdownActin31=-0.1,
                                                                         rdownRhiz31=-0.1,
                                                                         tLdownA21A=log(12),tLdownA21B=log(0.4),tLdownA21C=log(0.001),tLdownA21D=log(0.001),
                                                                         tLdownG21A=log(32),tLdownG21B=log(6),tLdownG21C=log(6),tLdownG21D=log(7),
                                                                         tLdownM21A=log(0.001),tLdownM21B=log(16),tLdownM21C=log(14),tLdownM21D=log(0.6),
                                                                         tLdownR21A=log(6),tLdownR21B=log(43),tLdownR21C=log(40),tLdownR21D=log(12),
                                                                         tLdownA31A=log(0.01),tLdownA31B=log(4),tLdownA31C=log(6),tLdownA31D=log(10),
                                                                         tLdownG31A=log(0.1),tLdownG31B=log(0.001),tLdownG31C=log(0.001),tLdownG31D=log(20),
                                                                         tLdownM31A=log(9),tLdownM31B=log(0.002),tLdownM31C=log(0.001),tLdownM31D=log(6),
                                                                         tLdownR31A=log(0.002),tLdownR31B=log(20),tLdownR31C=log(0.001),tLdownR31D=log(13)),
                                 data=list(t1=A21DA$Day,t2=A21DB$Day,
                                           t3=A21DC$Day,t4=A21DD$Day,
                                           t5=G21DA$Day,t6=G21DB$Day,
                                           t7=G21DC$Day,t8=G21DD$Day,
                                           t9=M21DA$Day,t10=M21DB$Day,
                                           t11=M21DC$Day,t12=M21DD$Day,
                                           t13=R21DA$Day,t14=R21DB$Day,
                                           t15=R21DC$Day,t16=R21DD$Day,
                                           t17=A31DA$Day,t18=A31DB$Day,
                                           t19=A31DC$Day,t20=A31DD$Day,
                                           t21=G31DA$Day,t22=G31DB$Day,
                                           t23=G31DC$Day,t24=G31DD$Day,
                                           t25=M31DA$Day,t26=M31DB$Day,
                                           t27=M31DC$Day,t28=M31DD$Day,
                                           t29=R31DA$Day,t30=R31DB$Day,
                                           t31=R31DC$Day,t32=R31DD$Day,
                                           Nasedat1=A21DA$SNF,Nasedat2=A21DB$SNF,
                                           Nasedat3=A21DC$SNF,Nasedat4=A21DD$SNF,
                                           Nasedat5=G21DA$SNF,Nasedat6=G21DB$SNF,
                                           Nasedat7=G21DC$SNF,Nasedat8=G21DD$SNF,
                                           Nasedat9=M21DA$SNF,Nasedat10=M21DB$SNF,
                                           Nasedat11=M21DC$SNF,Nasedat12=M21DD$SNF,
                                           Nasedat13=R21DA$SNF,Nasedat14=R21DB$SNF,
                                           Nasedat15=R21DC$SNF,Nasedat16=R21DD$SNF,
                                           Nasedat17=A31DA$SNF,Nasedat18=A31DB$SNF,
                                           Nasedat19=A31DC$SNF,Nasedat20=A31DD$SNF,
                                           Nasedat21=G31DA$SNF,Nasedat22=G31DB$SNF,
                                           Nasedat23=G31DC$SNF,Nasedat24=G31DD$SNF,
                                           Nasedat25=M31DA$SNF,Nasedat26=M31DB$SNF,
                                           Nasedat27=M31DC$SNF,Nasedat28=M31DD$SNF,
                                           Nasedat29=R31DA$SNF,Nasedat30=R31DB$SNF,
                                           Nasedat31=R31DC$SNF,Nasedat32=R31DD$SNF),
                                 control=list(maxit=20000))
summary(fit_sigmoid_down_symtemp)

## Species
sigmoid_down_sp_normNLL <- function(sdNase,
                                    rdownA,
                                    rdownG,
                                    rdownM,
                                    rdownR,
                                    tLdownA21A,tLdownA21B,tLdownA21C,tLdownA21D,
                                    tLdownG21A,tLdownG21B,tLdownG21C,tLdownG21D,
                                    tLdownM21A,tLdownM21B,tLdownM21C,tLdownM21D,
                                    tLdownR21A,tLdownR21B,tLdownR21C,tLdownR21D,
                                    tLdownA31A,tLdownA31B,tLdownA31C,tLdownA31D,
                                    tLdownG31A,tLdownG31B,tLdownG31C,tLdownG31D,
                                    tLdownM31A,tLdownM31B,tLdownM31C,tLdownM31D,
                                    tLdownR31A,tLdownR31B,tLdownR31C,tLdownR31D,
                                    K1,K2,K3,K4,K5,K6,K7,K8,
                                    K9,K10,K11,K12,K13,K14,K15,K16,
                                    K17,K18,K19,K20,K21,K22,K23,K24,
                                    K25,K26,K27,K28,K29,K30,K31,K32,
                                    t1,t2,t3,t4,t5,t6,t7,t8,
                                    t9,t10,t11,t12,t13,t14,t15,t16,
                                    t17,t18,t19,t20,t21,t22,t23,t24,
                                    t25,t26,t27,t28,t29,t30,t31,t32,
                                    Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                    Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                    Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                    Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32){
  Nasemean1 <- sigmoid(exp(K1),rdownA,exp(tLdownA21A),t1)
  Nasemean2 <- sigmoid(exp(K2),rdownA,exp(tLdownA21B),t2)
  Nasemean3 <- sigmoid(exp(K3),rdownA,exp(tLdownA21C),t3)
  Nasemean4 <- sigmoid(exp(K4),rdownA,exp(tLdownA21D),t4)
  Nasemean5 <- sigmoid(exp(K5),rdownG,exp(tLdownG21A),t5)
  Nasemean6 <- sigmoid(exp(K6),rdownG,exp(tLdownG21B),t6)
  Nasemean7 <- sigmoid(exp(K7),rdownG,exp(tLdownG21C),t7)
  Nasemean8 <- sigmoid(exp(K8),rdownG,exp(tLdownG21D),t8)
  Nasemean9 <- sigmoid(exp(K9),rdownM,exp(tLdownM21A),t9)
  Nasemean10 <- sigmoid(exp(K10),rdownM,exp(tLdownM21B),t10)
  Nasemean11 <- sigmoid(exp(K11),rdownM,exp(tLdownM21C),t11)
  Nasemean12 <- sigmoid(exp(K12),rdownM,exp(tLdownM21D),t12)
  Nasemean13 <- sigmoid(exp(K13),rdownR,exp(tLdownR21A),t13)
  Nasemean14 <- sigmoid(exp(K14),rdownR,exp(tLdownR21B),t14)
  Nasemean15 <- sigmoid(exp(K15),rdownR,exp(tLdownR21C),t15)
  Nasemean16 <- sigmoid(exp(K16),rdownR,exp(tLdownR21D),t16)
  Nasemean17 <- sigmoid(exp(K17),rdownA,exp(tLdownA31A),t17)
  Nasemean18 <- sigmoid(exp(K18),rdownA,exp(tLdownA31B),t18)
  Nasemean19 <- sigmoid(exp(K19),rdownA,exp(tLdownA31C),t19)
  Nasemean20 <- sigmoid(exp(K20),rdownA,exp(tLdownA31D),t20)
  Nasemean21 <- sigmoid(exp(K21),rdownG,exp(tLdownG31A),t21)
  Nasemean22 <- sigmoid(exp(K22),rdownG,exp(tLdownG31B),t22)
  Nasemean23 <- sigmoid(exp(K23),rdownG,exp(tLdownG31C),t23)
  Nasemean24 <- sigmoid(exp(K24),rdownG,exp(tLdownG31D),t24)
  Nasemean25 <- sigmoid(exp(K25),rdownM,exp(tLdownM31A),t25)
  Nasemean26 <- sigmoid(exp(K26),rdownM,exp(tLdownM31B),t26)
  Nasemean27 <- sigmoid(exp(K27),rdownM,exp(tLdownM31C),t27)
  Nasemean28 <- sigmoid(exp(K28),rdownM,exp(tLdownM31D),t28)
  Nasemean29 <- sigmoid(exp(K29),rdownR,exp(tLdownR31A),t29)
  Nasemean30 <- sigmoid(exp(K30),rdownR,exp(tLdownR31B),t30)
  Nasemean31 <- sigmoid(exp(K31),rdownR,exp(tLdownR31C),t31)
  Nasemean32 <- sigmoid(exp(K32),rdownR,exp(tLdownR31D),t32)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_down_sp <- mle2(sigmoid_down_sp_normNLL,start=list(sdNase=-1,
                                                               K1=log(0.014),K2=log(0.01),
                                                               K3=log(0.014),K4=log(0.017),
                                                               K5=log(0.007),K6=log(0.001),
                                                               K7=log(0.007),K8=log(0.008),
                                                               K9=log(0.053),K10=log(0.029),
                                                               K11=log(0.05),K12=log(0.078),
                                                               K13=log(0.055),K14=log(0.038),
                                                               K15=log(0.029),K16=log(0.032),
                                                               K17=log(0.029),K18=log(0.024),
                                                               K19=log(0.021),K20=log(0.012),
                                                               K21=log(0.057),K22=log(0.033),
                                                               K23=log(0.061),K24=log(0.016),
                                                               K25=log(0.033),K26=log(0.061),
                                                               K27=log(0.022),K28=log(0.036),
                                                               K29=log(0.043),K30=log(0.017),
                                                               K31=log(0.027),K32=log(0.021),
                                                               rdownA=-0.1,
                                                               rdownG=-0.1,
                                                               rdownM=-0.1,
                                                               rdownR=-0.1,
                                                               tLdownA21A=log(12),tLdownA21B=log(0.4),tLdownA21C=log(0.001),tLdownA21D=log(0.001),
                                                               tLdownG21A=log(32),tLdownG21B=log(6),tLdownG21C=log(6),tLdownG21D=log(7),
                                                               tLdownM21A=log(0.001),tLdownM21B=log(16),tLdownM21C=log(14),tLdownM21D=log(0.6),
                                                               tLdownR21A=log(6),tLdownR21B=log(43),tLdownR21C=log(40),tLdownR21D=log(12),
                                                               tLdownA31A=log(0.01),tLdownA31B=log(4),tLdownA31C=log(6),tLdownA31D=log(10),
                                                               tLdownG31A=log(0.1),tLdownG31B=log(0.001),tLdownG31C=log(0.001),tLdownG31D=log(20),
                                                               tLdownM31A=log(9),tLdownM31B=log(0.002),tLdownM31C=log(0.001),tLdownM31D=log(6),
                                                               tLdownR31A=log(0.002),tLdownR31B=log(20),tLdownR31C=log(0.001),tLdownR31D=log(13)),
                            data=list(t1=A21DA$Day,t2=A21DB$Day,
                                      t3=A21DC$Day,t4=A21DD$Day,
                                      t5=G21DA$Day,t6=G21DB$Day,
                                      t7=G21DC$Day,t8=G21DD$Day,
                                      t9=M21DA$Day,t10=M21DB$Day,
                                      t11=M21DC$Day,t12=M21DD$Day,
                                      t13=R21DA$Day,t14=R21DB$Day,
                                      t15=R21DC$Day,t16=R21DD$Day,
                                      t17=A31DA$Day,t18=A31DB$Day,
                                      t19=A31DC$Day,t20=A31DD$Day,
                                      t21=G31DA$Day,t22=G31DB$Day,
                                      t23=G31DC$Day,t24=G31DD$Day,
                                      t25=M31DA$Day,t26=M31DB$Day,
                                      t27=M31DC$Day,t28=M31DD$Day,
                                      t29=R31DA$Day,t30=R31DB$Day,
                                      t31=R31DC$Day,t32=R31DD$Day,
                                      Nasedat1=A21DA$SNF,Nasedat2=A21DB$SNF,
                                      Nasedat3=A21DC$SNF,Nasedat4=A21DD$SNF,
                                      Nasedat5=G21DA$SNF,Nasedat6=G21DB$SNF,
                                      Nasedat7=G21DC$SNF,Nasedat8=G21DD$SNF,
                                      Nasedat9=M21DA$SNF,Nasedat10=M21DB$SNF,
                                      Nasedat11=M21DC$SNF,Nasedat12=M21DD$SNF,
                                      Nasedat13=R21DA$SNF,Nasedat14=R21DB$SNF,
                                      Nasedat15=R21DC$SNF,Nasedat16=R21DD$SNF,
                                      Nasedat17=A31DA$SNF,Nasedat18=A31DB$SNF,
                                      Nasedat19=A31DC$SNF,Nasedat20=A31DD$SNF,
                                      Nasedat21=G31DA$SNF,Nasedat22=G31DB$SNF,
                                      Nasedat23=G31DC$SNF,Nasedat24=G31DD$SNF,
                                      Nasedat25=M31DA$SNF,Nasedat26=M31DB$SNF,
                                      Nasedat27=M31DC$SNF,Nasedat28=M31DD$SNF,
                                      Nasedat29=R31DA$SNF,Nasedat30=R31DB$SNF,
                                      Nasedat31=R31DC$SNF,Nasedat32=R31DD$SNF),
                            control=list(maxit=20000))
summary(fit_sigmoid_down_sp)

## Temperature
sigmoid_down_temp_normNLL <- function(sdNase,rdown1,rdown31,
                                      tLdownA21A,tLdownA21B,tLdownA21C,tLdownA21D,
                                      tLdownG21A,tLdownG21B,tLdownG21C,tLdownG21D,
                                      tLdownM21A,tLdownM21B,tLdownM21C,tLdownM21D,
                                      tLdownR21A,tLdownR21B,tLdownR21C,tLdownR21D,
                                      tLdownA31A,tLdownA31B,tLdownA31C,tLdownA31D,
                                      tLdownG31A,tLdownG31B,tLdownG31C,tLdownG31D,
                                      tLdownM31A,tLdownM31B,tLdownM31C,tLdownM31D,
                                      tLdownR31A,tLdownR31B,tLdownR31C,tLdownR31D,
                                      K1,K2,K3,K4,K5,K6,K7,K8,
                                      K9,K10,K11,K12,K13,K14,K15,K16,
                                      K17,K18,K19,K20,K21,K22,K23,K24,
                                      K25,K26,K27,K28,K29,K30,K31,K32,
                                      t1,t2,t3,t4,t5,t6,t7,t8,
                                      t9,t10,t11,t12,t13,t14,t15,t16,
                                      t17,t18,t19,t20,t21,t22,t23,t24,
                                      t25,t26,t27,t28,t29,t30,t31,t32,
                                      Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                      Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                      Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                      Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32){
  Nasemean1 <- sigmoid(exp(K1),rdown1,exp(tLdownA21A),t1)
  Nasemean2 <- sigmoid(exp(K2),rdown1,exp(tLdownA21B),t2)
  Nasemean3 <- sigmoid(exp(K3),rdown1,exp(tLdownA21C),t3)
  Nasemean4 <- sigmoid(exp(K4),rdown1,exp(tLdownA21D),t4)
  Nasemean5 <- sigmoid(exp(K5),rdown1,exp(tLdownG21A),t5)
  Nasemean6 <- sigmoid(exp(K6),rdown1,exp(tLdownG21B),t6)
  Nasemean7 <- sigmoid(exp(K7),rdown1,exp(tLdownG21C),t7)
  Nasemean8 <- sigmoid(exp(K8),rdown1,exp(tLdownG21D),t8)
  Nasemean9 <- sigmoid(exp(K9),rdown1,exp(tLdownM21A),t9)
  Nasemean10 <- sigmoid(exp(K10),rdown1,exp(tLdownM21B),t10)
  Nasemean11 <- sigmoid(exp(K11),rdown1,exp(tLdownM21C),t11)
  Nasemean12 <- sigmoid(exp(K12),rdown1,exp(tLdownM21D),t12)
  Nasemean13 <- sigmoid(exp(K13),rdown1,exp(tLdownR21A),t13)
  Nasemean14 <- sigmoid(exp(K14),rdown1,exp(tLdownR21B),t14)
  Nasemean15 <- sigmoid(exp(K15),rdown1,exp(tLdownR21C),t15)
  Nasemean16 <- sigmoid(exp(K16),rdown1,exp(tLdownR21D),t16)
  Nasemean17 <- sigmoid(exp(K17),rdown31,exp(tLdownA31A),t17)
  Nasemean18 <- sigmoid(exp(K18),rdown31,exp(tLdownA31B),t18)
  Nasemean19 <- sigmoid(exp(K19),rdown31,exp(tLdownA31C),t19)
  Nasemean20 <- sigmoid(exp(K20),rdown31,exp(tLdownA31D),t20)
  Nasemean21 <- sigmoid(exp(K21),rdown31,exp(tLdownG31A),t21)
  Nasemean22 <- sigmoid(exp(K22),rdown31,exp(tLdownG31B),t22)
  Nasemean23 <- sigmoid(exp(K23),rdown31,exp(tLdownG31C),t23)
  Nasemean24 <- sigmoid(exp(K24),rdown31,exp(tLdownG31D),t24)
  Nasemean25 <- sigmoid(exp(K25),rdown31,exp(tLdownM31A),t25)
  Nasemean26 <- sigmoid(exp(K26),rdown31,exp(tLdownM31B),t26)
  Nasemean27 <- sigmoid(exp(K27),rdown31,exp(tLdownM31C),t27)
  Nasemean28 <- sigmoid(exp(K28),rdown31,exp(tLdownM31D),t28)
  Nasemean29 <- sigmoid(exp(K29),rdown31,exp(tLdownR31A),t29)
  Nasemean30 <- sigmoid(exp(K30),rdown31,exp(tLdownR31B),t30)
  Nasemean31 <- sigmoid(exp(K31),rdown31,exp(tLdownR31C),t31)
  Nasemean32 <- sigmoid(exp(K32),rdown31,exp(tLdownR31D),t32)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_down_temp <- mle2(sigmoid_down_temp_normNLL,start=list(sdNase=-1,
                                                                   K1=log(0.014),K2=log(0.01),
                                                                   K3=log(0.014),K4=log(0.017),
                                                                   K5=log(0.007),K6=log(0.001),
                                                                   K7=log(0.007),K8=log(0.008),
                                                                   K9=log(0.053),K10=log(0.029),
                                                                   K11=log(0.05),K12=log(0.078),
                                                                   K13=log(0.055),K14=log(0.038),
                                                                   K15=log(0.029),K16=log(0.032),
                                                                   K17=log(0.029),K18=log(0.024),
                                                                   K19=log(0.021),K20=log(0.012),
                                                                   K21=log(0.057),K22=log(0.033),
                                                                   K23=log(0.061),K24=log(0.016),
                                                                   K25=log(0.033),K26=log(0.061),
                                                                   K27=log(0.022),K28=log(0.036),
                                                                   K29=log(0.043),K30=log(0.017),
                                                                   K31=log(0.027),K32=log(0.021),
                                                                   rdown1=-0.1,
                                                                   rdown31=-0.1,
                                                                   tLdownA21A=log(12),tLdownA21B=log(0.4),tLdownA21C=log(0.001),tLdownA21D=log(0.001),
                                                                   tLdownG21A=log(32),tLdownG21B=log(6),tLdownG21C=log(6),tLdownG21D=log(7),
                                                                   tLdownM21A=log(0.001),tLdownM21B=log(16),tLdownM21C=log(14),tLdownM21D=log(0.6),
                                                                   tLdownR21A=log(6),tLdownR21B=log(43),tLdownR21C=log(40),tLdownR21D=log(12),
                                                                   tLdownA31A=log(0.01),tLdownA31B=log(4),tLdownA31C=log(6),tLdownA31D=log(10),
                                                                   tLdownG31A=log(0.1),tLdownG31B=log(0.001),tLdownG31C=log(0.001),tLdownG31D=log(20),
                                                                   tLdownM31A=log(9),tLdownM31B=log(0.002),tLdownM31C=log(0.001),tLdownM31D=log(6),
                                                                   tLdownR31A=log(0.002),tLdownR31B=log(20),tLdownR31C=log(0.001),tLdownR31D=log(13)),
                              data=list(t1=A21DA$Day,t2=A21DB$Day,
                                        t3=A21DC$Day,t4=A21DD$Day,
                                        t5=G21DA$Day,t6=G21DB$Day,
                                        t7=G21DC$Day,t8=G21DD$Day,
                                        t9=M21DA$Day,t10=M21DB$Day,
                                        t11=M21DC$Day,t12=M21DD$Day,
                                        t13=R21DA$Day,t14=R21DB$Day,
                                        t15=R21DC$Day,t16=R21DD$Day,
                                        t17=A31DA$Day,t18=A31DB$Day,
                                        t19=A31DC$Day,t20=A31DD$Day,
                                        t21=G31DA$Day,t22=G31DB$Day,
                                        t23=G31DC$Day,t24=G31DD$Day,
                                        t25=M31DA$Day,t26=M31DB$Day,
                                        t27=M31DC$Day,t28=M31DD$Day,
                                        t29=R31DA$Day,t30=R31DB$Day,
                                        t31=R31DC$Day,t32=R31DD$Day,
                                        Nasedat1=A21DA$SNF,Nasedat2=A21DB$SNF,
                                        Nasedat3=A21DC$SNF,Nasedat4=A21DD$SNF,
                                        Nasedat5=G21DA$SNF,Nasedat6=G21DB$SNF,
                                        Nasedat7=G21DC$SNF,Nasedat8=G21DD$SNF,
                                        Nasedat9=M21DA$SNF,Nasedat10=M21DB$SNF,
                                        Nasedat11=M21DC$SNF,Nasedat12=M21DD$SNF,
                                        Nasedat13=R21DA$SNF,Nasedat14=R21DB$SNF,
                                        Nasedat15=R21DC$SNF,Nasedat16=R21DD$SNF,
                                        Nasedat17=A31DA$SNF,Nasedat18=A31DB$SNF,
                                        Nasedat19=A31DC$SNF,Nasedat20=A31DD$SNF,
                                        Nasedat21=G31DA$SNF,Nasedat22=G31DB$SNF,
                                        Nasedat23=G31DC$SNF,Nasedat24=G31DD$SNF,
                                        Nasedat25=M31DA$SNF,Nasedat26=M31DB$SNF,
                                        Nasedat27=M31DC$SNF,Nasedat28=M31DD$SNF,
                                        Nasedat29=R31DA$SNF,Nasedat30=R31DB$SNF,
                                        Nasedat31=R31DC$SNF,Nasedat32=R31DD$SNF),
                              control=list(maxit=20000))
summary(fit_sigmoid_down_temp)

## Biome
sigmoid_down_biome_normNLL <- function(sdNase,rdowntemp,rdowntrop,
                                       tLdownA21A,tLdownA21B,tLdownA21C,tLdownA21D,
                                       tLdownG21A,tLdownG21B,tLdownG21C,tLdownG21D,
                                       tLdownM21A,tLdownM21B,tLdownM21C,tLdownM21D,
                                       tLdownR21A,tLdownR21B,tLdownR21C,tLdownR21D,
                                       tLdownA31A,tLdownA31B,tLdownA31C,tLdownA31D,
                                       tLdownG31A,tLdownG31B,tLdownG31C,tLdownG31D,
                                       tLdownM31A,tLdownM31B,tLdownM31C,tLdownM31D,
                                       tLdownR31A,tLdownR31B,tLdownR31C,tLdownR31D,
                                       K1,K2,K3,K4,K5,K6,K7,K8,
                                       K9,K10,K11,K12,K13,K14,K15,K16,
                                       K17,K18,K19,K20,K21,K22,K23,K24,
                                       K25,K26,K27,K28,K29,K30,K31,K32,
                                       t1,t2,t3,t4,t5,t6,t7,t8,
                                       t9,t10,t11,t12,t13,t14,t15,t16,
                                       t17,t18,t19,t20,t21,t22,t23,t24,
                                       t25,t26,t27,t28,t29,t30,t31,t32,
                                       Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                       Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                       Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                       Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32){
  Nasemean1 <- sigmoid(exp(K1),rdowntemp,exp(tLdownA21A),t1)
  Nasemean2 <- sigmoid(exp(K2),rdowntemp,exp(tLdownA21B),t2)
  Nasemean3 <- sigmoid(exp(K3),rdowntemp,exp(tLdownA21C),t3)
  Nasemean4 <- sigmoid(exp(K4),rdowntemp,exp(tLdownA21D),t4)
  Nasemean5 <- sigmoid(exp(K5),rdowntemp,exp(tLdownG21A),t5)
  Nasemean6 <- sigmoid(exp(K6),rdowntrop,exp(tLdownG21B),t6)
  Nasemean7 <- sigmoid(exp(K7),rdowntrop,exp(tLdownG21C),t7)
  Nasemean8 <- sigmoid(exp(K8),rdowntrop,exp(tLdownG21D),t8)
  Nasemean9 <- sigmoid(exp(K9),rdowntrop,exp(tLdownM21A),t9)
  Nasemean10 <- sigmoid(exp(K10),rdowntrop,exp(tLdownM21B),t10)
  Nasemean11 <- sigmoid(exp(K11),rdowntrop,exp(tLdownM21C),t11)
  Nasemean12 <- sigmoid(exp(K12),rdowntrop,exp(tLdownM21D),t12)
  Nasemean13 <- sigmoid(exp(K13),rdowntemp,exp(tLdownR21A),t13)
  Nasemean14 <- sigmoid(exp(K14),rdowntemp,exp(tLdownR21B),t14)
  Nasemean15 <- sigmoid(exp(K15),rdowntemp,exp(tLdownR21C),t15)
  Nasemean16 <- sigmoid(exp(K16),rdowntemp,exp(tLdownR21D),t16)
  Nasemean17 <- sigmoid(exp(K17),rdowntemp,exp(tLdownA31A),t17)
  Nasemean18 <- sigmoid(exp(K18),rdowntemp,exp(tLdownA31B),t18)
  Nasemean19 <- sigmoid(exp(K19),rdowntemp,exp(tLdownA31C),t19)
  Nasemean20 <- sigmoid(exp(K20),rdowntemp,exp(tLdownA31D),t20)
  Nasemean21 <- sigmoid(exp(K21),rdowntrop,exp(tLdownG31A),t21)
  Nasemean22 <- sigmoid(exp(K22),rdowntrop,exp(tLdownG31B),t22)
  Nasemean23 <- sigmoid(exp(K23),rdowntrop,exp(tLdownG31C),t23)
  Nasemean24 <- sigmoid(exp(K24),rdowntrop,exp(tLdownG31D),t24)
  Nasemean25 <- sigmoid(exp(K25),rdowntrop,exp(tLdownM31A),t25)
  Nasemean26 <- sigmoid(exp(K26),rdowntrop,exp(tLdownM31B),t26)
  Nasemean27 <- sigmoid(exp(K27),rdowntrop,exp(tLdownM31C),t27)
  Nasemean28 <- sigmoid(exp(K28),rdowntrop,exp(tLdownM31D),t28)
  Nasemean29 <- sigmoid(exp(K29),rdowntemp,exp(tLdownR31A),t29)
  Nasemean30 <- sigmoid(exp(K30),rdowntemp,exp(tLdownR31B),t30)
  Nasemean31 <- sigmoid(exp(K31),rdowntemp,exp(tLdownR31C),t31)
  Nasemean32 <- sigmoid(exp(K32),rdowntemp,exp(tLdownR31D),t32)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_down_biome <- mle2(sigmoid_down_biome_normNLL,start=list(sdNase=-1,
                                                                     K1=log(0.014),K2=log(0.01),
                                                                     K3=log(0.014),K4=log(0.017),
                                                                     K5=log(0.007),K6=log(0.001),
                                                                     K7=log(0.007),K8=log(0.008),
                                                                     K9=log(0.053),K10=log(0.029),
                                                                     K11=log(0.05),K12=log(0.078),
                                                                     K13=log(0.055),K14=log(0.038),
                                                                     K15=log(0.029),K16=log(0.032),
                                                                     K17=log(0.029),K18=log(0.024),
                                                                     K19=log(0.021),K20=log(0.012),
                                                                     K21=log(0.057),K22=log(0.033),
                                                                     K23=log(0.061),K24=log(0.016),
                                                                     K25=log(0.033),K26=log(0.061),
                                                                     K27=log(0.022),K28=log(0.036),
                                                                     K29=log(0.043),K30=log(0.017),
                                                                     K31=log(0.027),K32=log(0.021),
                                                                     rdowntemp=-0.1,
                                                                     rdowntrop=-0.1,
                                                                     tLdownA21A=log(12),tLdownA21B=log(0.4),tLdownA21C=log(0.001),tLdownA21D=log(0.001),
                                                                     tLdownG21A=log(32),tLdownG21B=log(6),tLdownG21C=log(6),tLdownG21D=log(7),
                                                                     tLdownM21A=log(0.001),tLdownM21B=log(16),tLdownM21C=log(14),tLdownM21D=log(0.6),
                                                                     tLdownR21A=log(6),tLdownR21B=log(43),tLdownR21C=log(40),tLdownR21D=log(12),
                                                                     tLdownA31A=log(0.01),tLdownA31B=log(4),tLdownA31C=log(6),tLdownA31D=log(10),
                                                                     tLdownG31A=log(0.1),tLdownG31B=log(0.001),tLdownG31C=log(0.001),tLdownG31D=log(20),
                                                                     tLdownM31A=log(9),tLdownM31B=log(0.002),tLdownM31C=log(0.001),tLdownM31D=log(6),
                                                                     tLdownR31A=log(0.002),tLdownR31B=log(20),tLdownR31C=log(0.001),tLdownR31D=log(13)),
                               data=list(t1=A21DA$Day,t2=A21DB$Day,
                                         t3=A21DC$Day,t4=A21DD$Day,
                                         t5=G21DA$Day,t6=G21DB$Day,
                                         t7=G21DC$Day,t8=G21DD$Day,
                                         t9=M21DA$Day,t10=M21DB$Day,
                                         t11=M21DC$Day,t12=M21DD$Day,
                                         t13=R21DA$Day,t14=R21DB$Day,
                                         t15=R21DC$Day,t16=R21DD$Day,
                                         t17=A31DA$Day,t18=A31DB$Day,
                                         t19=A31DC$Day,t20=A31DD$Day,
                                         t21=G31DA$Day,t22=G31DB$Day,
                                         t23=G31DC$Day,t24=G31DD$Day,
                                         t25=M31DA$Day,t26=M31DB$Day,
                                         t27=M31DC$Day,t28=M31DD$Day,
                                         t29=R31DA$Day,t30=R31DB$Day,
                                         t31=R31DC$Day,t32=R31DD$Day,
                                         Nasedat1=A21DA$SNF,Nasedat2=A21DB$SNF,
                                         Nasedat3=A21DC$SNF,Nasedat4=A21DD$SNF,
                                         Nasedat5=G21DA$SNF,Nasedat6=G21DB$SNF,
                                         Nasedat7=G21DC$SNF,Nasedat8=G21DD$SNF,
                                         Nasedat9=M21DA$SNF,Nasedat10=M21DB$SNF,
                                         Nasedat11=M21DC$SNF,Nasedat12=M21DD$SNF,
                                         Nasedat13=R21DA$SNF,Nasedat14=R21DB$SNF,
                                         Nasedat15=R21DC$SNF,Nasedat16=R21DD$SNF,
                                         Nasedat17=A31DA$SNF,Nasedat18=A31DB$SNF,
                                         Nasedat19=A31DC$SNF,Nasedat20=A31DD$SNF,
                                         Nasedat21=G31DA$SNF,Nasedat22=G31DB$SNF,
                                         Nasedat23=G31DC$SNF,Nasedat24=G31DD$SNF,
                                         Nasedat25=M31DA$SNF,Nasedat26=M31DB$SNF,
                                         Nasedat27=M31DC$SNF,Nasedat28=M31DD$SNF,
                                         Nasedat29=R31DA$SNF,Nasedat30=R31DB$SNF,
                                         Nasedat31=R31DC$SNF,Nasedat32=R31DD$SNF),
                               control=list(maxit=20000))
summary(fit_sigmoid_down_biome)

## Biome x temperature interaction
sigmoid_down_biometemp_normNLL <- function(sdNase,rdowntemp21,rdowntrop21,rdowntemp31,rdowntrop31,
                                           tLdownA21A,tLdownA21B,tLdownA21C,tLdownA21D,
                                           tLdownG21A,tLdownG21B,tLdownG21C,tLdownG21D,
                                           tLdownM21A,tLdownM21B,tLdownM21C,tLdownM21D,
                                           tLdownR21A,tLdownR21B,tLdownR21C,tLdownR21D,
                                           tLdownA31A,tLdownA31B,tLdownA31C,tLdownA31D,
                                           tLdownG31A,tLdownG31B,tLdownG31C,tLdownG31D,
                                           tLdownM31A,tLdownM31B,tLdownM31C,tLdownM31D,
                                           tLdownR31A,tLdownR31B,tLdownR31C,tLdownR31D,
                                           K1,K2,K3,K4,K5,K6,K7,K8,
                                           K9,K10,K11,K12,K13,K14,K15,K16,
                                           K17,K18,K19,K20,K21,K22,K23,K24,
                                           K25,K26,K27,K28,K29,K30,K31,K32,
                                           t1,t2,t3,t4,t5,t6,t7,t8,
                                           t9,t10,t11,t12,t13,t14,t15,t16,
                                           t17,t18,t19,t20,t21,t22,t23,t24,
                                           t25,t26,t27,t28,t29,t30,t31,t32,
                                           Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                           Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                           Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                           Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32){
  Nasemean1 <- sigmoid(exp(K1),rdowntemp21,exp(tLdownA21A),t1)
  Nasemean2 <- sigmoid(exp(K2),rdowntemp21,exp(tLdownA21B),t2)
  Nasemean3 <- sigmoid(exp(K3),rdowntemp21,exp(tLdownA21C),t3)
  Nasemean4 <- sigmoid(exp(K4),rdowntemp21,exp(tLdownA21D),t4)
  Nasemean5 <- sigmoid(exp(K5),rdowntemp21,exp(tLdownG21A),t5)
  Nasemean6 <- sigmoid(exp(K6),rdowntrop21,exp(tLdownG21B),t6)
  Nasemean7 <- sigmoid(exp(K7),rdowntrop21,exp(tLdownG21C),t7)
  Nasemean8 <- sigmoid(exp(K8),rdowntrop21,exp(tLdownG21D),t8)
  Nasemean9 <- sigmoid(exp(K9),rdowntrop21,exp(tLdownM21A),t9)
  Nasemean10 <- sigmoid(exp(K10),rdowntrop21,exp(tLdownM21B),t10)
  Nasemean11 <- sigmoid(exp(K11),rdowntrop21,exp(tLdownM21C),t11)
  Nasemean12 <- sigmoid(exp(K12),rdowntrop21,exp(tLdownM21D),t12)
  Nasemean13 <- sigmoid(exp(K13),rdowntemp21,exp(tLdownR21A),t13)
  Nasemean14 <- sigmoid(exp(K14),rdowntemp21,exp(tLdownR21B),t14)
  Nasemean15 <- sigmoid(exp(K15),rdowntemp21,exp(tLdownR21C),t15)
  Nasemean16 <- sigmoid(exp(K16),rdowntemp21,exp(tLdownR21D),t16)
  Nasemean17 <- sigmoid(exp(K17),rdowntemp31,exp(tLdownA31A),t17)
  Nasemean18 <- sigmoid(exp(K18),rdowntemp31,exp(tLdownA31B),t18)
  Nasemean19 <- sigmoid(exp(K19),rdowntemp31,exp(tLdownA31C),t19)
  Nasemean20 <- sigmoid(exp(K20),rdowntemp31,exp(tLdownA31D),t20)
  Nasemean21 <- sigmoid(exp(K21),rdowntrop31,exp(tLdownG31A),t21)
  Nasemean22 <- sigmoid(exp(K22),rdowntrop31,exp(tLdownG31B),t22)
  Nasemean23 <- sigmoid(exp(K23),rdowntrop31,exp(tLdownG31C),t23)
  Nasemean24 <- sigmoid(exp(K24),rdowntrop31,exp(tLdownG31D),t24)
  Nasemean25 <- sigmoid(exp(K25),rdowntrop31,exp(tLdownM31A),t25)
  Nasemean26 <- sigmoid(exp(K26),rdowntrop31,exp(tLdownM31B),t26)
  Nasemean27 <- sigmoid(exp(K27),rdowntrop31,exp(tLdownM31C),t27)
  Nasemean28 <- sigmoid(exp(K28),rdowntrop31,exp(tLdownM31D),t28)
  Nasemean29 <- sigmoid(exp(K29),rdowntemp31,exp(tLdownR31A),t29)
  Nasemean30 <- sigmoid(exp(K30),rdowntemp31,exp(tLdownR31B),t30)
  Nasemean31 <- sigmoid(exp(K31),rdowntemp31,exp(tLdownR31C),t31)
  Nasemean32 <- sigmoid(exp(K32),rdowntemp31,exp(tLdownR31D),t32)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_down_biometemp <- mle2(sigmoid_down_biometemp_normNLL,start=list(sdNase=-1,
                                                                             K1=log(0.014),K2=log(0.01),
                                                                             K3=log(0.014),K4=log(0.017),
                                                                             K5=log(0.007),K6=log(0.001),
                                                                             K7=log(0.007),K8=log(0.008),
                                                                             K9=log(0.053),K10=log(0.029),
                                                                             K11=log(0.05),K12=log(0.078),
                                                                             K13=log(0.055),K14=log(0.038),
                                                                             K15=log(0.029),K16=log(0.032),
                                                                             K17=log(0.029),K18=log(0.024),
                                                                             K19=log(0.021),K20=log(0.012),
                                                                             K21=log(0.057),K22=log(0.033),
                                                                             K23=log(0.061),K24=log(0.016),
                                                                             K25=log(0.033),K26=log(0.061),
                                                                             K27=log(0.022),K28=log(0.036),
                                                                             K29=log(0.043),K30=log(0.017),
                                                                             K31=log(0.027),K32=log(0.021),
                                                                             rdowntemp21=-0.1,
                                                                             rdowntrop21=-0.1,
                                                                             rdowntemp31=-0.1,
                                                                             rdowntrop31=-0.1,
                                                                             tLdownA21A=log(12),tLdownA21B=log(0.4),tLdownA21C=log(0.001),tLdownA21D=log(0.001),
                                                                             tLdownG21A=log(32),tLdownG21B=log(6),tLdownG21C=log(6),tLdownG21D=log(7),
                                                                             tLdownM21A=log(0.001),tLdownM21B=log(16),tLdownM21C=log(14),tLdownM21D=log(0.6),
                                                                             tLdownR21A=log(6),tLdownR21B=log(43),tLdownR21C=log(40),tLdownR21D=log(12),
                                                                             tLdownA31A=log(0.01),tLdownA31B=log(4),tLdownA31C=log(6),tLdownA31D=log(10),
                                                                             tLdownG31A=log(0.1),tLdownG31B=log(0.001),tLdownG31C=log(0.001),tLdownG31D=log(20),
                                                                             tLdownM31A=log(9),tLdownM31B=log(0.002),tLdownM31C=log(0.001),tLdownM31D=log(6),
                                                                             tLdownR31A=log(0.002),tLdownR31B=log(20),tLdownR31C=log(0.001),tLdownR31D=log(13)),
                                   data=list(t1=A21DA$Day,t2=A21DB$Day,
                                             t3=A21DC$Day,t4=A21DD$Day,
                                             t5=G21DA$Day,t6=G21DB$Day,
                                             t7=G21DC$Day,t8=G21DD$Day,
                                             t9=M21DA$Day,t10=M21DB$Day,
                                             t11=M21DC$Day,t12=M21DD$Day,
                                             t13=R21DA$Day,t14=R21DB$Day,
                                             t15=R21DC$Day,t16=R21DD$Day,
                                             t17=A31DA$Day,t18=A31DB$Day,
                                             t19=A31DC$Day,t20=A31DD$Day,
                                             t21=G31DA$Day,t22=G31DB$Day,
                                             t23=G31DC$Day,t24=G31DD$Day,
                                             t25=M31DA$Day,t26=M31DB$Day,
                                             t27=M31DC$Day,t28=M31DD$Day,
                                             t29=R31DA$Day,t30=R31DB$Day,
                                             t31=R31DC$Day,t32=R31DD$Day,
                                             Nasedat1=A21DA$SNF,Nasedat2=A21DB$SNF,
                                             Nasedat3=A21DC$SNF,Nasedat4=A21DD$SNF,
                                             Nasedat5=G21DA$SNF,Nasedat6=G21DB$SNF,
                                             Nasedat7=G21DC$SNF,Nasedat8=G21DD$SNF,
                                             Nasedat9=M21DA$SNF,Nasedat10=M21DB$SNF,
                                             Nasedat11=M21DC$SNF,Nasedat12=M21DD$SNF,
                                             Nasedat13=R21DA$SNF,Nasedat14=R21DB$SNF,
                                             Nasedat15=R21DC$SNF,Nasedat16=R21DD$SNF,
                                             Nasedat17=A31DA$SNF,Nasedat18=A31DB$SNF,
                                             Nasedat19=A31DC$SNF,Nasedat20=A31DD$SNF,
                                             Nasedat21=G31DA$SNF,Nasedat22=G31DB$SNF,
                                             Nasedat23=G31DC$SNF,Nasedat24=G31DD$SNF,
                                             Nasedat25=M31DA$SNF,Nasedat26=M31DB$SNF,
                                             Nasedat27=M31DC$SNF,Nasedat28=M31DD$SNF,
                                             Nasedat29=R31DA$SNF,Nasedat30=R31DB$SNF,
                                             Nasedat31=R31DC$SNF,Nasedat32=R31DD$SNF),
                                   control=list(maxit=20000))
summary(fit_sigmoid_down_biometemp)

## Calculate delta AICc to find best model
AICctab(fit_sigmoid_down_same,fit_sigmoid_down_sptemp,fit_sigmoid_down_sym,
        fit_sigmoid_down_symtemp,fit_sigmoid_down_sp,fit_sigmoid_down_temp,
        fit_sigmoid_down_biome,fit_sigmoid_down_biometemp,nobs=272)

## Likelihood ratio test to see if temperature is a necessary predictor
anova(fit_sigmoid_down_biome,fit_sigmoid_down_biometemp) #p=0.01447; keep temperature effect
#biome x temp is best
#temperature increases rate
#temperate rate faster than tropical

### Compare tL values for best fit (biome x temperature interaction)
ft_biometemp_down<-coef(fit_sigmoid_down_biometemp)
tL.down<-exp(c(ft_biometemp_down[6],ft_biometemp_down[7],ft_biometemp_down[8],ft_biometemp_down[9],
               ft_biometemp_down[10],ft_biometemp_down[11],ft_biometemp_down[12],ft_biometemp_down[13],
               ft_biometemp_down[14],ft_biometemp_down[15],ft_biometemp_down[16],ft_biometemp_down[17],
               ft_biometemp_down[18],ft_biometemp_down[19],ft_biometemp_down[20],ft_biometemp_down[21],
               ft_biometemp_down[22],ft_biometemp_down[23],ft_biometemp_down[24],ft_biometemp_down[25],
               ft_biometemp_down[26],ft_biometemp_down[27],ft_biometemp_down[28],ft_biometemp_down[29],
               ft_biometemp_down[30],ft_biometemp_down[31],ft_biometemp_down[32],ft_biometemp_down[33],
               ft_biometemp_down[34],ft_biometemp_down[35],ft_biometemp_down[36],ft_biometemp_down[37]))
treat.down<-c("A21","A21","A21","A21","G21","G21","G21","G21","M21","M21","M21","M21","R21","R21","R21","R21",
              "A31","A31","A31","A31","G31","G31","G31","G31","M31","M31","M31","M31","R31","R31","R31","R31")
sym.down<-c("Actin","Actin","Actin","Actin","Rhiz","Rhiz","Rhiz","Rhiz","Actin","Actin","Actin","Actin","Rhiz","Rhiz","Rhiz","Rhiz",
            "Actin","Actin","Actin","Actin","Rhiz","Rhiz","Rhiz","Rhiz","Actin","Actin","Actin","Actin","Rhiz","Rhiz","Rhiz","Rhiz")
biome.down<-c("temp","temp","temp","temp","trop","trop","trop","trop","trop","trop","trop","trop","temp","temp","temp","temp",
              "temp","temp","temp","temp","trop","trop","trop","trop","trop","trop","trop","trop","temp","temp","temp","temp")
sp.down<-c("A","A","A","A","G","G","G","G","M","M","M","M","R","R","R","R",
           "A","A","A","A","G","G","G","G","M","M","M","M","R","R","R","R")
temp.down<-c("21","21","21","21","21","21","21","21","21","21","21","21","21","21","21","21",
             "31","31","31","31","31","31","31","31","31","31","31","31","31","31","31","31")

tL.aov<-aov(tL.down~treat.down)
summary(tL.aov) #NS
TukeyHSD(tL.aov)
#No difference between treatments

tL.down.null.lm<-lm(tL.down~1)
tL.down.tr.lm<-lm(tL.down~treat.down)
tL.down.sym.lm<-lm(tL.down~sym.down)
tL.down.temp.lm<-lm(tL.down~temp.down)
tL.down.biome.lm<-lm(tL.down~biome.down)
tL.down.sp.lm<-lm(tL.down~sp.down)
tL.down.sp.temp.lm<-lm(tL.down~sp.down*temp.down)
tL.down.biome.temp.lm<-lm(tL.down~biome.down*temp.down)
tL.down.sym.temp.lm<-lm(tL.down~sym.down*temp.down)

## Calculate delta AICc to find best model
AICctab(tL.down.null.lm,tL.down.tr.lm,tL.down.sym.lm,
        tL.down.temp.lm,tL.down.biome.lm,tL.down.sp.lm,tL.down.sp.temp.lm,
        tL.down.biome.temp.lm,tL.down.sym.temp.lm, nobs = 32)
#Symbiosis is best

## Likelihood ratio test to see if symbiosis is a necessary predictor
anova(tL.down.sym.lm, tL.down.null.lm) #p=0.129
#no difference between tL values is best model

### Fit up-regulation curves with Eq. 5 using maximum likelihood

## One r for all
sigmoid_up_same_normNLL <- function(sdNase,rup,tLupG21A,tLupG21B,tLupG21C,
                                       tLupM21A,tLupM21B,tLupM21C,tLupM21D,
                                       tLupR21A,tLupR21B,tLupR21C,
                                       tLupA31A,tLupA31B,tLupA31C,
                                       tLupG31A,tLupG31B,tLupG31C,tLupG31D,tLupG31E,tLupG31F,
                                       tLupM31A,tLupM31B,tLupM31C,tLupM31D,tLupM31E,tLupM31F,tLupM31G,tLupM31H,
                                       tLupR31A,tLupR31B,tLupR31C,tLupR31D,tLupR31E,tLupR31F,
                                       K1,K2,K3,K4,K5,K6,K7,K8,
                                       K9,K10,K11,K12,K13,K14,K15,K16,
                                       K17,K18,K19,K20,K21,K22,K23,K24,
                                       K25,K26,K27,K28,K29,K30,K31,K32,K33,
                                       t1,t2,t3,t4,t5,t6,t7,t8,
                                       t9,t10,t11,t12,t13,t14,t15,t16,
                                       t17,t18,t19,t20,t21,t22,t23,t24,
                                       t25,t26,t27,t28,t29,t30,t31,t32,t33,
                                       Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                       Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                       Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                       Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32,Nadedat33){
  Nasemean1 <- sigmoid(exp(K1),exp(rup),exp(tLupG21A),t1)
  Nasemean2 <- sigmoid(exp(K2),exp(rup),exp(tLupG21B),t2)
  Nasemean3 <- sigmoid(exp(K3),exp(rup),exp(tLupG21C),t3)
  Nasemean4 <- sigmoid(exp(K4),exp(rup),exp(tLupM21A),t4)
  Nasemean5 <- sigmoid(exp(K5),exp(rup),exp(tLupM21B),t5)
  Nasemean6 <- sigmoid(exp(K6),exp(rup),exp(tLupM21C),t6)
  Nasemean7 <- sigmoid(exp(K7),exp(rup),exp(tLupM21D),t7)
  Nasemean8 <- sigmoid(exp(K8),exp(rup),exp(tLupR21A),t8)
  Nasemean9 <- sigmoid(exp(K9),exp(rup),exp(tLupR21B),t9)
  Nasemean10 <- sigmoid(exp(K10),exp(rup),exp(tLupR21C),t10)
  Nasemean11 <- sigmoid(exp(K11),exp(rup),exp(tLupA31A),t11)
  Nasemean12 <- sigmoid(exp(K12),exp(rup),exp(tLupA31B),t12)
  Nasemean13 <- sigmoid(exp(K13),exp(rup),exp(tLupA31C),t13)
  Nasemean14 <- sigmoid(exp(K14),exp(rup),exp(tLupG31A),t14)
  Nasemean15 <- sigmoid(exp(K15),exp(rup),exp(tLupG31B),t15)
  Nasemean16 <- sigmoid(exp(K16),exp(rup),exp(tLupG31C),t16)
  Nasemean17 <- sigmoid(exp(K17),exp(rup),exp(tLupG31D),t17)
  Nasemean18 <- sigmoid(exp(K18),exp(rup),exp(tLupG31E),t18)
  Nasemean19 <- sigmoid(exp(K19),exp(rup),exp(tLupG31F),t19)
  Nasemean20 <- sigmoid(exp(K20),exp(rup),exp(tLupM31A),t20)
  Nasemean21 <- sigmoid(exp(K21),exp(rup),exp(tLupM31B),t21)
  Nasemean22 <- sigmoid(exp(K22),exp(rup),exp(tLupM31C),t22)
  Nasemean23 <- sigmoid(exp(K23),exp(rup),exp(tLupM31D),t23)
  Nasemean24 <- sigmoid(exp(K24),exp(rup),exp(tLupM31E),t24)
  Nasemean25 <- sigmoid(exp(K25),exp(rup),exp(tLupM31F),t25)
  Nasemean26 <- sigmoid(exp(K26),exp(rup),exp(tLupM31G),t26)
  Nasemean27 <- sigmoid(exp(K27),exp(rup),exp(tLupM31H),t27)
  Nasemean28 <- sigmoid(exp(K28),exp(rup),exp(tLupR31A),t28)
  Nasemean29 <- sigmoid(exp(K29),exp(rup),exp(tLupR31B),t29)
  Nasemean30 <- sigmoid(exp(K30),exp(rup),exp(tLupR31C),t30)
  Nasemean31 <- sigmoid(exp(K31),exp(rup),exp(tLupR31D),t31)
  Nasemean32 <- sigmoid(exp(K32),exp(rup),exp(tLupR31E),t32)
  Nasemean33 <- sigmoid(exp(K33),exp(rup),exp(tLupR31F),t33)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat33,mean=Nasemean33,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_up_same <- mle2(sigmoid_up_same_normNLL,start=list(sdNase=-1,
                                                                   K1=log(max(G21UA$SNF)),K2=log(max(G21UB$SNF)),
                                                                   K3=log(max(G21UC$SNF)),K4=log(max(M21UA$SNF)),
                                                                   K5=log(max(M21UB$SNF)),K6=log(max(M21UC$SNF)),
                                                                   K7=log(max(M21UD$SNF)),K8=log(max(R21UA$SNF)),
                                                                   K9=log(max(R21UB$SNF)),K10=log(max(R21UC$SNF)),
                                                                   K11=log(max(A31UA$SNF)),K12=log(max(A31UB$SNF)),
                                                                   K13=log(max(A31UC$SNF)),K14=log(max(G31UA$SNF)),
                                                                   K15=log(max(G31UB$SNF)),K16=log(max(G31UC$SNF)),
                                                                   K17=log(max(G31UD$SNF)),K18=log(max(G31UE$SNF)),
                                                                   K19=log(max(G31UF$SNF)),K20=log(max(M31UA$SNF)),
                                                                   K21=log(max(M31UB$SNF)),K22=log(max(M31UC$SNF)),
                                                                   K23=log(max(M31UD$SNF)),K24=log(max(M31UE$SNF)),
                                                                   K25=log(max(M31UF$SNF)),K26=log(max(M31UG$SNF)),
                                                                   K27=log(max(M31UH$SNF)),K28=log(max(R31UA$SNF)),
                                                                   K29=log(max(R31UB$SNF)),K30=log(max(R31UC$SNF)),
                                                                   K31=log(max(R31UD$SNF)),K32=log(max(R31UE$SNF)),
                                                                   K33=log(max(R31UF$SNF)),
                                                                   rup=log(0.1),tLupG21A=log(70),tLupG21B=log(110),tLupG21C=log(100),
                                                                   tLupM21A=log(150),tLupM21B=log(100),tLupM21C=log(140),tLupM21D=log(180),
                                                                   tLupR21A=log(60),tLupR21B=log(30),tLupR21C=log(40),
                                                                   tLupA31A=log(160),tLupA31B=log(110),tLupA31C=log(110),
                                                                   tLupG31A=log(100),tLupG31B=log(70),tLupG31C=log(110),tLupG31D=log(40),tLupG31E=log(70),tLupG31F=log(110),
                                                                   tLupM31A=log(110),tLupM31B=log(120),tLupM31C=log(150),tLupM31D=log(100),tLupM31E=log(130),tLupM31F=log(130),tLupM31G=log(130),tLupM31H=log(130),
                                                                   tLupR31A=log(70),tLupR31B=log(100),tLupR31C=log(120),tLupR31D=log(30),tLupR31E=log(120),tLupR31F=log(70)),
                             data=list(t1=G21UA$Day,t2=G21UB$Day,
                                       t3=G21UC$Day,t4=M21UA$Day,
                                       t5=M21UB$Day,t6=M21UC$Day,
                                       t7=M21UD$Day,t8=R21UA$Day,
                                       t9=R21UB$Day,t10=R21UC$Day,
                                       t11=A31UA$Day,t12=A31UB$Day,
                                       t13=A31UC$Day,t14=G31UA$Day,
                                       t15=G31UB$Day,t16=G31UC$Day,
                                       t17=G31UD$Day,t18=G31UE$Day,
                                       t19=G31UF$Day,t20=M31UA$Day,
                                       t21=M31UB$Day,t22=M31UC$Day,
                                       t23=M31UD$Day,t24=M31UE$Day,
                                       t25=M31UF$Day,t26=M31UG$Day,
                                       t27=M31UH$Day,t28=R31UA$Day,
                                       t29=R31UB$Day,t30=R31UC$Day,
                                       t31=R31UD$Day,t32=R31UE$Day,
                                       t33=R31UF$Day,
                                       Nasedat1=G21UA$SNF,Nasedat2=G21UB$SNF,
                                       Nasedat3=G21UC$SNF,Nasedat4=M21UA$SNF,
                                       Nasedat5=M21UB$SNF,Nasedat6=M21UC$SNF,
                                       Nasedat7=M21UD$SNF,Nasedat8=R21UA$SNF,
                                       Nasedat9=R21UB$SNF,Nasedat10=R21UC$SNF,
                                       Nasedat11=A31UA$SNF,Nasedat12=A31UB$SNF,
                                       Nasedat13=A31UC$SNF,Nasedat14=G31UA$SNF,
                                       Nasedat15=G31UB$SNF,Nasedat16=G31UC$SNF,
                                       Nasedat17=G31UD$SNF,Nasedat18=G31UE$SNF,
                                       Nasedat19=G31UF$SNF,Nasedat20=M31UA$SNF,
                                       Nasedat21=M31UB$SNF,Nasedat22=M31UC$SNF,
                                       Nasedat23=M31UD$SNF,Nasedat24=M31UE$SNF,
                                       Nasedat25=M31UF$SNF,Nasedat26=M31UG$SNF,
                                       Nasedat27=M31UH$SNF,Nasedat28=R31UA$SNF,
                                       Nasedat29=R31UB$SNF,Nasedat30=R31UC$SNF,
                                       Nasedat31=R31UD$SNF,Nasedat32=R31UE$SNF,
                                       Nasedat33=R31UF$SNF),
                             control=list(maxit=20000))
summary(fit_sigmoid_up_same)

## Species x temperature interaction (one r for each treatment)
sigmoid_up_sptemp_normNLL <- function(sdNase,rupG21,tLupG21A,tLupG21B,tLupG21C,
                                         rupM21,tLupM21A,tLupM21B,tLupM21C,tLupM21D,
                                         rupR21,tLupR21A,tLupR21B,tLupR21C,
                                         rupA31,tLupA31A,tLupA31B,tLupA31C,
                                         rupG31,tLupG31A,tLupG31B,tLupG31C,tLupG31D,tLupG31E,tLupG31F,
                                         rupM31,tLupM31A,tLupM31B,tLupM31C,tLupM31D,tLupM31E,tLupM31F,tLupM31G,tLupM31H,
                                         rupR31,tLupR31A,tLupR31B,tLupR31C,tLupR31D,tLupR31E,tLupR31F,
                                         K1,K2,K3,K4,K5,K6,K7,K8,
                                         K9,K10,K11,K12,K13,K14,K15,K16,
                                         K17,K18,K19,K20,K21,K22,K23,K24,
                                         K25,K26,K27,K28,K29,K30,K31,K32,K33,
                                         t1,t2,t3,t4,t5,t6,t7,t8,
                                         t9,t10,t11,t12,t13,t14,t15,t16,
                                         t17,t18,t19,t20,t21,t22,t23,t24,
                                         t25,t26,t27,t28,t29,t30,t31,t32,t33,
                                         Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                         Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                         Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                         Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32,Nadedat33){
  Nasemean1 <- sigmoid(exp(K1),exp(rupG21),exp(tLupG21A),t1)
  Nasemean2 <- sigmoid(exp(K2),exp(rupG21),exp(tLupG21B),t2)
  Nasemean3 <- sigmoid(exp(K3),exp(rupG21),exp(tLupG21C),t3)
  Nasemean4 <- sigmoid(exp(K4),exp(rupM21),exp(tLupM21A),t4)
  Nasemean5 <- sigmoid(exp(K5),exp(rupM21),exp(tLupM21B),t5)
  Nasemean6 <- sigmoid(exp(K6),exp(rupM21),exp(tLupM21C),t6)
  Nasemean7 <- sigmoid(exp(K7),exp(rupM21),exp(tLupM21D),t7)
  Nasemean8 <- sigmoid(exp(K8),exp(rupR21),exp(tLupR21A),t8)
  Nasemean9 <- sigmoid(exp(K9),exp(rupR21),exp(tLupR21B),t9)
  Nasemean10 <- sigmoid(exp(K10),exp(rupR21),exp(tLupR21C),t10)
  Nasemean11 <- sigmoid(exp(K11),exp(rupA31),exp(tLupA31A),t11)
  Nasemean12 <- sigmoid(exp(K12),exp(rupA31),exp(tLupA31B),t12)
  Nasemean13 <- sigmoid(exp(K13),exp(rupA31),exp(tLupA31C),t13)
  Nasemean14 <- sigmoid(exp(K14),exp(rupG31),exp(tLupG31A),t14)
  Nasemean15 <- sigmoid(exp(K15),exp(rupG31),exp(tLupG31B),t15)
  Nasemean16 <- sigmoid(exp(K16),exp(rupG31),exp(tLupG31C),t16)
  Nasemean17 <- sigmoid(exp(K17),exp(rupG31),exp(tLupG31D),t17)
  Nasemean18 <- sigmoid(exp(K18),exp(rupG31),exp(tLupG31E),t18)
  Nasemean19 <- sigmoid(exp(K19),exp(rupG31),exp(tLupG31F),t19)
  Nasemean20 <- sigmoid(exp(K20),exp(rupM31),exp(tLupM31A),t20)
  Nasemean21 <- sigmoid(exp(K21),exp(rupM31),exp(tLupM31B),t21)
  Nasemean22 <- sigmoid(exp(K22),exp(rupM31),exp(tLupM31C),t22)
  Nasemean23 <- sigmoid(exp(K23),exp(rupM31),exp(tLupM31D),t23)
  Nasemean24 <- sigmoid(exp(K24),exp(rupM31),exp(tLupM31E),t24)
  Nasemean25 <- sigmoid(exp(K25),exp(rupM31),exp(tLupM31F),t25)
  Nasemean26 <- sigmoid(exp(K26),exp(rupM31),exp(tLupM31G),t26)
  Nasemean27 <- sigmoid(exp(K27),exp(rupM31),exp(tLupM31H),t27)
  Nasemean28 <- sigmoid(exp(K28),exp(rupR31),exp(tLupR31A),t28)
  Nasemean29 <- sigmoid(exp(K29),exp(rupR31),exp(tLupR31B),t29)
  Nasemean30 <- sigmoid(exp(K30),exp(rupR31),exp(tLupR31C),t30)
  Nasemean31 <- sigmoid(exp(K31),exp(rupR31),exp(tLupR31D),t31)
  Nasemean32 <- sigmoid(exp(K32),exp(rupR31),exp(tLupR31E),t32)
  Nasemean33 <- sigmoid(exp(K33),exp(rupR31),exp(tLupR31F),t33)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat33,mean=Nasemean33,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_up_sptemp <- mle2(sigmoid_up_sptemp_normNLL,start=list(sdNase=-1,
                                                                       K1=log(max(G21UA$SNF)),K2=log(max(G21UB$SNF)),
                                                                       K3=log(max(G21UC$SNF)),K4=log(max(M21UA$SNF)),
                                                                       K5=log(max(M21UB$SNF)),K6=log(max(M21UC$SNF)),
                                                                       K7=log(max(M21UD$SNF)),K8=log(max(R21UA$SNF)),
                                                                       K9=log(max(R21UB$SNF)),K10=log(max(R21UC$SNF)),
                                                                       K11=log(max(A31UA$SNF)),K12=log(max(A31UB$SNF)),
                                                                       K13=log(max(A31UC$SNF)),K14=log(max(G31UA$SNF)),
                                                                       K15=log(max(G31UB$SNF)),K16=log(max(G31UC$SNF)),
                                                                       K17=log(max(G31UD$SNF)),K18=log(max(G31UE$SNF)),
                                                                       K19=log(max(G31UF$SNF)),K20=log(max(M31UA$SNF)),
                                                                       K21=log(max(M31UB$SNF)),K22=log(max(M31UC$SNF)),
                                                                       K23=log(max(M31UD$SNF)),K24=log(max(M31UE$SNF)),
                                                                       K25=log(max(M31UF$SNF)),K26=log(max(M31UG$SNF)),
                                                                       K27=log(max(M31UH$SNF)),K28=log(max(R31UA$SNF)),
                                                                       K29=log(max(R31UB$SNF)),K30=log(max(R31UC$SNF)),
                                                                       K31=log(max(R31UD$SNF)),K32=log(max(R31UE$SNF)),
                                                                       K33=log(max(R31UF$SNF)),
                                                                       rupG21=log(0.1),tLupG21A=log(70),tLupG21B=log(110),tLupG21C=log(100),
                                                                       rupM21=log(0.1),tLupM21A=log(150),tLupM21B=log(100),tLupM21C=log(140),tLupM21D=log(180),
                                                                       rupR21=log(0.05),tLupR21A=log(60),tLupR21B=log(30),tLupR21C=log(40),
                                                                       rupA31=log(0.2),tLupA31A=log(160),tLupA31B=log(110),tLupA31C=log(110),
                                                                       rupG31=log(0.05),tLupG31A=log(100),tLupG31B=log(70),tLupG31C=log(110),tLupG31D=log(40),tLupG31E=log(70),tLupG31F=log(110),
                                                                       rupM31=log(0.3),tLupM31A=log(110),tLupM31B=log(120),tLupM31C=log(150),tLupM31D=log(100),tLupM31E=log(130),tLupM31F=log(130),tLupM31G=log(130),tLupM31H=log(130),
                                                                       rupR31=log(0.06),tLupR31A=log(70),tLupR31B=log(100),tLupR31C=log(120),tLupR31D=log(30),tLupR31E=log(120),tLupR31F=log(70)),
                               data=list(t1=G21UA$Day,t2=G21UB$Day,
                                         t3=G21UC$Day,t4=M21UA$Day,
                                         t5=M21UB$Day,t6=M21UC$Day,
                                         t7=M21UD$Day,t8=R21UA$Day,
                                         t9=R21UB$Day,t10=R21UC$Day,
                                         t11=A31UA$Day,t12=A31UB$Day,
                                         t13=A31UC$Day,t14=G31UA$Day,
                                         t15=G31UB$Day,t16=G31UC$Day,
                                         t17=G31UD$Day,t18=G31UE$Day,
                                         t19=G31UF$Day,t20=M31UA$Day,
                                         t21=M31UB$Day,t22=M31UC$Day,
                                         t23=M31UD$Day,t24=M31UE$Day,
                                         t25=M31UF$Day,t26=M31UG$Day,
                                         t27=M31UH$Day,t28=R31UA$Day,
                                         t29=R31UB$Day,t30=R31UC$Day,
                                         t31=R31UD$Day,t32=R31UE$Day,
                                         t33=R31UF$Day,
                                         Nasedat1=G21UA$SNF,Nasedat2=G21UB$SNF,
                                         Nasedat3=G21UC$SNF,Nasedat4=M21UA$SNF,
                                         Nasedat5=M21UB$SNF,Nasedat6=M21UC$SNF,
                                         Nasedat7=M21UD$SNF,Nasedat8=R21UA$SNF,
                                         Nasedat9=R21UB$SNF,Nasedat10=R21UC$SNF,
                                         Nasedat11=A31UA$SNF,Nasedat12=A31UB$SNF,
                                         Nasedat13=A31UC$SNF,Nasedat14=G31UA$SNF,
                                         Nasedat15=G31UB$SNF,Nasedat16=G31UC$SNF,
                                         Nasedat17=G31UD$SNF,Nasedat18=G31UE$SNF,
                                         Nasedat19=G31UF$SNF,Nasedat20=M31UA$SNF,
                                         Nasedat21=M31UB$SNF,Nasedat22=M31UC$SNF,
                                         Nasedat23=M31UD$SNF,Nasedat24=M31UE$SNF,
                                         Nasedat25=M31UF$SNF,Nasedat26=M31UG$SNF,
                                         Nasedat27=M31UH$SNF,Nasedat28=R31UA$SNF,
                                         Nasedat29=R31UB$SNF,Nasedat30=R31UC$SNF,
                                         Nasedat31=R31UD$SNF,Nasedat32=R31UE$SNF,
                                         Nasedat33=R31UF$SNF),
                               control=list(maxit=20000))
summary(fit_sigmoid_up_sptemp)

## Symbiosis
sigmoid_up_sym_normNLL <- function(sdNase,rupRhiz,rupActin,
                                      tLupG21A,tLupG21B,tLupG21C,
                                      tLupM21A,tLupM21B,tLupM21C,tLupM21D,
                                      tLupR21A,tLupR21B,tLupR21C,
                                      tLupA31A,tLupA31B,tLupA31C,
                                      tLupG31A,tLupG31B,tLupG31C,tLupG31D,tLupG31E,tLupG31F,
                                      tLupM31A,tLupM31B,tLupM31C,tLupM31D,tLupM31E,tLupM31F,tLupM31G,tLupM31H,
                                      tLupR31A,tLupR31B,tLupR31C,tLupR31D,tLupR31E,tLupR31F,
                                      K1,K2,K3,K4,K5,K6,K7,K8,
                                      K9,K10,K11,K12,K13,K14,K15,K16,
                                      K17,K18,K19,K20,K21,K22,K23,K24,
                                      K25,K26,K27,K28,K29,K30,K31,K32,K33,
                                      t1,t2,t3,t4,t5,t6,t7,t8,
                                      t9,t10,t11,t12,t13,t14,t15,t16,
                                      t17,t18,t19,t20,t21,t22,t23,t24,
                                      t25,t26,t27,t28,t29,t30,t31,t32,t33,
                                      Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                      Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                      Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                      Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32,Nadedat33){
  Nasemean1 <- sigmoid(exp(K1),exp(rupRhiz),exp(tLupG21A),t1)
  Nasemean2 <- sigmoid(exp(K2),exp(rupRhiz),exp(tLupG21B),t2)
  Nasemean3 <- sigmoid(exp(K3),exp(rupRhiz),exp(tLupG21C),t3)
  Nasemean4 <- sigmoid(exp(K4),exp(rupActin),exp(tLupM21A),t4)
  Nasemean5 <- sigmoid(exp(K5),exp(rupActin),exp(tLupM21B),t5)
  Nasemean6 <- sigmoid(exp(K6),exp(rupActin),exp(tLupM21C),t6)
  Nasemean7 <- sigmoid(exp(K7),exp(rupActin),exp(tLupM21D),t7)
  Nasemean8 <- sigmoid(exp(K8),exp(rupRhiz),exp(tLupR21A),t8)
  Nasemean9 <- sigmoid(exp(K9),exp(rupRhiz),exp(tLupR21B),t9)
  Nasemean10 <- sigmoid(exp(K10),exp(rupRhiz),exp(tLupR21C),t10)
  Nasemean11 <- sigmoid(exp(K11),exp(rupActin),exp(tLupA31A),t11)
  Nasemean12 <- sigmoid(exp(K12),exp(rupActin),exp(tLupA31B),t12)
  Nasemean13 <- sigmoid(exp(K13),exp(rupActin),exp(tLupA31C),t13)
  Nasemean14 <- sigmoid(exp(K14),exp(rupRhiz),exp(tLupG31A),t14)
  Nasemean15 <- sigmoid(exp(K15),exp(rupRhiz),exp(tLupG31B),t15)
  Nasemean16 <- sigmoid(exp(K16),exp(rupRhiz),exp(tLupG31C),t16)
  Nasemean17 <- sigmoid(exp(K17),exp(rupRhiz),exp(tLupG31D),t17)
  Nasemean18 <- sigmoid(exp(K18),exp(rupRhiz),exp(tLupG31E),t18)
  Nasemean19 <- sigmoid(exp(K19),exp(rupRhiz),exp(tLupG31F),t19)
  Nasemean20 <- sigmoid(exp(K20),exp(rupActin),exp(tLupM31A),t20)
  Nasemean21 <- sigmoid(exp(K21),exp(rupActin),exp(tLupM31B),t21)
  Nasemean22 <- sigmoid(exp(K22),exp(rupActin),exp(tLupM31C),t22)
  Nasemean23 <- sigmoid(exp(K23),exp(rupActin),exp(tLupM31D),t23)
  Nasemean24 <- sigmoid(exp(K24),exp(rupActin),exp(tLupM31E),t24)
  Nasemean25 <- sigmoid(exp(K25),exp(rupActin),exp(tLupM31F),t25)
  Nasemean26 <- sigmoid(exp(K26),exp(rupActin),exp(tLupM31G),t26)
  Nasemean27 <- sigmoid(exp(K27),exp(rupActin),exp(tLupM31H),t27)
  Nasemean28 <- sigmoid(exp(K28),exp(rupRhiz),exp(tLupR31A),t28)
  Nasemean29 <- sigmoid(exp(K29),exp(rupRhiz),exp(tLupR31B),t29)
  Nasemean30 <- sigmoid(exp(K30),exp(rupRhiz),exp(tLupR31C),t30)
  Nasemean31 <- sigmoid(exp(K31),exp(rupRhiz),exp(tLupR31D),t31)
  Nasemean32 <- sigmoid(exp(K32),exp(rupRhiz),exp(tLupR31E),t32)
  Nasemean33 <- sigmoid(exp(K33),exp(rupRhiz),exp(tLupR31F),t33)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat33,mean=Nasemean33,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_up_sym <- mle2(sigmoid_up_sym_normNLL,start=list(sdNase=-1,
                                                                 K1=log(max(G21UA$SNF)),K2=log(max(G21UB$SNF)),
                                                                 K3=log(max(G21UC$SNF)),K4=log(max(M21UA$SNF)),
                                                                 K5=log(max(M21UB$SNF)),K6=log(max(M21UC$SNF)),
                                                                 K7=log(max(M21UD$SNF)),K8=log(max(R21UA$SNF)),
                                                                 K9=log(max(R21UB$SNF)),K10=log(max(R21UC$SNF)),
                                                                 K11=log(max(A31UA$SNF)),K12=log(max(A31UB$SNF)),
                                                                 K13=log(max(A31UC$SNF)),K14=log(max(G31UA$SNF)),
                                                                 K15=log(max(G31UB$SNF)),K16=log(max(G31UC$SNF)),
                                                                 K17=log(max(G31UD$SNF)),K18=log(max(G31UE$SNF)),
                                                                 K19=log(max(G31UF$SNF)),K20=log(max(M31UA$SNF)),
                                                                 K21=log(max(M31UB$SNF)),K22=log(max(M31UC$SNF)),
                                                                 K23=log(max(M31UD$SNF)),K24=log(max(M31UE$SNF)),
                                                                 K25=log(max(M31UF$SNF)),K26=log(max(M31UG$SNF)),
                                                                 K27=log(max(M31UH$SNF)),K28=log(max(R31UA$SNF)),
                                                                 K29=log(max(R31UB$SNF)),K30=log(max(R31UC$SNF)),
                                                                 K31=log(max(R31UD$SNF)),K32=log(max(R31UE$SNF)),
                                                                 K33=log(max(R31UF$SNF)),
                                                                 rupRhiz=log(0.07),rupActin=log(0.2),
                                                                 tLupG21A=log(70),tLupG21B=log(110),tLupG21C=log(100),
                                                                 tLupM21A=log(150),tLupM21B=log(100),tLupM21C=log(140),tLupM21D=log(180),
                                                                 tLupR21A=log(60),tLupR21B=log(30),tLupR21C=log(40),
                                                                 tLupA31A=log(160),tLupA31B=log(110),tLupA31C=log(110),
                                                                 tLupG31A=log(100),tLupG31B=log(70),tLupG31C=log(110),tLupG31D=log(40),tLupG31E=log(70),tLupG31F=log(110),
                                                                 tLupM31A=log(110),tLupM31B=log(120),tLupM31C=log(150),tLupM31D=log(100),tLupM31E=log(130),tLupM31F=log(130),tLupM31G=log(130),tLupM31H=log(130),
                                                                 tLupR31A=log(70),tLupR31B=log(100),tLupR31C=log(120),tLupR31D=log(30),tLupR31E=log(120),tLupR31F=log(70)),
                            data=list(t1=G21UA$Day,t2=G21UB$Day,
                                      t3=G21UC$Day,t4=M21UA$Day,
                                      t5=M21UB$Day,t6=M21UC$Day,
                                      t7=M21UD$Day,t8=R21UA$Day,
                                      t9=R21UB$Day,t10=R21UC$Day,
                                      t11=A31UA$Day,t12=A31UB$Day,
                                      t13=A31UC$Day,t14=G31UA$Day,
                                      t15=G31UB$Day,t16=G31UC$Day,
                                      t17=G31UD$Day,t18=G31UE$Day,
                                      t19=G31UF$Day,t20=M31UA$Day,
                                      t21=M31UB$Day,t22=M31UC$Day,
                                      t23=M31UD$Day,t24=M31UE$Day,
                                      t25=M31UF$Day,t26=M31UG$Day,
                                      t27=M31UH$Day,t28=R31UA$Day,
                                      t29=R31UB$Day,t30=R31UC$Day,
                                      t31=R31UD$Day,t32=R31UE$Day,
                                      t33=R31UF$Day,
                                      Nasedat1=G21UA$SNF,Nasedat2=G21UB$SNF,
                                      Nasedat3=G21UC$SNF,Nasedat4=M21UA$SNF,
                                      Nasedat5=M21UB$SNF,Nasedat6=M21UC$SNF,
                                      Nasedat7=M21UD$SNF,Nasedat8=R21UA$SNF,
                                      Nasedat9=R21UB$SNF,Nasedat10=R21UC$SNF,
                                      Nasedat11=A31UA$SNF,Nasedat12=A31UB$SNF,
                                      Nasedat13=A31UC$SNF,Nasedat14=G31UA$SNF,
                                      Nasedat15=G31UB$SNF,Nasedat16=G31UC$SNF,
                                      Nasedat17=G31UD$SNF,Nasedat18=G31UE$SNF,
                                      Nasedat19=G31UF$SNF,Nasedat20=M31UA$SNF,
                                      Nasedat21=M31UB$SNF,Nasedat22=M31UC$SNF,
                                      Nasedat23=M31UD$SNF,Nasedat24=M31UE$SNF,
                                      Nasedat25=M31UF$SNF,Nasedat26=M31UG$SNF,
                                      Nasedat27=M31UH$SNF,Nasedat28=R31UA$SNF,
                                      Nasedat29=R31UB$SNF,Nasedat30=R31UC$SNF,
                                      Nasedat31=R31UD$SNF,Nasedat32=R31UE$SNF,
                                      Nasedat33=R31UF$SNF),
                            control=list(maxit=20000))
summary(fit_sigmoid_up_sym)

## Biome
sigmoid_up_biome_normNLL <- function(sdNase,ruptemp,ruptrop,tLupG21A,tLupG21B,tLupG21C,
                                     tLupM21A,tLupM21B,tLupM21C,tLupM21D,
                                     tLupR21A,tLupR21B,tLupR21C,
                                     tLupA31A,tLupA31B,tLupA31C,
                                     tLupG31A,tLupG31B,tLupG31C,tLupG31D,tLupG31E,tLupG31F,
                                     tLupM31A,tLupM31B,tLupM31C,tLupM31D,tLupM31E,tLupM31F,tLupM31G,tLupM31H,
                                     tLupR31A,tLupR31B,tLupR31C,tLupR31D,tLupR31E,tLupR31F,
                                     K1,K2,K3,K4,K5,K6,K7,K8,
                                     K9,K10,K11,K12,K13,K14,K15,K16,
                                     K17,K18,K19,K20,K21,K22,K23,K24,
                                     K25,K26,K27,K28,K29,K30,K31,K32,K33,
                                     t1,t2,t3,t4,t5,t6,t7,t8,
                                     t9,t10,t11,t12,t13,t14,t15,t16,
                                     t17,t18,t19,t20,t21,t22,t23,t24,
                                     t25,t26,t27,t28,t29,t30,t31,t32,t33,
                                     Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                     Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                     Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                     Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32,Nadedat33){
  Nasemean1 <- sigmoid(exp(K1),exp(ruptrop),exp(tLupG21A),t1)
  Nasemean2 <- sigmoid(exp(K2),exp(ruptrop),exp(tLupG21B),t2)
  Nasemean3 <- sigmoid(exp(K3),exp(ruptrop),exp(tLupG21C),t3)
  Nasemean4 <- sigmoid(exp(K4),exp(ruptrop),exp(tLupM21A),t4)
  Nasemean5 <- sigmoid(exp(K5),exp(ruptrop),exp(tLupM21B),t5)
  Nasemean6 <- sigmoid(exp(K6),exp(ruptrop),exp(tLupM21C),t6)
  Nasemean7 <- sigmoid(exp(K7),exp(ruptrop),exp(tLupM21D),t7)
  Nasemean8 <- sigmoid(exp(K8),exp(ruptemp),exp(tLupR21A),t8)
  Nasemean9 <- sigmoid(exp(K9),exp(ruptemp),exp(tLupR21B),t9)
  Nasemean10 <- sigmoid(exp(K10),exp(ruptemp),exp(tLupR21C),t10)
  Nasemean11 <- sigmoid(exp(K11),exp(ruptemp),exp(tLupA31A),t11)
  Nasemean12 <- sigmoid(exp(K12),exp(ruptemp),exp(tLupA31B),t12)
  Nasemean13 <- sigmoid(exp(K13),exp(ruptemp),exp(tLupA31C),t13)
  Nasemean14 <- sigmoid(exp(K14),exp(ruptrop),exp(tLupG31A),t14)
  Nasemean15 <- sigmoid(exp(K15),exp(ruptrop),exp(tLupG31B),t15)
  Nasemean16 <- sigmoid(exp(K16),exp(ruptrop),exp(tLupG31C),t16)
  Nasemean17 <- sigmoid(exp(K17),exp(ruptrop),exp(tLupG31D),t17)
  Nasemean18 <- sigmoid(exp(K18),exp(ruptrop),exp(tLupG31E),t18)
  Nasemean19 <- sigmoid(exp(K19),exp(ruptrop),exp(tLupG31F),t19)
  Nasemean20 <- sigmoid(exp(K20),exp(ruptrop),exp(tLupM31A),t20)
  Nasemean21 <- sigmoid(exp(K21),exp(ruptrop),exp(tLupM31B),t21)
  Nasemean22 <- sigmoid(exp(K22),exp(ruptrop),exp(tLupM31C),t22)
  Nasemean23 <- sigmoid(exp(K23),exp(ruptrop),exp(tLupM31D),t23)
  Nasemean24 <- sigmoid(exp(K24),exp(ruptrop),exp(tLupM31E),t24)
  Nasemean25 <- sigmoid(exp(K25),exp(ruptrop),exp(tLupM31F),t25)
  Nasemean26 <- sigmoid(exp(K26),exp(ruptrop),exp(tLupM31G),t26)
  Nasemean27 <- sigmoid(exp(K27),exp(ruptrop),exp(tLupM31H),t27)
  Nasemean28 <- sigmoid(exp(K28),exp(ruptemp),exp(tLupR31A),t28)
  Nasemean29 <- sigmoid(exp(K29),exp(ruptemp),exp(tLupR31B),t29)
  Nasemean30 <- sigmoid(exp(K30),exp(ruptemp),exp(tLupR31C),t30)
  Nasemean31 <- sigmoid(exp(K31),exp(ruptemp),exp(tLupR31D),t31)
  Nasemean32 <- sigmoid(exp(K32),exp(ruptemp),exp(tLupR31E),t32)
  Nasemean33 <- sigmoid(exp(K33),exp(ruptemp),exp(tLupR31F),t33)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat33,mean=Nasemean33,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_up_biome <- mle2(sigmoid_up_biome_normNLL,start=list(sdNase=-1,
                                                                 K1=log(max(G21UA$SNF)),K2=log(max(G21UB$SNF)),
                                                                 K3=log(max(G21UC$SNF)),K4=log(max(M21UA$SNF)),
                                                                 K5=log(max(M21UB$SNF)),K6=log(max(M21UC$SNF)),
                                                                 K7=log(max(M21UD$SNF)),K8=log(max(R21UA$SNF)),
                                                                 K9=log(max(R21UB$SNF)),K10=log(max(R21UC$SNF)),
                                                                 K11=log(max(A31UA$SNF)),K12=log(max(A31UB$SNF)),
                                                                 K13=log(max(A31UC$SNF)),K14=log(max(G31UA$SNF)),
                                                                 K15=log(max(G31UB$SNF)),K16=log(max(G31UC$SNF)),
                                                                 K17=log(max(G31UD$SNF)),K18=log(max(G31UE$SNF)),
                                                                 K19=log(max(G31UF$SNF)),K20=log(max(M31UA$SNF)),
                                                                 K21=log(max(M31UB$SNF)),K22=log(max(M31UC$SNF)),
                                                                 K23=log(max(M31UD$SNF)),K24=log(max(M31UE$SNF)),
                                                                 K25=log(max(M31UF$SNF)),K26=log(max(M31UG$SNF)),
                                                                 K27=log(max(M31UH$SNF)),K28=log(max(R31UA$SNF)),
                                                                 K29=log(max(R31UB$SNF)),K30=log(max(R31UC$SNF)),
                                                                 K31=log(max(R31UD$SNF)),K32=log(max(R31UE$SNF)),
                                                                 K33=log(max(R31UF$SNF)),
                                                                 ruptrop=log(0.1),ruptemp=log(0.1),tLupG21A=log(70),tLupG21B=log(110),tLupG21C=log(100),
                                                                 tLupM21A=log(150),tLupM21B=log(100),tLupM21C=log(140),tLupM21D=log(180),
                                                                 tLupR21A=log(60),tLupR21B=log(30),tLupR21C=log(40),
                                                                 tLupA31A=log(160),tLupA31B=log(110),tLupA31C=log(110),
                                                                 tLupG31A=log(100),tLupG31B=log(70),tLupG31C=log(110),tLupG31D=log(40),tLupG31E=log(70),tLupG31F=log(110),
                                                                 tLupM31A=log(110),tLupM31B=log(120),tLupM31C=log(150),tLupM31D=log(100),tLupM31E=log(130),tLupM31F=log(130),tLupM31G=log(130),tLupM31H=log(130),
                                                                 tLupR31A=log(70),tLupR31B=log(100),tLupR31C=log(120),tLupR31D=log(30),tLupR31E=log(120),tLupR31F=log(70)),
                             data=list(t1=G21UA$Day,t2=G21UB$Day,
                                       t3=G21UC$Day,t4=M21UA$Day,
                                       t5=M21UB$Day,t6=M21UC$Day,
                                       t7=M21UD$Day,t8=R21UA$Day,
                                       t9=R21UB$Day,t10=R21UC$Day,
                                       t11=A31UA$Day,t12=A31UB$Day,
                                       t13=A31UC$Day,t14=G31UA$Day,
                                       t15=G31UB$Day,t16=G31UC$Day,
                                       t17=G31UD$Day,t18=G31UE$Day,
                                       t19=G31UF$Day,t20=M31UA$Day,
                                       t21=M31UB$Day,t22=M31UC$Day,
                                       t23=M31UD$Day,t24=M31UE$Day,
                                       t25=M31UF$Day,t26=M31UG$Day,
                                       t27=M31UH$Day,t28=R31UA$Day,
                                       t29=R31UB$Day,t30=R31UC$Day,
                                       t31=R31UD$Day,t32=R31UE$Day,
                                       t33=R31UF$Day,
                                       Nasedat1=G21UA$SNF,Nasedat2=G21UB$SNF,
                                       Nasedat3=G21UC$SNF,Nasedat4=M21UA$SNF,
                                       Nasedat5=M21UB$SNF,Nasedat6=M21UC$SNF,
                                       Nasedat7=M21UD$SNF,Nasedat8=R21UA$SNF,
                                       Nasedat9=R21UB$SNF,Nasedat10=R21UC$SNF,
                                       Nasedat11=A31UA$SNF,Nasedat12=A31UB$SNF,
                                       Nasedat13=A31UC$SNF,Nasedat14=G31UA$SNF,
                                       Nasedat15=G31UB$SNF,Nasedat16=G31UC$SNF,
                                       Nasedat17=G31UD$SNF,Nasedat18=G31UE$SNF,
                                       Nasedat19=G31UF$SNF,Nasedat20=M31UA$SNF,
                                       Nasedat21=M31UB$SNF,Nasedat22=M31UC$SNF,
                                       Nasedat23=M31UD$SNF,Nasedat24=M31UE$SNF,
                                       Nasedat25=M31UF$SNF,Nasedat26=M31UG$SNF,
                                       Nasedat27=M31UH$SNF,Nasedat28=R31UA$SNF,
                                       Nasedat29=R31UB$SNF,Nasedat30=R31UC$SNF,
                                       Nasedat31=R31UD$SNF,Nasedat32=R31UE$SNF,
                                       Nasedat33=R31UF$SNF),
                             control=list(maxit=20000))
summary(fit_sigmoid_up_biome)

## Species
sigmoid_up_sp_normNLL <- function(sdNase,rupG,tLupG21A,tLupG21B,tLupG21C,
                                  rupM,tLupM21A,tLupM21B,tLupM21C,tLupM21D,
                                  rupR,tLupR21A,tLupR21B,tLupR21C,
                                  rupA,tLupA31A,tLupA31B,tLupA31C,
                                  tLupG31A,tLupG31B,tLupG31C,tLupG31D,tLupG31E,tLupG31F,
                                  tLupM31A,tLupM31B,tLupM31C,tLupM31D,tLupM31E,tLupM31F,tLupM31G,tLupM31H,
                                  tLupR31A,tLupR31B,tLupR31C,tLupR31D,tLupR31E,tLupR31F,
                                  K1,K2,K3,K4,K5,K6,K7,K8,
                                  K9,K10,K11,K12,K13,K14,K15,K16,
                                  K17,K18,K19,K20,K21,K22,K23,K24,
                                  K25,K26,K27,K28,K29,K30,K31,K32,K33,
                                  t1,t2,t3,t4,t5,t6,t7,t8,
                                  t9,t10,t11,t12,t13,t14,t15,t16,
                                  t17,t18,t19,t20,t21,t22,t23,t24,
                                  t25,t26,t27,t28,t29,t30,t31,t32,t33,
                                  Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                  Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                  Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                  Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32,Nadedat33){
  Nasemean1 <- sigmoid(exp(K1),exp(rupG),exp(tLupG21A),t1)
  Nasemean2 <- sigmoid(exp(K2),exp(rupG),exp(tLupG21B),t2)
  Nasemean3 <- sigmoid(exp(K3),exp(rupG),exp(tLupG21C),t3)
  Nasemean4 <- sigmoid(exp(K4),exp(rupM),exp(tLupM21A),t4)
  Nasemean5 <- sigmoid(exp(K5),exp(rupM),exp(tLupM21B),t5)
  Nasemean6 <- sigmoid(exp(K6),exp(rupM),exp(tLupM21C),t6)
  Nasemean7 <- sigmoid(exp(K7),exp(rupM),exp(tLupM21D),t7)
  Nasemean8 <- sigmoid(exp(K8),exp(rupR),exp(tLupR21A),t8)
  Nasemean9 <- sigmoid(exp(K9),exp(rupR),exp(tLupR21B),t9)
  Nasemean10 <- sigmoid(exp(K10),exp(rupR),exp(tLupR21C),t10)
  Nasemean11 <- sigmoid(exp(K11),exp(rupA),exp(tLupA31A),t11)
  Nasemean12 <- sigmoid(exp(K12),exp(rupA),exp(tLupA31B),t12)
  Nasemean13 <- sigmoid(exp(K13),exp(rupA),exp(tLupA31C),t13)
  Nasemean14 <- sigmoid(exp(K14),exp(rupG),exp(tLupG31A),t14)
  Nasemean15 <- sigmoid(exp(K15),exp(rupG),exp(tLupG31B),t15)
  Nasemean16 <- sigmoid(exp(K16),exp(rupG),exp(tLupG31C),t16)
  Nasemean17 <- sigmoid(exp(K17),exp(rupG),exp(tLupG31D),t17)
  Nasemean18 <- sigmoid(exp(K18),exp(rupG),exp(tLupG31E),t18)
  Nasemean19 <- sigmoid(exp(K19),exp(rupG),exp(tLupG31F),t19)
  Nasemean20 <- sigmoid(exp(K20),exp(rupM),exp(tLupM31A),t20)
  Nasemean21 <- sigmoid(exp(K21),exp(rupM),exp(tLupM31B),t21)
  Nasemean22 <- sigmoid(exp(K22),exp(rupM),exp(tLupM31C),t22)
  Nasemean23 <- sigmoid(exp(K23),exp(rupM),exp(tLupM31D),t23)
  Nasemean24 <- sigmoid(exp(K24),exp(rupM),exp(tLupM31E),t24)
  Nasemean25 <- sigmoid(exp(K25),exp(rupM),exp(tLupM31F),t25)
  Nasemean26 <- sigmoid(exp(K26),exp(rupM),exp(tLupM31G),t26)
  Nasemean27 <- sigmoid(exp(K27),exp(rupM),exp(tLupM31H),t27)
  Nasemean28 <- sigmoid(exp(K28),exp(rupR),exp(tLupR31A),t28)
  Nasemean29 <- sigmoid(exp(K29),exp(rupR),exp(tLupR31B),t29)
  Nasemean30 <- sigmoid(exp(K30),exp(rupR),exp(tLupR31C),t30)
  Nasemean31 <- sigmoid(exp(K31),exp(rupR),exp(tLupR31D),t31)
  Nasemean32 <- sigmoid(exp(K32),exp(rupR),exp(tLupR31E),t32)
  Nasemean33 <- sigmoid(exp(K33),exp(rupR),exp(tLupR31F),t33)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat33,mean=Nasemean33,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_up_sp <- mle2(sigmoid_up_sp_normNLL,start=list(sdNase=-1,
                                                           K1=log(max(G21UA$SNF)),K2=log(max(G21UB$SNF)),
                                                           K3=log(max(G21UC$SNF)),K4=log(max(M21UA$SNF)),
                                                           K5=log(max(M21UB$SNF)),K6=log(max(M21UC$SNF)),
                                                           K7=log(max(M21UD$SNF)),K8=log(max(R21UA$SNF)),
                                                           K9=log(max(R21UB$SNF)),K10=log(max(R21UC$SNF)),
                                                           K11=log(max(A31UA$SNF)),K12=log(max(A31UB$SNF)),
                                                           K13=log(max(A31UC$SNF)),K14=log(max(G31UA$SNF)),
                                                           K15=log(max(G31UB$SNF)),K16=log(max(G31UC$SNF)),
                                                           K17=log(max(G31UD$SNF)),K18=log(max(G31UE$SNF)),
                                                           K19=log(max(G31UF$SNF)),K20=log(max(M31UA$SNF)),
                                                           K21=log(max(M31UB$SNF)),K22=log(max(M31UC$SNF)),
                                                           K23=log(max(M31UD$SNF)),K24=log(max(M31UE$SNF)),
                                                           K25=log(max(M31UF$SNF)),K26=log(max(M31UG$SNF)),
                                                           K27=log(max(M31UH$SNF)),K28=log(max(R31UA$SNF)),
                                                           K29=log(max(R31UB$SNF)),K30=log(max(R31UC$SNF)),
                                                           K31=log(max(R31UD$SNF)),K32=log(max(R31UE$SNF)),
                                                           K33=log(max(R31UF$SNF)),
                                                           rupG=log(0.1),tLupG21A=log(70),tLupG21B=log(110),tLupG21C=log(100),
                                                           rupM=log(0.1),tLupM21A=log(150),tLupM21B=log(100),tLupM21C=log(140),tLupM21D=log(180),
                                                           rupR=log(0.1),tLupR21A=log(60),tLupR21B=log(30),tLupR21C=log(40),
                                                           rupA=log(0.1),tLupA31A=log(160),tLupA31B=log(110),tLupA31C=log(110),
                                                           tLupG31A=log(100),tLupG31B=log(70),tLupG31C=log(110),tLupG31D=log(40),tLupG31E=log(70),tLupG31F=log(110),
                                                           tLupM31A=log(110),tLupM31B=log(120),tLupM31C=log(150),tLupM31D=log(100),tLupM31E=log(130),tLupM31F=log(130),tLupM31G=log(130),tLupM31H=log(130),
                                                           tLupR31A=log(70),tLupR31B=log(100),tLupR31C=log(120),tLupR31D=log(30),tLupR31E=log(120),tLupR31F=log(70)),
                          data=list(t1=G21UA$Day,t2=G21UB$Day,
                                    t3=G21UC$Day,t4=M21UA$Day,
                                    t5=M21UB$Day,t6=M21UC$Day,
                                    t7=M21UD$Day,t8=R21UA$Day,
                                    t9=R21UB$Day,t10=R21UC$Day,
                                    t11=A31UA$Day,t12=A31UB$Day,
                                    t13=A31UC$Day,t14=G31UA$Day,
                                    t15=G31UB$Day,t16=G31UC$Day,
                                    t17=G31UD$Day,t18=G31UE$Day,
                                    t19=G31UF$Day,t20=M31UA$Day,
                                    t21=M31UB$Day,t22=M31UC$Day,
                                    t23=M31UD$Day,t24=M31UE$Day,
                                    t25=M31UF$Day,t26=M31UG$Day,
                                    t27=M31UH$Day,t28=R31UA$Day,
                                    t29=R31UB$Day,t30=R31UC$Day,
                                    t31=R31UD$Day,t32=R31UE$Day,
                                    t33=R31UF$Day,
                                    Nasedat1=G21UA$SNF,Nasedat2=G21UB$SNF,
                                    Nasedat3=G21UC$SNF,Nasedat4=M21UA$SNF,
                                    Nasedat5=M21UB$SNF,Nasedat6=M21UC$SNF,
                                    Nasedat7=M21UD$SNF,Nasedat8=R21UA$SNF,
                                    Nasedat9=R21UB$SNF,Nasedat10=R21UC$SNF,
                                    Nasedat11=A31UA$SNF,Nasedat12=A31UB$SNF,
                                    Nasedat13=A31UC$SNF,Nasedat14=G31UA$SNF,
                                    Nasedat15=G31UB$SNF,Nasedat16=G31UC$SNF,
                                    Nasedat17=G31UD$SNF,Nasedat18=G31UE$SNF,
                                    Nasedat19=G31UF$SNF,Nasedat20=M31UA$SNF,
                                    Nasedat21=M31UB$SNF,Nasedat22=M31UC$SNF,
                                    Nasedat23=M31UD$SNF,Nasedat24=M31UE$SNF,
                                    Nasedat25=M31UF$SNF,Nasedat26=M31UG$SNF,
                                    Nasedat27=M31UH$SNF,Nasedat28=R31UA$SNF,
                                    Nasedat29=R31UB$SNF,Nasedat30=R31UC$SNF,
                                    Nasedat31=R31UD$SNF,Nasedat32=R31UE$SNF,
                                    Nasedat33=R31UF$SNF),
                          control=list(maxit=20000))
summary(fit_sigmoid_up_sp)

## Temperature
sigmoid_up_temp_normNLL <- function(sdNase,rup1,rup31,tLupG21A,tLupG21B,tLupG21C,
                                    tLupM21A,tLupM21B,tLupM21C,tLupM21D,
                                    tLupR21A,tLupR21B,tLupR21C,
                                    tLupA31A,tLupA31B,tLupA31C,
                                    tLupG31A,tLupG31B,tLupG31C,tLupG31D,tLupG31E,tLupG31F,
                                    tLupM31A,tLupM31B,tLupM31C,tLupM31D,tLupM31E,tLupM31F,tLupM31G,tLupM31H,
                                    tLupR31A,tLupR31B,tLupR31C,tLupR31D,tLupR31E,tLupR31F,
                                    K1,K2,K3,K4,K5,K6,K7,K8,
                                    K9,K10,K11,K12,K13,K14,K15,K16,
                                    K17,K18,K19,K20,K21,K22,K23,K24,
                                    K25,K26,K27,K28,K29,K30,K31,K32,K33,
                                    t1,t2,t3,t4,t5,t6,t7,t8,
                                    t9,t10,t11,t12,t13,t14,t15,t16,
                                    t17,t18,t19,t20,t21,t22,t23,t24,
                                    t25,t26,t27,t28,t29,t30,t31,t32,t33,
                                    Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                    Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                    Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                    Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32,Nadedat33){
  Nasemean1 <- sigmoid(exp(K1),exp(rup1),exp(tLupG21A),t1)
  Nasemean2 <- sigmoid(exp(K2),exp(rup1),exp(tLupG21B),t2)
  Nasemean3 <- sigmoid(exp(K3),exp(rup1),exp(tLupG21C),t3)
  Nasemean4 <- sigmoid(exp(K4),exp(rup1),exp(tLupM21A),t4)
  Nasemean5 <- sigmoid(exp(K5),exp(rup1),exp(tLupM21B),t5)
  Nasemean6 <- sigmoid(exp(K6),exp(rup1),exp(tLupM21C),t6)
  Nasemean7 <- sigmoid(exp(K7),exp(rup1),exp(tLupM21D),t7)
  Nasemean8 <- sigmoid(exp(K8),exp(rup1),exp(tLupR21A),t8)
  Nasemean9 <- sigmoid(exp(K9),exp(rup1),exp(tLupR21B),t9)
  Nasemean10 <- sigmoid(exp(K10),exp(rup1),exp(tLupR21C),t10)
  Nasemean11 <- sigmoid(exp(K11),exp(rup31),exp(tLupA31A),t11)
  Nasemean12 <- sigmoid(exp(K12),exp(rup31),exp(tLupA31B),t12)
  Nasemean13 <- sigmoid(exp(K13),exp(rup31),exp(tLupA31C),t13)
  Nasemean14 <- sigmoid(exp(K14),exp(rup31),exp(tLupG31A),t14)
  Nasemean15 <- sigmoid(exp(K15),exp(rup31),exp(tLupG31B),t15)
  Nasemean16 <- sigmoid(exp(K16),exp(rup31),exp(tLupG31C),t16)
  Nasemean17 <- sigmoid(exp(K17),exp(rup31),exp(tLupG31D),t17)
  Nasemean18 <- sigmoid(exp(K18),exp(rup31),exp(tLupG31E),t18)
  Nasemean19 <- sigmoid(exp(K19),exp(rup31),exp(tLupG31F),t19)
  Nasemean20 <- sigmoid(exp(K20),exp(rup31),exp(tLupM31A),t20)
  Nasemean21 <- sigmoid(exp(K21),exp(rup31),exp(tLupM31B),t21)
  Nasemean22 <- sigmoid(exp(K22),exp(rup31),exp(tLupM31C),t22)
  Nasemean23 <- sigmoid(exp(K23),exp(rup31),exp(tLupM31D),t23)
  Nasemean24 <- sigmoid(exp(K24),exp(rup31),exp(tLupM31E),t24)
  Nasemean25 <- sigmoid(exp(K25),exp(rup31),exp(tLupM31F),t25)
  Nasemean26 <- sigmoid(exp(K26),exp(rup31),exp(tLupM31G),t26)
  Nasemean27 <- sigmoid(exp(K27),exp(rup31),exp(tLupM31H),t27)
  Nasemean28 <- sigmoid(exp(K28),exp(rup31),exp(tLupR31A),t28)
  Nasemean29 <- sigmoid(exp(K29),exp(rup31),exp(tLupR31B),t29)
  Nasemean30 <- sigmoid(exp(K30),exp(rup31),exp(tLupR31C),t30)
  Nasemean31 <- sigmoid(exp(K31),exp(rup31),exp(tLupR31D),t31)
  Nasemean32 <- sigmoid(exp(K32),exp(rup31),exp(tLupR31E),t32)
  Nasemean33 <- sigmoid(exp(K33),exp(rup31),exp(tLupR31F),t33)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat33,mean=Nasemean33,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_up_temp <- mle2(sigmoid_up_temp_normNLL,start=list(sdNase=-1,
                                                               K1=log(max(G21UA$SNF)),K2=log(max(G21UB$SNF)),
                                                               K3=log(max(G21UC$SNF)),K4=log(max(M21UA$SNF)),
                                                               K5=log(max(M21UB$SNF)),K6=log(max(M21UC$SNF)),
                                                               K7=log(max(M21UD$SNF)),K8=log(max(R21UA$SNF)),
                                                               K9=log(max(R21UB$SNF)),K10=log(max(R21UC$SNF)),
                                                               K11=log(max(A31UA$SNF)),K12=log(max(A31UB$SNF)),
                                                               K13=log(max(A31UC$SNF)),K14=log(max(G31UA$SNF)),
                                                               K15=log(max(G31UB$SNF)),K16=log(max(G31UC$SNF)),
                                                               K17=log(max(G31UD$SNF)),K18=log(max(G31UE$SNF)),
                                                               K19=log(max(G31UF$SNF)),K20=log(max(M31UA$SNF)),
                                                               K21=log(max(M31UB$SNF)),K22=log(max(M31UC$SNF)),
                                                               K23=log(max(M31UD$SNF)),K24=log(max(M31UE$SNF)),
                                                               K25=log(max(M31UF$SNF)),K26=log(max(M31UG$SNF)),
                                                               K27=log(max(M31UH$SNF)),K28=log(max(R31UA$SNF)),
                                                               K29=log(max(R31UB$SNF)),K30=log(max(R31UC$SNF)),
                                                               K31=log(max(R31UD$SNF)),K32=log(max(R31UE$SNF)),
                                                               K33=log(max(R31UF$SNF)),
                                                               rup1=log(0.1),rup31=log(0.1),tLupG21A=log(70),tLupG21B=log(110),tLupG21C=log(100),
                                                               tLupM21A=log(150),tLupM21B=log(100),tLupM21C=log(140),tLupM21D=log(180),
                                                               tLupR21A=log(60),tLupR21B=log(30),tLupR21C=log(40),
                                                               tLupA31A=log(160),tLupA31B=log(110),tLupA31C=log(110),
                                                               tLupG31A=log(100),tLupG31B=log(70),tLupG31C=log(110),tLupG31D=log(40),tLupG31E=log(70),tLupG31F=log(110),
                                                               tLupM31A=log(110),tLupM31B=log(120),tLupM31C=log(150),tLupM31D=log(100),tLupM31E=log(130),tLupM31F=log(130),tLupM31G=log(130),tLupM31H=log(130),
                                                               tLupR31A=log(70),tLupR31B=log(100),tLupR31C=log(120),tLupR31D=log(30),tLupR31E=log(120),tLupR31F=log(70)),
                            data=list(t1=G21UA$Day,t2=G21UB$Day,
                                      t3=G21UC$Day,t4=M21UA$Day,
                                      t5=M21UB$Day,t6=M21UC$Day,
                                      t7=M21UD$Day,t8=R21UA$Day,
                                      t9=R21UB$Day,t10=R21UC$Day,
                                      t11=A31UA$Day,t12=A31UB$Day,
                                      t13=A31UC$Day,t14=G31UA$Day,
                                      t15=G31UB$Day,t16=G31UC$Day,
                                      t17=G31UD$Day,t18=G31UE$Day,
                                      t19=G31UF$Day,t20=M31UA$Day,
                                      t21=M31UB$Day,t22=M31UC$Day,
                                      t23=M31UD$Day,t24=M31UE$Day,
                                      t25=M31UF$Day,t26=M31UG$Day,
                                      t27=M31UH$Day,t28=R31UA$Day,
                                      t29=R31UB$Day,t30=R31UC$Day,
                                      t31=R31UD$Day,t32=R31UE$Day,
                                      t33=R31UF$Day,
                                      Nasedat1=G21UA$SNF,Nasedat2=G21UB$SNF,
                                      Nasedat3=G21UC$SNF,Nasedat4=M21UA$SNF,
                                      Nasedat5=M21UB$SNF,Nasedat6=M21UC$SNF,
                                      Nasedat7=M21UD$SNF,Nasedat8=R21UA$SNF,
                                      Nasedat9=R21UB$SNF,Nasedat10=R21UC$SNF,
                                      Nasedat11=A31UA$SNF,Nasedat12=A31UB$SNF,
                                      Nasedat13=A31UC$SNF,Nasedat14=G31UA$SNF,
                                      Nasedat15=G31UB$SNF,Nasedat16=G31UC$SNF,
                                      Nasedat17=G31UD$SNF,Nasedat18=G31UE$SNF,
                                      Nasedat19=G31UF$SNF,Nasedat20=M31UA$SNF,
                                      Nasedat21=M31UB$SNF,Nasedat22=M31UC$SNF,
                                      Nasedat23=M31UD$SNF,Nasedat24=M31UE$SNF,
                                      Nasedat25=M31UF$SNF,Nasedat26=M31UG$SNF,
                                      Nasedat27=M31UH$SNF,Nasedat28=R31UA$SNF,
                                      Nasedat29=R31UB$SNF,Nasedat30=R31UC$SNF,
                                      Nasedat31=R31UD$SNF,Nasedat32=R31UE$SNF,
                                      Nasedat33=R31UF$SNF),
                            control=list(maxit=20000))
summary(fit_sigmoid_up_temp)

## Biome x temperature interaction
sigmoid_up_biometemp_normNLL <- function(sdNase,ruptrop21,ruptrop31,ruptemp21,ruptemp31,
                                         tLupG21A,tLupG21B,tLupG21C,
                                         tLupM21A,tLupM21B,tLupM21C,tLupM21D,
                                         tLupR21A,tLupR21B,tLupR21C,
                                         tLupA31A,tLupA31B,tLupA31C,
                                         tLupG31A,tLupG31B,tLupG31C,tLupG31D,tLupG31E,tLupG31F,
                                         tLupM31A,tLupM31B,tLupM31C,tLupM31D,tLupM31E,tLupM31F,tLupM31G,tLupM31H,
                                         tLupR31A,tLupR31B,tLupR31C,tLupR31D,tLupR31E,tLupR31F,
                                         K1,K2,K3,K4,K5,K6,K7,K8,
                                         K9,K10,K11,K12,K13,K14,K15,K16,
                                         K17,K18,K19,K20,K21,K22,K23,K24,
                                         K25,K26,K27,K28,K29,K30,K31,K32,K33,
                                         t1,t2,t3,t4,t5,t6,t7,t8,
                                         t9,t10,t11,t12,t13,t14,t15,t16,
                                         t17,t18,t19,t20,t21,t22,t23,t24,
                                         t25,t26,t27,t28,t29,t30,t31,t32,t33,
                                         Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                         Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                         Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                         Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32,Nadedat33){
  Nasemean1 <- sigmoid(exp(K1),exp(ruptrop21),exp(tLupG21A),t1)
  Nasemean2 <- sigmoid(exp(K2),exp(ruptrop21),exp(tLupG21B),t2)
  Nasemean3 <- sigmoid(exp(K3),exp(ruptrop21),exp(tLupG21C),t3)
  Nasemean4 <- sigmoid(exp(K4),exp(ruptrop21),exp(tLupM21A),t4)
  Nasemean5 <- sigmoid(exp(K5),exp(ruptrop21),exp(tLupM21B),t5)
  Nasemean6 <- sigmoid(exp(K6),exp(ruptrop21),exp(tLupM21C),t6)
  Nasemean7 <- sigmoid(exp(K7),exp(ruptrop21),exp(tLupM21D),t7)
  Nasemean8 <- sigmoid(exp(K8),exp(ruptemp21),exp(tLupR21A),t8)
  Nasemean9 <- sigmoid(exp(K9),exp(ruptemp21),exp(tLupR21B),t9)
  Nasemean10 <- sigmoid(exp(K10),exp(ruptemp21),exp(tLupR21C),t10)
  Nasemean11 <- sigmoid(exp(K11),exp(ruptemp31),exp(tLupA31A),t11)
  Nasemean12 <- sigmoid(exp(K12),exp(ruptemp31),exp(tLupA31B),t12)
  Nasemean13 <- sigmoid(exp(K13),exp(ruptemp31),exp(tLupA31C),t13)
  Nasemean14 <- sigmoid(exp(K14),exp(ruptrop31),exp(tLupG31A),t14)
  Nasemean15 <- sigmoid(exp(K15),exp(ruptrop31),exp(tLupG31B),t15)
  Nasemean16 <- sigmoid(exp(K16),exp(ruptrop31),exp(tLupG31C),t16)
  Nasemean17 <- sigmoid(exp(K17),exp(ruptrop31),exp(tLupG31D),t17)
  Nasemean18 <- sigmoid(exp(K18),exp(ruptrop31),exp(tLupG31E),t18)
  Nasemean19 <- sigmoid(exp(K19),exp(ruptrop31),exp(tLupG31F),t19)
  Nasemean20 <- sigmoid(exp(K20),exp(ruptrop31),exp(tLupM31A),t20)
  Nasemean21 <- sigmoid(exp(K21),exp(ruptrop31),exp(tLupM31B),t21)
  Nasemean22 <- sigmoid(exp(K22),exp(ruptrop31),exp(tLupM31C),t22)
  Nasemean23 <- sigmoid(exp(K23),exp(ruptrop31),exp(tLupM31D),t23)
  Nasemean24 <- sigmoid(exp(K24),exp(ruptrop31),exp(tLupM31E),t24)
  Nasemean25 <- sigmoid(exp(K25),exp(ruptrop31),exp(tLupM31F),t25)
  Nasemean26 <- sigmoid(exp(K26),exp(ruptrop31),exp(tLupM31G),t26)
  Nasemean27 <- sigmoid(exp(K27),exp(ruptrop31),exp(tLupM31H),t27)
  Nasemean28 <- sigmoid(exp(K28),exp(ruptemp31),exp(tLupR31A),t28)
  Nasemean29 <- sigmoid(exp(K29),exp(ruptemp31),exp(tLupR31B),t29)
  Nasemean30 <- sigmoid(exp(K30),exp(ruptemp31),exp(tLupR31C),t30)
  Nasemean31 <- sigmoid(exp(K31),exp(ruptemp31),exp(tLupR31D),t31)
  Nasemean32 <- sigmoid(exp(K32),exp(ruptemp31),exp(tLupR31E),t32)
  Nasemean33 <- sigmoid(exp(K33),exp(ruptemp31),exp(tLupR31F),t33)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat33,mean=Nasemean33,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_up_biometemp <- mle2(sigmoid_up_biometemp_normNLL,start=list(sdNase=-1,
                                                                         K1=log(max(G21UA$SNF)),K2=log(max(G21UB$SNF)),
                                                                         K3=log(max(G21UC$SNF)),K4=log(max(M21UA$SNF)),
                                                                         K5=log(max(M21UB$SNF)),K6=log(max(M21UC$SNF)),
                                                                         K7=log(max(M21UD$SNF)),K8=log(max(R21UA$SNF)),
                                                                         K9=log(max(R21UB$SNF)),K10=log(max(R21UC$SNF)),
                                                                         K11=log(max(A31UA$SNF)),K12=log(max(A31UB$SNF)),
                                                                         K13=log(max(A31UC$SNF)),K14=log(max(G31UA$SNF)),
                                                                         K15=log(max(G31UB$SNF)),K16=log(max(G31UC$SNF)),
                                                                         K17=log(max(G31UD$SNF)),K18=log(max(G31UE$SNF)),
                                                                         K19=log(max(G31UF$SNF)),K20=log(max(M31UA$SNF)),
                                                                         K21=log(max(M31UB$SNF)),K22=log(max(M31UC$SNF)),
                                                                         K23=log(max(M31UD$SNF)),K24=log(max(M31UE$SNF)),
                                                                         K25=log(max(M31UF$SNF)),K26=log(max(M31UG$SNF)),
                                                                         K27=log(max(M31UH$SNF)),K28=log(max(R31UA$SNF)),
                                                                         K29=log(max(R31UB$SNF)),K30=log(max(R31UC$SNF)),
                                                                         K31=log(max(R31UD$SNF)),K32=log(max(R31UE$SNF)),
                                                                         K33=log(max(R31UF$SNF)),
                                                                         ruptrop21=log(0.1),ruptrop31=log(0.1),ruptemp21=log(0.1),ruptemp31=log(0.1),
                                                                         tLupG21A=log(70),tLupG21B=log(110),tLupG21C=log(100),
                                                                         tLupM21A=log(150),tLupM21B=log(100),tLupM21C=log(140),tLupM21D=log(180),
                                                                         tLupR21A=log(60),tLupR21B=log(30),tLupR21C=log(40),
                                                                         tLupA31A=log(160),tLupA31B=log(110),tLupA31C=log(110),
                                                                         tLupG31A=log(100),tLupG31B=log(70),tLupG31C=log(110),tLupG31D=log(40),tLupG31E=log(70),tLupG31F=log(110),
                                                                         tLupM31A=log(110),tLupM31B=log(120),tLupM31C=log(150),tLupM31D=log(100),tLupM31E=log(130),tLupM31F=log(130),tLupM31G=log(130),tLupM31H=log(130),
                                                                         tLupR31A=log(70),tLupR31B=log(100),tLupR31C=log(120),tLupR31D=log(30),tLupR31E=log(120),tLupR31F=log(70)),
                                 data=list(t1=G21UA$Day,t2=G21UB$Day,
                                           t3=G21UC$Day,t4=M21UA$Day,
                                           t5=M21UB$Day,t6=M21UC$Day,
                                           t7=M21UD$Day,t8=R21UA$Day,
                                           t9=R21UB$Day,t10=R21UC$Day,
                                           t11=A31UA$Day,t12=A31UB$Day,
                                           t13=A31UC$Day,t14=G31UA$Day,
                                           t15=G31UB$Day,t16=G31UC$Day,
                                           t17=G31UD$Day,t18=G31UE$Day,
                                           t19=G31UF$Day,t20=M31UA$Day,
                                           t21=M31UB$Day,t22=M31UC$Day,
                                           t23=M31UD$Day,t24=M31UE$Day,
                                           t25=M31UF$Day,t26=M31UG$Day,
                                           t27=M31UH$Day,t28=R31UA$Day,
                                           t29=R31UB$Day,t30=R31UC$Day,
                                           t31=R31UD$Day,t32=R31UE$Day,
                                           t33=R31UF$Day,
                                           Nasedat1=G21UA$SNF,Nasedat2=G21UB$SNF,
                                           Nasedat3=G21UC$SNF,Nasedat4=M21UA$SNF,
                                           Nasedat5=M21UB$SNF,Nasedat6=M21UC$SNF,
                                           Nasedat7=M21UD$SNF,Nasedat8=R21UA$SNF,
                                           Nasedat9=R21UB$SNF,Nasedat10=R21UC$SNF,
                                           Nasedat11=A31UA$SNF,Nasedat12=A31UB$SNF,
                                           Nasedat13=A31UC$SNF,Nasedat14=G31UA$SNF,
                                           Nasedat15=G31UB$SNF,Nasedat16=G31UC$SNF,
                                           Nasedat17=G31UD$SNF,Nasedat18=G31UE$SNF,
                                           Nasedat19=G31UF$SNF,Nasedat20=M31UA$SNF,
                                           Nasedat21=M31UB$SNF,Nasedat22=M31UC$SNF,
                                           Nasedat23=M31UD$SNF,Nasedat24=M31UE$SNF,
                                           Nasedat25=M31UF$SNF,Nasedat26=M31UG$SNF,
                                           Nasedat27=M31UH$SNF,Nasedat28=R31UA$SNF,
                                           Nasedat29=R31UB$SNF,Nasedat30=R31UC$SNF,
                                           Nasedat31=R31UD$SNF,Nasedat32=R31UE$SNF,
                                           Nasedat33=R31UF$SNF),
                                 control=list(maxit=20000))
summary(fit_sigmoid_up_biometemp)

## Symbiosis x temperature interaction
sigmoid_up_symtemp_normNLL <- function(sdNase,ruprhiz21,ruprhiz31,rupactin21,rupactin31,
                                       tLupG21A,tLupG21B,tLupG21C,
                                       tLupM21A,tLupM21B,tLupM21C,tLupM21D,
                                       tLupR21A,tLupR21B,tLupR21C,
                                       tLupA31A,tLupA31B,tLupA31C,
                                       tLupG31A,tLupG31B,tLupG31C,tLupG31D,tLupG31E,tLupG31F,
                                       tLupM31A,tLupM31B,tLupM31C,tLupM31D,tLupM31E,tLupM31F,tLupM31G,tLupM31H,
                                       tLupR31A,tLupR31B,tLupR31C,tLupR31D,tLupR31E,tLupR31F,
                                       K1,K2,K3,K4,K5,K6,K7,K8,
                                       K9,K10,K11,K12,K13,K14,K15,K16,
                                       K17,K18,K19,K20,K21,K22,K23,K24,
                                       K25,K26,K27,K28,K29,K30,K31,K32,K33,
                                       t1,t2,t3,t4,t5,t6,t7,t8,
                                       t9,t10,t11,t12,t13,t14,t15,t16,
                                       t17,t18,t19,t20,t21,t22,t23,t24,
                                       t25,t26,t27,t28,t29,t30,t31,t32,t33,
                                       Nasedat1,Nasedat2,Nasedat3,Nasedat4,Nasedat5,Nasedat6,Nasedat7,Nasedat8,
                                       Nasedat9,Nasedat10,Nasedat11,Nasedat12,Nasedat13,Nasedat14,Nasedat15,Nasedat16,
                                       Nadedat17,Nadedat18,Nadedat19,Nadedat20,Nadedat21,Nadedat22,Nadedat23,Nadedat24,
                                       Nadedat25,Nadedat26,Nadedat27,Nadedat28,Nadedat29,Nadedat30,Nadedat31,Nadedat32,Nadedat33){
  Nasemean1 <- sigmoid(exp(K1),exp(ruprhiz21),exp(tLupG21A),t1)
  Nasemean2 <- sigmoid(exp(K2),exp(ruprhiz21),exp(tLupG21B),t2)
  Nasemean3 <- sigmoid(exp(K3),exp(ruprhiz21),exp(tLupG21C),t3)
  Nasemean4 <- sigmoid(exp(K4),exp(rupactin21),exp(tLupM21A),t4)
  Nasemean5 <- sigmoid(exp(K5),exp(rupactin21),exp(tLupM21B),t5)
  Nasemean6 <- sigmoid(exp(K6),exp(rupactin21),exp(tLupM21C),t6)
  Nasemean7 <- sigmoid(exp(K7),exp(rupactin21),exp(tLupM21D),t7)
  Nasemean8 <- sigmoid(exp(K8),exp(ruprhiz21),exp(tLupR21A),t8)
  Nasemean9 <- sigmoid(exp(K9),exp(ruprhiz21),exp(tLupR21B),t9)
  Nasemean10 <- sigmoid(exp(K10),exp(ruprhiz21),exp(tLupR21C),t10)
  Nasemean11 <- sigmoid(exp(K11),exp(rupactin31),exp(tLupA31A),t11)
  Nasemean12 <- sigmoid(exp(K12),exp(rupactin31),exp(tLupA31B),t12)
  Nasemean13 <- sigmoid(exp(K13),exp(rupactin31),exp(tLupA31C),t13)
  Nasemean14 <- sigmoid(exp(K14),exp(ruprhiz31),exp(tLupG31A),t14)
  Nasemean15 <- sigmoid(exp(K15),exp(ruprhiz31),exp(tLupG31B),t15)
  Nasemean16 <- sigmoid(exp(K16),exp(ruprhiz31),exp(tLupG31C),t16)
  Nasemean17 <- sigmoid(exp(K17),exp(ruprhiz31),exp(tLupG31D),t17)
  Nasemean18 <- sigmoid(exp(K18),exp(ruprhiz31),exp(tLupG31E),t18)
  Nasemean19 <- sigmoid(exp(K19),exp(ruprhiz31),exp(tLupG31F),t19)
  Nasemean20 <- sigmoid(exp(K20),exp(rupactin31),exp(tLupM31A),t20)
  Nasemean21 <- sigmoid(exp(K21),exp(rupactin31),exp(tLupM31B),t21)
  Nasemean22 <- sigmoid(exp(K22),exp(rupactin31),exp(tLupM31C),t22)
  Nasemean23 <- sigmoid(exp(K23),exp(rupactin31),exp(tLupM31D),t23)
  Nasemean24 <- sigmoid(exp(K24),exp(rupactin31),exp(tLupM31E),t24)
  Nasemean25 <- sigmoid(exp(K25),exp(rupactin31),exp(tLupM31F),t25)
  Nasemean26 <- sigmoid(exp(K26),exp(rupactin31),exp(tLupM31G),t26)
  Nasemean27 <- sigmoid(exp(K27),exp(rupactin31),exp(tLupM31H),t27)
  Nasemean28 <- sigmoid(exp(K28),exp(ruprhiz31),exp(tLupR31A),t28)
  Nasemean29 <- sigmoid(exp(K29),exp(ruprhiz31),exp(tLupR31B),t29)
  Nasemean30 <- sigmoid(exp(K30),exp(ruprhiz31),exp(tLupR31C),t30)
  Nasemean31 <- sigmoid(exp(K31),exp(ruprhiz31),exp(tLupR31D),t31)
  Nasemean32 <- sigmoid(exp(K32),exp(ruprhiz31),exp(tLupR31E),t32)
  Nasemean33 <- sigmoid(exp(K33),exp(ruprhiz31),exp(tLupR31F),t33)
  -(sum(dnorm(Nasedat1,mean=Nasemean1,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat2,mean=Nasemean2,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat3,mean=Nasemean3,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat4,mean=Nasemean4,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat5,mean=Nasemean5,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat6,mean=Nasemean6,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat7,mean=Nasemean7,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat8,mean=Nasemean8,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat9,mean=Nasemean9,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat10,mean=Nasemean10,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat11,mean=Nasemean11,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat12,mean=Nasemean12,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat13,mean=Nasemean13,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat14,mean=Nasemean14,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat15,mean=Nasemean15,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat16,mean=Nasemean16,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat17,mean=Nasemean17,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat18,mean=Nasemean18,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat19,mean=Nasemean19,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat20,mean=Nasemean20,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat21,mean=Nasemean21,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat22,mean=Nasemean22,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat23,mean=Nasemean23,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat24,mean=Nasemean24,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat25,mean=Nasemean25,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26,mean=Nasemean26,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat27,mean=Nasemean27,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat28,mean=Nasemean28,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat29,mean=Nasemean29,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat30,mean=Nasemean30,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31,mean=Nasemean31,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat32,mean=Nasemean32,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat33,mean=Nasemean33,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

fit_sigmoid_up_symtemp <- mle2(sigmoid_up_symtemp_normNLL,start=list(sdNase=-1,
                                                                     K1=log(max(G21UA$SNF)),K2=log(max(G21UB$SNF)),
                                                                     K3=log(max(G21UC$SNF)),K4=log(max(M21UA$SNF)),
                                                                     K5=log(max(M21UB$SNF)),K6=log(max(M21UC$SNF)),
                                                                     K7=log(max(M21UD$SNF)),K8=log(max(R21UA$SNF)),
                                                                     K9=log(max(R21UB$SNF)),K10=log(max(R21UC$SNF)),
                                                                     K11=log(max(A31UA$SNF)),K12=log(max(A31UB$SNF)),
                                                                     K13=log(max(A31UC$SNF)),K14=log(max(G31UA$SNF)),
                                                                     K15=log(max(G31UB$SNF)),K16=log(max(G31UC$SNF)),
                                                                     K17=log(max(G31UD$SNF)),K18=log(max(G31UE$SNF)),
                                                                     K19=log(max(G31UF$SNF)),K20=log(max(M31UA$SNF)),
                                                                     K21=log(max(M31UB$SNF)),K22=log(max(M31UC$SNF)),
                                                                     K23=log(max(M31UD$SNF)),K24=log(max(M31UE$SNF)),
                                                                     K25=log(max(M31UF$SNF)),K26=log(max(M31UG$SNF)),
                                                                     K27=log(max(M31UH$SNF)),K28=log(max(R31UA$SNF)),
                                                                     K29=log(max(R31UB$SNF)),K30=log(max(R31UC$SNF)),
                                                                     K31=log(max(R31UD$SNF)),K32=log(max(R31UE$SNF)),
                                                                     K33=log(max(R31UF$SNF)),
                                                                     ruprhiz21=log(0.1),ruprhiz31=log(0.1),rupactin21=log(0.2),rupactin31=log(0.2),
                                                                     tLupG21A=log(70),tLupG21B=log(110),tLupG21C=log(100),
                                                                     tLupM21A=log(150),tLupM21B=log(100),tLupM21C=log(140),tLupM21D=log(180),
                                                                     tLupR21A=log(60),tLupR21B=log(30),tLupR21C=log(40),
                                                                     tLupA31A=log(160),tLupA31B=log(110),tLupA31C=log(110),
                                                                     tLupG31A=log(100),tLupG31B=log(70),tLupG31C=log(110),tLupG31D=log(40),tLupG31E=log(70),tLupG31F=log(110),
                                                                     tLupM31A=log(110),tLupM31B=log(120),tLupM31C=log(150),tLupM31D=log(100),tLupM31E=log(130),tLupM31F=log(130),tLupM31G=log(130),tLupM31H=log(130),
                                                                     tLupR31A=log(70),tLupR31B=log(100),tLupR31C=log(120),tLupR31D=log(30),tLupR31E=log(120),tLupR31F=log(70)),
                               data=list(t1=G21UA$Day,t2=G21UB$Day,
                                         t3=G21UC$Day,t4=M21UA$Day,
                                         t5=M21UB$Day,t6=M21UC$Day,
                                         t7=M21UD$Day,t8=R21UA$Day,
                                         t9=R21UB$Day,t10=R21UC$Day,
                                         t11=A31UA$Day,t12=A31UB$Day,
                                         t13=A31UC$Day,t14=G31UA$Day,
                                         t15=G31UB$Day,t16=G31UC$Day,
                                         t17=G31UD$Day,t18=G31UE$Day,
                                         t19=G31UF$Day,t20=M31UA$Day,
                                         t21=M31UB$Day,t22=M31UC$Day,
                                         t23=M31UD$Day,t24=M31UE$Day,
                                         t25=M31UF$Day,t26=M31UG$Day,
                                         t27=M31UH$Day,t28=R31UA$Day,
                                         t29=R31UB$Day,t30=R31UC$Day,
                                         t31=R31UD$Day,t32=R31UE$Day,
                                         t33=R31UF$Day,
                                         Nasedat1=G21UA$SNF,Nasedat2=G21UB$SNF,
                                         Nasedat3=G21UC$SNF,Nasedat4=M21UA$SNF,
                                         Nasedat5=M21UB$SNF,Nasedat6=M21UC$SNF,
                                         Nasedat7=M21UD$SNF,Nasedat8=R21UA$SNF,
                                         Nasedat9=R21UB$SNF,Nasedat10=R21UC$SNF,
                                         Nasedat11=A31UA$SNF,Nasedat12=A31UB$SNF,
                                         Nasedat13=A31UC$SNF,Nasedat14=G31UA$SNF,
                                         Nasedat15=G31UB$SNF,Nasedat16=G31UC$SNF,
                                         Nasedat17=G31UD$SNF,Nasedat18=G31UE$SNF,
                                         Nasedat19=G31UF$SNF,Nasedat20=M31UA$SNF,
                                         Nasedat21=M31UB$SNF,Nasedat22=M31UC$SNF,
                                         Nasedat23=M31UD$SNF,Nasedat24=M31UE$SNF,
                                         Nasedat25=M31UF$SNF,Nasedat26=M31UG$SNF,
                                         Nasedat27=M31UH$SNF,Nasedat28=R31UA$SNF,
                                         Nasedat29=R31UB$SNF,Nasedat30=R31UC$SNF,
                                         Nasedat31=R31UD$SNF,Nasedat32=R31UE$SNF,
                                         Nasedat33=R31UF$SNF),
                               control=list(maxit=20000))
summary(fit_sigmoid_up_symtemp)

# Calculate delta AICc to find best model for r
AICctab(fit_sigmoid_up_same,fit_sigmoid_up_sptemp,fit_sigmoid_up_sym,fit_sigmoid_up_biome,
        fit_sigmoid_up_sp,fit_sigmoid_up_temp,fit_sigmoid_up_biometemp,fit_sigmoid_up_symtemp,nobs=515)
#best is symbiosis x temperature interaction
#symbiosis is 0.6 deltaAICc worse

# Likelihood ratio test to see if temperature is a necessary predictor
anova(fit_sigmoid_up_sym,fit_sigmoid_up_symtemp) #p=0.051 # marginally significant difference
#actinorhizal has higher r than rhizobial
#temperature increases rate

### Compare tL values for best fit (symbiosis type)
ft_sym_up<-coef(fit_sigmoid_up_sym)
tL.up<-exp(c(ft_sym_up[4],ft_sym_up[5],ft_sym_up[6],
             ft_sym_up[7],ft_sym_up[8],ft_sym_up[9],ft_sym_up[10],
             ft_sym_up[11],ft_sym_up[12],ft_sym_up[13],
             ft_sym_up[14],ft_sym_up[15],ft_sym_up[16],
             ft_sym_up[17],ft_sym_up[18],ft_sym_up[19],ft_sym_up[20],ft_sym_up[21],ft_sym_up[22],
             ft_sym_up[23],ft_sym_up[24],ft_sym_up[25],ft_sym_up[26],ft_sym_up[27],ft_sym_up[28],ft_sym_up[29],ft_sym_up[30],
             ft_sym_up[31],ft_sym_up[32],ft_sym_up[33],ft_sym_up[34],ft_sym_up[35],ft_sym_up[36]))
treat.up<-c("G21","G21","G21","M21","M21","M21","M21","R21","R21","R21",
         "A31","A31","A31","G31","G31","G31","G31","G31","G31",
         "M31","M31","M31","M31","M31","M31","M31","M31",
         "R31","R31","R31","R31","R31","R31")
sym.up<-c("Rhiz","Rhiz","Rhiz","Actin","Actin","Actin","Actin","Rhiz","Rhiz","Rhiz",
       "Actin","Actin","Actin","Rhiz","Rhiz","Rhiz","Rhiz","Rhiz","Rhiz",
       "Actin","Actin","Actin","Actin","Actin","Actin","Actin","Actin",
       "Rhiz","Rhiz","Rhiz","Rhiz","Rhiz","Rhiz")
sp.up<-c("G","G","G","M","M","M","M","R","R","R",
      "A","A","A","G","G","G","G","G","G",
      "M","M","M","M","M","M","M","M",
      "R","R","R","R","R","R")
biome.up<-c("trop","trop","trop","trop","trop","trop","trop","temp","temp","temp",
         "temp","temp","temp","trop","trop","trop","trop","trop","trop",
         "trop","trop","trop","trop","trop","trop","trop","trop",
         "temp","temp","temp","temp","temp","temp")
temp.up<-c("21","21","21","21","21","21","21","21","21","21",
        "31","31","31","31","31","31","31","31","31",
        "31","31","31","31","31","31","31","31",
        "31","31","31","31","31","31")

tL.aov<-aov(tL.up~treat.up)
summary(tL.aov) #p=0.001
TukeyHSD(tL.aov)
#sig:
#M21-G31
#R21-A31
#R21-M21
#R21-M31

tL.up.null.lm<-lm(tL.up~1)
tL.up.tr.lm<-lm(tL.up~treat.up)
tL.up.sym.lm<-lm(tL.up~sym.up)
tL.up.temp.lm<-lm(tL.up~temp.up)
tL.up.biome.lm<-lm(tL.up~biome.up)
tL.up.sp.lm<-lm(tL.up~sp.up)
tL.up.sp.temp.lm<-lm(tL.up~sp.up*temp.up)
tL.up.biome.temp.lm<-lm(tL.up~biome.up*temp.up)
tL.up.sym.temp.lm<-lm(tL.up~sym.up*temp.up)

AICctab(tL.up.null.lm,tL.up.tr.lm,tL.up.sym.lm,
        tL.up.temp.lm,tL.up.biome.lm,tL.up.sp.lm,tL.up.sp.temp.lm,
        tL.up.biome.temp.lm,tL.up.sym.temp.lm, nobs = 33)

anova(tL.up.null.lm,tL.up.sym.lm) #p=5.843e-05

### Coefficients and 95%CI for supplementary tables

## Down-regulation
# For r
rdowntemp21.SE<-summary(fit_sigmoid_down_biometemp)@coef["rdowntemp21", "Std. Error"]
rdowntemp21.mu<-coef(fit_sigmoid_down_biometemp)[2]
quantile(rnorm(10000,rdowntemp21.mu,rdowntemp21.SE),c(0.025,0.975))

rdowntrop21.SE<-summary(fit_sigmoid_down_biometemp)@coef["rdowntrop21", "Std. Error"]
rdowntrop21.mu<-coef(fit_sigmoid_down_biometemp)[3]
quantile(rnorm(10000,rdowntrop21.mu,rdowntrop21.SE),c(0.025,0.975))

rdowntemp31.SE<-summary(fit_sigmoid_down_biometemp)@coef["rdowntemp31", "Std. Error"]
rdowntemp31.mu<-coef(fit_sigmoid_down_biometemp)[4]
quantile(rnorm(10000,rdowntemp31.mu,rdowntemp31.SE),c(0.025,0.975))

rdowntrop31.SE<-summary(fit_sigmoid_down_biometemp)@coef["rdowntrop31", "Std. Error"]
rdowntrop31.mu<-coef(fit_sigmoid_down_biometemp)[5]
quantile(rnorm(10000,rdowntrop31.mu,rdowntrop31.SE),c(0.025,0.975))

# For tL
summary(tL.down.null.lm)$coefficients[1,1]
quantile(rnorm(10000,summary(tL.down.null.lm)$coefficients[1,1],summary(tL.down.null.lm)$coefficients[1,2]),c(0.025,0.975))

## Up-regulation
# For r
rupactin.SE<-summary(fit_sigmoid_up_sym)@coef["rupActin", "Std. Error"]
rupactin.mu<-coef(fit_sigmoid_up_sym)[3]
exp(rupactin.mu)
exp(quantile(rnorm(10000,rupactin.mu,rupactin.SE),c(0.025,0.975)))

ruprhiz.SE<-summary(fit_sigmoid_up_sym)@coef["rupRhiz", "Std. Error"]
ruprhiz.mu<-coef(fit_sigmoid_up_sym)[2]
exp(ruprhiz.mu)
exp(quantile(rnorm(10000,ruprhiz.mu,ruprhiz.SE),c(0.025,0.975)))

## For tL
# Get the variance-covariance matrix
summary(tL.up.sym.lm)
vcov_matrix <- vcov(tL.up.sym.lm)

# Extract variances and covariance
var_intercept <- vcov_matrix["(Intercept)", "(Intercept)"]
var_symRhiz <- vcov_matrix["sym.upRhiz", "sym.upRhiz"]
cov_intercept_symRhiz <- vcov_matrix["(Intercept)", "sym.upRhiz"]

# Calculate the SEM for symRhiz
sem_symRhiz <- sqrt(var_intercept + var_symRhiz + 2 * cov_intercept_symRhiz)

tL.up.actin.est<-summary(tL.up.sym.lm)$coefficients[1,1]
tL.up.actin.se<-summary(tL.up.sym.lm)$coefficients[1,2]
tL.up.rhiz.est<-summary(tL.up.sym.lm)$coefficients[2,1]+tL.up.actin.est
tL.up.rhiz.se<-sem_symRhiz

quantile(rnorm(10000,tL.up.actin.est,tL.up.actin.se),c(0.025,0.975))
quantile(rnorm(10000,tL.up.rhiz.est,tL.up.rhiz.se),c(0.025,0.975))

################################################################################
### Figures
################################################################################

### Figure 1
par(pty="s")
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,5,4,4))
par(mar=c(0,0,0,0))

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=c(0,0.5,1,1.5,2),las=1,cex.axis=1.2)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G21DA$Day,rev(G21DA$Day)),y=c((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[10])))*G21DA$SNF.low/exp(ft_biometemp_down[42]),rev((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[10])))*G21DA$SNF.high/exp(ft_biometemp_down[42]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DB$Day,rev(G21DB$Day)),y=c((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[11])))*G21DB$SNF.low/exp(ft_biometemp_down[43]),rev((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[11])))*G21DB$SNF.high/exp(ft_biometemp_down[43]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DC$Day,rev(G21DC$Day)),y=c((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[12])))*G21DC$SNF.low/exp(ft_biometemp_down[44]),rev((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[12])))*G21DC$SNF.high/exp(ft_biometemp_down[44]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DD$Day,rev(G21DD$Day)),y=c((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[13])))*G21DD$SNF.low/exp(ft_biometemp_down[45]),rev((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[13])))*G21DD$SNF.high/exp(ft_biometemp_down[45]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[10])))*SNF/exp(ft_biometemp_down[42])~Day,data=G21DA,lwd=1,type="p",col="black",pch=0)
points((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[11])))*SNF/exp(ft_biometemp_down[43])~Day,data=G21DB,lwd=1,type="p",col="black",pch=1)
points((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[12])))*SNF/exp(ft_biometemp_down[44])~Day,data=G21DC,lwd=1,type="p",col="black",pch=2)
points((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[13])))*SNF/exp(ft_biometemp_down[45])~Day,data=G21DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[10]))),ft_biometemp_down[3],exp(ft_biometemp_down[10]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[11]))),ft_biometemp_down[3],exp(ft_biometemp_down[11]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[12]))),ft_biometemp_down[3],exp(ft_biometemp_down[12]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[13]))),ft_biometemp_down[3],exp(ft_biometemp_down[13]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)
curve(sigmoid((1+exp(ft_biometemp_down[3]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[3],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=2,col="blue",lwd=4,add=T) #tropical 21

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R21DA$Day,rev(R21DA$Day)),y=c((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[18])))*R21DA$SNF.low/exp(ft_biometemp_down[50]),rev((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[18])))*R21DA$SNF.high/exp(ft_biometemp_down[50]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DB$Day,rev(R21DB$Day)),y=c((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[19])))*R21DB$SNF.low/exp(ft_biometemp_down[51]),rev((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[19])))*R21DB$SNF.high/exp(ft_biometemp_down[51]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DC$Day,rev(R21DC$Day)),y=c((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[20])))*R21DC$SNF.low/exp(ft_biometemp_down[52]),rev((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[20])))*R21DC$SNF.high/exp(ft_biometemp_down[52]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DD$Day,rev(R21DD$Day)),y=c((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[21])))*R21DD$SNF.low/exp(ft_biometemp_down[53]),rev((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[21])))*R21DD$SNF.high/exp(ft_biometemp_down[53]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[18])))*SNF/exp(ft_biometemp_down[50])~Day,data=R21DA,lwd=1,type="p",col="black",pch=0)
points((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[19])))*SNF/exp(ft_biometemp_down[51])~Day,data=R21DB,lwd=1,type="p",col="black",pch=1)
points((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[20])))*SNF/exp(ft_biometemp_down[52])~Day,data=R21DC,lwd=1,type="p",col="black",pch=2)
points((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[21])))*SNF/exp(ft_biometemp_down[53])~Day,data=R21DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[18]))),ft_biometemp_down[2],exp(ft_biometemp_down[18]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[19]))),ft_biometemp_down[2],exp(ft_biometemp_down[19]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[20]))),ft_biometemp_down[2],exp(ft_biometemp_down[20]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[21]))),ft_biometemp_down[2],exp(ft_biometemp_down[21]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)
curve(sigmoid((1+exp(ft_biometemp_down[2]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[2],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=1,col="blue",lwd=4,add=T) #temperate 21

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M21DA$Day,rev(M21DA$Day)),y=c((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[14])))*M21DA$SNF.low/exp(ft_biometemp_down[46]),rev((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[14])))*M21DA$SNF.high/exp(ft_biometemp_down[46]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DB$Day,rev(M21DB$Day)),y=c((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[15])))*M21DB$SNF.low/exp(ft_biometemp_down[47]),rev((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[15])))*M21DB$SNF.high/exp(ft_biometemp_down[47]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DC$Day,rev(M21DC$Day)),y=c((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[16])))*M21DC$SNF.low/exp(ft_biometemp_down[48]),rev((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[16])))*M21DC$SNF.high/exp(ft_biometemp_down[48]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DD$Day,rev(M21DD$Day)),y=c((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[17])))*M21DD$SNF.low/exp(ft_biometemp_down[49]),rev((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[17])))*M21DD$SNF.high/exp(ft_biometemp_down[49]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[14])))*SNF/exp(ft_biometemp_down[46])~Day,data=M21DA,lwd=1,type="p",col="black",lty=2,pch=0)
points((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[15])))*SNF/exp(ft_biometemp_down[47])~Day,data=M21DB,lwd=1,type="p",col="black",lty=2,pch=1)
points((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[16])))*SNF/exp(ft_biometemp_down[48])~Day,data=M21DC,lwd=1,type="p",col="black",lty=2,pch=2)
points((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[17])))*SNF/exp(ft_biometemp_down[49])~Day,data=M21DD,lwd=1,type="p",col="black",lty=2,pch=3)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[14]))),ft_biometemp_down[3],exp(ft_biometemp_down[14]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[15]))),ft_biometemp_down[3],exp(ft_biometemp_down[15]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[16]))),ft_biometemp_down[3],exp(ft_biometemp_down[16]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[17]))),ft_biometemp_down[3],exp(ft_biometemp_down[17]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)
curve(sigmoid((1+exp(ft_biometemp_down[3]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[3],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=2,col="blue",lwd=4,add=T) #tropical 21

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A21DA$Day,rev(A21DA$Day)),y=c((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[6])))*A21DA$SNF.low/exp(ft_biometemp_down[38]),rev((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[6])))*A21DA$SNF.high/exp(ft_biometemp_down[38]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DB$Day,rev(A21DB$Day)),y=c((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[7])))*A21DB$SNF.low/exp(ft_biometemp_down[39]),rev((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[7])))*A21DB$SNF.high/exp(ft_biometemp_down[39]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DC$Day,rev(A21DC$Day)),y=c((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[8])))*A21DC$SNF.low/exp(ft_biometemp_down[40]),rev((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[8])))*A21DC$SNF.high/exp(ft_biometemp_down[40]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DD$Day,rev(A21DD$Day)),y=c((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[9])))*A21DD$SNF.low/exp(ft_biometemp_down[41]),rev((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[9])))*A21DD$SNF.high/exp(ft_biometemp_down[41]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[6])))*SNF/exp(ft_biometemp_down[38])~Day,data=A21DA,lwd=1,type="p",col="black",pch=0)
points((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[7])))*SNF/exp(ft_biometemp_down[39])~Day,data=A21DB,lwd=1,type="p",col="black",pch=1)
points((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[8])))*SNF/exp(ft_biometemp_down[40])~Day,data=A21DC,lwd=1,type="p",col="black",pch=2)
points((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[9])))*SNF/exp(ft_biometemp_down[41])~Day,data=A21DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[6]))),ft_biometemp_down[2],exp(ft_biometemp_down[6]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[7]))),ft_biometemp_down[2],exp(ft_biometemp_down[7]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[8]))),ft_biometemp_down[2],exp(ft_biometemp_down[8]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[9]))),ft_biometemp_down[2],exp(ft_biometemp_down[9]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
curve(sigmoid((1+exp(ft_biometemp_down[2]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[2],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=1,col="blue",lwd=4,add=T) #temperate 21
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=c(0,.5,1,1.5,2),las=1,cex.axis=1.2)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G31DA$Day,rev(G31DA$Day)),y=c((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[26])))*G31DA$SNF.low/exp(ft_biometemp_down[58]),rev((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[26])))*G31DA$SNF.high/exp(ft_biometemp_down[58]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DB$Day,rev(G31DB$Day)),y=c((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[27])))*G31DB$SNF.low/exp(ft_biometemp_down[59]),rev((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[27])))*G31DB$SNF.high/exp(ft_biometemp_down[59]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DC$Day,rev(G31DC$Day)),y=c((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[28])))*G31DC$SNF.low/exp(ft_biometemp_down[60]),rev((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[28])))*G31DC$SNF.high/exp(ft_biometemp_down[60]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DD$Day,rev(G31DD$Day)),y=c((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[29])))*G31DD$SNF.low/exp(ft_biometemp_down[61]),rev((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[29])))*G31DD$SNF.high/exp(ft_biometemp_down[61]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[26])))*SNF/exp(ft_biometemp_down[58])~Day,data=G31DA,lwd=1,type="p",col="black",pch=0)
points((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[27])))*SNF/exp(ft_biometemp_down[59])~Day,data=G31DB,lwd=1,type="p",col="black",pch=1)
points((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[28])))*SNF/exp(ft_biometemp_down[60])~Day,data=G31DC,lwd=1,type="p",col="black",pch=2)
points((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[29])))*SNF/exp(ft_biometemp_down[61])~Day,data=G31DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[26]))),ft_biometemp_down[5],exp(ft_biometemp_down[26]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[27]))),ft_biometemp_down[5],exp(ft_biometemp_down[27]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[28]))),ft_biometemp_down[5],exp(ft_biometemp_down[28]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[29]))),ft_biometemp_down[5],exp(ft_biometemp_down[29]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[5],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=2,col="red",lwd=4,add=T) #tropical 31

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R31DA$Day,rev(R31DA$Day)),y=c((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[34])))*R31DA$SNF.low/exp(ft_biometemp_down[66]),rev((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[34])))*R31DA$SNF.high/exp(ft_biometemp_down[66]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DB$Day,rev(R31DB$Day)),y=c((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[35])))*R31DB$SNF.low/exp(ft_biometemp_down[67]),rev((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[35])))*R31DB$SNF.high/exp(ft_biometemp_down[67]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DC$Day,rev(R31DC$Day)),y=c((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[36])))*R31DC$SNF.low/exp(ft_biometemp_down[68]),rev((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[36])))*R31DC$SNF.high/exp(ft_biometemp_down[68]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DD$Day,rev(R31DD$Day)),y=c((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[37])))*R31DD$SNF.low/exp(ft_biometemp_down[69]),rev((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[37])))*R31DD$SNF.high/exp(ft_biometemp_down[69]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[34])))*SNF/exp(ft_biometemp_down[66])~Day,data=R31DA,lwd=1,type="p",col="black",pch=0)
points((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[35])))*SNF/exp(ft_biometemp_down[67])~Day,data=R31DB,lwd=1,type="p",col="black",pch=1)
points((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[36])))*SNF/exp(ft_biometemp_down[68])~Day,data=R31DC,lwd=1,type="p",col="black",pch=2)
points((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[37])))*SNF/exp(ft_biometemp_down[69])~Day,data=R31DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[34]))),ft_biometemp_down[4],exp(ft_biometemp_down[34]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[35]))),ft_biometemp_down[4],exp(ft_biometemp_down[35]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[36]))),ft_biometemp_down[4],exp(ft_biometemp_down[36]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[37]))),ft_biometemp_down[4],exp(ft_biometemp_down[37]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[4],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=1,col="red",lwd=4,add=T) #temperate 31

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M31DA$Day,rev(M31DA$Day)),y=c((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[30])))*M31DA$SNF.low/exp(ft_biometemp_down[62]),rev((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[30])))*M31DA$SNF.high/exp(ft_biometemp_down[62]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DB$Day,rev(M31DB$Day)),y=c((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[31])))*M31DB$SNF.low/exp(ft_biometemp_down[63]),rev((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[31])))*M31DB$SNF.high/exp(ft_biometemp_down[63]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DC$Day,rev(M31DC$Day)),y=c((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[32])))*M31DC$SNF.low/exp(ft_biometemp_down[64]),rev((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[32])))*M31DC$SNF.high/exp(ft_biometemp_down[64]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DD$Day,rev(M31DD$Day)),y=c((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[33])))*M31DD$SNF.low/exp(ft_biometemp_down[65]),rev((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[33])))*M31DD$SNF.high/exp(ft_biometemp_down[65]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[30])))*SNF/exp(ft_biometemp_down[62])~Day,data=M31DA,lwd=1,type="p",col="black",pch=0)
points((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[31])))*SNF/exp(ft_biometemp_down[63])~Day,data=M31DB,lwd=1,type="p",col="black",pch=1)
points((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[32])))*SNF/exp(ft_biometemp_down[64])~Day,data=M31DC,lwd=1,type="p",col="black",pch=2)
points((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[33])))*SNF/exp(ft_biometemp_down[65])~Day,data=M31DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[30]))),ft_biometemp_down[5],exp(ft_biometemp_down[30]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[31]))),ft_biometemp_down[5],exp(ft_biometemp_down[31]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[32]))),ft_biometemp_down[5],exp(ft_biometemp_down[32]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[33]))),ft_biometemp_down[5],exp(ft_biometemp_down[33]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[5],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=2,col="red",lwd=4,add=T) #tropical 31

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  h'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A31DA$Day,rev(A31DA$Day)),y=c((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[22])))*A31DA$SNF.low/exp(ft_biometemp_down[54]),rev((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[22])))*A31DA$SNF.high/exp(ft_biometemp_down[54]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DB$Day,rev(A31DB$Day)),y=c((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[23])))*A31DB$SNF.low/exp(ft_biometemp_down[55]),rev((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[23])))*A31DB$SNF.high/exp(ft_biometemp_down[55]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DC$Day,rev(A31DC$Day)),y=c((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[24])))*A31DC$SNF.low/exp(ft_biometemp_down[56]),rev((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[24])))*A31DC$SNF.high/exp(ft_biometemp_down[56]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DD$Day,rev(A31DD$Day)),y=c((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[25])))*A31DD$SNF.low/exp(ft_biometemp_down[57]),rev((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[25])))*A31DD$SNF.high/exp(ft_biometemp_down[57]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[22])))*SNF/exp(ft_biometemp_down[54])~Day,data=A31DA,lwd=1,type="p",col="black",pch=0)
points((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[23])))*SNF/exp(ft_biometemp_down[55])~Day,data=A31DB,lwd=1,type="p",col="black",pch=1)
points((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[24])))*SNF/exp(ft_biometemp_down[56])~Day,data=A31DC,lwd=1,type="p",col="black",pch=2)
points((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[25])))*SNF/exp(ft_biometemp_down[57])~Day,data=A31DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[22]))),ft_biometemp_down[4],exp(ft_biometemp_down[22]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[23]))),ft_biometemp_down[4],exp(ft_biometemp_down[23]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[24]))),ft_biometemp_down[4],exp(ft_biometemp_down[24]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[25]))),ft_biometemp_down[4],exp(ft_biometemp_down[25]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[4],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=1,col="red",lwd=4,add=T) #temperate 31
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('N fixation (normalized to 1 at day = 0)'),side=2,line=2.9,cex=1,outer=T)
mtext(expression('Time since switching to high N (days)'),side=1,line=3,cex=1,outer=T)

### Figure 2
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=c(0,.5,1,1.5,2),cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G21UA$Day,rev(G21UA$Day)),y=c(G21UA$SNF.low/exp(ft_sym_up[37]),rev(G21UA$SNF.high/exp(ft_sym_up[37]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UB$Day,rev(G21UB$Day)),y=c(G21UB$SNF.low/exp(ft_sym_up[38]),rev(G21UB$SNF.high/exp(ft_sym_up[38]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UC$Day,rev(G21UC$Day)),y=c(G21UC$SNF.low/exp(ft_sym_up[39]),rev(G21UC$SNF.high/exp(ft_sym_up[39]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/exp(ft_sym_up[37])~Day,data=G21UA,lwd=1,type="p",col="black",pch=0)
points(SNF/exp(ft_sym_up[38])~Day,data=G21UB,lwd=1,type="p",col="black",pch=1)
points(SNF/exp(ft_sym_up[39])~Day,data=G21UC,lwd=1,type="p",col="black",pch=2)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[4]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[5]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[6]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),tL.up.rhiz.est,x),from=0,to=200,lty=1,col="orchid3",lwd=5,add=T)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R21UA$Day,rev(R21UA$Day)),y=c(R21UA$SNF.low/exp(ft_sym_up[44]),rev(R21UA$SNF.high/exp(ft_sym_up[44]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UB$Day,rev(R21UB$Day)),y=c(R21UB$SNF.low/exp(ft_sym_up[45]),rev(R21UB$SNF.high/exp(ft_sym_up[45]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UC$Day,rev(R21UC$Day)),y=c(R21UC$SNF.low/exp(ft_sym_up[46]),rev(R21UC$SNF.high/exp(ft_sym_up[46]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/exp(ft_sym_up[44])~Day,data=R21UA,lwd=1,type="p",col="black",pch=0)
points(SNF/exp(ft_sym_up[45])~Day,data=R21UB,lwd=1,type="p",col="black",pch=1)
points(SNF/exp(ft_sym_up[46])~Day,data=R21UC,lwd=1,type="p",col="black",pch=2)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[11]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[12]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[13]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),tL.up.rhiz.est,x),from=0,to=200,lty=1,col="orchid3",lwd=5,add=T)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M21UA$Day,rev(M21UA$Day)),y=c(M21UA$SNF.low/exp(ft_sym_up[40]),rev(M21UA$SNF.high/exp(ft_sym_up[40]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UB$Day,rev(M21UB$Day)),y=c(M21UB$SNF.low/exp(ft_sym_up[41]),rev(M21UB$SNF.high/exp(ft_sym_up[41]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UC$Day,rev(M21UC$Day)),y=c(M21UC$SNF.low/exp(ft_sym_up[42]),rev(M21UC$SNF.high/exp(ft_sym_up[42]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UD$Day,rev(M21UD$Day)),y=c(M21UD$SNF.low/exp(ft_sym_up[43]),rev(M21UD$SNF.high/exp(ft_sym_up[43]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/exp(ft_sym_up[40])~Day,data=M21UA,lwd=1,type="p",col="black",pch=0)
points(SNF/exp(ft_sym_up[41])~Day,data=M21UB,lwd=1,type="p",col="black",pch=1)
points(SNF/exp(ft_sym_up[42])~Day,data=M21UC,lwd=1,type="p",col="black",pch=2)
points(SNF/exp(ft_sym_up[43])~Day,data=M21UD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[7]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[8]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[9]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[10]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),tL.up.actin.est,x),from=0,to=200,lty=1,col="goldenrod1",lwd=5,add=T)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n",bty="n")

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=c(0,0.5,1,1.5,2),cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G31UA$Day,rev(G31UA$Day)),y=c(G31UA$SNF.low/exp(ft_sym_up[50]),rev(G31UA$SNF.high/exp(ft_sym_up[50]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UB$Day,rev(G31UB$Day)),y=c(G31UB$SNF.low/exp(ft_sym_up[51]),rev(G31UB$SNF.high/exp(ft_sym_up[51]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UC$Day,rev(G31UC$Day)),y=c(G31UC$SNF.low/exp(ft_sym_up[52]),rev(G31UC$SNF.high/exp(ft_sym_up[52]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UD$Day,rev(G31UD$Day)),y=c(G31UD$SNF.low/exp(ft_sym_up[53]),rev(G31UD$SNF.high/exp(ft_sym_up[53]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UE$Day,rev(G31UE$Day)),y=c(G31UE$SNF.low/exp(ft_sym_up[54]),rev(G31UE$SNF.high/exp(ft_sym_up[54]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UF$Day,rev(G31UF$Day)),y=c(G31UF$SNF.low/exp(ft_sym_up[55]),rev(G31UF$SNF.high/exp(ft_sym_up[55]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/exp(ft_sym_up[50])~Day,data=G31UA,lwd=1,type="p",col="black",pch=0)
points(SNF/exp(ft_sym_up[51])~Day,data=G31UB,lwd=1,type="p",col="black",pch=1)
points(SNF/exp(ft_sym_up[52])~Day,data=G31UC,lwd=1,type="p",col="black",pch=2)
points(SNF/exp(ft_sym_up[53])~Day,data=G31UD,lwd=1,type="p",col="black",pch=3)
points(SNF/exp(ft_sym_up[54])~Day,data=G31UE,lwd=1,type="p",col="black",pch=4)
points(SNF/exp(ft_sym_up[55])~Day,data=G31UF,lwd=1,type="p",col="black",pch=5)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[17]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[18]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[19]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[20]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[21]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[22]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),tL.up.rhiz.est,x),from=0,to=200,lty=1,col="orchid3",lwd=5,add=T)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R31UA$Day,rev(R31UA$Day)),y=c(R31UA$SNF.low/exp(ft_sym_up[64]),rev(R31UA$SNF.high/exp(ft_sym_up[64]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UB$Day,rev(R31UB$Day)),y=c(R31UB$SNF.low/exp(ft_sym_up[65]),rev(R31UB$SNF.high/exp(ft_sym_up[65]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UC$Day,rev(R31UC$Day)),y=c(R31UC$SNF.low/exp(ft_sym_up[66]),rev(R31UC$SNF.high/exp(ft_sym_up[66]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UD$Day,rev(R31UD$Day)),y=c(R31UD$SNF.low/exp(ft_sym_up[67]),rev(R31UD$SNF.high/exp(ft_sym_up[67]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UE$Day,rev(R31UE$Day)),y=c(R31UE$SNF.low/exp(ft_sym_up[68]),rev(R31UE$SNF.high/exp(ft_sym_up[68]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UF$Day,rev(R31UF$Day)),y=c(R31UF$SNF.low/exp(ft_sym_up[69]),rev(R31UF$SNF.high/exp(ft_sym_up[69]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/exp(ft_sym_up[64])~Day,data=R31UA,lwd=1,type="p",col="black",pch=0)
points(SNF/exp(ft_sym_up[65])~Day,data=R31UB,lwd=1,type="p",col="black",pch=1)
points(SNF/exp(ft_sym_up[66])~Day,data=R31UC,lwd=1,type="p",col="black",pch=2)
points(SNF/exp(ft_sym_up[67])~Day,data=R31UD,lwd=1,type="p",col="black",pch=3)
points(SNF/exp(ft_sym_up[68])~Day,data=R31UE,lwd=1,type="p",col="black",pch=4)
points(SNF/exp(ft_sym_up[69])~Day,data=R31UF,lwd=1,type="p",col="black",pch=5)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[31]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[32]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[33]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[34]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[35]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[36]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),tL.up.rhiz.est,x),from=0,to=200,lty=1,col="orchid3",lwd=5,add=T)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M31UA$Day,rev(M31UA$Day)),y=c(M31UA$SNF.low/exp(ft_sym_up[56]),rev(M31UA$SNF.high/exp(ft_sym_up[56]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UB$Day,rev(M31UB$Day)),y=c(M31UB$SNF.low/exp(ft_sym_up[57]),rev(M31UB$SNF.high/exp(ft_sym_up[57]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UC$Day,rev(M31UC$Day)),y=c(M31UC$SNF.low/exp(ft_sym_up[58]),rev(M31UC$SNF.high/exp(ft_sym_up[58]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UD$Day,rev(M31UD$Day)),y=c(M31UD$SNF.low/exp(ft_sym_up[59]),rev(M31UD$SNF.high/exp(ft_sym_up[59]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UE$Day,rev(M31UE$Day)),y=c(M31UE$SNF.low/exp(ft_sym_up[60]),rev(M31UE$SNF.high/exp(ft_sym_up[60]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UF$Day,rev(M31UF$Day)),y=c(M31UF$SNF.low/exp(ft_sym_up[61]),rev(M31UF$SNF.high/exp(ft_sym_up[61]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UG$Day,rev(M31UG$Day)),y=c(M31UG$SNF.low/exp(ft_sym_up[62]),rev(M31UG$SNF.high/exp(ft_sym_up[62]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UH$Day,rev(M31UH$Day)),y=c(M31UH$SNF.low/exp(ft_sym_up[63]),rev(M31UH$SNF.high/exp(ft_sym_up[63]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/exp(ft_sym_up[56])~Day,data=M31UA,lwd=1,type="p",col="black",pch=0)
points(SNF/exp(ft_sym_up[57])~Day,data=M31UB,lwd=1,type="p",col="black",pch=1)
points(SNF/exp(ft_sym_up[58])~Day,data=M31UC,lwd=1,type="p",col="black",pch=2)
points(SNF/exp(ft_sym_up[59])~Day,data=M31UD,lwd=1,type="p",col="black",pch=3)
points(SNF/exp(ft_sym_up[60])~Day,data=M31UE,lwd=1,type="p",col="black",pch=4)
points(SNF/exp(ft_sym_up[61])~Day,data=M31UF,lwd=1,type="p",col="black",pch=5)
points(SNF/exp(ft_sym_up[62])~Day,data=M31UG,lwd=1,type="p",col="black",pch=6)
points(SNF/exp(ft_sym_up[63])~Day,data=M31UH,lwd=1,type="p",col="black",pch=7)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[23]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[24]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[25]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[26]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[27]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[28]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[29]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[30]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),tL.up.actin.est,x),from=0,to=200,lty=1,col="goldenrod1",lwd=5,add=T)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,2),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,.5,1,1.5,2),labels=F,cex.axis=1.2)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A31UA$Day,rev(A31UA$Day)),y=c(A31UA$SNF.low/exp(ft_sym_up[47]),rev(A31UA$SNF.high/exp(ft_sym_up[47]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UB$Day,rev(A31UB$Day)),y=c(A31UB$SNF.low/exp(ft_sym_up[48]),rev(A31UB$SNF.high/exp(ft_sym_up[48]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UC$Day,rev(A31UC$Day)),y=c(A31UC$SNF.low/exp(ft_sym_up[49]),rev(A31UC$SNF.high/exp(ft_sym_up[49]))),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/exp(ft_sym_up[47])~Day,data=A31UA,lwd=1,type="p",col="black",pch=0)
points(SNF/exp(ft_sym_up[48])~Day,data=A31UB,lwd=1,type="p",col="black",pch=1)
points(SNF/exp(ft_sym_up[49])~Day,data=A31UC,lwd=1,type="p",col="black",pch=2)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[14]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[15]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[16]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),tL.up.actin.est,x),from=0,to=200,lty=1,col="goldenrod1",lwd=5,add=T)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('N fixation (normalized to 1 at maximum)'),side=2,line=2.9,cex=1,outer=T)
mtext(expression('Time since switching to low N (days)'),side=1,line=3,cex=1,outer=T)

### Figure 3
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))
par(oma=c(0,0,0,0))

plot(0:100,(0:100)/100,col="white",xlab="Time (days)",
     ylab="N fixation (normalized to 1 at day = 0)",las=1,ylim=c(0,1),main=NA)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[6]))),ft_biometemp_down[2],exp(ft_biometemp_down[6]),x),from=0,to=100,lty=1,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[7]))),ft_biometemp_down[2],exp(ft_biometemp_down[7]),x),from=0,to=100,lty=1,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[8]))),ft_biometemp_down[2],exp(ft_biometemp_down[8]),x),from=0,to=100,lty=1,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[9]))),ft_biometemp_down[2],exp(ft_biometemp_down[9]),x),from=0,to=100,lty=1,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[10]))),ft_biometemp_down[3],exp(ft_biometemp_down[10]),x),from=0,to=100,lty=2,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[11]))),ft_biometemp_down[3],exp(ft_biometemp_down[11]),x),from=0,to=100,lty=2,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[12]))),ft_biometemp_down[3],exp(ft_biometemp_down[12]),x),from=0,to=100,lty=2,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[13]))),ft_biometemp_down[3],exp(ft_biometemp_down[13]),x),from=0,to=100,lty=2,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[14]))),ft_biometemp_down[3],exp(ft_biometemp_down[14]),x),from=0,to=100,lty=2,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[15]))),ft_biometemp_down[3],exp(ft_biometemp_down[15]),x),from=0,to=100,lty=2,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[16]))),ft_biometemp_down[3],exp(ft_biometemp_down[16]),x),from=0,to=100,lty=2,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[17]))),ft_biometemp_down[3],exp(ft_biometemp_down[17]),x),from=0,to=100,lty=2,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[18]))),ft_biometemp_down[2],exp(ft_biometemp_down[18]),x),from=0,to=100,lty=1,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[19]))),ft_biometemp_down[2],exp(ft_biometemp_down[19]),x),from=0,to=100,lty=1,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[20]))),ft_biometemp_down[2],exp(ft_biometemp_down[20]),x),from=0,to=100,lty=1,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[21]))),ft_biometemp_down[2],exp(ft_biometemp_down[21]),x),from=0,to=100,lty=1,col="blue",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[22]))),ft_biometemp_down[4],exp(ft_biometemp_down[22]),x),from=0,to=100,lty=1,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[23]))),ft_biometemp_down[4],exp(ft_biometemp_down[23]),x),from=0,to=100,lty=1,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[24]))),ft_biometemp_down[4],exp(ft_biometemp_down[24]),x),from=0,to=100,lty=1,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[25]))),ft_biometemp_down[4],exp(ft_biometemp_down[25]),x),from=0,to=100,lty=1,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[26]))),ft_biometemp_down[5],exp(ft_biometemp_down[26]),x),from=0,to=100,lty=2,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[27]))),ft_biometemp_down[5],exp(ft_biometemp_down[27]),x),from=0,to=100,lty=2,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[28]))),ft_biometemp_down[5],exp(ft_biometemp_down[28]),x),from=0,to=100,lty=2,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[29]))),ft_biometemp_down[5],exp(ft_biometemp_down[29]),x),from=0,to=100,lty=2,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[30]))),ft_biometemp_down[5],exp(ft_biometemp_down[30]),x),from=0,to=100,lty=2,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[31]))),ft_biometemp_down[5],exp(ft_biometemp_down[31]),x),from=0,to=100,lty=2,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[32]))),ft_biometemp_down[5],exp(ft_biometemp_down[32]),x),from=0,to=100,lty=2,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[33]))),ft_biometemp_down[5],exp(ft_biometemp_down[33]),x),from=0,to=100,lty=2,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[34]))),ft_biometemp_down[4],exp(ft_biometemp_down[34]),x),from=0,to=100,lty=1,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[35]))),ft_biometemp_down[4],exp(ft_biometemp_down[35]),x),from=0,to=100,lty=1,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[36]))),ft_biometemp_down[4],exp(ft_biometemp_down[36]),x),from=0,to=100,lty=1,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[37]))),ft_biometemp_down[4],exp(ft_biometemp_down[37]),x),from=0,to=100,lty=1,col="red",lwd=0.5,add=T)
curve(sigmoid((1+exp(ft_biometemp_down[2]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[2],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=1,col="blue",lwd=4,add=T) #temperate 21
curve(sigmoid((1+exp(ft_biometemp_down[3]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[3],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=2,col="blue",lwd=4,add=T) #tropical 21
curve(sigmoid((1+exp(ft_biometemp_down[4]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[4],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=1,col="red",lwd=4,add=T) #temperate 31
curve(sigmoid((1+exp(ft_biometemp_down[5]*mean(exp(ft_biometemp_down[6:37])))),ft_biometemp_down[5],mean(exp(ft_biometemp_down[6:37])),x),from=0,to=100,lty=2,col="red",lwd=4,add=T) #tropical 31
legend("topright", legend = c(expression('Temperate; 21/15 '*degree*'C'), expression('Temperate; 31/25 '*degree*'C'),
                              expression('Tropical; 21/15 '*degree*'C'), expression('Tropical; 31/25 '*degree*'C')),
       lty=c(1,1,2,2), lwd=c(3,3,3,3), col = c("blue", "red", "blue", "red"), bty="n")

### Calculations for figure 4
biometemp.down.summary <- summary(fit_sigmoid_down_biometemp)
coef_down <- as.data.frame(biometemp.down.summary@coef)

## Parametric bootstrap of Eq. 5 parameters
down.tL.v<-rnorm(10000,summary(tL.down.null.lm)$coefficients[1],summary(tL.down.null.lm)$coefficients[2])
temp21.r.v<-rnorm(10000,coef(fit_sigmoid_down_biometemp)[2],coef_down["rdowntemp21", "Std. Error"])
trop21.r.v<-rnorm(10000,coef(fit_sigmoid_down_biometemp)[3],coef_down["rdowntrop21", "Std. Error"])
temp31.r.v<-rnorm(10000,coef(fit_sigmoid_down_biometemp)[4],coef_down["rdowntemp31", "Std. Error"])
trop31.r.v<-rnorm(10000,coef(fit_sigmoid_down_biometemp)[5],coef_down["rdowntrop31", "Std. Error"])

## Time to 5% down-regulation
temp21.t05<-down.tL.v-log((1+exp(temp21.r.v*down.tL.v))/0.95-1)/temp21.r.v
t05.temp21.est<-mean(exp(ft_biometemp_down[6:37]))-log((1+exp(ft_biometemp_down[2]*mean(exp(ft_biometemp_down[6:37]))))/0.95-1)/ft_biometemp_down[2]
temp21.d05.plot<-quantile(temp21.t05,c(0.025,0.975))

trop21.t05<-down.tL.v-log((1+exp(trop21.r.v*down.tL.v))/0.95-1)/trop21.r.v
t05.trop21.est<-mean(exp(ft_biometemp_down[6:37]))-log((1+exp(ft_biometemp_down[3]*mean(exp(ft_biometemp_down[6:37]))))/0.95-1)/ft_biometemp_down[3]
trop21.d05.plot<-quantile(trop21.t05,c(0.025,0.975))

temp31.t05<-down.tL.v-log((1+exp(temp31.r.v*down.tL.v))/0.95-1)/temp31.r.v
t05.temp31.est<-mean(exp(ft_biometemp_down[6:37]))-log((1+exp(ft_biometemp_down[4]*mean(exp(ft_biometemp_down[6:37]))))/0.95-1)/ft_biometemp_down[4]
temp31.d05.plot<-quantile(temp31.t05,c(0.025,0.975))

trop31.t05<-down.tL.v-log((1+exp(trop31.r.v*down.tL.v))/0.95-1)/trop31.r.v
t05.trop31.est<-mean(exp(ft_biometemp_down[6:37]))-log((1+exp(ft_biometemp_down[5]*mean(exp(ft_biometemp_down[6:37]))))/0.95-1)/ft_biometemp_down[5]
trop31.d05.plot<-quantile(trop31.t05,c(0.025,0.975))

d05.plot<-data.frame("Treatment"=c("temp21","temp21","temp31","temp31","trop21","trop21","trop31","trop31"),
                     "d05"=as.numeric(c(temp21.d05.plot,temp31.d05.plot,trop21.d05.plot,trop31.d05.plot)))

d05.all<-c(exp(ft_biometemp_down[6])-log((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[6])))/0.95-1)/ft_biometemp_down[2],
           exp(ft_biometemp_down[7])-log((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[7])))/0.95-1)/ft_biometemp_down[2],
           exp(ft_biometemp_down[8])-log((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[8])))/0.95-1)/ft_biometemp_down[2],
           exp(ft_biometemp_down[9])-log((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[9])))/0.95-1)/ft_biometemp_down[2],
           exp(ft_biometemp_down[10])-log((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[10])))/0.95-1)/ft_biometemp_down[3],
           exp(ft_biometemp_down[11])-log((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[11])))/0.95-1)/ft_biometemp_down[3],
           exp(ft_biometemp_down[12])-log((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[12])))/0.95-1)/ft_biometemp_down[3],
           exp(ft_biometemp_down[13])-log((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[13])))/0.95-1)/ft_biometemp_down[3],
           exp(ft_biometemp_down[14])-log((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[14])))/0.95-1)/ft_biometemp_down[3],
           exp(ft_biometemp_down[15])-log((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[15])))/0.95-1)/ft_biometemp_down[3],
           exp(ft_biometemp_down[16])-log((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[16])))/0.95-1)/ft_biometemp_down[3],
           exp(ft_biometemp_down[17])-log((1+exp(ft_biometemp_down[3]*exp(ft_biometemp_down[17])))/0.95-1)/ft_biometemp_down[3],
           exp(ft_biometemp_down[18])-log((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[18])))/0.95-1)/ft_biometemp_down[2],
           exp(ft_biometemp_down[19])-log((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[19])))/0.95-1)/ft_biometemp_down[2],
           exp(ft_biometemp_down[20])-log((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[20])))/0.95-1)/ft_biometemp_down[2],
           exp(ft_biometemp_down[21])-log((1+exp(ft_biometemp_down[2]*exp(ft_biometemp_down[21])))/0.95-1)/ft_biometemp_down[2],
           exp(ft_biometemp_down[22])-log((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[22])))/0.95-1)/ft_biometemp_down[4],
           exp(ft_biometemp_down[23])-log((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[23])))/0.95-1)/ft_biometemp_down[4],
           exp(ft_biometemp_down[24])-log((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[24])))/0.95-1)/ft_biometemp_down[4],
           exp(ft_biometemp_down[25])-log((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[25])))/0.95-1)/ft_biometemp_down[4],
           exp(ft_biometemp_down[26])-log((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[26])))/0.95-1)/ft_biometemp_down[5],
           exp(ft_biometemp_down[27])-log((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[27])))/0.95-1)/ft_biometemp_down[5],
           exp(ft_biometemp_down[28])-log((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[28])))/0.95-1)/ft_biometemp_down[5],
           exp(ft_biometemp_down[29])-log((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[29])))/0.95-1)/ft_biometemp_down[5],
           exp(ft_biometemp_down[30])-log((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[30])))/0.95-1)/ft_biometemp_down[5],
           exp(ft_biometemp_down[31])-log((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[31])))/0.95-1)/ft_biometemp_down[5],
           exp(ft_biometemp_down[32])-log((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[32])))/0.95-1)/ft_biometemp_down[5],
           exp(ft_biometemp_down[33])-log((1+exp(ft_biometemp_down[5]*exp(ft_biometemp_down[33])))/0.95-1)/ft_biometemp_down[5],
           exp(ft_biometemp_down[34])-log((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[34])))/0.95-1)/ft_biometemp_down[4],
           exp(ft_biometemp_down[35])-log((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[35])))/0.95-1)/ft_biometemp_down[4],
           exp(ft_biometemp_down[36])-log((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[36])))/0.95-1)/ft_biometemp_down[4],
           exp(ft_biometemp_down[37])-log((1+exp(ft_biometemp_down[4]*exp(ft_biometemp_down[37])))/0.95-1)/ft_biometemp_down[4])

d05.all.plot<-data.frame("Treatment"=treat.down,"Biome.Temp"=c(1,1,1,1,3,3,3,3,3,3,3,3,1,1,1,1,2,2,2,2,4,4,4,4,4,4,4,4,2,2,2,2),
                         "d05"=as.numeric(d05.all))

## Time to 95% down-regulation
temp21.t95<-down.tL.v-log(19+20*exp(temp21.r.v*down.tL.v))/temp21.r.v
t95.temp21.est<-mean(exp(ft_biometemp_down[6:37]))-log(19+20*exp(ft_biometemp_down[2]*mean(exp(ft_biometemp_down[6:37]))))/ft_biometemp_down[2]
temp21.d95.plot<-quantile(temp21.t95,c(0.025,0.975))

trop21.t95<-down.tL.v-log(19+20*exp(trop21.r.v*down.tL.v))/trop21.r.v
t95.trop21.est<-mean(exp(ft_biometemp_down[6:37]))-log(19+20*exp(ft_biometemp_down[3]*mean(exp(ft_biometemp_down[6:37]))))/ft_biometemp_down[3]
trop21.d95.plot<-quantile(trop21.t95,c(0.025,0.975))

temp31.t95<-down.tL.v-log(19+20*exp(temp31.r.v*down.tL.v))/temp31.r.v
t95.temp31.est<-mean(exp(ft_biometemp_down[6:37]))-log(19+20*exp(ft_biometemp_down[4]*mean(exp(ft_biometemp_down[6:37]))))/ft_biometemp_down[4]
temp31.d95.plot<-quantile(temp31.t95,c(0.025,0.975))

trop31.t95<-down.tL.v-log(19+20*exp(trop31.r.v*down.tL.v))/trop31.r.v
t95.trop31.est<-mean(exp(ft_biometemp_down[6:37]))-log(19+20*exp(ft_biometemp_down[5]*mean(exp(ft_biometemp_down[6:37]))))/ft_biometemp_down[5]
trop31.d95.plot<-quantile(trop31.t95,c(0.025,0.975))

d95.plot<-data.frame("Treatment"=c("temp21","temp21","temp31","temp31","trop21","trop21","trop31","trop31"),
                     "d95"=as.numeric(c(temp21.d95.plot,temp31.d95.plot,trop21.d95.plot,trop31.d95.plot)))

d95.all<-c(exp(ft_biometemp_down[6])-log(19+20*exp(ft_biometemp_down[2]*exp(ft_biometemp_down[6])))/ft_biometemp_down[2],
           exp(ft_biometemp_down[7])-log(19+20*exp(ft_biometemp_down[2]*exp(ft_biometemp_down[7])))/ft_biometemp_down[2],
           exp(ft_biometemp_down[8])-log(19+20*exp(ft_biometemp_down[2]*exp(ft_biometemp_down[8])))/ft_biometemp_down[2],
           exp(ft_biometemp_down[9])-log(19+20*exp(ft_biometemp_down[2]*exp(ft_biometemp_down[9])))/ft_biometemp_down[2],
           exp(ft_biometemp_down[10])-log(19+20*exp(ft_biometemp_down[3]*exp(ft_biometemp_down[10])))/ft_biometemp_down[3],
           exp(ft_biometemp_down[11])-log(19+20*exp(ft_biometemp_down[3]*exp(ft_biometemp_down[11])))/ft_biometemp_down[3],
           exp(ft_biometemp_down[12])-log(19+20*exp(ft_biometemp_down[3]*exp(ft_biometemp_down[12])))/ft_biometemp_down[3],
           exp(ft_biometemp_down[13])-log(19+20*exp(ft_biometemp_down[3]*exp(ft_biometemp_down[13])))/ft_biometemp_down[3],
           exp(ft_biometemp_down[14])-log(19+20*exp(ft_biometemp_down[3]*exp(ft_biometemp_down[14])))/ft_biometemp_down[3],
           exp(ft_biometemp_down[15])-log(19+20*exp(ft_biometemp_down[3]*exp(ft_biometemp_down[15])))/ft_biometemp_down[3],
           exp(ft_biometemp_down[16])-log(19+20*exp(ft_biometemp_down[3]*exp(ft_biometemp_down[16])))/ft_biometemp_down[3],
           exp(ft_biometemp_down[17])-log(19+20*exp(ft_biometemp_down[3]*exp(ft_biometemp_down[17])))/ft_biometemp_down[3],
           exp(ft_biometemp_down[18])-log(19+20*exp(ft_biometemp_down[2]*exp(ft_biometemp_down[18])))/ft_biometemp_down[2],
           exp(ft_biometemp_down[19])-log(19+20*exp(ft_biometemp_down[2]*exp(ft_biometemp_down[19])))/ft_biometemp_down[2],
           exp(ft_biometemp_down[20])-log(19+20*exp(ft_biometemp_down[2]*exp(ft_biometemp_down[20])))/ft_biometemp_down[2],
           exp(ft_biometemp_down[21])-log(19+20*exp(ft_biometemp_down[2]*exp(ft_biometemp_down[21])))/ft_biometemp_down[2],
           exp(ft_biometemp_down[22])-log(19+20*exp(ft_biometemp_down[4]*exp(ft_biometemp_down[22])))/ft_biometemp_down[4],
           exp(ft_biometemp_down[23])-log(19+20*exp(ft_biometemp_down[4]*exp(ft_biometemp_down[23])))/ft_biometemp_down[4],
           exp(ft_biometemp_down[24])-log(19+20*exp(ft_biometemp_down[4]*exp(ft_biometemp_down[24])))/ft_biometemp_down[4],
           exp(ft_biometemp_down[25])-log(19+20*exp(ft_biometemp_down[4]*exp(ft_biometemp_down[25])))/ft_biometemp_down[4],
           exp(ft_biometemp_down[26])-log(19+20*exp(ft_biometemp_down[5]*exp(ft_biometemp_down[26])))/ft_biometemp_down[5],
           exp(ft_biometemp_down[27])-log(19+20*exp(ft_biometemp_down[5]*exp(ft_biometemp_down[27])))/ft_biometemp_down[5],
           exp(ft_biometemp_down[28])-log(19+20*exp(ft_biometemp_down[5]*exp(ft_biometemp_down[28])))/ft_biometemp_down[5],
           exp(ft_biometemp_down[29])-log(19+20*exp(ft_biometemp_down[5]*exp(ft_biometemp_down[29])))/ft_biometemp_down[5],
           exp(ft_biometemp_down[30])-log(19+20*exp(ft_biometemp_down[5]*exp(ft_biometemp_down[30])))/ft_biometemp_down[5],
           exp(ft_biometemp_down[31])-log(19+20*exp(ft_biometemp_down[5]*exp(ft_biometemp_down[31])))/ft_biometemp_down[5],
           exp(ft_biometemp_down[32])-log(19+20*exp(ft_biometemp_down[5]*exp(ft_biometemp_down[32])))/ft_biometemp_down[5],
           exp(ft_biometemp_down[33])-log(19+20*exp(ft_biometemp_down[5]*exp(ft_biometemp_down[33])))/ft_biometemp_down[5],
           exp(ft_biometemp_down[34])-log(19+20*exp(ft_biometemp_down[4]*exp(ft_biometemp_down[34])))/ft_biometemp_down[4],
           exp(ft_biometemp_down[35])-log(19+20*exp(ft_biometemp_down[4]*exp(ft_biometemp_down[35])))/ft_biometemp_down[4],
           exp(ft_biometemp_down[36])-log(19+20*exp(ft_biometemp_down[4]*exp(ft_biometemp_down[36])))/ft_biometemp_down[4],
           exp(ft_biometemp_down[37])-log(19+20*exp(ft_biometemp_down[4]*exp(ft_biometemp_down[37])))/ft_biometemp_down[4])

d95.all.plot<-data.frame("Treatment"=treat.down,"Biome.Temp"=c(1,1,1,1,3,3,3,3,3,3,3,3,1,1,1,1,2,2,2,2,4,4,4,4,4,4,4,4,2,2,2,2),
                         "d95"=as.numeric(d95.all))

## Time from 5 to 95% down-regulation
t0595.temp21.v<-(log(0.05+exp(temp21.r.v*down.tL.v))-log(0.95)-log(19+20*exp(temp21.r.v*down.tL.v)))/temp21.r.v
t0595.trop21.v<-(log(0.05+exp(trop21.r.v*down.tL.v))-log(0.95)-log(19+20*exp(trop21.r.v*down.tL.v)))/trop21.r.v
t0595.temp31.v<-(log(0.05+exp(temp31.r.v*down.tL.v))-log(0.95)-log(19+20*exp(temp31.r.v*down.tL.v)))/temp31.r.v
t0595.trop31.v<-(log(0.05+exp(trop31.r.v*down.tL.v))-log(0.95)-log(19+20*exp(trop31.r.v*down.tL.v)))/trop31.r.v

temp21.d0595.plot<-quantile(t0595.temp21.v,c(0.025,0.975))
temp31.d0595.plot<-quantile(t0595.temp31.v,c(0.025,0.975))
trop21.d0595.plot<-quantile(t0595.trop21.v,c(0.025,0.975))
trop31.d0595.plot<-quantile(t0595.trop31.v,c(0.025,0.975))

d0595.plot<-data.frame("Treatment"=c("temp21","temp21","temp31","temp31","trop21","trop21","trop31","trop31"),
                       "d0595"=as.numeric(c(temp21.d0595.plot,temp31.d0595.plot,trop21.d0595.plot,trop31.d0595.plot)))

d0595.all<-d95.all-d05.all

d0595.all.plot<-data.frame("Treatment"=treat.down,"Biome.Temp"=c(1,1,1,1,3,3,3,3,3,3,3,3,1,1,1,1,2,2,2,2,4,4,4,4,4,4,4,4,2,2,2,2),
                           "d0595"=as.numeric(d0595.all))

## r down-regulation
temp21.r.plot<-quantile(temp21.r.v,c(0.025,0.975))
temp31.r.plot<-quantile(temp31.r.v,c(0.025,0.975))
trop21.r.plot<-quantile(trop21.r.v,c(0.025,0.975))
trop31.r.plot<-quantile(trop31.r.v,c(0.025,0.975))

r.plot<-data.frame("Treatment"=c("temp21","temp21","temp31","temp31","trop21","trop21","trop31","trop31"),
                   "r"=as.numeric(c(temp21.r.plot,temp31.r.plot,trop21.r.plot,trop31.r.plot)))

### Figure 4
nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),widths=c(2,2,2,2),heights=c(2,2,2,2),T)
layout.show(nf)
par(oma=c(2,2,2,2))
par(mar=c(7,4,1,0))
par(pty="s")

boxplot(d05~Treatment, data=d05.plot,yaxt="n",xaxt="n",ylim=c(0,20),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("Temp. Cold","Temp. Warm","Trop. Cold","Trop. Warm"),cex.axis=1,las=2)
axis(2,at=c(0,5,10,15,20),labels=c(0,5,10,15,20),las=1,cex.axis=1)
mtext(expression('Time to 5% down-regulation (days)'),side=2,line=3.5,cex=1)
arrows(1,temp21.d05.plot[1],1,temp21.d05.plot[2],angle=90,length=0.1,code=3,lwd=2,col="blue")
arrows(2,temp31.d05.plot[1],2,temp31.d05.plot[2],angle=90,length=0.1,code=3,lwd=2,col="red")
arrows(3,trop21.d05.plot[1],3,trop21.d05.plot[2],angle=90,length=0.1,code=3,lwd=2,col="blue")
arrows(4,trop31.d05.plot[1],4,trop31.d05.plot[2],angle=90,length=0.1,code=3,lwd=2,col="red")
points(1,t05.temp21.est,pch=16,cex=2,lwd=2,col="blue")
points(2,t05.temp31.est,pch=16,cex=2,lwd=2,col="red")
points(3,t05.trop21.est,pch=17,cex=2,lwd=2,col="blue")
points(4,t05.trop31.est,pch=17,cex=2,lwd=2,col="red")
stripchart(d05~Biome.Temp, vertical = TRUE, data = d05.all.plot, 
           method = "jitter", add = TRUE, pch=1,col = 'black',cex=1)
mtext(text="a",side=3,cex=1,adj=0)

boxplot(d95~Treatment, data=d95.plot,yaxt="n",xaxt="n",ylim=c(0,80),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("Temp. Cold","Temp. Warm","Trop. Cold","Trop. Warm"),cex.axis=1,las=2)
axis(2,at=c(0,20,40,60,80),labels=c(0,20,40,60,80),las=1,cex.axis=1)
mtext(expression('Time to 95% down-regulation (days)'),side=2,line=3.5,cex=1)
arrows(1,temp21.d95.plot[1],1,temp21.d95.plot[2],angle=90,length=0.1,code=3,lwd=2,col="blue")
arrows(2,temp31.d95.plot[1],2,temp31.d95.plot[2],angle=90,length=0.1,code=3,lwd=2,col="red")
arrows(3,trop21.d95.plot[1],3,trop21.d95.plot[2],angle=90,length=0.1,code=3,lwd=2,col="blue")
arrows(4,trop31.d95.plot[1],4,trop31.d95.plot[2],angle=90,length=0.1,code=3,lwd=2,col="red")
points(1,t95.temp21.est,pch=16,cex=2,lwd=2,col="blue")
points(2,t95.temp31.est,pch=16,cex=2,lwd=2,col="red")
points(3,t95.trop21.est,pch=17,cex=2,lwd=2,col="blue")
points(4,t95.trop31.est,pch=17,cex=2,lwd=2,col="red")
stripchart(d95~Biome.Temp, vertical = TRUE, data = d95.all.plot, 
           method = "jitter", add = TRUE, pch=1,col = 'black',cex=1)
mtext(text="b",side=3,cex=1,adj=0)

boxplot(d0595~Treatment, data=d0595.plot,yaxt="n",xaxt="n",ylim=c(0,80),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("Temp. Cold","Temp. Warm","Trop. Cold","Trop. Warm"),cex.axis=1,las=2)
axis(2,at=c(0,20,40,60,80),labels=c(0,20,40,60,80),las=1,cex.axis=1.1)
mtext(expression('Time from 5% to 95% down-regulation (days)'),side=2,line=3.5,cex=1)
arrows(1,temp21.d0595.plot[1],1,temp21.d0595.plot[2],angle=90,length=0.1,code=3,lwd=2,col="blue")
arrows(2,temp31.d0595.plot[1],2,temp31.d0595.plot[2],angle=90,length=0.1,code=3,lwd=2,col="red")
arrows(3,trop21.d0595.plot[1],3,trop21.d0595.plot[2],angle=90,length=0.1,code=3,lwd=2,col="blue")
arrows(4,trop31.d0595.plot[1],4,trop31.d0595.plot[2],angle=90,length=0.1,code=3,lwd=2,col="red")
points(1,c(t95.temp21.est - t05.temp21.est),pch=16,cex=2,lwd=2,col="blue")
points(2,c(t95.temp31.est - t05.temp31.est),pch=16,cex=2,lwd=2,col="red")
points(3,c(t95.trop21.est - t05.trop21.est),pch=17,cex=2,lwd=2,col="blue")
points(4,c(t95.trop31.est - t05.trop31.est),pch=17,cex=2,lwd=2,col="red")
stripchart(d0595~Biome.Temp, vertical = TRUE, data = d0595.all.plot, 
           method = "jitter", add = TRUE, pch=1,col = 'black',cex=1)
mtext(text="c",side=3,cex=1,adj=0)

boxplot(r~Treatment, data=r.plot,yaxt="n",xaxt="n",ylim=c(-0.2,0),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("Temp. Cold","Temp. Warm","Trop. Cold","Trop. Warm"),cex.axis=1,las=2)
axis(2,at=c(-0.2,-.15,-.1,-.05,0),labels=c(-0.20,-.15,-.10,-.05,0),las=1,cex.axis=1.1)
mtext(expression('Rate of change; '*italic(r)*' (day'^-1*')'),side=2,line=3.5,cex=1)
arrows(1,temp21.r.plot[1],1,temp21.r.plot[2],angle=90,length=0.1,code=3,lwd=2,col="blue")
arrows(2,temp31.r.plot[1],2,temp31.r.plot[2],angle=90,length=0.1,code=3,lwd=2,col="red")
arrows(3,trop21.r.plot[1],3,trop21.r.plot[2],angle=90,length=0.1,code=3,lwd=2,col="blue")
arrows(4,trop31.r.plot[1],4,trop31.r.plot[2],angle=90,length=0.1,code=3,lwd=2,col="red")
points(1,ft_biometemp_down[2],pch=16,cex=2,lwd=2,col="blue")
points(2,ft_biometemp_down[4],pch=16,cex=2,lwd=2,col="red")
points(3,ft_biometemp_down[3],pch=17,cex=2,lwd=2,col="blue")
points(4,ft_biometemp_down[5],pch=17,cex=2,lwd=2,col="red")
mtext(text="d",side=3,cex=1,adj=0)

### Figure 5
par(oma=c(1,2,1,0))
par(mar=c(3,4,1,1))
par(mfrow = c(1,2))

plot(0:200,(0:200)/200,col="white",xlab="Time (days)",
     ylab="N fixation (normalized to 1 at maximum)",las=1,ylim=c(0,1.2),main=NA)
mtext(text="a",side=3,cex=1.2,adj=0)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[4]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[5]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[6]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[7]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[8]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[9]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[10]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[11]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[12]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[13]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[14]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[15]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[16]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[17]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[18]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[19]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[20]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[21]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[22]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[23]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[24]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[25]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[26]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[27]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[28]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[29]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),exp(ft_sym_up[30]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[31]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[32]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[33]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[34]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[35]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),exp(ft_sym_up[36]),x),from=0,to=200,lty=1,col="orchid3",lwd=0.5,add=T)
curve(sigmoid(1,exp(ft_sym_up[3]),tL.up.actin.est,x),from=0,to=200,lty=1,col="goldenrod1",lwd=5,add=T)
curve(sigmoid(1,exp(ft_sym_up[2]),tL.up.rhiz.est,x),from=0,to=200,lty=1,col="orchid3",lwd=5,add=T)

plot(0:100,(0:100)/100,col="white",xlab="Time from 5% maximum N fixation (days)",
     ylab=NA,las=1,ylim=c(0,1.2),main=NA)
mtext(text="b",side=3,cex=1.2,adj=0)
curve(up.rate(1,exp(ft_sym_up[2]),x),from=0,to=100,lty=1,col="orchid3",lwd=5,add=T)
curve(up.rate(1,exp(ft_sym_up[3]),x),from=0,to=100,lty=1,col="goldenrod1",lwd=5,add=T)
legend("topleft", legend = c("Actinorhizal", "Rhizobial"), 
       lwd = c(3,3), lty = c(1, 1), bty = "n", col = c("goldenrod1","orchid3"))

### Calculations for Figure 6
sym.up.summary <- summary(fit_sigmoid_up_sym)
coef_up <- as.data.frame(sym.up.summary@coef)

## Parametric bootstrap of Eq. 5 parameters
actin.tL.v<-rnorm(10000,tL.up.actin.est,tL.up.actin.se)
actin.r.v<-exp(rnorm(10000,coef(fit_sigmoid_up_sym)[3],coef_up["rupActin", "Std. Error"]))
rhiz.tL.v<-rnorm(10000,tL.up.rhiz.est,tL.up.rhiz.se)
rhiz.r.v<-exp(rnorm(10000,coef(fit_sigmoid_up_sym)[2],coef_up["rupRhiz", "Std. Error"]))

## Time for 5% up-regulation
actin.t05<-actin.tL.v-log(19)/actin.r.v
t05.actin.est<-mean(exp(ft_sym_up[c(7:10,14:16,23:30)]))-log(19)/exp(ft_sym_up[3])
actin.u05.plot<-quantile(actin.t05,c(0.025,0.975))

rhiz.t05<-rhiz.tL.v-log(19)/rhiz.r.v
t05.rhiz.est<-mean(exp(ft_sym_up[c(4:6,11:13,17:22,31:36)]))-log(19)/exp(ft_sym_up[2])
rhiz.u05.plot<-quantile(rhiz.t05,c(0.025,0.975))

u05.plot<-data.frame("Treatment"=c("actin","actin","rhiz","rhiz"),
                     "u05"=as.numeric(c(actin.u05.plot,rhiz.u05.plot)))

u05.all<-c(exp(ft_sym_up[4])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[5])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[6])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[7])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[8])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[9])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[10])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[11])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[12])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[13])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[14])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[15])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[16])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[17])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[18])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[19])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[20])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[21])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[22])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[23])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[24])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[25])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[26])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[27])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[28])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[29])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[30])-log(19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[31])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[32])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[33])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[34])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[35])-log(19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[36])-log(19)/exp(ft_sym_up[2]))

u05.all.plot<-data.frame("Treatment"=treat.up,"sym"=c(2,2,2,1,1,1,1,2,2,2,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,2,2,2,2,2,2),
                         "u05"=as.numeric(u05.all))

## Time to 95% up-regulation
actin.t95<-actin.tL.v-log(1/19)/actin.r.v
t95.actin.est<-mean(exp(ft_sym_up[c(7:10,14:16,23:30)]))-log(1/19)/exp(ft_sym_up[3])
actin.u95.plot<-quantile(actin.t95,c(0.025,0.975))

rhiz.t95<-rhiz.tL.v-log(1/19)/rhiz.r.v
t95.rhiz.est<-mean(exp(ft_sym_up[c(4:6,11:13,17:22,31:36)]))-log(1/19)/exp(ft_sym_up[2])
rhiz.u95.plot<-quantile(rhiz.t95,c(0.025,0.975))

u95.plot<-data.frame("Treatment"=c("actin","actin","rhiz","rhiz"),
                     "u95"=as.numeric(c(actin.u95.plot,rhiz.u95.plot)))

u95.all<-c(exp(ft_sym_up[4])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[5])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[6])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[7])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[8])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[9])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[10])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[11])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[12])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[13])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[14])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[15])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[16])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[17])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[18])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[19])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[20])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[21])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[22])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[23])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[24])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[25])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[26])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[27])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[28])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[29])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[30])-log(1/19)/exp(ft_sym_up[3]),
           exp(ft_sym_up[31])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[32])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[33])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[34])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[35])-log(1/19)/exp(ft_sym_up[2]),
           exp(ft_sym_up[36])-log(1/19)/exp(ft_sym_up[2]))

u95.all.plot<-data.frame("Treatment"=treat.up,"sym"=c(2,2,2,1,1,1,1,2,2,2,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,2,2,2,2,2,2),
                         "u95"=as.numeric(u95.all))

## Time from 5 to 95% up-regulation
actin.u0595.plot<-quantile(c(2*log(19)/actin.r.v),c(0.025,0.975))
rhiz.u0595.plot<-quantile(c(2*log(19)/rhiz.r.v),c(0.025,0.975))

u0595.plot<-data.frame("Treatment"=c("actin","actin","rhiz","rhiz"),
                       "Symbiosis"=c("actin","actin","rhiz","rhiz"),
                       "u0595"=as.numeric(c(actin.u0595.plot,rhiz.u0595.plot)))

u0595.all<-u95.all-u05.all

u0595.all.plot<-data.frame("Treatment"=treat.up,"sym.temp"=c(2,2,2,1,1,1,1,2,2,2,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,2,2,2,2,2,2),
                           "u0595"=as.numeric(u0595.all))

## r up-regulation
actin.r.plot<-quantile(actin.r.v,c(0.025,0.975))
rhiz.r.plot<-quantile(rhiz.r.v,c(0.025,0.975))

r.plot<-data.frame("Treatment"=c("actin","actin","rhiz","rhiz"),
                   "r"=as.numeric(c(actin.r.plot,rhiz.r.plot)))

### FIGURE 6
nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),widths=c(2,2,2,2),heights=c(2,2,2,2),T)
layout.show(nf)
par(oma=c(2,2,2,2))
par(mar=c(7,4,1,0))
par(pty="s")

boxplot(u05~Treatment, data=u05.plot,yaxt="n",xaxt="n",ylim=c(0,200),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2),labels=c("Actinorhizal","Rhizobial"),cex.axis=1,las=2)
axis(2,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,cex.axis=1)
mtext(expression('Time to 5% up-regulation (days)'),side=2,line=3.5,cex=1)
arrows(1,actin.u05.plot[1],1,actin.u05.plot[2],angle=90,length=0.1,code=3,lwd=2,col="goldenrod1")
arrows(2,rhiz.u05.plot[1],2,rhiz.u05.plot[2],angle=90,length=0.1,code=3,lwd=2,col="orchid3")
points(1,t05.actin.est,pch=15,cex=2,lwd=2,col="goldenrod1")
points(2,t05.rhiz.est,pch=18,cex=2,lwd=2,col="orchid3")
stripchart(u05~sym, vertical = TRUE, data = u05.all.plot, 
           method = "jitter", add = TRUE, pch=1,col = 'black',cex=1)
mtext(text="a",side=3,cex=1,adj=0)

boxplot(u95~Treatment, data=u95.plot,yaxt="n",xaxt="n",ylim=c(0,200),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2),labels=c("Actinorhizal","Rhizobial"),cex.axis=1,las=2)
axis(2,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,cex.axis=1)
mtext(expression('Time to 95% up-regulation (days)'),side=2,line=3.5,cex=1)
arrows(1,actin.u95.plot[1],1,actin.u95.plot[2],angle=90,length=0.1,code=3,lwd=2,col="goldenrod1")
arrows(2,rhiz.u95.plot[1],2,rhiz.u95.plot[2],angle=90,length=0.1,code=3,lwd=2,col="orchid3")
points(1,t95.actin.est,pch=15,cex=2,lwd=2,col="goldenrod1")
points(2,t95.rhiz.est,pch=18,cex=2,lwd=2,col="orchid3")
stripchart(u95~sym, vertical = TRUE, data = u95.all.plot, 
           method = "jitter", add = TRUE, pch=1,col = 'black',cex=1)
mtext(text="b",side=3,cex=1,adj=0)

boxplot(u0595~Symbiosis, data=u0595.plot,yaxt="n",xaxt="n",ylim=c(0,80),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2),labels=c("Actinorhizal","Rhizobial"),cex.axis=1,las=2)
axis(2,at=c(0,20,40,60,80),labels=c(0,20,40,60,80),las=1,cex.axis=1.1)
mtext(expression('Time from 5% to 95% up-regulation (days)'),side=2,line=3.5,cex=1)
arrows(1,actin.u0595.plot[1],1,actin.u0595.plot[2],angle=90,length=0.1,code=3,lwd=2,col="goldenrod1")
arrows(2,rhiz.u0595.plot[1],2,rhiz.u0595.plot[2],angle=90,length=0.1,code=3,lwd=2,col="orchid3")
points(1,c(t95.actin.est - t05.actin.est),pch=15,cex=2,lwd=2,col="goldenrod1")
points(2,c(t95.rhiz.est - t05.rhiz.est),pch=18,cex=2,lwd=2,col="orchid3")
mtext(text="c",side=3,cex=1,adj=0)

boxplot(r~Treatment, data=r.plot,yaxt="n",xaxt="n",ylim=c(0,0.4),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2),labels=c("Actinorhizal","Rhizobial"),cex.axis=1,las=2)
axis(2,at=c(0,0.1,0.2,0.3,0.4),labels=c(0,0.1,0.2,0.3,0.4),las=1,cex.axis=1.1)
mtext(expression('Rate of change; '*italic(r)*' (day'^-1*')'),side=2,line=3.5,cex=1)
arrows(1,actin.r.plot[1],1,actin.r.plot[2],angle=90,length=0.1,code=3,lwd=2,col="goldenrod1")
arrows(2,rhiz.r.plot[1],2,rhiz.r.plot[2],angle=90,length=0.1,code=3,lwd=2,col="orchid3")
points(1,exp(ft_sym_up[3]),pch=15,cex=2,lwd=2,col="goldenrod1")
points(2,exp(ft_sym_up[2]),pch=18,cex=2,lwd=2,col="orchid3")
mtext(text="d",side=3,cex=1,adj=0)

### Figure S1
par(mfrow=c(1,1))
par(mar=c(5,5,2,1))
par(oma=c(1,1,1,1))
boxplot(N_ug_per_L~Factor,data=rinse,xaxt="n",ylim=c(0,5000),las=1,cex.axis=1.2,
        border=c("turquoise","orange","turquoise","orange","turquoise","orange"),xlab=NA,ylab=NA,col="white")
mtext(expression('N in KCl extract ( '*mu*'g N L'^-1*')'),side=2,line=4,cex=1.2)
text(cex=1.2,x=c(1.5,3.5,5.5),y=-250,c("Control","Fertilized","3L Rinse"),
     adj=1,xpd=T,srt=60)

### Figure S2
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G21DA$Day,rev(G21DA$Day)),y=c(G21DA$SNF.low/3.84,rev(G21DA$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DB$Day,rev(G21DB$Day)),y=c(G21DB$SNF.low/3.84,rev(G21DB$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DC$Day,rev(G21DC$Day)),y=c(G21DC$SNF.low/3.84,rev(G21DC$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DD$Day,rev(G21DD$Day)),y=c(G21DD$SNF.low/3.84,rev(G21DD$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/3.84~Day,data=G21DA,lwd=1,type="p",col="black",pch=0)
points(SNF/3.84~Day,data=G21DB,lwd=1,type="p",col="black",pch=1)
points(SNF/3.84~Day,data=G21DC,lwd=1,type="p",col="black",pch=2)
points(SNF/3.84~Day,data=G21DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid(exp(ft_biometemp_down[42])/3.84,ft_biometemp_down[3],exp(ft_biometemp_down[10]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[43])/3.84,ft_biometemp_down[3],exp(ft_biometemp_down[11]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[44])/3.84,ft_biometemp_down[3],exp(ft_biometemp_down[12]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[45])/3.84,ft_biometemp_down[3],exp(ft_biometemp_down[13]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R21DA$Day,rev(R21DA$Day)),y=c(R21DA$SNF.low/4.27,rev(R21DA$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DB$Day,rev(R21DB$Day)),y=c(R21DB$SNF.low/4.27,rev(R21DB$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DC$Day,rev(R21DC$Day)),y=c(R21DC$SNF.low/4.27,rev(R21DC$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DD$Day,rev(R21DD$Day)),y=c(R21DD$SNF.low/4.27,rev(R21DD$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/4.27~Day,data=R21DA,lwd=1,type="p",col="black",pch=0)
points(SNF/4.27~Day,data=R21DB,lwd=1,type="p",col="black",pch=1)
points(SNF/4.27~Day,data=R21DC,lwd=1,type="p",col="black",pch=2)
points(SNF/4.27~Day,data=R21DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid(exp(ft_biometemp_down[50])/4.27,ft_biometemp_down[2],exp(ft_biometemp_down[18]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[51])/4.27,ft_biometemp_down[2],exp(ft_biometemp_down[19]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[52])/4.27,ft_biometemp_down[2],exp(ft_biometemp_down[20]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[53])/4.27,ft_biometemp_down[2],exp(ft_biometemp_down[21]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M21DA$Day,rev(M21DA$Day)),y=c(M21DA$SNF.low/4.98,rev(M21DA$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DB$Day,rev(M21DB$Day)),y=c(M21DB$SNF.low/4.98,rev(M21DB$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DC$Day,rev(M21DC$Day)),y=c(M21DC$SNF.low/4.98,rev(M21DC$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DD$Day,rev(M21DD$Day)),y=c(M21DD$SNF.low/4.98,rev(M21DD$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/4.98~Day,data=M21DA,lwd=1,type="p",col="black",lty=2,pch=0)
points(SNF/4.98~Day,data=M21DB,lwd=1,type="p",col="black",lty=2,pch=1)
points(SNF/4.98~Day,data=M21DC,lwd=1,type="p",col="black",lty=2,pch=2)
points(SNF/4.98~Day,data=M21DD,lwd=1,type="p",col="black",lty=2,pch=3)
curve(sigmoid(exp(ft_biometemp_down[46])/4.98,ft_biometemp_down[3],exp(ft_biometemp_down[14]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[47])/4.98,ft_biometemp_down[3],exp(ft_biometemp_down[15]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[48])/4.98,ft_biometemp_down[3],exp(ft_biometemp_down[16]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[49])/4.98,ft_biometemp_down[3],exp(ft_biometemp_down[17]),x),from=0,to=100,lty=2,col="blue",lwd=1,add=T)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A21DA$Day,rev(A21DA$Day)),y=c(A21DA$SNF.low/3.15,rev(A21DA$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DB$Day,rev(A21DB$Day)),y=c(A21DB$SNF.low/3.15,rev(A21DB$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DC$Day,rev(A21DC$Day)),y=c(A21DC$SNF.low/3.15,rev(A21DC$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DD$Day,rev(A21DD$Day)),y=c(A21DD$SNF.low/3.15,rev(A21DD$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/3.15~Day,data=A21DA,lwd=1,type="p",col="black",pch=0)
points(SNF/3.15~Day,data=A21DB,lwd=1,type="p",col="black",pch=1)
points(SNF/3.15~Day,data=A21DC,lwd=1,type="p",col="black",pch=2)
points(SNF/3.15~Day,data=A21DD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid(exp(ft_biometemp_down[38])/3.15,ft_biometemp_down[2],exp(ft_biometemp_down[6]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[39])/3.15,ft_biometemp_down[2],exp(ft_biometemp_down[7]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[40])/3.15,ft_biometemp_down[2],exp(ft_biometemp_down[8]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[41])/3.15,ft_biometemp_down[2],exp(ft_biometemp_down[9]),x),from=0,to=100,lty=1,col="blue",lwd=1,add=T)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=T,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G31DA$Day,rev(G31DA$Day)),y=c(G31DA$SNF.low/3.84,rev(G31DA$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DB$Day,rev(G31DB$Day)),y=c(G31DB$SNF.low/3.84,rev(G31DB$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DC$Day,rev(G31DC$Day)),y=c(G31DC$SNF.low/3.84,rev(G31DC$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DD$Day,rev(G31DD$Day)),y=c(G31DD$SNF.low/3.84,rev(G31DD$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/3.84~Day,data=G31DA,lwd=1,type="p",col="black",pch=0)
points(SNF/3.84~Day,data=G31DB,lwd=1,type="p",col="black",pch=0)
points(SNF/3.84~Day,data=G31DC,lwd=1,type="p",col="black",pch=0)
points(SNF/3.84~Day,data=G31DD,lwd=1,type="p",col="black",pch=0)
curve(sigmoid(exp(ft_biometemp_down[58])/3.84,ft_biometemp_down[5],exp(ft_biometemp_down[26]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[59])/3.84,ft_biometemp_down[5],exp(ft_biometemp_down[27]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[60])/3.84,ft_biometemp_down[5],exp(ft_biometemp_down[28]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[61])/3.84,ft_biometemp_down[5],exp(ft_biometemp_down[29]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R31DA$Day,rev(R31DA$Day)),y=c(R31DA$SNF.low/4.27,rev(R31DA$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DB$Day,rev(R31DB$Day)),y=c(R31DB$SNF.low/4.27,rev(R31DB$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DC$Day,rev(R31DC$Day)),y=c(R31DC$SNF.low/4.27,rev(R31DC$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DD$Day,rev(R31DD$Day)),y=c(R31DD$SNF.low/4.27,rev(R31DD$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/4.27~Day,data=R31DA,lwd=1,type="p",col="black",pch=0)
points(SNF/4.27~Day,data=R31DB,lwd=1,type="p",col="black",pch=0)
points(SNF/4.27~Day,data=R31DC,lwd=1,type="p",col="black",pch=0)
points(SNF/4.27~Day,data=R31DD,lwd=1,type="p",col="black",pch=0)
curve(sigmoid(exp(ft_biometemp_down[66])/4.27,ft_biometemp_down[4],exp(ft_biometemp_down[34]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[67])/4.27,ft_biometemp_down[4],exp(ft_biometemp_down[35]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[68])/4.27,ft_biometemp_down[4],exp(ft_biometemp_down[36]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[69])/4.27,ft_biometemp_down[4],exp(ft_biometemp_down[37]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M31DA$Day,rev(M31DA$Day)),y=c(M31DA$SNF.low/4.98,rev(M31DA$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DB$Day,rev(M31DB$Day)),y=c(M31DB$SNF.low/4.98,rev(M31DB$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DC$Day,rev(M31DC$Day)),y=c(M31DC$SNF.low/4.98,rev(M31DC$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DD$Day,rev(M31DD$Day)),y=c(M31DD$SNF.low/4.98,rev(M31DD$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/4.98~Day,data=M31DA,lwd=1,type="p",col="black",pch=0)
points(SNF/4.98~Day,data=M31DB,lwd=1,type="p",col="black",pch=0)
points(SNF/4.98~Day,data=M31DC,lwd=1,type="p",col="black",pch=0)
points(SNF/4.98~Day,data=M31DD,lwd=1,type="p",col="black",pch=0)
curve(sigmoid(exp(ft_biometemp_down[62])/4.98,ft_biometemp_down[5],exp(ft_biometemp_down[30]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[63])/4.98,ft_biometemp_down[5],exp(ft_biometemp_down[31]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[64])/4.98,ft_biometemp_down[5],exp(ft_biometemp_down[32]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[65])/4.98,ft_biometemp_down[5],exp(ft_biometemp_down[33]),x),from=0,to=100,lty=2,col="red",lwd=1,add=T)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  h'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A31DA$Day,rev(A31DA$Day)),y=c(A31DA$SNF.low/3.15,rev(A31DA$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DB$Day,rev(A31DB$Day)),y=c(A31DB$SNF.low/3.15,rev(A31DB$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DC$Day,rev(A31DC$Day)),y=c(A31DC$SNF.low/3.15,rev(A31DC$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DD$Day,rev(A31DD$Day)),y=c(A31DD$SNF.low/3.15,rev(A31DD$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/3.15~Day,data=A31DA,lwd=1,type="p",col="black",pch=0)
points(SNF/3.15~Day,data=A31DB,lwd=1,type="p",col="black",pch=0)
points(SNF/3.15~Day,data=A31DC,lwd=1,type="p",col="black",pch=0)
points(SNF/3.15~Day,data=A31DD,lwd=1,type="p",col="black",pch=0)
curve(sigmoid(exp(ft_biometemp_down[54])/3.15,ft_biometemp_down[4],exp(ft_biometemp_down[22]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[55])/3.15,ft_biometemp_down[4],exp(ft_biometemp_down[23]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[56])/3.15,ft_biometemp_down[4],exp(ft_biometemp_down[24]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
curve(sigmoid(exp(ft_biometemp_down[57])/3.15,ft_biometemp_down[4],exp(ft_biometemp_down[25]),x),from=0,to=100,lty=1,col="red",lwd=1,add=T)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('N fixation per whole-symbiosis respiration'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(mol N'[2]*' mol CO'[2]^-1*')'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to high N (days)'),side=1,line=3,cex=1,outer=T)

### Figure S3
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G21UA$Day,rev(G21UA$Day)),y=c(G21UA$SNF.low/3.84,rev(G21UA$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UB$Day,rev(G21UB$Day)),y=c(G21UB$SNF.low/3.84,rev(G21UB$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UC$Day,rev(G21UC$Day)),y=c(G21UC$SNF.low/3.84,rev(G21UC$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/3.84~Day,data=G21UA,lwd=1,type="p",col="black",pch=0)
points(SNF/3.84~Day,data=G21UB,lwd=1,type="p",col="black",pch=1)
points(SNF/3.84~Day,data=G21UC,lwd=1,type="p",col="black",pch=2)
curve(sigmoid(exp(ft_sym_up[37])/3.84,exp(ft_sym_up[2]),exp(ft_sym_up[4]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[38])/3.84,exp(ft_sym_up[2]),exp(ft_sym_up[5]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[39])/3.84,exp(ft_sym_up[2]),exp(ft_sym_up[6]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R21UA$Day,rev(R21UA$Day)),y=c(R21UA$SNF.low/4.27,rev(R21UA$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UB$Day,rev(R21UB$Day)),y=c(R21UB$SNF.low/4.27,rev(R21UB$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UC$Day,rev(R21UC$Day)),y=c(R21UC$SNF.low/4.27,rev(R21UC$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/4.27~Day,data=R21UA,lwd=1,type="p",col="black",pch=0)
points(SNF/4.27~Day,data=R21UB,lwd=1,type="p",col="black",pch=1)
points(SNF/4.27~Day,data=R21UC,lwd=1,type="p",col="black",pch=2)
curve(sigmoid(exp(ft_sym_up[44])/4.27,exp(ft_sym_up[2]),exp(ft_sym_up[11]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[45])/4.27,exp(ft_sym_up[2]),exp(ft_sym_up[12]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[46])/4.27,exp(ft_sym_up[2]),exp(ft_sym_up[13]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M21UA$Day,rev(M21UA$Day)),y=c(M21UA$SNF.low/4.98,rev(M21UA$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UB$Day,rev(M21UB$Day)),y=c(M21UB$SNF.low/4.98,rev(M21UB$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UC$Day,rev(M21UC$Day)),y=c(M21UC$SNF.low/4.98,rev(M21UC$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UD$Day,rev(M21UD$Day)),y=c(M21UD$SNF.low/4.98,rev(M21UD$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/4.98~Day,data=M21UA,lwd=1,type="p",col="black",pch=0)
points(SNF/4.98~Day,data=M21UB,lwd=1,type="p",col="black",pch=1)
points(SNF/4.98~Day,data=M21UC,lwd=1,type="p",col="black",pch=2)
points(SNF/4.98~Day,data=M21UD,lwd=1,type="p",col="black",pch=3)
curve(sigmoid(exp(ft_sym_up[40])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[7]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[41])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[8]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[42])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[9]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[43])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[10]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.015),xaxt="n",yaxt="n",bty="n")

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=T,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G31UA$Day,rev(G31UA$Day)),y=c(G31UA$SNF.low/3.84,rev(G31UA$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UB$Day,rev(G31UB$Day)),y=c(G31UB$SNF.low/3.84,rev(G31UB$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UC$Day,rev(G31UC$Day)),y=c(G31UC$SNF.low/3.84,rev(G31UC$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UD$Day,rev(G31UD$Day)),y=c(G31UD$SNF.low/3.84,rev(G31UD$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UE$Day,rev(G31UE$Day)),y=c(G31UE$SNF.low/3.84,rev(G31UE$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UF$Day,rev(G31UF$Day)),y=c(G31UF$SNF.low/3.84,rev(G31UF$SNF.high/3.84)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/3.84~Day,data=G31UA,lwd=1,type="p",col="black",pch=0)
points(SNF/3.84~Day,data=G31UB,lwd=1,type="p",col="black",pch=1)
points(SNF/3.84~Day,data=G31UC,lwd=1,type="p",col="black",pch=2)
points(SNF/3.84~Day,data=G31UD,lwd=1,type="p",col="black",pch=3)
points(SNF/3.84~Day,data=G31UE,lwd=1,type="p",col="black",pch=4)
points(SNF/3.84~Day,data=G31UF,lwd=1,type="p",col="black",pch=5)
curve(sigmoid(exp(ft_sym_up[50])/3.84,exp(ft_sym_up[2]),exp(ft_sym_up[17]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[51])/3.84,exp(ft_sym_up[2]),exp(ft_sym_up[18]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[52])/3.84,exp(ft_sym_up[2]),exp(ft_sym_up[19]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[53])/3.84,exp(ft_sym_up[2]),exp(ft_sym_up[20]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[54])/3.84,exp(ft_sym_up[2]),exp(ft_sym_up[21]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[55])/3.84,exp(ft_sym_up[2]),exp(ft_sym_up[22]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R31UA$Day,rev(R31UA$Day)),y=c(R31UA$SNF.low/4.27,rev(R31UA$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UB$Day,rev(R31UB$Day)),y=c(R31UB$SNF.low/4.27,rev(R31UB$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UC$Day,rev(R31UC$Day)),y=c(R31UC$SNF.low/4.27,rev(R31UC$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UD$Day,rev(R31UD$Day)),y=c(R31UD$SNF.low/4.27,rev(R31UD$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UE$Day,rev(R31UE$Day)),y=c(R31UE$SNF.low/4.27,rev(R31UE$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UF$Day,rev(R31UF$Day)),y=c(R31UF$SNF.low/4.27,rev(R31UF$SNF.high/4.27)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/4.27~Day,data=R31UA,lwd=1,type="p",col="black",pch=0)
points(SNF/4.27~Day,data=R31UB,lwd=1,type="p",col="black",pch=1)
points(SNF/4.27~Day,data=R31UC,lwd=1,type="p",col="black",pch=2)
points(SNF/4.27~Day,data=R31UD,lwd=1,type="p",col="black",pch=3)
points(SNF/4.27~Day,data=R31UE,lwd=1,type="p",col="black",pch=4)
points(SNF/4.27~Day,data=R31UF,lwd=1,type="p",col="black",pch=5)
curve(sigmoid(exp(ft_sym_up[64])/4.27,exp(ft_sym_up[2]),exp(ft_sym_up[31]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[65])/4.27,exp(ft_sym_up[2]),exp(ft_sym_up[32]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[66])/4.27,exp(ft_sym_up[2]),exp(ft_sym_up[33]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[67])/4.27,exp(ft_sym_up[2]),exp(ft_sym_up[34]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[68])/4.27,exp(ft_sym_up[2]),exp(ft_sym_up[35]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[69])/4.27,exp(ft_sym_up[2]),exp(ft_sym_up[36]),x),from=0,to=200,lty=1,col="orchid3",lwd=1,add=T)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M31UA$Day,rev(M31UA$Day)),y=c(M31UA$SNF.low/4.98,rev(M31UA$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UB$Day,rev(M31UB$Day)),y=c(M31UB$SNF.low/4.98,rev(M31UB$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UC$Day,rev(M31UC$Day)),y=c(M31UC$SNF.low/4.98,rev(M31UC$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UD$Day,rev(M31UD$Day)),y=c(M31UD$SNF.low/4.98,rev(M31UD$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UE$Day,rev(M31UE$Day)),y=c(M31UE$SNF.low/4.98,rev(M31UE$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UF$Day,rev(M31UF$Day)),y=c(M31UF$SNF.low/4.98,rev(M31UF$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UG$Day,rev(M31UG$Day)),y=c(M31UG$SNF.low/4.98,rev(M31UG$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UH$Day,rev(M31UH$Day)),y=c(M31UH$SNF.low/4.98,rev(M31UH$SNF.high/4.98)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/4.98~Day,data=M31UA,lwd=1,type="p",col="black",pch=0)
points(SNF/4.98~Day,data=M31UB,lwd=1,type="p",col="black",pch=1)
points(SNF/4.98~Day,data=M31UC,lwd=1,type="p",col="black",pch=2)
points(SNF/4.98~Day,data=M31UD,lwd=1,type="p",col="black",pch=3)
points(SNF/4.98~Day,data=M31UE,lwd=1,type="p",col="black",pch=4)
points(SNF/4.98~Day,data=M31UF,lwd=1,type="p",col="black",pch=5)
points(SNF/4.98~Day,data=M31UG,lwd=1,type="p",col="black",pch=6)
points(SNF/4.98~Day,data=M31UH,lwd=1,type="p",col="black",pch=7)
curve(sigmoid(exp(ft_sym_up[56])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[23]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[57])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[24]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[58])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[25]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[59])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[26]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[60])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[27]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[61])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[28]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[62])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[29]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[63])/4.98,exp(ft_sym_up[3]),exp(ft_sym_up[30]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,0.012),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.002,0.004,0.006,0.008,0.010,0.012),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A31UA$Day,rev(A31UA$Day)),y=c(A31UA$SNF.low/3.15,rev(A31UA$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UB$Day,rev(A31UB$Day)),y=c(A31UB$SNF.low/3.15,rev(A31UB$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UC$Day,rev(A31UC$Day)),y=c(A31UC$SNF.low/3.15,rev(A31UC$SNF.high/3.15)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(SNF/3.15~Day,data=A31UA,lwd=1,type="p",col="black",pch=0)
points(SNF/3.15~Day,data=A31UB,lwd=1,type="p",col="black",pch=1)
points(SNF/3.15~Day,data=A31UC,lwd=1,type="p",col="black",pch=2)
curve(sigmoid(exp(ft_sym_up[47])/3.15,exp(ft_sym_up[3]),exp(ft_sym_up[14]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[48])/3.15,exp(ft_sym_up[3]),exp(ft_sym_up[15]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
curve(sigmoid(exp(ft_sym_up[49])/3.15,exp(ft_sym_up[3]),exp(ft_sym_up[16]),x),from=0,to=200,lty=1,col="goldenrod1",lwd=1,add=T)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('N fixation per whole-symbiosis respiration'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(mol N'[2]*' mol CO'[2]^-1*')'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to low N (days)'),side=1,line=3,cex=1,outer=T)

### Figure S4
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,5),xaxt="n",yaxt="n")
axis(1,at=c(0,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(G21DA.co2$Day,rev(G21DA.co2$Day)),y=c(G21DA.co2$Photo.low.post/G21DA.co2$Photo.mean.post[1],rev(G21DA.co2$Photo.high.post/G21DA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DB.co2$Day,rev(G21DB.co2$Day)),y=c(G21DB.co2$Photo.low.post/G21DB.co2$Photo.mean.post[1],rev(G21DB.co2$Photo.high.post/G21DB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DC.co2$Day,rev(G21DC.co2$Day)),y=c(G21DC.co2$Photo.low.post/G21DC.co2$Photo.mean.post[1],rev(G21DC.co2$Photo.high.post/G21DC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DD.co2$Day,rev(G21DD.co2$Day)),y=c(G21DD.co2$Photo.low.post/G21DD.co2$Photo.mean.post[1],rev(G21DD.co2$Photo.high.post/G21DD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,5),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(R21DA.co2$Day,rev(R21DA.co2$Day)),y=c(R21DA.co2$Photo.low.post/R21DA.co2$Photo.mean.post[1],rev(R21DA.co2$Photo.high.post/R21DA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DB.co2$Day,rev(R21DB.co2$Day)),y=c(R21DB.co2$Photo.low.post/R21DB.co2$Photo.mean.post[1],rev(R21DB.co2$Photo.high.post/R21DB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DC.co2$Day,rev(R21DC.co2$Day)),y=c(R21DC.co2$Photo.low.post/R21DC.co2$Photo.mean.post[1],rev(R21DC.co2$Photo.high.post/R21DC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DD.co2$Day,rev(R21DD.co2$Day)),y=c(R21DD.co2$Photo.low.post/R21DD.co2$Photo.mean.post[1],rev(R21DD.co2$Photo.high.post/R21DD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,5),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(M21DA.co2$Day,rev(M21DA.co2$Day)),y=c(M21DA.co2$Photo.low.post/M21DA.co2$Photo.mean.post[1],rev(M21DA.co2$Photo.high.post/M21DA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DB.co2$Day,rev(M21DB.co2$Day)),y=c(M21DB.co2$Photo.low.post/M21DB.co2$Photo.mean.post[1],rev(M21DB.co2$Photo.high.post/M21DB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DC.co2$Day,rev(M21DC.co2$Day)),y=c(M21DC.co2$Photo.low.post/M21DC.co2$Photo.mean.post[1],rev(M21DC.co2$Photo.high.post/M21DC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DD.co2$Day,rev(M21DD.co2$Day)),y=c(M21DD.co2$Photo.low.post/M21DD.co2$Photo.mean.post[1],rev(M21DD.co2$Photo.high.post/M21DD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M21DA.co2,lwd=1,type="b",col="black",lty=1,pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M21DB.co2,lwd=1,type="b",col="black",lty=1,pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M21DC.co2,lwd=1,type="b",col="black",lty=1,pch=2)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M21DD.co2,lwd=1,type="b",col="black",lty=1,pch=3)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,5),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5),labels=F,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(A21DA.co2$Day,rev(A21DA.co2$Day)),y=c(A21DA.co2$Photo.low.post/A21DA.co2$Photo.mean.post[1],rev(A21DA.co2$Photo.high.post/A21DA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DB.co2$Day,rev(A21DB.co2$Day)),y=c(A21DB.co2$Photo.low.post/A21DB.co2$Photo.mean.post[1],rev(A21DB.co2$Photo.high.post/A21DB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DC.co2$Day,rev(A21DC.co2$Day)),y=c(A21DC.co2$Photo.low.post/A21DC.co2$Photo.mean.post[1],rev(A21DC.co2$Photo.high.post/A21DC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DD.co2$Day,rev(A21DD.co2$Day)),y=c(A21DD.co2$Photo.low.post/A21DD.co2$Photo.mean.post[1],rev(A21DD.co2$Photo.high.post/A21DD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,5),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5),labels=T,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(G31DA.co2$Day,rev(G31DA.co2$Day)),y=c(G31DA.co2$Photo.low.post/G31DA.co2$Photo.mean.post[1],rev(G31DA.co2$Photo.high.post/G31DA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DB.co2$Day,rev(G31DB.co2$Day)),y=c(G31DB.co2$Photo.low.post/G31DB.co2$Photo.mean.post[1],rev(G31DB.co2$Photo.high.post/G31DB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DC.co2$Day,rev(G31DC.co2$Day)),y=c(G31DC.co2$Photo.low.post/G31DC.co2$Photo.mean.post[1],rev(G31DC.co2$Photo.high.post/G31DC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DD.co2$Day,rev(G31DD.co2$Day)),y=c(G31DD.co2$Photo.low.post/G31DD.co2$Photo.mean.post[1],rev(G31DD.co2$Photo.high.post/G31DD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,5),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(R31DA.co2$Day,rev(R31DA.co2$Day)),y=c(R31DA.co2$Photo.low.post/R31DA.co2$Photo.mean.post[1],rev(R31DA.co2$Photo.high.post/R31DA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DB.co2$Day,rev(R31DB.co2$Day)),y=c(R31DB.co2$Photo.low.post/R31DB.co2$Photo.mean.post[1],rev(R31DB.co2$Photo.high.post/R31DB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DC.co2$Day,rev(R31DC.co2$Day)),y=c(R31DC.co2$Photo.low.post/R31DC.co2$Photo.mean.post[1],rev(R31DC.co2$Photo.high.post/R31DC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DD.co2$Day,rev(R31DD.co2$Day)),y=c(R31DD.co2$Photo.low.post/R31DD.co2$Photo.mean.post[1],rev(R31DD.co2$Photo.high.post/R31DD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,5),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(M31DA.co2$Day,rev(M31DA.co2$Day)),y=c(M31DA.co2$Photo.low.post/M31DA.co2$Photo.mean.post[1],rev(M31DA.co2$Photo.high.post/M31DA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DB.co2$Day,rev(M31DB.co2$Day)),y=c(M31DB.co2$Photo.low.post/M31DB.co2$Photo.mean.post[1],rev(M31DB.co2$Photo.high.post/M31DB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DC.co2$Day,rev(M31DC.co2$Day)),y=c(M31DC.co2$Photo.low.post/M31DC.co2$Photo.mean.post[1],rev(M31DC.co2$Photo.high.post/M31DC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DD.co2$Day,rev(M31DD.co2$Day)),y=c(M31DD.co2$Photo.low.post/M31DD.co2$Photo.mean.post[1],rev(M31DD.co2$Photo.high.post/M31DD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,5),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5),labels=F,cex.axis=1.2,las=1)
title(main=expression('  h'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(A31DA.co2$Day,rev(A31DA.co2$Day)),y=c(A31DA.co2$Photo.low.post/A31DA.co2$Photo.mean.post[1],rev(A31DA.co2$Photo.high.post/A31DA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DB.co2$Day,rev(A31DB.co2$Day)),y=c(A31DB.co2$Photo.low.post/A31DB.co2$Photo.mean.post[1],rev(A31DB.co2$Photo.high.post/A31DB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DC.co2$Day,rev(A31DC.co2$Day)),y=c(A31DC.co2$Photo.low.post/A31DC.co2$Photo.mean.post[1],rev(A31DC.co2$Photo.high.post/A31DC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DD.co2$Day,rev(A31DD.co2$Day)),y=c(A31DD.co2$Photo.low.post/A31DD.co2$Photo.mean.post[1],rev(A31DD.co2$Photo.high.post/A31DD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A31DD.co2,lwd=1,type="b",col="black",pch=0)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('Apparent whole-plant photosynthesis'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(normalized to 1 at day = 0)'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to high N (days)'),side=1,line=3,cex=1,outer=T)

### Figure S5
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G21DA.co2$Day,rev(G21DA.co2$Day)),y=c(G21DA.co2$Photo.low.post,rev(G21DA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DB.co2$Day,rev(G21DB.co2$Day)),y=c(G21DB.co2$Photo.low.post,rev(G21DB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DC.co2$Day,rev(G21DC.co2$Day)),y=c(G21DC.co2$Photo.low.post,rev(G21DC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DD.co2$Day,rev(G21DD.co2$Day)),y=c(G21DD.co2$Photo.low.post,rev(G21DD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=G21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=G21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=G21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post~Day,data=G21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R21DA.co2$Day,rev(R21DA.co2$Day)),y=c(R21DA.co2$Photo.low.post,rev(R21DA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DB.co2$Day,rev(R21DB.co2$Day)),y=c(R21DB.co2$Photo.low.post,rev(R21DB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DC.co2$Day,rev(R21DC.co2$Day)),y=c(R21DC.co2$Photo.low.post,rev(R21DC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DD.co2$Day,rev(R21DD.co2$Day)),y=c(R21DD.co2$Photo.low.post,rev(R21DD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=R21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=R21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=R21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post~Day,data=R21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M21DA.co2$Day,rev(M21DA.co2$Day)),y=c(M21DA.co2$Photo.low.post,rev(M21DA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DB.co2$Day,rev(M21DB.co2$Day)),y=c(M21DB.co2$Photo.low.post,rev(M21DB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DC.co2$Day,rev(M21DC.co2$Day)),y=c(M21DC.co2$Photo.low.post,rev(M21DC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DD.co2$Day,rev(M21DD.co2$Day)),y=c(M21DD.co2$Photo.low.post,rev(M21DD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=M21DA.co2,lwd=1,type="b",col="black",lty=1,pch=0)
points(Photo.mean.post~Day,data=M21DB.co2,lwd=1,type="b",col="black",lty=1,pch=1)
points(Photo.mean.post~Day,data=M21DC.co2,lwd=1,type="b",col="black",lty=1,pch=2)
points(Photo.mean.post~Day,data=M21DD.co2,lwd=1,type="b",col="black",lty=1,pch=3)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A21DA.co2$Day,rev(A21DA.co2$Day)),y=c(A21DA.co2$Photo.low.post,rev(A21DA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DB.co2$Day,rev(A21DB.co2$Day)),y=c(A21DB.co2$Photo.low.post,rev(A21DB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DC.co2$Day,rev(A21DC.co2$Day)),y=c(A21DC.co2$Photo.low.post,rev(A21DC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DD.co2$Day,rev(A21DD.co2$Day)),y=c(A21DD.co2$Photo.low.post,rev(A21DD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=A21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=A21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=A21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post~Day,data=A21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=T,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G31DA.co2$Day,rev(G31DA.co2$Day)),y=c(G31DA.co2$Photo.low.post,rev(G31DA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DB.co2$Day,rev(G31DB.co2$Day)),y=c(G31DB.co2$Photo.low.post,rev(G31DB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DC.co2$Day,rev(G31DC.co2$Day)),y=c(G31DC.co2$Photo.low.post,rev(G31DC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DD.co2$Day,rev(G31DD.co2$Day)),y=c(G31DD.co2$Photo.low.post,rev(G31DD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=G31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=G31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=G31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=G31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R31DA.co2$Day,rev(R31DA.co2$Day)),y=c(R31DA.co2$Photo.low.post,rev(R31DA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DB.co2$Day,rev(R31DB.co2$Day)),y=c(R31DB.co2$Photo.low.post,rev(R31DB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DC.co2$Day,rev(R31DC.co2$Day)),y=c(R31DC.co2$Photo.low.post,rev(R31DC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DD.co2$Day,rev(R31DD.co2$Day)),y=c(R31DD.co2$Photo.low.post,rev(R31DD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=R31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=R31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=R31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=R31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M31DA.co2$Day,rev(M31DA.co2$Day)),y=c(M31DA.co2$Photo.low.post,rev(M31DA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DB.co2$Day,rev(M31DB.co2$Day)),y=c(M31DB.co2$Photo.low.post,rev(M31DB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DC.co2$Day,rev(M31DC.co2$Day)),y=c(M31DC.co2$Photo.low.post,rev(M31DC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DD.co2$Day,rev(M31DD.co2$Day)),y=c(M31DD.co2$Photo.low.post,rev(M31DD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=M31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=M31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=M31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=M31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  h'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A31DA.co2$Day,rev(A31DA.co2$Day)),y=c(A31DA.co2$Photo.low.post,rev(A31DA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DB.co2$Day,rev(A31DB.co2$Day)),y=c(A31DB.co2$Photo.low.post,rev(A31DB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DC.co2$Day,rev(A31DC.co2$Day)),y=c(A31DC.co2$Photo.low.post,rev(A31DC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DD.co2$Day,rev(A31DD.co2$Day)),y=c(A31DD.co2$Photo.low.post,rev(A31DD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=A31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=A31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=A31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=A31DD.co2,lwd=1,type="b",col="black",pch=0)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('Apparent whole-plant photosynthesis'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(nmol CO'[2]*' s'^-1*')'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to high N (days)'),side=1,line=3,cex=1,outer=T)

### Figure S6
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,7),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5,6,7),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(G21UA.co2$Day,rev(G21UA.co2$Day)),y=c(G21UA.co2$Photo.low.post/G21UA.co2$Photo.mean.post[1],rev(G21UA.co2$Photo.high.post/G21UA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UB.co2$Day,rev(G21UB.co2$Day)),y=c(G21UB.co2$Photo.low.post/G21UB.co2$Photo.mean.post[1],rev(G21UB.co2$Photo.high.post/G21UB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UC.co2$Day,rev(G21UC.co2$Day)),y=c(G21UC.co2$Photo.low.post/G21UC.co2$Photo.mean.post[1],rev(G21UC.co2$Photo.high.post/G21UC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G21UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,7),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5,6,7),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(R21UA.co2$Day,rev(R21UA.co2$Day)),y=c(R21UA.co2$Photo.low.post/R21UA.co2$Photo.mean.post[1],rev(R21UA.co2$Photo.high.post/R21UA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UB.co2$Day,rev(R21UB.co2$Day)),y=c(R21UB.co2$Photo.low.post/R21UB.co2$Photo.mean.post[1],rev(R21UB.co2$Photo.high.post/R21UB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UC.co2$Day,rev(R21UC.co2$Day)),y=c(R21UC.co2$Photo.low.post/R21UC.co2$Photo.mean.post[1],rev(R21UC.co2$Photo.high.post/R21UC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R21UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,7),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5,6,7),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(M21UA.co2$Day,rev(M21UA.co2$Day)),y=c(M21UA.co2$Photo.low.post/M21UA.co2$Photo.mean.post[1],rev(M21UA.co2$Photo.high.post/M21UA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UB.co2$Day,rev(M21UB.co2$Day)),y=c(M21UB.co2$Photo.low.post/M21UB.co2$Photo.mean.post[1],rev(M21UB.co2$Photo.high.post/M21UB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UC.co2$Day,rev(M21UC.co2$Day)),y=c(M21UC.co2$Photo.low.post/M21UC.co2$Photo.mean.post[1],rev(M21UC.co2$Photo.high.post/M21UC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UD.co2$Day,rev(M21UD.co2$Day)),y=c(M21UD.co2$Photo.low.post/M21UD.co2$Photo.mean.post[1],rev(M21UD.co2$Photo.high.post/M21UD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M21UC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M21UD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,700),xaxt="n",yaxt="n",bty="n")

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,7),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5,6,7),labels=T,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(G31UA.co2$Day,rev(G31UA.co2$Day)),y=c(G31UA.co2$Photo.low.post/G31UA.co2$Photo.mean.post[1],rev(G31UA.co2$Photo.high.post/G31UA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UB.co2$Day,rev(G31UB.co2$Day)),y=c(G31UB.co2$Photo.low.post/G31UB.co2$Photo.mean.post[1],rev(G31UB.co2$Photo.high.post/G31UB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UC.co2$Day,rev(G31UC.co2$Day)),y=c(G31UC.co2$Photo.low.post/G31UC.co2$Photo.mean.post[1],rev(G31UC.co2$Photo.high.post/G31UC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UD.co2$Day,rev(G31UD.co2$Day)),y=c(G31UD.co2$Photo.low.post/G31UD.co2$Photo.mean.post[1],rev(G31UD.co2$Photo.high.post/G31UD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UE.co2$Day,rev(G31UE.co2$Day)),y=c(G31UE.co2$Photo.low.post/G31UE.co2$Photo.mean.post[1],rev(G31UE.co2$Photo.high.post/G31UE.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UF.co2$Day,rev(G31UF.co2$Day)),y=c(G31UF.co2$Photo.low.post/G31UF.co2$Photo.mean.post[1],rev(G31UF.co2$Photo.high.post/G31UF.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=G31UF.co2,lwd=1,type="b",col="black",pch=5)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,7),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5,6,7),labels=F,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(R31UA.co2$Day,rev(R31UA.co2$Day)),y=c(R31UA.co2$Photo.low.post/R31UA.co2$Photo.mean.post[1],rev(R31UA.co2$Photo.high.post/R31UA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UB.co2$Day,rev(R31UB.co2$Day)),y=c(R31UB.co2$Photo.low.post/R31UB.co2$Photo.mean.post[1],rev(R31UB.co2$Photo.high.post/R31UB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UC.co2$Day,rev(R31UC.co2$Day)),y=c(R31UC.co2$Photo.low.post/R31UC.co2$Photo.mean.post[1],rev(R31UC.co2$Photo.high.post/R31UC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UD.co2$Day,rev(R31UD.co2$Day)),y=c(R31UD.co2$Photo.low.post/R31UD.co2$Photo.mean.post[1],rev(R31UD.co2$Photo.high.post/R31UD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UE.co2$Day,rev(R31UE.co2$Day)),y=c(R31UE.co2$Photo.low.post/R31UE.co2$Photo.mean.post[1],rev(R31UE.co2$Photo.high.post/R31UE.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UF.co2$Day,rev(R31UF.co2$Day)),y=c(R31UF.co2$Photo.low.post/R31UF.co2$Photo.mean.post[1],rev(R31UF.co2$Photo.high.post/R31UF.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=R31UF.co2,lwd=1,type="b",col="black",pch=5)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,7),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5,6,7),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(M31UA.co2$Day,rev(M31UA.co2$Day)),y=c(M31UA.co2$Photo.low.post/M31UA.co2$Photo.mean.post[1],rev(M31UA.co2$Photo.high.post/M31UA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UB.co2$Day,rev(M31UB.co2$Day)),y=c(M31UB.co2$Photo.low.post/M31UB.co2$Photo.mean.post[1],rev(M31UB.co2$Photo.high.post/M31UB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UC.co2$Day,rev(M31UC.co2$Day)),y=c(M31UC.co2$Photo.low.post/M31UC.co2$Photo.mean.post[1],rev(M31UC.co2$Photo.high.post/M31UC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UD.co2$Day,rev(M31UD.co2$Day)),y=c(M31UD.co2$Photo.low.post/M31UD.co2$Photo.mean.post[1],rev(M31UD.co2$Photo.high.post/M31UD.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UE.co2$Day,rev(M31UE.co2$Day)),y=c(M31UE.co2$Photo.low.post/M31UE.co2$Photo.mean.post[1],rev(M31UE.co2$Photo.high.post/M31UE.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UF.co2$Day,rev(M31UF.co2$Day)),y=c(M31UF.co2$Photo.low.post/M31UF.co2$Photo.mean.post[1],rev(M31UF.co2$Photo.high.post/M31UF.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UG.co2$Day,rev(M31UG.co2$Day)),y=c(M31UG.co2$Photo.low.post/M31UG.co2$Photo.mean.post[1],rev(M31UG.co2$Photo.high.post/M31UG.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UH.co2$Day,rev(M31UH.co2$Day)),y=c(M31UH.co2$Photo.low.post/M31UH.co2$Photo.mean.post[1],rev(M31UH.co2$Photo.high.post/M31UH.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31UF.co2,lwd=1,type="b",col="black",pch=5)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31UG.co2,lwd=1,type="b",col="black",pch=6)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=M31UH.co2,lwd=1,type="b",col="black",pch=7)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,7),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,1,2,3,4,5,6,7),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(A31UA.co2$Day,rev(A31UA.co2$Day)),y=c(A31UA.co2$Photo.low.post/A31UA.co2$Photo.mean.post[1],rev(A31UA.co2$Photo.high.post/A31UA.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UB.co2$Day,rev(A31UB.co2$Day)),y=c(A31UB.co2$Photo.low.post/A31UB.co2$Photo.mean.post[1],rev(A31UB.co2$Photo.high.post/A31UB.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UC.co2$Day,rev(A31UC.co2$Day)),y=c(A31UC.co2$Photo.low.post/A31UC.co2$Photo.mean.post[1],rev(A31UC.co2$Photo.high.post/A31UC.co2$Photo.mean.post[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post/Photo.mean.post[1]~Day,data=A31UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('Apparent whole-plant photosynthesis'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(normalized to 1 at day = 0)'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to low N (days)'),side=1,line=3,cex=1,outer=T)

### Figure S7
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G21UA.co2$Day,rev(G21UA.co2$Day)),y=c(G21UA.co2$Photo.low.post,rev(G21UA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UB.co2$Day,rev(G21UB.co2$Day)),y=c(G21UB.co2$Photo.low.post,rev(G21UB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UC.co2$Day,rev(G21UC.co2$Day)),y=c(G21UC.co2$Photo.low.post,rev(G21UC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=G21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=G21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=G21UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R21UA.co2$Day,rev(R21UA.co2$Day)),y=c(R21UA.co2$Photo.low.post,rev(R21UA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UB.co2$Day,rev(R21UB.co2$Day)),y=c(R21UB.co2$Photo.low.post,rev(R21UB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UC.co2$Day,rev(R21UC.co2$Day)),y=c(R21UC.co2$Photo.low.post,rev(R21UC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=R21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=R21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=R21UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M21UA.co2$Day,rev(M21UA.co2$Day)),y=c(M21UA.co2$Photo.low.post,rev(M21UA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UB.co2$Day,rev(M21UB.co2$Day)),y=c(M21UB.co2$Photo.low.post,rev(M21UB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UC.co2$Day,rev(M21UC.co2$Day)),y=c(M21UC.co2$Photo.low.post,rev(M21UC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UD.co2$Day,rev(M21UD.co2$Day)),y=c(M21UD.co2$Photo.low.post,rev(M21UD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=M21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=M21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=M21UC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post~Day,data=M21UD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n",bty="n")

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=T,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G31UA.co2$Day,rev(G31UA.co2$Day)),y=c(G31UA.co2$Photo.low.post,rev(G31UA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UB.co2$Day,rev(G31UB.co2$Day)),y=c(G31UB.co2$Photo.low.post,rev(G31UB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UC.co2$Day,rev(G31UC.co2$Day)),y=c(G31UC.co2$Photo.low.post,rev(G31UC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UD.co2$Day,rev(G31UD.co2$Day)),y=c(G31UD.co2$Photo.low.post,rev(G31UD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UE.co2$Day,rev(G31UE.co2$Day)),y=c(G31UE.co2$Photo.low.post,rev(G31UE.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UF.co2$Day,rev(G31UF.co2$Day)),y=c(G31UF.co2$Photo.low.post,rev(G31UF.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=G31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=G31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=G31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post~Day,data=G31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Photo.mean.post~Day,data=G31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Photo.mean.post~Day,data=G31UF.co2,lwd=1,type="b",col="black",pch=5)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R31UA.co2$Day,rev(R31UA.co2$Day)),y=c(R31UA.co2$Photo.low.post,rev(R31UA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UB.co2$Day,rev(R31UB.co2$Day)),y=c(R31UB.co2$Photo.low.post,rev(R31UB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UC.co2$Day,rev(R31UC.co2$Day)),y=c(R31UC.co2$Photo.low.post,rev(R31UC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UD.co2$Day,rev(R31UD.co2$Day)),y=c(R31UD.co2$Photo.low.post,rev(R31UD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UE.co2$Day,rev(R31UE.co2$Day)),y=c(R31UE.co2$Photo.low.post,rev(R31UE.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UF.co2$Day,rev(R31UF.co2$Day)),y=c(R31UF.co2$Photo.low.post,rev(R31UF.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=R31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=R31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=R31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post~Day,data=R31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Photo.mean.post~Day,data=R31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Photo.mean.post~Day,data=R31UF.co2,lwd=1,type="b",col="black",pch=5)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M31UA.co2$Day,rev(M31UA.co2$Day)),y=c(M31UA.co2$Photo.low.post,rev(M31UA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UB.co2$Day,rev(M31UB.co2$Day)),y=c(M31UB.co2$Photo.low.post,rev(M31UB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UC.co2$Day,rev(M31UC.co2$Day)),y=c(M31UC.co2$Photo.low.post,rev(M31UC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UD.co2$Day,rev(M31UD.co2$Day)),y=c(M31UD.co2$Photo.low.post,rev(M31UD.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UE.co2$Day,rev(M31UE.co2$Day)),y=c(M31UE.co2$Photo.low.post,rev(M31UE.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UF.co2$Day,rev(M31UF.co2$Day)),y=c(M31UF.co2$Photo.low.post,rev(M31UF.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UG.co2$Day,rev(M31UG.co2$Day)),y=c(M31UG.co2$Photo.low.post,rev(M31UG.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UH.co2$Day,rev(M31UH.co2$Day)),y=c(M31UH.co2$Photo.low.post,rev(M31UH.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=M31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=M31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=M31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Photo.mean.post~Day,data=M31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Photo.mean.post~Day,data=M31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Photo.mean.post~Day,data=M31UF.co2,lwd=1,type="b",col="black",pch=5)
points(Photo.mean.post~Day,data=M31UG.co2,lwd=1,type="b",col="black",pch=6)
points(Photo.mean.post~Day,data=M31UH.co2,lwd=1,type="b",col="black",pch=7)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,600),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,100,200,300,400,500,600),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A31UA.co2$Day,rev(A31UA.co2$Day)),y=c(A31UA.co2$Photo.low.post,rev(A31UA.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UB.co2$Day,rev(A31UB.co2$Day)),y=c(A31UB.co2$Photo.low.post,rev(A31UB.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UC.co2$Day,rev(A31UC.co2$Day)),y=c(A31UC.co2$Photo.low.post,rev(A31UC.co2$Photo.high.post)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Photo.mean.post~Day,data=A31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Photo.mean.post~Day,data=A31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Photo.mean.post~Day,data=A31UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('Apparent whole-plant photosynthesis'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(nmol CO'[2]*' s'^-1*')'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to low N (days)'),side=1,line=3,cex=1,outer=T)

### Figure S8
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(G21DA.co2$Day,rev(G21DA.co2$Day)),y=c(G21DA.co2$Resp.low/G21DA.co2$Resp.mean[1],rev(G21DA.co2$Resp.high/G21DA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DB.co2$Day,rev(G21DB.co2$Day)),y=c(G21DB.co2$Resp.low/G21DB.co2$Resp.mean[1],rev(G21DB.co2$Resp.high/G21DB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DC.co2$Day,rev(G21DC.co2$Day)),y=c(G21DC.co2$Resp.low/G21DC.co2$Resp.mean[1],rev(G21DC.co2$Resp.high/G21DC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DD.co2$Day,rev(G21DD.co2$Day)),y=c(G21DD.co2$Resp.low/G21DD.co2$Resp.mean[1],rev(G21DD.co2$Resp.high/G21DD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=G21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=G21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=G21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean/Resp.mean[1]~Day,data=G21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(R21DA.co2$Day,rev(R21DA.co2$Day)),y=c(R21DA.co2$Resp.low/R21DA.co2$Resp.mean[1],rev(R21DA.co2$Resp.high/R21DA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DB.co2$Day,rev(R21DB.co2$Day)),y=c(R21DB.co2$Resp.low/R21DB.co2$Resp.mean[1],rev(R21DB.co2$Resp.high/R21DB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DC.co2$Day,rev(R21DC.co2$Day)),y=c(R21DC.co2$Resp.low/R21DC.co2$Resp.mean[1],rev(R21DC.co2$Resp.high/R21DC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DD.co2$Day,rev(R21DD.co2$Day)),y=c(R21DD.co2$Resp.low/R21DD.co2$Resp.mean[1],rev(R21DD.co2$Resp.high/R21DD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=R21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=R21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=R21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean/Resp.mean[1]~Day,data=R21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(M21DA.co2$Day,rev(M21DA.co2$Day)),y=c(M21DA.co2$Resp.low/M21DA.co2$Resp.mean[1],rev(M21DA.co2$Resp.high/M21DA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DB.co2$Day,rev(M21DB.co2$Day)),y=c(M21DB.co2$Resp.low/M21DB.co2$Resp.mean[1],rev(M21DB.co2$Resp.high/M21DB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DC.co2$Day,rev(M21DC.co2$Day)),y=c(M21DC.co2$Resp.low/M21DC.co2$Resp.mean[1],rev(M21DC.co2$Resp.high/M21DC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DD.co2$Day,rev(M21DD.co2$Day)),y=c(M21DD.co2$Resp.low/M21DD.co2$Resp.mean[1],rev(M21DD.co2$Resp.high/M21DD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=M21DA.co2,lwd=1,type="b",col="black",lty=1,pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=M21DB.co2,lwd=1,type="b",col="black",lty=1,pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=M21DC.co2,lwd=1,type="b",col="black",lty=1,pch=2)
points(Resp.mean/Resp.mean[1]~Day,data=M21DD.co2,lwd=1,type="b",col="black",lty=1,pch=3)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(A21DA.co2$Day,rev(A21DA.co2$Day)),y=c(A21DA.co2$Resp.low/A21DA.co2$Resp.mean[1],rev(A21DA.co2$Resp.high/A21DA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DB.co2$Day,rev(A21DB.co2$Day)),y=c(A21DB.co2$Resp.low/A21DB.co2$Resp.mean[1],rev(A21DB.co2$Resp.high/A21DB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DC.co2$Day,rev(A21DC.co2$Day)),y=c(A21DC.co2$Resp.low/A21DC.co2$Resp.mean[1],rev(A21DC.co2$Resp.high/A21DC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DD.co2$Day,rev(A21DD.co2$Day)),y=c(A21DD.co2$Resp.low/A21DD.co2$Resp.mean[1],rev(A21DD.co2$Resp.high/A21DD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=A21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=A21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=A21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean/Resp.mean[1]~Day,data=A21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=T,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(G31DA.co2$Day,rev(G31DA.co2$Day)),y=c(G31DA.co2$Resp.low/G31DA.co2$Resp.mean[1],rev(G31DA.co2$Resp.high/G31DA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DB.co2$Day,rev(G31DB.co2$Day)),y=c(G31DB.co2$Resp.low/G31DB.co2$Resp.mean[1],rev(G31DB.co2$Resp.high/G31DB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DC.co2$Day,rev(G31DC.co2$Day)),y=c(G31DC.co2$Resp.low/G31DC.co2$Resp.mean[1],rev(G31DC.co2$Resp.high/G31DC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DD.co2$Day,rev(G31DD.co2$Day)),y=c(G31DD.co2$Resp.low/G31DD.co2$Resp.mean[1],rev(G31DD.co2$Resp.high/G31DD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=G31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=G31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=G31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=G31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(R31DA.co2$Day,rev(R31DA.co2$Day)),y=c(R31DA.co2$Resp.low/R31DA.co2$Resp.mean[1],rev(R31DA.co2$Resp.high/R31DA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DB.co2$Day,rev(R31DB.co2$Day)),y=c(R31DB.co2$Resp.low/R31DB.co2$Resp.mean[1],rev(R31DB.co2$Resp.high/R31DB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DC.co2$Day,rev(R31DC.co2$Day)),y=c(R31DC.co2$Resp.low/R31DC.co2$Resp.mean[1],rev(R31DC.co2$Resp.high/R31DC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DD.co2$Day,rev(R31DD.co2$Day)),y=c(R31DD.co2$Resp.low/R31DD.co2$Resp.mean[1],rev(R31DD.co2$Resp.high/R31DD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=R31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=R31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=R31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=R31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(M31DA.co2$Day,rev(M31DA.co2$Day)),y=c(M31DA.co2$Resp.low/M31DA.co2$Resp.mean[1],rev(M31DA.co2$Resp.high/M31DA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DB.co2$Day,rev(M31DB.co2$Day)),y=c(M31DB.co2$Resp.low/M31DB.co2$Resp.mean[1],rev(M31DB.co2$Resp.high/M31DB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DC.co2$Day,rev(M31DC.co2$Day)),y=c(M31DC.co2$Resp.low/M31DC.co2$Resp.mean[1],rev(M31DC.co2$Resp.high/M31DC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DD.co2$Day,rev(M31DD.co2$Day)),y=c(M31DD.co2$Resp.low/M31DD.co2$Resp.mean[1],rev(M31DD.co2$Resp.high/M31DD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=M31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=M31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=M31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=M31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  h'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(A31DA.co2$Day,rev(A31DA.co2$Day)),y=c(A31DA.co2$Resp.low/A31DA.co2$Resp.mean[1],rev(A31DA.co2$Resp.high/A31DA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DB.co2$Day,rev(A31DB.co2$Day)),y=c(A31DB.co2$Resp.low/A31DB.co2$Resp.mean[1],rev(A31DB.co2$Resp.high/A31DB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DC.co2$Day,rev(A31DC.co2$Day)),y=c(A31DC.co2$Resp.low/A31DC.co2$Resp.mean[1],rev(A31DC.co2$Resp.high/A31DC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DD.co2$Day,rev(A31DD.co2$Day)),y=c(A31DD.co2$Resp.low/A31DD.co2$Resp.mean[1],rev(A31DD.co2$Resp.high/A31DD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=A31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=A31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=A31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=A31DD.co2,lwd=1,type="b",col="black",pch=0)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('Whole-symbiosis respiration'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(normalized to 1 at day = 0)'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to high N (days)'),side=1,line=3,cex=1,outer=T)

### Figure S9
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G21DA.co2$Day,rev(G21DA.co2$Day)),y=c(G21DA.co2$Resp.low,rev(G21DA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DB.co2$Day,rev(G21DB.co2$Day)),y=c(G21DB.co2$Resp.low,rev(G21DB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DC.co2$Day,rev(G21DC.co2$Day)),y=c(G21DC.co2$Resp.low,rev(G21DC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21DD.co2$Day,rev(G21DD.co2$Day)),y=c(G21DD.co2$Resp.low,rev(G21DD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=G21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=G21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=G21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean~Day,data=G21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R21DA.co2$Day,rev(R21DA.co2$Day)),y=c(R21DA.co2$Resp.low,rev(R21DA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DB.co2$Day,rev(R21DB.co2$Day)),y=c(R21DB.co2$Resp.low,rev(R21DB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DC.co2$Day,rev(R21DC.co2$Day)),y=c(R21DC.co2$Resp.low,rev(R21DC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21DD.co2$Day,rev(R21DD.co2$Day)),y=c(R21DD.co2$Resp.low,rev(R21DD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=R21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=R21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=R21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean~Day,data=R21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M21DA.co2$Day,rev(M21DA.co2$Day)),y=c(M21DA.co2$Resp.low,rev(M21DA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DB.co2$Day,rev(M21DB.co2$Day)),y=c(M21DB.co2$Resp.low,rev(M21DB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DC.co2$Day,rev(M21DC.co2$Day)),y=c(M21DC.co2$Resp.low,rev(M21DC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21DD.co2$Day,rev(M21DD.co2$Day)),y=c(M21DD.co2$Resp.low,rev(M21DD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=M21DA.co2,lwd=1,type="b",col="black",lty=1,pch=0)
points(Resp.mean~Day,data=M21DB.co2,lwd=1,type="b",col="black",lty=1,pch=1)
points(Resp.mean~Day,data=M21DC.co2,lwd=1,type="b",col="black",lty=1,pch=2)
points(Resp.mean~Day,data=M21DD.co2,lwd=1,type="b",col="black",lty=1,pch=3)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A21DA.co2$Day,rev(A21DA.co2$Day)),y=c(A21DA.co2$Resp.low,rev(A21DA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DB.co2$Day,rev(A21DB.co2$Day)),y=c(A21DB.co2$Resp.low,rev(A21DB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DC.co2$Day,rev(A21DC.co2$Day)),y=c(A21DC.co2$Resp.low,rev(A21DC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A21DD.co2$Day,rev(A21DD.co2$Day)),y=c(A21DD.co2$Resp.low,rev(A21DD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=A21DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=A21DB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=A21DC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean~Day,data=A21DD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G31DA.co2$Day,rev(G31DA.co2$Day)),y=c(G31DA.co2$Resp.low,rev(G31DA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DB.co2$Day,rev(G31DB.co2$Day)),y=c(G31DB.co2$Resp.low,rev(G31DB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DC.co2$Day,rev(G31DC.co2$Day)),y=c(G31DC.co2$Resp.low,rev(G31DC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31DD.co2$Day,rev(G31DD.co2$Day)),y=c(G31DD.co2$Resp.low,rev(G31DD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=G31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=G31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=G31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=G31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R31DA.co2$Day,rev(R31DA.co2$Day)),y=c(R31DA.co2$Resp.low,rev(R31DA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DB.co2$Day,rev(R31DB.co2$Day)),y=c(R31DB.co2$Resp.low,rev(R31DB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DC.co2$Day,rev(R31DC.co2$Day)),y=c(R31DC.co2$Resp.low,rev(R31DC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31DD.co2$Day,rev(R31DD.co2$Day)),y=c(R31DD.co2$Resp.low,rev(R31DD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=R31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=R31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=R31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=R31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M31DA.co2$Day,rev(M31DA.co2$Day)),y=c(M31DA.co2$Resp.low,rev(M31DA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DB.co2$Day,rev(M31DB.co2$Day)),y=c(M31DB.co2$Resp.low,rev(M31DB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DC.co2$Day,rev(M31DC.co2$Day)),y=c(M31DC.co2$Resp.low,rev(M31DC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31DD.co2$Day,rev(M31DD.co2$Day)),y=c(M31DD.co2$Resp.low,rev(M31DD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=M31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=M31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=M31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=M31DD.co2,lwd=1,type="b",col="black",pch=0)

plot(0:100,(0:100)/100,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  h'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A31DA.co2$Day,rev(A31DA.co2$Day)),y=c(A31DA.co2$Resp.low,rev(A31DA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DB.co2$Day,rev(A31DB.co2$Day)),y=c(A31DB.co2$Resp.low,rev(A31DB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DC.co2$Day,rev(A31DC.co2$Day)),y=c(A31DC.co2$Resp.low,rev(A31DC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31DD.co2$Day,rev(A31DD.co2$Day)),y=c(A31DD.co2$Resp.low,rev(A31DD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=A31DA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=A31DB.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=A31DC.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=A31DD.co2,lwd=1,type="b",col="black",pch=0)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('Whole-symbiosis respiration'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(nmol CO'[2]*' s'^-1*')'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to high N (days)'),side=1,line=3,cex=1,outer=T)

### Figure S10
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(G21UA.co2$Day,rev(G21UA.co2$Day)),y=c(G21UA.co2$Resp.low/G21UA.co2$Resp.mean[1],rev(G21UA.co2$Resp.high/G21UA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UB.co2$Day,rev(G21UB.co2$Day)),y=c(G21UB.co2$Resp.low/G21UB.co2$Resp.mean[1],rev(G21UB.co2$Resp.high/G21UB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UC.co2$Day,rev(G21UC.co2$Day)),y=c(G21UC.co2$Resp.low/G21UC.co2$Resp.mean[1],rev(G21UC.co2$Resp.high/G21UC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=G21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=G21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=G21UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(R21UA.co2$Day,rev(R21UA.co2$Day)),y=c(R21UA.co2$Resp.low/R21UA.co2$Resp.mean[1],rev(R21UA.co2$Resp.high/R21UA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UB.co2$Day,rev(R21UB.co2$Day)),y=c(R21UB.co2$Resp.low/R21UB.co2$Resp.mean[1],rev(R21UB.co2$Resp.high/R21UB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UC.co2$Day,rev(R21UC.co2$Day)),y=c(R21UC.co2$Resp.low/R21UC.co2$Resp.mean[1],rev(R21UC.co2$Resp.high/R21UC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=R21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=R21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=R21UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(M21UA.co2$Day,rev(M21UA.co2$Day)),y=c(M21UA.co2$Resp.low/M21UA.co2$Resp.mean[1],rev(M21UA.co2$Resp.high/M21UA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UB.co2$Day,rev(M21UB.co2$Day)),y=c(M21UB.co2$Resp.low/M21UB.co2$Resp.mean[1],rev(M21UB.co2$Resp.high/M21UB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UC.co2$Day,rev(M21UC.co2$Day)),y=c(M21UC.co2$Resp.low/M21UC.co2$Resp.mean[1],rev(M21UC.co2$Resp.high/M21UC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UD.co2$Day,rev(M21UD.co2$Day)),y=c(M21UD.co2$Resp.low/M21UD.co2$Resp.mean[1],rev(M21UD.co2$Resp.high/M21UD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=M21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=M21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=M21UC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean/Resp.mean[1]~Day,data=M21UD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,300),xaxt="n",yaxt="n",bty="n")

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=T,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(G31UA.co2$Day,rev(G31UA.co2$Day)),y=c(G31UA.co2$Resp.low/G31UA.co2$Resp.mean[1],rev(G31UA.co2$Resp.high/G31UA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UB.co2$Day,rev(G31UB.co2$Day)),y=c(G31UB.co2$Resp.low/G31UB.co2$Resp.mean[1],rev(G31UB.co2$Resp.high/G31UB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UC.co2$Day,rev(G31UC.co2$Day)),y=c(G31UC.co2$Resp.low/G31UC.co2$Resp.mean[1],rev(G31UC.co2$Resp.high/G31UC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UD.co2$Day,rev(G31UD.co2$Day)),y=c(G31UD.co2$Resp.low/G31UD.co2$Resp.mean[1],rev(G31UD.co2$Resp.high/G31UD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UE.co2$Day,rev(G31UE.co2$Day)),y=c(G31UE.co2$Resp.low/G31UE.co2$Resp.mean[1],rev(G31UE.co2$Resp.high/G31UE.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UF.co2$Day,rev(G31UF.co2$Day)),y=c(G31UF.co2$Resp.low/G31UF.co2$Resp.mean[1],rev(G31UF.co2$Resp.high/G31UF.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=G31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=G31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=G31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean/Resp.mean[1]~Day,data=G31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Resp.mean/Resp.mean[1]~Day,data=G31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Resp.mean/Resp.mean[1]~Day,data=G31UF.co2,lwd=1,type="b",col="black",pch=5)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(R31UA.co2$Day,rev(R31UA.co2$Day)),y=c(R31UA.co2$Resp.low/R31UA.co2$Resp.mean[1],rev(R31UA.co2$Resp.high/R31UA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UB.co2$Day,rev(R31UB.co2$Day)),y=c(R31UB.co2$Resp.low/R31UB.co2$Resp.mean[1],rev(R31UB.co2$Resp.high/R31UB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UC.co2$Day,rev(R31UC.co2$Day)),y=c(R31UC.co2$Resp.low/R31UC.co2$Resp.mean[1],rev(R31UC.co2$Resp.high/R31UC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UD.co2$Day,rev(R31UD.co2$Day)),y=c(R31UD.co2$Resp.low/R31UD.co2$Resp.mean[1],rev(R31UD.co2$Resp.high/R31UD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UE.co2$Day,rev(R31UE.co2$Day)),y=c(R31UE.co2$Resp.low/R31UE.co2$Resp.mean[1],rev(R31UE.co2$Resp.high/R31UE.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UF.co2$Day,rev(R31UF.co2$Day)),y=c(R31UF.co2$Resp.low/R31UF.co2$Resp.mean[1],rev(R31UF.co2$Resp.high/R31UF.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=R31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=R31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=R31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean/Resp.mean[1]~Day,data=R31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Resp.mean/Resp.mean[1]~Day,data=R31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Resp.mean/Resp.mean[1]~Day,data=R31UF.co2,lwd=1,type="b",col="black",pch=5)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(M31UA.co2$Day,rev(M31UA.co2$Day)),y=c(M31UA.co2$Resp.low/M31UA.co2$Resp.mean[1],rev(M31UA.co2$Resp.high/M31UA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UB.co2$Day,rev(M31UB.co2$Day)),y=c(M31UB.co2$Resp.low/M31UB.co2$Resp.mean[1],rev(M31UB.co2$Resp.high/M31UB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UC.co2$Day,rev(M31UC.co2$Day)),y=c(M31UC.co2$Resp.low/M31UC.co2$Resp.mean[1],rev(M31UC.co2$Resp.high/M31UC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UD.co2$Day,rev(M31UD.co2$Day)),y=c(M31UD.co2$Resp.low/M31UD.co2$Resp.mean[1],rev(M31UD.co2$Resp.high/M31UD.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UE.co2$Day,rev(M31UE.co2$Day)),y=c(M31UE.co2$Resp.low/M31UE.co2$Resp.mean[1],rev(M31UE.co2$Resp.high/M31UE.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UF.co2$Day,rev(M31UF.co2$Day)),y=c(M31UF.co2$Resp.low/M31UF.co2$Resp.mean[1],rev(M31UF.co2$Resp.high/M31UF.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UG.co2$Day,rev(M31UG.co2$Day)),y=c(M31UG.co2$Resp.low/M31UG.co2$Resp.mean[1],rev(M31UG.co2$Resp.high/M31UG.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UH.co2$Day,rev(M31UH.co2$Day)),y=c(M31UH.co2$Resp.low/M31UH.co2$Resp.mean[1],rev(M31UH.co2$Resp.high/M31UH.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=M31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=M31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=M31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean/Resp.mean[1]~Day,data=M31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Resp.mean/Resp.mean[1]~Day,data=M31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Resp.mean/Resp.mean[1]~Day,data=M31UF.co2,lwd=1,type="b",col="black",pch=5)
points(Resp.mean/Resp.mean[1]~Day,data=M31UG.co2,lwd=1,type="b",col="black",pch=6)
points(Resp.mean/Resp.mean[1]~Day,data=M31UH.co2,lwd=1,type="b",col="black",pch=7)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,0.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
abline(h=1,col="red")
polygon(x=c(A31UA.co2$Day,rev(A31UA.co2$Day)),y=c(A31UA.co2$Resp.low/A31UA.co2$Resp.mean[1],rev(A31UA.co2$Resp.high/A31UA.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UB.co2$Day,rev(A31UB.co2$Day)),y=c(A31UB.co2$Resp.low/A31UB.co2$Resp.mean[1],rev(A31UB.co2$Resp.high/A31UB.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UC.co2$Day,rev(A31UC.co2$Day)),y=c(A31UC.co2$Resp.low/A31UC.co2$Resp.mean[1],rev(A31UC.co2$Resp.high/A31UC.co2$Resp.mean[1])),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean/Resp.mean[1]~Day,data=A31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean/Resp.mean[1]~Day,data=A31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean/Resp.mean[1]~Day,data=A31UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('Whole-symbiosis respiration'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(normalized to 1 at day = 0)'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to low N (days)'),side=1,line=3,cex=1,outer=T)

### Figure S11
nf<-layout(matrix(seq(1,8,1),2,4,byrow=T),rep(3,8),rep(3,8),T)
layout.show(nf)
par(oma=c(5,7,4,4))
par(mar=c(0,0,0,0))

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2,las=1)
title(main=expression('  a'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G21UA.co2$Day,rev(G21UA.co2$Day)),y=c(G21UA.co2$Resp.low,rev(G21UA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UB.co2$Day,rev(G21UB.co2$Day)),y=c(G21UB.co2$Resp.low,rev(G21UB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G21UC.co2$Day,rev(G21UC.co2$Day)),y=c(G21UC.co2$Resp.low,rev(G21UC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=G21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=G21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=G21UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  b'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R21UA.co2$Day,rev(R21UA.co2$Day)),y=c(R21UA.co2$Resp.low,rev(R21UA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UB.co2$Day,rev(R21UB.co2$Day)),y=c(R21UB.co2$Resp.low,rev(R21UB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R21UC.co2$Day,rev(R21UC.co2$Day)),y=c(R21UC.co2$Resp.low,rev(R21UC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=R21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=R21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=R21UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1)


plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=F,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  c'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M21UA.co2$Day,rev(M21UA.co2$Day)),y=c(M21UA.co2$Resp.low,rev(M21UA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UB.co2$Day,rev(M21UB.co2$Day)),y=c(M21UB.co2$Resp.low,rev(M21UB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UC.co2$Day,rev(M21UC.co2$Day)),y=c(M21UC.co2$Resp.low,rev(M21UC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M21UD.co2$Day,rev(M21UD.co2$Day)),y=c(M21UD.co2$Resp.low,rev(M21UD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=M21UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=M21UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=M21UC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean~Day,data=M21UD.co2,lwd=1,type="b",col="black",pch=3)
mtext(expression(italic(Morella)),side=3,line=1,cex=1)
mtext(expression('21/15 '*degree*'C'),side=4,line=1,cex=1)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n",bty="n")

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=T,cex.axis=1.2,las=1)
title(main=expression('  d'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(G31UA.co2$Day,rev(G31UA.co2$Day)),y=c(G31UA.co2$Resp.low,rev(G31UA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UB.co2$Day,rev(G31UB.co2$Day)),y=c(G31UB.co2$Resp.low,rev(G31UB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UC.co2$Day,rev(G31UC.co2$Day)),y=c(G31UC.co2$Resp.low,rev(G31UC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UD.co2$Day,rev(G31UD.co2$Day)),y=c(G31UD.co2$Resp.low,rev(G31UD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UE.co2$Day,rev(G31UE.co2$Day)),y=c(G31UE.co2$Resp.low,rev(G31UE.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(G31UF.co2$Day,rev(G31UF.co2$Day)),y=c(G31UF.co2$Resp.low,rev(G31UF.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=G31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=G31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=G31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean~Day,data=G31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Resp.mean~Day,data=G31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Resp.mean~Day,data=G31UF.co2,lwd=1,type="b",col="black",pch=5)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  e'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(R31UA.co2$Day,rev(R31UA.co2$Day)),y=c(R31UA.co2$Resp.low,rev(R31UA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UB.co2$Day,rev(R31UB.co2$Day)),y=c(R31UB.co2$Resp.low,rev(R31UB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UC.co2$Day,rev(R31UC.co2$Day)),y=c(R31UC.co2$Resp.low,rev(R31UC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UD.co2$Day,rev(R31UD.co2$Day)),y=c(R31UD.co2$Resp.low,rev(R31UD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UE.co2$Day,rev(R31UE.co2$Day)),y=c(R31UE.co2$Resp.low,rev(R31UE.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(R31UF.co2$Day,rev(R31UF.co2$Day)),y=c(R31UF.co2$Resp.low,rev(R31UF.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=R31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=R31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=R31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean~Day,data=R31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Resp.mean~Day,data=R31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Resp.mean~Day,data=R31UF.co2,lwd=1,type="b",col="black",pch=5)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  f'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(M31UA.co2$Day,rev(M31UA.co2$Day)),y=c(M31UA.co2$Resp.low,rev(M31UA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UB.co2$Day,rev(M31UB.co2$Day)),y=c(M31UB.co2$Resp.low,rev(M31UB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UC.co2$Day,rev(M31UC.co2$Day)),y=c(M31UC.co2$Resp.low,rev(M31UC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UD.co2$Day,rev(M31UD.co2$Day)),y=c(M31UD.co2$Resp.low,rev(M31UD.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UE.co2$Day,rev(M31UE.co2$Day)),y=c(M31UE.co2$Resp.low,rev(M31UE.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UF.co2$Day,rev(M31UF.co2$Day)),y=c(M31UF.co2$Resp.low,rev(M31UF.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UG.co2$Day,rev(M31UG.co2$Day)),y=c(M31UG.co2$Resp.low,rev(M31UG.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(M31UH.co2$Day,rev(M31UH.co2$Day)),y=c(M31UH.co2$Resp.low,rev(M31UH.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=M31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=M31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=M31UC.co2,lwd=1,type="b",col="black",pch=2)
points(Resp.mean~Day,data=M31UD.co2,lwd=1,type="b",col="black",pch=3)
points(Resp.mean~Day,data=M31UE.co2,lwd=1,type="b",col="black",pch=4)
points(Resp.mean~Day,data=M31UF.co2,lwd=1,type="b",col="black",pch=5)
points(Resp.mean~Day,data=M31UG.co2,lwd=1,type="b",col="black",pch=6)
points(Resp.mean~Day,data=M31UH.co2,lwd=1,type="b",col="black",pch=7)

plot(0:200,(0:200)/200,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,100),xaxt="n",yaxt="n")
axis(1,at=c(0,40,80,120,160,200),labels=T,cex.axis=1.2)
axis(2,at=c(0,20,40,60,80,100),labels=F,cex.axis=1.2,las=1)
title(main=expression('  g'),cex.main=1.2,adj=0,line=-1)
polygon(x=c(A31UA.co2$Day,rev(A31UA.co2$Day)),y=c(A31UA.co2$Resp.low,rev(A31UA.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UB.co2$Day,rev(A31UB.co2$Day)),y=c(A31UB.co2$Resp.low,rev(A31UB.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
polygon(x=c(A31UC.co2$Day,rev(A31UC.co2$Day)),y=c(A31UC.co2$Resp.low,rev(A31UC.co2$Resp.high)),
        col=adjustcolor("black",alpha.f = 0.2),border=NA)
points(Resp.mean~Day,data=A31UA.co2,lwd=1,type="b",col="black",pch=0)
points(Resp.mean~Day,data=A31UB.co2,lwd=1,type="b",col="black",pch=1)
points(Resp.mean~Day,data=A31UC.co2,lwd=1,type="b",col="black",pch=2)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1)
mtext(expression('31/25 '*degree*'C'),side=4,line=1,cex=1)

mtext(expression('Whole-symbiosis respiration'),side=2,line=5.5,cex=1,outer=T)
mtext(expression('(nmol CO'[2]*' s'^-1*')'),side=2,line=3.5,cex=1,outer=T)
mtext(expression('Time since switching to low N (days)'),side=1,line=3,cex=1,outer=T)