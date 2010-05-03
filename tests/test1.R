library(mrdrc)

## Dataset: offspring

offspring <- read.table("offspring.txt", header = TRUE)

offspring.m1 <- pnsdrm(offspring$conc, offspring$number, type = "continuous",
model = "semi-parametric", fct = list(L.3(), LL.3()))

offspring.m1

offspring.m1$xtables

plot(offspring.m1)

offspring.m2 <- pnsdrm(offspring$conc, offspring$number, type = "continuous",
model = "semi-parametric", fct = list(L.3(), LL.3()), reference=0, respLev=c(5,25,35,85,95))

offspring.m2

offspring.m2$xtables

plot(offspring.m2)


## TerMet

TerMet <- read.table("TerMet.txt", header = TRUE)

TerMet.m1 <- pnsdrm(TerMet$dose, TerMet$rgr, type = "continuous",
model = "semi-parametric", fct = list(L.4(), LL.4(), W1.4(), W2.4()))

TerMet.m1

TerMet.m1$xtables

plot(TerMet.m1)

# For comparison
TerMet.m2 <- drm(rgr~dose, data=TerMet, fct=LL.4())
drc:::ED(TerMet.m2, c(10,20,50))


## termec2

termec2 <- read.table("termec2.txt", header = TRUE)

termec2.m1 <- pnsdrm(termec2$dose1, termec2$rgr1, type = "continuous",
model = "semi-parametric", fct = list(LL.3(), W1.3(), W2.3()), respLev = c(10, 20, 50, 80, 90))

termec2.m1

termec2.m1$xtables

plot(termec2.m1)

# For comparison
termec2.m2<-drm(rgr1~dose1, data=termec2, fct=LL.3())
drc:::ED(termec2.m2, c(10,20,50,80,90))


## secalonic (no replicates)

secalonic.m1 <- pnsdrm(secalonic$dose, secalonic$rootl, type = "continuous",
model = "semi-parametric", fct = list(LL.3(), W1.3(), W2.3()), 
respLev = c(1, 3, 5, 7, 10, 15, 20, 50, 80, 90, 95, 97, 99))

secalonic.m1

secalonic.m1$xtables

plot(secalonic.m1)

secalonic.m2 <- pnsdrm(secalonic$dose, secalonic$rootl, type = "continuous",
model = "semi-parametric", fct = list(LL.3(), W1.3(), W2.3()), 
respLev = c(1, 3, 5, 7, 10, 15, 20, 50, 80, 90, 95, 97, 99), reference = 0)

secalonic.m2$xtables


# For comparison
secalonic.m3<-drm(rootl~dose, data=secalonic, fct=LL.3())
drc:::ED(secalonic.m3, c(1, 3, 5, 7, 10, 15, 20, 50, 80, 90, 95, 97, 99))


## S.alba2

S.alba2 <- read.table("S.alba2.txt", header = TRUE)

S.alba2.m1 <- pnsdrm(S.alba2$Dose, S.alba2$DryMatter, type = "continuous",
model = "semi-parametric", fct = list(LL.4(), LL.5(), W1.4(), W2.4()), respLev = c(10, 20, 50, 80, 90))

S.alba2.m1

S.alba2.m1$xtables

plot(S.alba2.m1)

# For comparison
S.alba2.m2<-drm(DryMatter~Dose, data=S.alba2, fct=LL.4())
drc:::ED(S.alba2.m2, c(10,20,50,80,90))

 
## heartrate

heartrate.m1 <- pnsdrm(heartrate$pressure, heartrate$rate, type = "continuous",
model = "semi-parametric", fct = list(LL.4(), LL.5(), W1.4(), W2.4()), respLev = c(10, 20, 50, 80, 90))

heartrate.m1

heartrate.m1$xtables

plot(heartrate.m1)

# For comparison
heartrate.m2<-drm(rate~pressure, data=heartrate, fct=LL.4())
drc:::ED(heartrate.m2, c(10,20,50,80,90))


## finney71

finney71.m1 <- pnsdrm(finney71$dose, finney71$affected, finney71$total, type = "binomial",
model = "semi-parametric", fct = list(LL.2(), W1.2(), W2.2()), respLev = c(10, 20, 50, 80, 90))

finney71.m1

finney71.m1$xtables

plot(finney71.m1)

# For comparison
finney71.m2<-drm(affected/total~dose, weights = total, data=finney71, fct=LL.2())
drc:::ED(finney71.m2, c(10,20,50,80,90))
 

## earthworms ... not working in the GUI (the argument logex = TRUE is needed)

earthworms.m1 <- pnsdrm(earthworms$dose, earthworms$number, earthworms$total, type = "binomial",
model = "semi-parametric", fct = list(LL.2(), W2.2()), respLev = c(5, 7, 10, 20, 50, 80, 90), logex = TRUE)

earthworms.m1

earthworms.m1$xtables

plot(earthworms.m1)

# For comparison
earthworms.m2<-drm(number/total~dose, weights = total, data=earthworms, fct=LL.2(), type = "binomial")
drc:::ED(earthworms.m2, c(5,7,10,20,50,80,90))
 

## drymatter12

drymatter12 <- read.table("drymatter12.txt", header = TRUE)

drymatter12.m1 <- pnsdrm(drymatter12$DOSE, drymatter12$Drymatter, type = "continuous",
model = "semi-parametric", fct = list(LL.4(), LL.5(), W1.4(), W2.4()), respLev = c(10, 20, 50, 80, 90))

drymatter12.m1

drymatter12.m1$xtables

plot(drymatter12.m1)

## not working in the GUI (the argument logex = TRUE is needed)
drymatter12.m2 <- pnsdrm(drymatter12$DOSE, drymatter12$Drymatter, type = "continuous",
model = "semi-parametric", fct = list(LL.4(), LL.5(), W1.4(), W2.4()), respLev = c(10, 20, 50, 80, 90),
logex = TRUE)

drymatter12.m2$xtables

# For comparison
drymatter12.m3<-drm(Drymatter~DOSE, data=drymatter12, fct=LL.4())
drc:::ED(drymatter12.m3, c(10,20,50,80,90))


## daphnids2

daphnids2 <- read.table("daphnids2.txt", header = TRUE)

daphnids2.m1 <- pnsdrm(daphnids2$dose, daphnids2$no, daphnids2$total, type = "binomial",
model = "semi-parametric", fct = list(LL.2(), W2.2()), respLev = c(10, 20, 50, 80, 90))

daphnids2.m1

daphnids2.m1$xtables

plot(daphnids2.m1)

# For comparison
daphnids2.m2<-drm(no/total~dose, weights = total, data=daphnids2, fct=LL.2())
drc:::ED(daphnids2.m2, c(10,20,50,80,90))
 

## BinPro0.001EsfDaf48_2

binpro <- read.csv("BinPro0.001EsfDaf48_2.csv")

binpro.m1 <- pnsdrm(binpro$dose, binpro$no, binpro$total, type = "binomial",
model = "semi-parametric", fct = list(LL.2(), W2.2()), respLev = c(1, 5, 10, 20, 50, 80, 90),
reference = 0)

binpro.m1

binpro.m1$xtables

plot(binpro.m1)

# For comparison
binpro.m2<-drm(no/total~dose, weights = total, data=binpro, fct=LL.2())
drc:::ED(binpro.m2, c(1, 5, 10,20,50,80,90))
 

## barley

barley <- read.csv("barley.csv")

barley.m1 <- pnsdrm(barley$Dose, barley$DW, type = "continuous",
model = "semi-parametric", fct = list(LL.3(), W1.3(), W2.3()), respLev = c(10, 20, 25, 50, 75, 80, 90))
## used to be  (L.3() not working in general?)
#model = "semi-parametric", fct = list(L.3(), LL.3(), W1.3(), W2.3()), respLev = c(10, 20, 25, 50, 75, 80, 90))

barley.m1

barley.m1$xtables

plot(barley.m1)

# For comparison
barley.m2<-drm(DW~Dose, data=barley, fct=LL.4())
drc:::ED(barley.m2, c(10,20,25,50,75,80,90))
