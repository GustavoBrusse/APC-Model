################################################################################
############  GUSTAVO BRUSSE- EDSD 2018/2019 ###################################
############  APC model - Assignment         ###################################
################################################################################

install.packages("HMDHFDplus")
library(HMDHFDplus) 
library(Epi)

## Question 1
## Reading the data from HMD
irl.D <- readHMDweb(CNTRY = "IRL", item="Deaths_1x1", username = "gustavo.brusse@gmail.com", password = "1536426637")
irl.Y <- readHMDweb( CNTRY = "IRL", item="Exposures_1x1", username = "gustavo.brusse@gmail.com", password = "1536426637")

head(irl.D)
head(irl.Y)

# choosing the period and age
t1=1990
t2=2010
a1=30
a2=65

event <- irl.D[which(irl.D$Year>=t1 & irl.D$Year<=t2 & irl.D$Age<=a2 & irl.D$Age>=a1 ),]
exposed <- irl.Y[which(irl.Y$Year>=t1 & irl.Y$Year<=t2 & irl.Y$Age<=a2 & irl.Y$Age>=a1 ),]

## Question 2
## age-period-cohort analysis of male and female mortality rates
# Males
male <- data.frame(A=event$Age, P=event$Year, D=event$Male, Y=exposed$Male)
head(male)

# Males model
m.cp <- apc.fit( male, model = "factor", parm = "ACP", ref.c=2000, scale=10^5)
par(mai=c(1,1,1,1))

# Females
female <- data.frame(A=event$Age, P=event$Year, D=event$Female, Y=exposed$Female)
head(female)

# Females model
f.cp <- apc.fit(female, model = "factor", parm = "ACP", ref.c=2000, scale=10^5)

## show the results as curves in the same display
plot(m.cp,lwd=2)
title(main="APC Mortality Rates, Ireland (1990-2010): Males vs. Females ")
legend("topright", legend=c("Males", "Females"), lwd = 2, col=c(1,2)) 
lines(f.cp,col=2,lwd=2)

# Question 4: 
# Graph the M/F mortality rate-ratio in an apc.frame. 
# You may want to consult the function ci.ratio from the Epi package.
# Comparison of rates
# Males vs. Females
Age.comp <- cbind(Age=m.cp$Age[,1],Male=m.cp$Age[,2],Female=f.cp$Age[,2])
Age.comp  
Per.comp <- cbind(Period=m.cp$Per[,1],Male=m.cp$Per[,2],Female=f.cp$Per[,2]);Per.comp # Period
Coh.comp <- cbind(Cohort=m.cp$Coh[,1],Male=m.cp$Coh[,2],Female=f.cp$Coh[,2]);Coh.comp # Cohort

# Plot ratio between males and females according to age and calendar time
par(mfrow=c(1,1))
par(mai=c(1,1,1,1))
plot(f.cp,"Male vs. Female RR", col="transparent")
title(main="APC Mortality Rates: Ratio between Males and Females ")
##### 4. APC.FRAME PLOT #####
par(mfrow=c(1,1))
par(mai=c(1,1,1,1))
par( mar=c(4,4,1,4), mgp=c(3,1,0)/1.6, las=1 )
apc.frame( a.lab = seq(0,30,5),
           cp.lab = seq(1955,2010,10),
           r.lab = c(c(0.5, 0.75, 1.5),c(0.1,0.2)*10),
           rr.ref = 1,
           a.tic = seq(0,30,5),
           cp.tic = seq(1955,2010,10),
           r.tic = c(c(0.5, 0.75, 1.5),c(0.1,0.2)*10),
           a.txt = "Age",
           cp.txt = "Calendar time",
           r.txt = "Male vs. Female Rate Ratio",
           rr.txt = " ",
           ref.line=TRUE,
           gap = 8,
           col.grid = gray(0.85),
           sides = c(1,2,4))
matshade(m.cp$Age[,1], ci.ratio(m.cp$Age[,-1],f.cp$Age[,-1]), col=4)
pc.matshade(m.cp$Per[,1], ci.ratio(m.cp$Per[,-1],f.cp$Per[,-1]), col=3)
pc.matshade(m.cp$Coh[,1], ci.ratio(m.cp$Coh[,-1],f.cp$Coh[,-1]), col=6)
abline(h=1)





##### 2. APC ANALYSIS #####

### FEMALES ### 

# Age-, age-period-, age-cohort-, age-drift- and age-period-cohort-model
A.F <- glm( D ~ factor(A), offset = log(Y),
            family = poisson, data = male )
AP.F <- glm( D ~ factor(A) + factor(P) + offset( log(Y) ),
             family=poisson, data=male )
AC.F <- glm( D ~ factor(A) + factor(P-A) + offset( log(Y) ),
             family=poisson, data=male )
Ad.F <- glm( D ~ factor(A) + P + offset( log(Y) ),
             family=poisson, data=male )
APC.F <- glm( D ~ factor(A) + factor(P) + factor(P-A), offset = log(Y),
              family = poisson, data = male )


# test all the models in a sequence
anova(A.F, Ad.F, AP.F, APC.F, AC.F, Ad.F, test = "Chisq")



