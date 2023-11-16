# uniform layer effect model 
# Yamaguchi, Kazuo. 1987. “Models for Comparing Mobility Tables: Toward Parsimony and Substance.” American Sociological Review 52(4):482-494.

# packages
library(tidyverse)
library(MASS)
library(magrittr)
library(knitr)
library(vcdExtra)
library(car)

# Data
tab <- matrix(c(1275,364,274,272,17,
                1055,597,394,443,31,
                1043,587,1045,951,47,
                1159,791,1323,2046,52,
                666,496,1031,1632,646,
                
                474,129,87,124,11,
                300,218,171,220,8,
                438,254,669,703,16,
                601,388,932,1789,37,
                76,56,125,295,191,
                
                127,101,24,30,12,
                86,207,64,61,13,
                43,73,122,60,13,
                35,51,62,66,11,
                109,206,184,253,325), 
              ncol = 5, 
              byrow = TRUE, 
              dimnames = list(1:15,1:5)) %>% 
  as.table() 
df <- as.data.frame(tab)
names(df) <- c("O","D","Freq")
df$O <- 1:5
df$C <- rep(c(1,2,3), each = 5)
df$O <- as.factor(df$O)
df$D <- as.factor(df$D)
df$C <- as.factor(df$C)

# Diagonal parameters
df %<>% mutate(Diag = ifelse(O == D, O, 0) %>% factor())

# Row and column scores
df$U <- as.numeric(df$O)
df$V <- as.numeric(df$D)

# Off diagonal cells (Diag*C)
# No association between R and C, given L
fit_1 <- glm(Freq ~ O*C + D*C + Diag*C, data = df, family = poisson)
summary(fit_1)

# U0-Cross-nationally homogeneous uniform association
fit_2 <- glm(Freq ~ O*C + D*C + Diag*C + U:V, data = df, family = poisson)
summary(fit_2)
b_2 <- fit_2$coefficients[grep("U:V", names(fit_2$coefficients))]
b_2

# U1-Cross-nationally uniform uniform association
fit_3 <- glm(Freq ~ O*C + D*C + Diag*C + U:V*C, data = df, family = poisson)
summary(fit_3)
b_3 <- fit_3$coefficients[grep("U:V", names(fit_3$coefficients))]
b_3

# U2-Cross-nationally uniform uniform association
df$C12 <- car::recode(df$C,"1=1;2=1;3=2",as.factor = TRUE)
fit_4 <- glm(Freq ~ O*C + D*C + Diag*C + U:V*C12, data = df, family = poisson)
summary(fit_4)
b_4 <- fit_4$coefficients[grep("U:V", names(fit_4$coefficients))]
b_4
      
# R0-Cross-nationally homogeneous row effect association
fit_5 <- glm(Freq ~ O*C + D*C + Diag*C + O:V, data = df, family = poisson)
summary(fit_5)

# R1-Cross-nationally uniform row effect association
fit_6 <- glm(Freq ~ O*C + D*C + Diag*C + O:V + U:V*C, data = df, family = poisson)
summary(fit_6)
b_6 <- fit_6$coefficients[grep("V:U", names(fit_6$coefficients))]
b_6

# R2-Cross-nationally uniform row effect association
fit_7 <- glm(Freq ~ O*C + D*C + Diag*C + O:V + U:V*C12, data = df, family = poisson)
summary(fit_7)
b_7 <- fit_7$coefficients[grep("V:U", names(fit_7$coefficients))]
b_7

# C0-Cross-nationallyhomogeneous columneffect association
fit_8 <- glm(Freq ~ O*C + D*C + Diag*C + U:D, data = df, family = poisson)
summary(fit_8)

# C1-Cross-nationally uniform columneffect association
fit_9 <- glm(Freq ~ O*C + D*C + Diag*C + U:D + U:V*C, data = df, family = poisson)
summary(fit_9)
b_9 <- fit_9$coefficients[grep("U:V", names(fit_9$coefficients))]
b_9

# C2-Cross-nationally uniform columneffect association
fit_10 <- glm(Freq ~ O*C + D*C + Diag*C + U:D + U:V*C12, data = df, family = poisson)
summary(fit_10)
b_10 <- fit_10$coefficients[grep("U:V", names(fit_10$coefficients))]
b_10

# H0-Cross-nationally homogeneous homogeneous row & column effect association
#fit_11 <- gnm(Freq ~ O*C + D*C + Diag*C + MultHomog(O,D), data = df, family = poisson)
#summary(fit_11)

# H1-Cross-nationally uniform homogeneousrow & column effect association
#fit_12 <- gnm(Freq ~ O*C + D*C + Diag*C + MultHomog(O,D) + U:V*C, data = df, family = poisson)
#summary(fit_12)

# H2-Cross-nationally uniform homogeneousrow & column effect association
#fit_13 <- gnm(Freq ~ O*C + D*C + Diag*C + MultHomog(O,D) + U:V:C12, data = df, family = poisson)
#summary(fit_13)

# R+C0 -Cross-nationally uniform row & column effect association
fit_14 <- glm(Freq ~ O*C + D*C + Diag*C + U:D + O:V, data = df, family = poisson)
summary(fit_14)

# R+C1 -Cross-nationally uniform row & column effect association
fit_15 <- glm(Freq ~ O*C + D*C + Diag*C + U:D + O:V + U:V*C, data = df, family = poisson)
summary(fit_15)
b_15 <- fit_15$coefficients[grep("U:V", names(fit_15$coefficients))]
b_15

# R+C2-Cross-nationally uniform row & column effect association
fit_16 <- glm(Freq ~ O*C + D*C + Diag*C + U:D + O:V + U:V*C12, data = df, family = poisson)
summary(fit_16)
b_16 <- fit_16$coefficients[grep("U:V", names(fit_16$coefficients))]
b_16

df$SYM <- c(rep(c(0,1,2,3,4),3),
            rep(c(1,0,5,6,7),3),
            rep(c(2,5,0,8,9),3),
            rep(c(3,6,8,0,10),3),
            rep(c(4,7,9,10,0),3)) %>% 
  factor()

# QS0-Cross-nationally homogeneous quasi-symmetry
fit_17 <- glm(Freq ~ O*C + D*C + Diag*C + SYM, data = df, family = poisson)
summary(fit_17)


# QS1-Cross-nationally uniform quasi-symmetr
fit_18 <- glm(Freq ~ O*C + D*C + Diag*C + SYM + U:V*C, data = df, family = poisson)
summary(fit_18)
b_18 <- fit_18$coefficients[grep("U:V", names(fit_18$coefficients))]
b_18

# QS2-Cross-nationally uniform quasi-symmetr
fit_19 <- glm(Freq ~ O*C + D*C + Diag*C + SYM + U:V*C12, data = df, family = poisson)
summary(fit_19)
b_19 <- fit_19$coefficients[grep("U:V", names(fit_19$coefficients))]
b_19

# FI0-Cross-nationallyhomogeneous full two-way interaction
fit_20 <- glm(Freq ~ O*C + D*C + Diag*C + O:D, data = df, family = poisson)
summary(fit_20)

# FI1-Cross-nationally uniform full two-way interaction
fit_21 <- glm(Freq ~ O*C + D*C + Diag*C + O:D + U:V*C, data = df, family = poisson)
summary(fit_21)
b_21 <- fit_21$coefficients[grep("U:V", names(fit_21$coefficients))]
b_21

# FI2-Cross-nationally uniform full two-way interaction
fit_22 <- glm(Freq ~ O*C + D*C + Diag*C + O:D + U:V*C12, data = df, family = poisson)
summary(fit_22)
b_22 <- fit_22$coefficients[grep("U:V", names(fit_22$coefficients))]
b_22
