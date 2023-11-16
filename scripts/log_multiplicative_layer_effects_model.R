# Log-multiplicative layer effects model
# Xie, Yu. 1992. “The Log-Multiplicative Layer Effect Model for Comparing Mobility Tables.” American Sociological Review 57(3):380-395.

# packages
library(tidyverse)
library(magrittr)
library(gnm)

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
# Null association between R and C, given L
fit_1 <- glm(Freq ~ O*C + D*C + Diag*C, data = df, family = poisson)
summary(fit_1)

# R0 Cross-nationally homogeneous row effect association
fit_2 <- glm(Freq ~ O*C + D*C + Diag*C + O:V, data = df, family = poisson)
summary(fit_2)

# R1 Cross-nationally uniform row effect association
fit_3 <- glm(Freq ~ O*C + D*C + Diag*C + O:V + U:V*C, data = df, family = poisson)
summary(fit_3)

# Rx Cross-nationally log-multiplicative row effect association
fit_4 <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(C,O,V), data = df, family = poisson)
summary(fit_4)

# C0-Cross-nationally homogeneous column effect association
fit_5 <- glm(Freq ~ O*C + D*C + Diag*C + U:D, data = df, family = poisson)
summary(fit_5)

# C1-Cross-nationally uniform column effect association
fit_6 <- glm(Freq ~ O*C + D*C + Diag*C + U:D + U:V*C, data = df, family = poisson)
summary(fit_6)

# Cx-Cross-nationally log-multiplicative column effect association
fit_7 <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(C,U,D), data = df, family = poisson)
summary(fit_7)

# (R+C) 0 - Cross-nationally homogeneous row and column effects association I
fit_8 <- gnm(Freq ~ O*C + D*C + Diag*C  + O:V + U:D, data = df, family = poisson)
summary(fit_8)

# (R+C) u - Cross-nationally uniform row and column effects association I
fit_9 <- gnm(Freq ~ O*C + D*C + Diag*C  + O:V + U:D + U:V*C, data = df, family = poisson)
summary(fit_9)

# (R+C) x - Cross-nationally log-multiplicative row and column effects association I
fit_10 <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(C,O:V+U:D), data = df, family = poisson)
summary(fit_10)

# RC0 Cross-nationally homogeneous row and column effects association
fit_11 <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(O,D), data = df, family = poisson)
summary(fit_11)

# RCx Cross-nationally log-multiplicative row and column effects association
fit_12 <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(Exp(C),O,D), data = df, family = poisson)
summary(fit_12)

# FI0 Cross-nationally homogeneous full two-way R and C interaction
fit_13 <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(1,O*D), data = df, family = poisson)
summary(fit_13)

# FIu Cross-nationally uniform full two-way R and C interaction
fit_14 <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(1,O*D) + U:V*C, data = df, family = poisson)
summary(fit_14)

# FIx Cross-nationally log-multiplicative full two-way R and C interaction
fit_15 <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(Exp(C),O*D), data = df, family = poisson)
summary(fit_15)