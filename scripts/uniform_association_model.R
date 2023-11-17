# Uniform association model
# Duncan, Otis Dudley. 1979. “How Destination Depends on Origin in the Occupational Mobility Table.” American Journal of Sociology 84(4):793–803.

# packages
library(tidyverse)
library(MASS)
library(vcdExtra)
library(knitr)
library(broom)


# Data (p. 795)
occ <- c("I","II","III","IV","Va","Vb","VI","VII")

tab <- matrix(
  c(50, 19, 26, 8, 7, 11, 6, 2,
    16, 40, 34, 18, 11, 20, 8, 3, 
    12, 35, 65, 66, 35, 88, 23, 21, 
    11, 20, 58, 110, 40, 183, 64, 32, 
    2, 8, 12, 23, 25, 46, 28, 12,
    12, 28, 102, 162, 90, 554, 230, 177,
    0, 6, 19, 40, 21, 158, 143, 71, 
    0, 3, 14, 32, 15, 126, 91, 106
    ),
  nrow = 8, 
  byrow = TRUE, 
  dimnames = list(occ,occ)) |> 
  as.table() 
df <- as.data.frame(tab)
names(df) <- c("O","D","Freq")

# Row and column scores
df$U <- as.numeric(df$O)
df$V <- as.numeric(df$D)

# Diagonal parameters
df <- df |> mutate(Diag = ifelse(U == V, U, 0) |> factor())

# Weight 
df <- df |> mutate(W = ifelse(Diag != 0, 0, 1))

# Data frame
df

# (1) Independence
fit_1 <- glm(Freq ~ O + D, data = df, family = poisson)
summary(fit_1)
tidy(fit_1)
glance(fit_1)

# (2) Row effects
fit_2 <- glm(Freq ~ O + D + O:V, data = df, family = poisson)
summary(fit_2)
tidy(fit_2)
glance(fit_2)

# (3) Quasi independence, diagonal omitted
fit_3 <- glm(Freq ~ O + D, data = df, weights = W, family = poisson)
summary(fit_3)
tidy(fit_3)
glance(fit_3)

fit_3b <- glm(Freq ~ O + D + Diag, data = df, family = poisson)
summary(fit_3b)
tidy(fit_3b)
glance(fit_3b)

# (4) Uniform association, diagonal omitted
fit_4 <- glm(Freq ~  O + D + U:V, data = df, weights = W, family = poisson)
summary(fit_4)
tidy(fit_4)
glance(fit_4)

# log-odds
fit_4$coefficients["U:V"] |> exp()

fit_4b <- glm(Freq ~  O + D + U:V + Diag, data = df, family = poisson)
summary(fit_4b)
tidy(fit_4b)
glance(fit_4b)

# log-odds
fit_4b$coefficients["U:V"] |> exp()

# (5) Row effects, diagonal omitted*
fit_5 <- glm(Freq ~ O + D + O:V, data = df, weights = W, family = poisson)
summary(fit_5)
tidy(fit_5)
glance(fit_5)

fit_5b <- glm(Freq ~ O + D + O:V + Diag, data = df, family = poisson)
summary(fit_5b)
tidy(fit_5b)
glance(fit_5b)

# Tabel 1
fitted_1 <- augment(fit_1, type.predict = "response", newdata = df)
fitted_1 <- fitted_1 |> dplyr::select(O, D, Freq, Fitted_1 = .fitted)
fitted_2 <- augment(fit_2, type.predict = "response", newdata = df)
fitted_2 <- fitted_2 |> dplyr::select(O, D, Freq, Fitted_2 = .fitted)
fitted_3 <- augment(fit_3, type.predict = "response", newdata = df)
fitted_3 <- fitted_3 |> dplyr::select(O, D, Freq, Fitted_3 = .fitted)
fitted_4 <- augment(fit_4, type.predict = "response", newdata = df)
fitted_4 <- fitted_4 |> dplyr::select(O, D, Freq, Fitted_4 = .fitted)
fitted_5 <- augment(fit_5, type.predict = "response", newdata = df)
fitted_5 <- fitted_5 |> dplyr::select(O, D, Freq, Fitted_5 = .fitted)

fit_df <- df |>
  left_join(fitted_1) |>
  left_join(fitted_2) |>
  left_join(fitted_3) |>
  left_join(fitted_4) |>
  left_join(fitted_5) |>
  dplyr::select(O,D,Freq,Fitted_1,Fitted_2,Fitted_3,Fitted_4,Fitted_5)

obs <- fit_df |> xtabs(Freq ~ O + D, data = _)
t_fitted_1 <- fit_df |> xtabs(Fitted_1 ~ O + D, data = _)
t_fitted_2 <- fit_df |> xtabs(Fitted_2 ~ O + D, data = _)
t_fitted_3 <- fit_df |> xtabs(Fitted_3 ~ O + D, data = _)
t_fitted_4 <- fit_df |> xtabs(Fitted_4 ~ O + D, data = _)
t_fitted_5 <- fit_df |> xtabs(Fitted_5 ~ O + D, data = _)

# Observed Counts
obs
# Fitted Counts, Model (4)
t_fitted_4 |> round(1)
# Fitted Counts, Model (5)
t_fitted_5 |> round(1)

d_1 <- (obs/t_fitted_1) |> diag()
d_2 <- (obs/t_fitted_2) |> diag()
d_3 <- (obs/t_fitted_3) |> diag()
d_4 <- (obs/t_fitted_4) |> diag()
d_5 <- (obs/t_fitted_5) |> diag()

# Diagonal Ratios
tibble(occ,d_1,d_2,d_3,d_4,d_5)

# Table 2
Model <- c("(1) Independence",
           "(2) Row effects",
           "(3) Quasi independence, diagonal omitted",
           "(4) Uniform association, diagonal omitted",
           "(5) Row effects, diagonal omitted")
df <- c(fit_1$df.residual,fit_2$df.residual,fit_3$df.residual,fit_4$df.residual,fit_5$df.residual)
X2 <- c(fit_1$deviance,fit_2$deviance,fit_3$deviance,fit_4$deviance,fit_5$deviance)
tibble(Model, df, X2) |> 
  kable(digits = 1)

# Figure 2
b_4 <- exp(fit_4$coefficients["U:V"] * (7:0) )

b_5 <- fit_5$coefficients[grep(":V", names(fit_5$coefficients))]
b_5 <- exp(b_5 * (-1)) |> dplyr::recode(.missing = 1)

b_4 <- data.frame(y = b_4, model = 4, x = 1:8)
b_5 <- data.frame(y = b_5, model = 5, x = 1:8)
d <- bind_rows(b_4,b_5)
d$x <- factor(d$x)
d$model <- factor(d$model)

d |> 
  mutate(x = factor(x, levels = 8:1)) |>
  ggplot(aes(x = x, y = y, group = model, linetype = model)) + 
  geom_line() + 
  ylim(1,3) + 
  scale_linetype_manual(values = c("longdash","solid")) +
  labs(x = "FATHER'S OCCUPATION (i)", y = "b_i/b_8") +
  theme_classic() 

# mosaic plot
mosaic(fit_1, ~O+D, residuals_type = "rstandard")
mosaic(fit_2, ~O+D, residuals_type = "rstandard")
mosaic(fit_3b, ~O+D, residuals_type = "rstandard")
mosaic(fit_4b, ~O+D, residuals_type = "rstandard")
mosaic(fit_5b, ~O+D, residuals_type = "rstandard")

