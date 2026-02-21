# Check gnm Model 7 parameter structure
library(tidyverse)
library(gnm)

# Data
Freq <- c(1096,  1847,  1255,  925,    3321,  6123,  6830,  5524,
          1541,  3134,  3145, 3300,    1915,  4429,  7035,  9421,
          4183,  5139,  1857, 1272,    8080,  8586,  4788,  4294,
          6033,  9211,  5046, 1058,   28130, 44589, 20074,  3408,
          4354, 13430, 18670, 9821,    2250,  9075, 18286, 14358,
         14587, 31470, 16390, 3751,    8242, 17532, 12825,  3956,
          1517,  5820,  6197, 2372,     721,  2909,  4141,  2070,
          3581,  9268,  5463, 1007,    1341,  3808,  3163,   815,
          1454,  3109,  1055,  888,     563,  1909,  1018,  1051,
          3237,  3851,   377,  102,     731,   858,   247,    84,
         14882, 22182,  5363, 1136,   11650, 15818,  5524,  2122,
          6033,  3475,    63,   18,    1603,  1005,    30,    16,

          5968,  8783,  7701, 6483,    8733, 14329, 19386, 28143,
          1011,  2162,  3536, 6649,     697,  1299,  2362, 10796,
          3214,  3621,  2485, 3177,     793,  1134,  1292,  3597,
         11532, 16837,  6975, 1839,    2563,  2995,  2060,  1600,
          1009,  2719,  3521, 3409,     296,   503,   626,  1273,
          1586,  3025,  1726,  668,     245,   415,   238,   218,
           387,   941,   564,  316,      86,   138,    79,    48,
           994,  1988,   542,  145,     158,   259,   101,    56,
           171,   409,   223,  245,      65,   172,    99,   174,
           293,   290,    67,   31,      32,    62,    18,    30,
          4288,  4916,  1452,  766,     616,   794,   347,   300,
           370,   186,     3,    4,      67,    37,     5,     2)

income <- gl(4, 1, 192)
occ <- gl(12, 8, 192)
edu1 <- rep(c(1,2), each = 4, length.out = 96)
edu2 <- rep(c(3,4), each = 4, length.out = 96)
educ <- c(edu1, edu2) |> factor()

tab5.4 <- data.frame(Freq = Freq, occ = occ, educ = educ, income = income)

# Model 7
set.seed(123)
m7 <- gnm(Freq ~ educ + occ + income +
            Mult(1, occ, educ + income) +
            Mult(1, educ, income),
          data = tab5.4,
          family = poisson,
          tolerance = 1e-12)

cat("========== gnm Model 7 ==========\n")
cat("L²:", deviance(m7), "df:", df.residual(m7), "\n\n")

# All coefficients
cat("All coefficients:\n")
coef_names <- names(coef(m7))
print(coef_names)

cat("\n\nNumber of parameters:", length(coef(m7)), "\n")

# Mult term coefficients
cat("\n========== Mult term structure ==========\n")
mult1_coefs <- coef(m7)[grep("Mult.*inst = 1", names(coef(m7)))]
cat("\nMult(1, occ, educ + income) coefficients:\n")
print(mult1_coefs)

mult2_coefs <- coef(m7)[grep("Mult.*inst = 2", names(coef(m7)))]
cat("\nMult(1, educ, income) coefficients:\n")
print(mult2_coefs)

# Get scores using getContrasts
cat("\n========== Occupation scores ==========\n")
mu_occ <- getContrasts(model = m7,
                       set = pickCoef(m7, "[.]occ"),
                       ref = "mean",
                       scaleRef = "mean",
                       scaleWeights = "unit")
print(mu_occ$qvframe)

# Education scores from Mult term 1
cat("\n========== Education scores (from Mult 1) ==========\n")
edu_coefs_mult1 <- coef(m7)[grep("Mult.*inst = 1.*educ", names(coef(m7)))]
print(edu_coefs_mult1)

# Education scores from Mult term 2
cat("\n========== Education scores (from Mult 2) ==========\n")
edu_coefs_mult2 <- coef(m7)[grep("Mult.*inst = 2.*educ", names(coef(m7)))]
print(edu_coefs_mult2)

# Income scores from Mult term 1
cat("\n========== Income scores (from Mult 1) ==========\n")
inc_coefs_mult1 <- coef(m7)[grep("Mult.*inst = 1.*income", names(coef(m7)))]
print(inc_coefs_mult1)

# Income scores from Mult term 2
cat("\n========== Income scores (from Mult 2) ==========\n")
inc_coefs_mult2 <- coef(m7)[grep("Mult.*inst = 2.*income", names(coef(m7)))]
print(inc_coefs_mult2)

# Check if education scores are the same between Mult terms
cat("\n========== Score comparison ==========\n")
if (length(edu_coefs_mult1) > 0 && length(edu_coefs_mult2) > 0) {
  cat("Education scores correlation between Mult terms:",
      cor(edu_coefs_mult1, edu_coefs_mult2), "\n")
}
if (length(inc_coefs_mult1) > 0 && length(inc_coefs_mult2) > 0) {
  cat("Income scores correlation between Mult terms:",
      cor(inc_coefs_mult1, inc_coefs_mult2), "\n")
}
