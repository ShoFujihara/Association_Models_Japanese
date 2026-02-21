# Detailed check of gnm Model 7 parameterization
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

cat("========== All Raw Coefficients ==========\n")
all_coefs <- coef(m7)
for (i in seq_along(all_coefs)) {
  cat(sprintf("%2d: %-45s = %10.6f\n", i, names(all_coefs)[i], all_coefs[i]))
}

cat("\n========== Counts ==========\n")
cat("Total parameters:", length(all_coefs), "\n")
cat("NA parameters:", sum(is.na(all_coefs)), "\n")

# Group by type
cat("\n========== Main effects ==========\n")
main_idx <- grep("^(\\(Intercept\\)|educ[0-9]|occ[0-9]|income[0-9])", names(all_coefs))
print(all_coefs[main_idx])

cat("\n========== Mult term 1 (occ, educ+income) ==========\n")
mult1_phi <- all_coefs[grep("Mult\\(\\., occ, educ", names(all_coefs))]
mult1_occ <- all_coefs[grep("Mult\\(1, \\., educ \\+ income\\)", names(all_coefs))]
mult1_educ <- all_coefs[grep("Mult\\(1, occ, \\. \\+ income\\)", names(all_coefs))]
mult1_inc <- all_coefs[grep("Mult\\(1, occ, educ \\+ \\.\\)", names(all_coefs))]

cat("phi1:", mult1_phi, "\n")
cat("occ scores:", mult1_occ, "\n")
cat("educ scores:", mult1_educ, "\n")
cat("income scores:", mult1_inc, "\n")

cat("\n========== Mult term 2 (educ, income) ==========\n")
mult2_phi <- all_coefs[grep("Mult\\(\\., educ, income\\)", names(all_coefs))]
mult2_educ <- all_coefs[grep("Mult\\(1, \\., income\\)", names(all_coefs))]
mult2_inc <- all_coefs[grep("Mult\\(1, educ, \\.\\)", names(all_coefs))]

cat("phi2:", mult2_phi, "\n")
cat("educ scores:", mult2_educ, "\n")
cat("income scores:", mult2_inc, "\n")

# Check if educ scores are correlated between terms
cat("\n========== Education scores comparison ==========\n")
if (length(mult1_educ) == length(mult2_educ) && length(mult1_educ) > 0) {
  cat("Term 1 educ:", mult1_educ, "\n")
  cat("Term 2 educ:", mult2_educ, "\n")
  cat("Correlation:", cor(mult1_educ, mult2_educ), "\n")
}

# Calculate expected frequency for first cell manually
cat("\n========== Manual calculation for cell [occ=1, educ=1, income=1] ==========\n")
intercept <- all_coefs["(Intercept)"]
# For reference categories (occ1, educ1, income1), main effects are 0

# Mult term 1: phi1 * occ_score * (educ_score + income_score)
# For occ1: Mult(1, ., educ + income).occ1
# For educ1: Mult(1, occ, . + income).educ1
# For income1: income1 doesn't appear in Mult1, so score = 0?

occ1_score <- all_coefs["Mult(1, ., educ + income).occ1"]
educ1_score_m1 <- all_coefs["Mult(1, occ, . + income).educ1"]

cat("intercept:", intercept, "\n")
cat("phi1:", mult1_phi, "\n")
cat("occ1 score:", occ1_score, "\n")
cat("educ1 score (term1):", educ1_score_m1, "\n")

# Check what income1 score is in term 1
cat("\nNote: income1 in term 1 is missing (reference = 0)\n")

# Calculate log mu for [1,1,1]
log_mu_111 <- intercept + mult1_phi * occ1_score * educ1_score_m1
cat("Expected (term 1 only):", exp(log_mu_111), "\n")
cat("Observed:", Freq[1], "\n")
cat("gnm fitted:", fitted(m7)[1], "\n")
