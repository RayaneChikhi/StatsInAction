library(readr)
NHANES <- read_csv("NHANES_age_prediction 3.csv", show_col_types = FALSE)

data = NHANES[,c("DIQ010","age_group","RIDAGEYR","RIAGENDR","PAQ605","BMXBMI","LBXGLU","LBXGLT","LBXIN")]
colnames(data) = c("Diabete","age_group","Age","Sex","Phys_activ","BMI","Glu","Glu2h","BIL")
head(data)

# Question 1)

data_comparison <- subset(data, Diabete == 1 | Diabete == 2)
boxplot(Glu ~ Diabete, data = data_comparison)

# We notice that the variances for each diabete category are very different (indeed the
# boxes have varying lengths). Furthermore, the distribution of Glu within each category does not appear to be normal:
# for instance, in group 2, there are many outliers. In group 1, the distribution does not appear to be symmetrical.
# Thus, we use a Wilcoxon test, which does not assume normality of the data.

head(data_comparison)
wilcox.test(Glu ~ Diabete, data = data_comparison)

# We obtain a p-value of 7.9e-6 < 0.05: the test is significant. We conclude that
# there is a significant difference between these two categories.

# We also use a Welch test, since the question is about comparing means
# and this test does not assume equal variances.

t.test(Glu ~ Diabete, data = data_comparison)

# Question 2)

diabetic <- subset(data, Diabete == 1 & (age_group == "Adult" | age_group == "Senior"))
boxplot(Glu ~ age_group, data = diabetic)

# Again the data distribution for the adults is highly asymmetrical, meaning
# it would not be coherent to assume its distribution is normal. Furthermore,
# the variance between both groups is highly different. Hence, we use
# a Wilcoxon test once again.

wilcox.test(Glu ~ age_group, data = diabetic)

# We obtain a p-value of 0.08857 > 0.05. Hence, the test does not provide
# significant evidence of a difference between the two groups at the 5% level.
# We also cannot conclude that the two groups are the same.

# We also use a Welch test, since the question is about comparing means.

t.test(Glu ~ age_group, data = diabetic)

# Question 3)

adults_diabetic <- subset(data, Diabete == 1 & age_group == "Adult" & (Phys_activ == 1 | Phys_activ == 2))
boxplot(Glu ~ Phys_activ, data = adults_diabetic)
table(adults_diabetic$Phys_activ)

# We notice that there is a very small amount of diabetic adults that answered
# "yes" or "no" for their physical activity. 4 people answered yes and 9 answered no.
# Hence, the comparison is based on very small subgroup sizes, so the inference
# is not very reliable.

# If we use a Wilcoxon test, we will most likely obtain a large p-value for that reason.

wilcox.test(Glu ~ Phys_activ, data = adults_diabetic)

# Indeed, the p-value is very large: 0.9384. The test does not provide significant
# evidence of a difference between the two groups, but the sample sizes are too small
# for this comparison to be very informative.

# We also use a Welch test, since the question is about comparing means.

t.test(Glu ~ Phys_activ, data = adults_diabetic)

# Question 4)

table_sexe <- table(data$Sex, data$Diabete)
print(table_sexe)

table_binaire <- cbind(table_sexe[,1], table_sexe[,2] + table_sexe[,3])
colnames(table_binaire) <- c("Yes", "No_or_Borderline")
print(table_binaire)

chisq.test(table_binaire)

fisher.test(table_binaire)

# Exercice 2

liver_data <- load("liver_data.rda")
liver[1:5, 1:5]

# Question 1

correlations <- cor(liver[, 1], liver[, -1])
indice_max <- which.max(abs(correlations))
noms_genes <- colnames(liver)[-1]

gene_top <- noms_genes[indice_max]
value_top <- correlations[indice_max]

print(paste("The most correlated gene is :", gene_top))
print(paste("Its correlation with cholesterol is :", value_top))

plot(liver[[gene_top]], liver$cholesterol,
     main = "Cholesterol as a function of the most correlated gene",
     xlab = gene_top, ylab = "cholesterol",
     pch = 19, col = "blue")

model = lm(cholesterol ~ liver[[gene_top]], data = liver)
summary(model)
abline(model, col = "red", lwd = 2)

# In the summary, we obtain a p-value of 6.31e-12 for the slope (i.e.
# for the coefficient beta_1 associated with the most correlated gene).
# Hence we conclude that there is a significant relation between this gene
# expression and cholesterol.

# Question 2

genes <- colnames(liver)[-1]
p_values <- numeric(length(genes))
names(p_values) <- genes 

for (i in seq_along(genes)) {
    current_gene <- genes[i]
    model <- lm(cholesterol ~ liver[[current_gene]], data = liver)
    p_values[i] <- summary(model)$coefficients[2,4]
}

print(p_values)

# Question 3

p_sorted <- sort(p_values)
i <- 1:length(p_sorted)

plot(i, p_sorted, pch = 20,
     xlab = "Index", ylab = "Ordered p-values",
     main = "Ordered p-values and BH threshold")

lines(i, 0.05 * i / length(p_sorted), col = "red", lwd = 2)

# We compare the ordered p-values to the Benjamini-Hochberg threshold.
# The smallest p-values are clearly below this threshold, which suggests
# that many genes are significantly associated with cholesterol after
# controlling the false discovery rate.

# Question 4
# We use the Benjamini-Hochberg method

p_adj <- p.adjust(p_values, method = "BH")

discoveries <- p_adj[p_adj < 0.05]

length(discoveries)
# Thus, we discover 941 genes.

# Question 5

# We use Bonferonni

seuil_bonf <- 0.05 / 3116
discoveries_bonf <- p_values[p_values < seuil_bonf]

length(discoveries_bonf)

# We thus discover 220 genes. As expected, we discover fewer genes,
# since this method is more conservative: it may increase the
# number of false negatives.

# Exercice 3

# Question 1

library(readr)
data_3 <- read.table("data_exo3.csv", sep = "t", header = TRUE, dec = ",")
head(data_3)
plot(data_3, pch = 21, cex = 3, lwd = 5, col = "darkblue", bg = "skyblue")

# Question 2

# We start by creating a grid to plot our fits

x_grid <- seq(min(data_3$x), max(data_3$x), length.out = 100)
num_deg <- 12
colors <- rainbow(num_deg)
adjr2 <- numeric(num_deg)
aics <- numeric(num_deg)

for (i in 1:num_deg) {
    fit <- lm(y ~ poly(x, i), data = data_3)
    s <- summary(fit)
    adjr2[i] <- s$adj.r.squared
    aics[i] <- AIC(fit)

    y_pred = predict(fit, newdata = data.frame(x = x_grid))
    lines(x_grid, y_pred, lwd = 2, col = colors[i])
}
legend("topright", legend = paste("Degree", 1:num_deg), 
       col = colors, lwd = 2, cex = 0.7, bg = "white")

best_deg_r2_adj <- which.max(adjr2)
best_deg_AIC <- which.min(aics)

print(paste("The degree chosen by the R2-adj criterion is given by:", best_deg_r2_adj))
print(paste("The degree chosen by the AIC criterion is given by:", best_deg_AIC))

fit_r2_adj <- lm(y ~ poly(x, best_deg_r2_adj), data = data_3)
pred_r2_adj <- predict(fit_r2_adj, newdata = data.frame(x = x_grid))
lines(x_grid, pred_r2_adj, lwd = 5, col = "firebrick")

summary(fit_r2_adj)

# The adjusted R^2 and AIC criteria both suggest that degree 12 should be chosen.
# However, it seems as though the associated polynomial is overfitting our
# relatively small dataset. Graphically, we decide to keep the degree 5 polynomial,
# which seems to represent the data appropriately without overfitting.

chosen_polynomial <- lm(y ~ poly(x, 5), data = data_3)
plot(data_3, pch = 21, cex = 3, lwd = 5, col = "darkblue", bg = "skyblue")
y_pred <- predict(chosen_polynomial, newdata = data.frame(x = x_grid))
lines(x_grid, y_pred, lwd = 5, col = "firebrick")
summary(chosen_polynomial)

# Question 3

fit_nls <- nls(y ~ a / (1 + exp(-b * (x - c))), 
               data = data_3, 
               start = list(a = max(data_3$y), 
                            b = 0.5, 
                            c = mean(data_3$x)))

y_nls <- predict(fit_nls, newdata = data.frame(x = x_grid))
lines(x_grid, y_nls, col = "purple", lwd = 4)

summary(fit_nls)

# Two distinct modeling approaches were evaluated for this dataset:
# a 5th-degree polynomial regression and a non-linear sigmoidal model.
# From a purely statistical performance standpoint, the polynomial model
# gives a very good fit, with a Residual Standard Error of 3.71 and a high
# R^2 of 0.955. However, the non-linear model is also interesting,
# since it seems to capture the general shape of the data well.

# Question 4

x_new_data <- data.frame(x = 1:10)

conf_int <- predict(chosen_polynomial, newdata = x_new_data, interval = "confidence")
pred_int <- predict(chosen_polynomial, newdata = x_new_data, interval = "prediction")

print("Confidence Intervals (Average):")
print(conf_int)

print("Prediction Intervals (Individual):")
print(pred_int)

# Graph for the fitted values, confidence intervals and prediction intervals

plot(data_3$x, data_3$y, pch = 21, cex = 2, lwd = 3, col = "darkblue", bg = "skyblue",
     xlim = c(min(data_3$x), max(x_new_data$x)),
     ylim = range(c(data_3$y, conf_int, pred_int)),
     xlab = "x", ylab = "y",
     main = "Fitted values with confidence and prediction intervals")

x_grid2 <- seq(min(data_3$x), max(x_new_data$x), length.out = 200)
y_grid2 <- predict(chosen_polynomial, newdata = data.frame(x = x_grid2))
lines(x_grid2, y_grid2, lwd = 3, col = "firebrick")

lines(x_new_data$x, conf_int[, "lwr"], col = "darkgreen", lwd = 2, lty = 2)
lines(x_new_data$x, conf_int[, "upr"], col = "darkgreen", lwd = 2, lty = 2)

lines(x_new_data$x, pred_int[, "lwr"], col = "purple", lwd = 2, lty = 3)
lines(x_new_data$x, pred_int[, "upr"], col = "purple", lwd = 2, lty = 3)

points(x_new_data$x, conf_int[, "fit"], pch = 19, col = "black")

legend("topright",
       legend = c("Observed data", "Polynomial fit", "Confidence interval", "Prediction interval"),
       col = c("darkblue", "firebrick", "darkgreen", "purple"),
       pch = c(21, NA, NA, NA),
       pt.bg = c("skyblue", NA, NA, NA),
       lty = c(NA, 1, 2, 3),
       lwd = c(NA, 3, 2, 2),
       bty = "n")


# Exercice 4

# Question 1)

library(readr)
sp <- read_csv("sp500_history.csv", show_col_types = FALSE)
head(sp)
sp["daily_return"] <- sp["Close"] - sp["Open"]
head(sp)

# Question 2)

# First we fit the normal model by using the MLE estimators of the
# mean and standard deviation.

returns <- sp$daily_return
mu_hat <- mean(returns)
sigma_hat <- sd(returns)
cat("Estimated Mean (mu):", mu_hat, "\n")
cat("Estimated Std Dev (sigma):", sigma_hat, "\n")

# Let us first plot the fitted normal distribution along with the samples.

hist(returns, probability = TRUE, main = "Normal Fit to Daily Returns", col = "lightblue")
x_seq <- seq(min(returns), max(returns), length.out = 100)
y_norm <- dnorm(x_seq, mean = mu_hat, sd = sigma_hat)
lines(x_seq, y_norm, col = "red", lwd = 2)

# We can go further by using a QQ plot:

qqnorm(returns)
qqline(returns, col = "red")

# We see that the normal distribution is not a good fit. Indeed, the
# empirical distribution has much more extreme outliers than a normal distribution
# typically would.

# Question 3

dcomponents_hand <- function(theta, x) {
  p <- length(theta$pi)
  n <- length(x)
  probs <- matrix(0, nrow = n, ncol = p)
  for (j in 1:length(theta$pi)) {
    probs[, j] <- theta$pi[j] * dnorm(x, mean = theta$mu[j], sd = theta$sigma[j])
  }
  return(probs)
}

mixture_gaussian_hand <- function(x, p, max_iter = 500, threshold = 1e-8) {
  n <- length(x)
  
  # Smart initialization using quantiles
  initial_mu <- as.numeric(quantile(x, probs = seq(0.1, 0.9, length.out = p)))
  theta0 <- list(
    pi = rep(1/p, p), 
    mu = initial_mu, 
    sigma = rep(sd(x), p)
  )
  
  deviance <- numeric(max_iter)
  theta <- vector("list", max_iter)
  
  likelihoods <- dcomponents_hand(theta0, x)
  row_sums <- rowSums(likelihoods)
  deviance[1] <- -2 * sum(log(row_sums))
  theta[[1]] <- theta0
  
  for (t in 1:(max_iter - 1)) {
    # E STEP: Calculate posterior probabilities
    tau <- likelihoods / row_sums
    
    # M STEP: Update parameters
    pi_new <- colMeans(tau)
    mu_new <- colSums(tau * x) / colSums(tau)
    var_new <- colSums(tau * (matrix(x, n, p) - matrix(mu_new, n, p, byrow = TRUE))^2) / colSums(tau)
    sigma_new <- sqrt(var_new)
    
    theta[[t + 1]] <- list(pi = pi_new, mu = mu_new, sigma = sigma_new)
    
    # Convergence check
    likelihoods <- dcomponents_hand(theta[[t + 1]], x)
    row_sums <- rowSums(likelihoods)
    deviance[t + 1] <- -2 * sum(log(row_sums))
    
    if (abs(deviance[t + 1] - deviance[t]) < threshold) break
  }
  
  list(parameters = theta[[t + 1]], loglik = -deviance[t + 1] / 2, deviance = deviance[t + 1])
}

calc_mixture_density <- function(x_grid, params) {
  y_dens <- rep(0, length(x_grid))
  for (j in 1:length(params$pi)) {
    y_dens <- y_dens + params$pi[j] * dnorm(x_grid, mean = params$mu[j], sd = params$sigma[j])
  }
  return(y_dens)
}

# Generating the 5 plots (p = 2 to 6)

par(mfrow = c(2, 3))
x_grid <- seq(min(returns), max(returns), length.out = 1000)

for (p in 2:6) {
  fit <- mixture_gaussian_hand(returns, p = p)
  
  hist(returns, breaks = 50, probability = TRUE, 
       main = paste("GMM: p =", p), 
       col = "lightgray", border = "white",
       xlab = "Returns", ylab = "Density")
  
  y_values <- calc_mixture_density(x_grid, fit$parameters)
  lines(x_grid, y_values, col = "darkgreen", lwd = 2)
  
  legend("topright", legend = paste("LL:", round(fit$loglik, 2)), bty = "n", cex = 0.8)
}

par(mfrow = c(1, 1))

# Visually, the Mixture of Normals is a vastly superior fit compared to the single Normal model.
# While the single Normal (p = 1) fails to capture the high peak at the center and the extreme
# outliers in the tails, the Mixture models provide a much better description of the empirical density.
# In particular, p = 2 or p = 3 already seem to fit the data much better than the Normal model.

# Question 4.a)

# The shift from a Gaussian (Normal) framework to a Student-t distribution is
# primarily driven by the empirical failure of the Normal distribution to
# account for extreme events (outliers) in financial data.
#
# The key advantage is the 'nu' (degrees of freedom) parameter.
# Unlike the Normal distribution, the Student-t has heavier tails.
# This allows it to assign a more realistic probability to large returns.
#
# The Student-t model also remains relatively simple, with only 3 parameters:
# m (location), a (scale), and nu (degrees of freedom).

# Question 4.b)

log_lik_student <- function(theta, x) {
  nu <- theta[1]
  m <- theta[2]
  a <- theta[3]
  
  if (nu <= 0 || a <= 0) return(-Inf)
  
  n <- length(x)
  log_c <- lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(pi * nu * a^2)
  log_pdf <- log_c - ((nu + 1) / 2) * log(1 + (x - m)^2 / (nu * a^2))
  
  return(sum(log_pdf))
}

fit_t <- optim(par = c(5, mean(returns), sd(returns)), 
               fn = function(theta) -log_lik_student(theta, returns),
               method = "Nelder-Mead")

theta_mle <- fit_t$par
cat("MLE Estimates: nu =", theta_mle[1], "m =", theta_mle[2], "a =", theta_mle[3], "\n")

n <- length(returns)
logL_t <- -fit_t$value
bic_t <- -2 * logL_t + 3 * log(n)
cat("BIC for Student-t Model:", bic_t, "\n")

dt_ls <- function(x, nu, m, a) {
  term1 <- gamma((nu + 1) / 2) / (gamma(nu / 2) * sqrt(pi * nu * a^2))
  term2 <- (1 + (x - m)^2 / (nu * a^2))^(-(nu + 1) / 2)
  return(term1 * term2)
}

hist(returns, breaks = 60, probability = TRUE, 
     col = "lightgray", border = "white",
     main = "Student-t Fit vs. Empirical Distribution",
     xlab = "Daily Returns", ylab = "Density")

x_grid <- seq(min(returns), max(returns), length.out = 1000)
y_t <- dt_ls(x_grid, theta_mle[1], theta_mle[2], theta_mle[3])

lines(x_grid, y_t, col = "purple", lwd = 3)

legend("topright", 
       legend = c(paste("nu (df) =", round(theta_mle[1], 2)),
                  paste("m (loc) =", round(theta_mle[2], 2)),
                  paste("a (scale) =", round(theta_mle[3], 2))),
       col = "purple", lwd = 2, bty = "n")

y_sorted <- sort(returns)
n <- length(returns)

probs <- (1:n - 0.5) / n

theoretical_quantiles <- theta_mle[2] + theta_mle[3] * qt(probs, df = theta_mle[1])

plot(theoretical_quantiles, y_sorted, 
     main = "Q-Q Plot: Location-Scale Student-t",
     xlab = "Theoretical Student-t Quantiles",
     ylab = "Empirical Quantiles (Returns)",
     pch = 20, col = "darkorchid")

abline(0, 1, col = "red", lwd = 2)
grid()

# The Student-t fit still struggles with some outliers, although it seems better
# than the normal distribution fit.

# Question 5)

n_obs <- length(returns)

# BIC for Normal Model
logL_norm <- sum(dnorm(returns, mean = mean(returns), sd = sd(returns), log = TRUE))
bic_norm <- -2 * logL_norm + 2 * log(n_obs)

# BIC for GMM models (p = 2 to 6)
bic_gmm <- numeric(5)
loglik_gmm <- numeric(5)

for (p in 2:6) {
  fit_temp <- mixture_gaussian_hand(returns, p = p)
  bic_gmm[p - 1] <- -2 * fit_temp$loglik + (3 * p - 1) * log(n_obs)
  loglik_gmm[p - 1] <- fit_temp$loglik
}

# BIC for Student-t Model
bic_t <- -2 * logL_t + 3 * log(n_obs)

bic_comparison <- data.frame(
  Model = c("Normal (p = 1)", paste("GMM (p =", 2:6, ")"), "Student-t"),
  LogLikelihood = c(logL_norm, loglik_gmm, logL_t),
  BIC = c(bic_norm, bic_gmm, bic_t)
)

print(bic_comparison)

winner <- bic_comparison$Model[which.min(bic_comparison$BIC)]
cat("\nAccording to the BIC criterion, the best model is:", winner, "\n")

# --- Theoretical Quantiles Calculation ---

# 1. Normal (Q2)
theoretical_q_norm <- mu_hat + sigma_hat * qnorm(probs)

# 2. GMM p=3 (Q3)
# We need an additional function for GMM quantiles using interpolation
get_gmm_quantiles <- function(probs, params) {
  # Create a fine grid for interpolation
  x_grid_q <- seq(min(returns) * 1.5, max(returns) * 1.5, length.out = 10000)
  cdf_gmm <- rep(0, length(x_grid_q))
  for (j in 1:length(params$pi)) {
    cdf_gmm <- cdf_gmm + params$pi[j] * pnorm(x_grid_q, mean = params$mu[j], sd = params$sigma[j])
  }
  # Interpolate to find quantiles corresponding to probs
  return(approx(x = cdf_gmm, y = x_grid_q, xout = probs)$y)
}

fit_gmm3 <- mixture_gaussian_hand(returns, p = 3)
theoretical_q_gmm3 <- get_gmm_quantiles(probs, fit_gmm3$parameters)

# 3. Student-t (Q4b)
theoretical_q_t <- theta_mle[2] + theta_mle[3] * qt(probs, df = theta_mle[1])

# --- Side-by-Side Plotting ---

par(mfrow = c(1, 3), mar = c(5, 4, 4, 1))

# Plot 1: Normal Model
plot(theoretical_q_norm, y_sorted, 
     main = "Q-Q Plot: Normal",
     xlab = "Theoretical Normal Quantiles",
     ylab = "Empirical Quantiles",
     pch = 20, col = "red")
abline(0, 1, col = "black", lwd = 2)
grid()

# Plot 2: GMM p=3 Model
plot(theoretical_q_gmm3, y_sorted, 
     main = "Q-Q Plot: GMM p=3",
     xlab = "Theoretical GMM Quantiles",
     ylab = "Empirical Quantiles",
     pch = 20, col = "darkgreen")
abline(0, 1, col = "black", lwd = 2)
grid()

# Plot 3: Student-t Model
plot(theoretical_q_t, y_sorted, 
     main = "Q-Q Plot: Student-t",
     xlab = "Theoretical Student-t Quantiles",
     ylab = "Empirical Quantiles",
     pch = 20, col = "darkorchid")
abline(0, 1, col = "black", lwd = 2)
grid()

# Reset layout
par(mfrow = c(1, 1))

# The Q-Q plots confirm that the GMM provides the best fit among the models compared,
# which is consistent with the BIC results.

# Question 6)

# Due to the large amount of data, the central limit theorem allows us to use a t-test.

test_mu <- t.test(sp$daily_return, mu = 0)
print(test_mu)

print("We do not reject H0: the data do not provide significant evidence that the mean daily return differs from 0.")

