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

# Question 1) 

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

plot(data_3, 
     pch = 21, cex = 3, lwd = 5, 
     col = "darkblue", bg = "skyblue",
     main = "Données data_exo3",
     xlab = "x", ylab = "y")

# Question 2

x_grid <- seq(min(data_3$x), max(data_3$x), length.out = 100)
num_deg <- 12
colors  <- rainbow(num_deg)
adjr2   <- numeric(num_deg)
aics    <- numeric(num_deg)
bics    <- numeric(num_deg)

plot(data_3, 
     pch = 21, cex = 3, lwd = 5,
     col = "darkblue", bg = "skyblue",
     main = "Ajustements polynomiaux (degrés 1 à 12)",
     xlab = "x", ylab = "y")

for (i in 1:num_deg) {
  fit      <- lm(y ~ poly(x, i), data = data_3)
  s        <- summary(fit)
  adjr2[i] <- s$adj.r.squared
  aics[i]  <- AIC(fit)
  bics[i]  <- BIC(fit)
  
  y_pred <- predict(fit, newdata = data.frame(x = x_grid))
  lines(x_grid, y_pred, lwd = 2, col = colors[i])
}

legend("topright", legend = paste("Degré", 1:num_deg),
       col = colors, lwd = 2, cex = 0.7, bg = "white")

best_deg_r2  <- which.max(adjr2)
best_deg_aic <- which.min(aics)
best_deg_bic <- which.min(bics)

cat("Degré choisi par R²-ajusté :", best_deg_r2,  "\n")
cat("Degré choisi par AIC        :", best_deg_aic, "\n")
cat("Degré choisi par BIC        :", best_deg_bic, "\n")

par(mfrow = c(1, 3))

plot(1:num_deg, adjr2, type = "b", pch = 19, col = "steelblue",
     main = "R² ajusté selon le degré",
     xlab = "Degré", ylab = "R² ajusté")
abline(v = best_deg_r2, col = "red", lty = 2)

plot(1:num_deg, aics, type = "b", pch = 19, col = "darkorange",
     main = "AIC selon le degré",
     xlab = "Degré", ylab = "AIC")
abline(v = best_deg_aic, col = "red", lty = 2)

plot(1:num_deg, bics, type = "b", pch = 19, col = "darkgreen",
     main = "BIC selon le degré",
     xlab = "Degré", ylab = "BIC")
abline(v = best_deg_bic, col = "red", lty = 2)

par(mfrow = c(1, 1))

# Both AIC and R²-adjusted pick degree 12, but with only 30 observations
# that's clearly overfitting — p-values for terms above degree 6 are huge.
# BIC penalizes complexity more heavily and tends to pick a lower degree.
# Degree 5 looks like the sweet spot: captures the trend without going overboard.

chosen_deg <- 5
fit_chosen <- lm(y ~ poly(x, chosen_deg), data = data_3)
summary(fit_chosen)

plot(data_3, 
     pch = 21, cex = 3, lwd = 5,
     col = "darkblue", bg = "skyblue",
     main = "Polynôme de degré 5 sélectionné",
     xlab = "x", ylab = "y")
y_pred_chosen <- predict(fit_chosen, newdata = data.frame(x = x_grid))
lines(x_grid, y_pred_chosen, lwd = 3, col = "firebrick")

# Question 3

# The data follows an S-shaped decay toward a lower plateau, so a sigmoid fits naturally.
# We use: y = a / (1 + exp(b * (x - c))) + d
# where a = amplitude, b = steepness, c = inflection point, d = lower asymptote

fit_nls <- nls(y ~ a / (1 + exp(b * (x - cc))) + d,
               data  = data_3,
               start = list(a = 33, b = 1, cc = 3, d = 2))

summary(fit_nls)

y_nls <- predict(fit_nls, newdata = data.frame(x = x_grid))

plot(data_3,
     pch = 21, cex = 3, lwd = 5,
     col = "darkblue", bg = "skyblue",
     main = "Comparaison : Polynôme degré 5 vs NLS",
     xlab = "x", ylab = "y")
lines(x_grid, y_pred_chosen, lwd = 3, col = "firebrick", lty = 1)
lines(x_grid, y_nls,         lwd = 3, col = "purple",    lty = 2)
legend("topright",
       legend = c("Polynôme degré 5", "NLS : a / (1 + exp(b*(x-c))) + d"),
       col = c("firebrick", "purple"), lwd = 3, lty = c(1, 2), bty = "n")

# The polynomial fits better statistically, but the sigmoid only needs 4 parameters
# vs 6 and is much easier to interpret. It also extrapolates sensibly —
# converging to d as x grows — while the polynomial just diverges.

# Question 4

x_new <- data.frame(x = 1:10)

conf_int <- predict(fit_chosen, newdata = x_new, interval = "confidence", level = 0.95)
pred_int <- predict(fit_chosen, newdata = x_new, interval = "prediction", level = 0.95)

cat("\nIntervalles de confiance pour E(Y_new) :\n")
print(round(data.frame(x = 1:10, conf_int), 3))

cat("\nIntervalles de prédiction pour Y_new :\n")
print(round(data.frame(x = 1:10, pred_int), 3))

x_grid_full <- seq(0, 12, length.out = 300)

conf_full <- predict(fit_chosen, 
                     newdata = data.frame(x = x_grid_full), 
                     interval = "confidence", level = 0.95)
pred_full <- predict(fit_chosen, 
                     newdata = data.frame(x = x_grid_full), 
                     interval = "prediction", level = 0.95)

plot(data_3,
     pch = 21, cex = 3, lwd = 5,
     col = "darkblue", bg = "skyblue",
     xlim = c(0, 12), ylim = range(pred_full),
     main = "Intervalles de confiance et de prédiction (degré 5)",
     xlab = "x", ylab = "y")

polygon(c(x_grid_full, rev(x_grid_full)),
        c(pred_full[, "lwr"], rev(pred_full[, "upr"])),
        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)

polygon(c(x_grid_full, rev(x_grid_full)),
        c(conf_full[, "lwr"], rev(conf_full[, "upr"])),
        col = rgb(0.4, 0.6, 1, 0.4), border = NA)

lines(x_grid_full, conf_full[, "fit"], lwd = 3, col = "firebrick")

lines(x_grid_full, conf_full[, "lwr"], lwd = 1.5, col = "blue",  lty = 2)
lines(x_grid_full, conf_full[, "upr"], lwd = 1.5, col = "blue",  lty = 2)
lines(x_grid_full, pred_full[, "lwr"], lwd = 1.5, col = "black", lty = 3)
lines(x_grid_full, pred_full[, "upr"], lwd = 1.5, col = "black", lty = 3)

points(1:10, predict(fit_chosen, newdata = x_new),
       pch = 23, cex = 2, col = "black", bg = "yellow", lwd = 2)

legend("topright",
       legend = c("Ajustement (degré 5)",
                  "IC 95% pour E(Y_new)",
                  "IP 95% pour Y_new",
                  "x_new = 1,...,10"),
       col    = c("firebrick", "blue", "black", "black"),
       lwd    = c(3, 1.5, 1.5, NA),
       lty    = c(1, 2, 3, NA),
       pch    = c(NA, NA, NA, 23),
       pt.bg  = c(NA, NA, NA, "yellow"),
       bty    = "n", cex = 0.85)

# The CI is narrow where data is dense and widens outside that range.
# The PI is always wider since it also accounts for individual observation variability.
# For x_new in {1,...,10} both intervals stay reasonable, but beyond x = 11
# the polynomial extrapolation gets unreliable fast.

# Exercice 4

# Question 1)

library(readr)
sp <- read_csv("sp500_history.csv", show_col_types = FALSE)
head(sp)
sp["daily_return"] <- sp["Close"] - sp["Open"]
head(sp)

# Question 2)

# We fit the normal model using MLE estimators for the mean and standard deviation.

returns <- sp$daily_return
mu_hat <- mean(returns)
sigma_hat <- sd(returns)
cat("Estimated Mean (mu):", mu_hat, "\n")
cat("Estimated Std Dev (sigma):", sigma_hat, "\n")

# Let's plot the fitted normal distribution over the data.

hist(returns, probability = TRUE, main = "Normal Fit to Daily Returns", col = "lightblue")
x_seq <- seq(min(returns), max(returns), length.out = 100)
y_norm <- dnorm(x_seq, mean = mu_hat, sd = sigma_hat)
lines(x_seq, y_norm, col = "red", lwd = 2)

# We can go further with a QQ plot:

qqnorm(returns)
qqline(returns, col = "red")

# The normal distribution clearly doesn't fit well — the data has way more
# extreme values than a normal would ever produce.

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

# The GMM is a much better fit than a single Normal. The Normal misses the sharp
# peak in the center and the heavy tails, while even p = 2 or p = 3 already
# captures the shape of the data quite well.

# Question 4.a)

# The main reason we move from Normal to Student-t is that the Normal simply
# fails to account for extreme returns in financial data.
#
# The key is the 'nu' (degrees of freedom) parameter — it controls how heavy
# the tails are, making large returns much more plausible than under a Normal.
#
# It's also still a simple model with just 3 parameters:
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

# The Student-t still struggles with some outliers, but it's noticeably
# better than the Normal.

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

# The Q-Q plots confirm that the GMM fits best, which lines up with what the BIC told us.

# Question 6)

# With this much data, the CLT kicks in and we can safely use a t-test.

test_mu <- t.test(sp$daily_return, mu = 0)
print(test_mu)

print("We do not reject H0: the data do not provide significant evidence that the mean daily return differs from 0.")