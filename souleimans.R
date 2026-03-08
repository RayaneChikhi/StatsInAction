library(readr)
NHANES <- read_csv("NHANES_age_prediction 3.csv")

data = NHANES[,c("DIQ010","age_group","RIDAGEYR","RIAGENDR","PAQ605","BMXBMI","LBXGLU","LBXGLT","LBXIN")]
colnames(data) = c("Diabete","age_group","Age","Sex","Phys_activ","BMI","Glu","Glu2h","BIL")
head(data)

# Question 1)

boxplot(Glu ~ Diabete, data = data)

# We notice that the variances for each diabete category is very different (indeed the
# boxes have varying lengths). Furthermore, the distribution of Glu within each category does not appear to be normal:
# for instance, in group 2, there are many outliers. In group 1, the distribution does not appear to be symmetrical.
# Thus, instead of a t-test, we use a Wilcoxon test which does not assume any distribution on the data

data_comparison <- subset(data, Diabete==1 | Diabete==2)
head(data_comparison)
wilcox.test(Glu ~ Diabete, data = data_comparison)

# We obtain a p-value of 7.9e-6 < 0.05 : the test is significant. We conclude that
# the mean values for these two categories are not the same.

# Question 2)

diabetic <- subset(data, Diabete==1)
boxplot(Glu ~ age_group, data=diabetic)

# Again the data distribution for the adults is highly asymmetrical meaning
# it would not be coherent to assume its distribution is normal. Furthermore
# the variance between both groups is highly different. Hence, instead of a t-test
# we aim for a Wilcoxon test once again.

wilcox.test(Glu ~ age_group, data=diabetic)

# We obtain a p-value of 0.08857 > 0.05. Hence, we cannot conclude that the means are
# different. (Neither can we conclude that they are the same: we do not know at this point).

# Question 3)

adults_diabetic <- subset(data, Diabete==1 & age_group == "Adult" & (Phys_activ == 1 | Phys_activ == 2))
boxplot(Glu ~ Phys_activ, data=adults_diabetic)
table(adults_diabetic$Phys_activ)

# We notice that there is a very little amount of diabetic adults that answered
# "yes" or "no" for their physical activity. 4 people answered yes and 9 answered no.
# Hence it is already quite obvious that we cannot test for the adults diabetic respondents
# if the mean level of blood glucose after fasting is the same for those who have a vigorous work activity and for those who have not 

# If we use a Wilcoxon test, we will most likely obtain a large p-value for that reason

wilcox.test(Glu ~ Phys_activ, data=adults_diabetic)

# Indeed the p-value is extremely large: 0.9384. We cannot test for this.

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
indice_max <- which.max(correlations)
noms_genes <- colnames(liver)[-1]

gene_top <- noms_genes[indice_max]
value_top <- correlations[indice_max]

print(paste("The most correlated gene is :", gene_top))
print(paste("Its correlation with cholesterol is :", value_top))

plot(cholesterol ~ A_42_P631473, data=liver, main="Cholesterol as a function of the most correlated gene", pch=19,col="blue")

model = lm(cholesterol ~ A_42_P631473, data=liver)
summary(model)
abline(lm(cholesterol ~ A_42_P631473, data=liver), col="red", lwd=2)

# In the summary, we obtain a p-value of 6.31e-12 for the slope (i.e
# for the coefficient beta_1 associated with the most correlated gene)
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

p_values <- sort(p_values)
plot(p_values)
curve(x/3116, from=1, to=3116, col="red", lwd=2, 
      main="Graph of f(x) = x/3116",
      xlab="Index", ylab="Theoretical threshold")
abline(a=0, b=1/3116, col="red", lwd=2)

# The p_values curve is clearly below x->x/3116 meaning the p-values
# are significant. (Il y a un concept sous-jacent du cours mais je l'ai pas en tête là)


# Question 4
# We use the Benjamini-Hochberg method

p_adj <- p.adjust(p_values, method = "BH")

discoveries <- p_adj[p_adj < 0.05]

length(discoveries)
# Thus we discover 941 genes

# Question 5

# We use Bonferonni

seuil_bonf <- 0.05 / 3116
discoveries_bonf <- p_values[p_values < seuil_bonf]

length(discoveries_bonf)

# We thus discover 220 genes. As expected, we discover less genes
# since this method is more conservative: it might increase the
# number of false negatives.

# Exercice 3

# Question 1

library(readr)
data_3 <- read.table("data_exo3.csv", sep = "t", header = TRUE, dec = ",")
head(data_3)
plot(data_3, pch=21,cex=3, lwd=5, col="darkblue", bg="skyblue")


# Question 2

# We start by creating a grid to plot our fits

x_grid <- seq(min(data_3$x), max(data_3$x), length.out = 100)
num_deg <- 12
colors <- rainbow(num_deg)
adjr2 <- numeric(num_deg)
aics <- numeric(num_deg)

for (i in 1:num_deg) {
    fit <- lm(y~poly(x,i), data=data_3)
    s <- summary(fit)
    adjr2[i] <- s$adj.r.squared
    aics[i] <- AIC(fit)

    y_pred = predict(fit, newdata = data.frame(x=x_grid))
    lines(x_grid,y_pred,lwd=2, col = colors[i])
}
legend("topright", legend=paste("Degree", 1:num_deg), 
       col=colors, lwd=2, cex=0.7, bg="white")

best_deg_r2_adj<- which.max(adjr2)
best_deg_AIC <- which.min(aics)


print(paste("The degree chosen by the R2-adj criterion is given by:", best_deg_r2_adj))
print(paste("The degree chosen by the AIC criterion is given by:", best_deg_AIC))

fit_r2_adj <- lm(y~poly(x, best_deg_r2_adj), data=data_3)
pred_r2_adj <- predict(fit_r2_adj, newdata=data.frame(x=x_grid))
lines(x_grid, pred_r2_adj, lwd=5, col="firebrick")

summary(fit_r2_adj)

# The adjusted R^2 and AIC criterions both claim that degree 12 is the one
# that must be chosen. However, it seems as though the associated
# polynomial is overfitting on our (low amount of) data. This is
# furthermore confirmed by the study of the p-values that are
# drastically high for high degree monomials (> degree 6). 
# Graphically, we decide to conserve the degree 5 polynomial
# which seems to appropriately represent the data without overfitting.

chosen_polynomial <- lm(y~poly(x, 5), data=data_3)
plot(data_3, pch=21,cex=3, lwd=5, col="darkblue", bg="skyblue")
y_pred <- predict(chosen_polynomial, newdata=data.frame(x=x_grid))
lines(x_grid, y_pred, lwd=5, col="firebrick")
summary(chosen_polynomial)

# Question 3

fit_nls <- nls(y ~ a*(1-b*exp(-c*x)), 
               data=data_3, 
               start = list(a = 1, b = 15, c = 0.5))

y_nls <- predict(fit_nls, newdata = data.frame(x = x_grid))
lines(x_grid, y_nls, col = "purple", lwd = 4)

summary(fit_nls)

#Two distinct modeling approaches were evaluated for this dataset: 
#a 5th-degree polynomial regression and a non-linear exponential 
#model.From a purely statistical performance standpoint,
# the polynomial model is superior, yielding a Residual Standard Error
# of 3.71 and a high Adjusted $R^2$ of 0.955. These metrics indicate
# that the polynomial fit captures the local fluctuations of the 30 
# observations with high precision.In contrast, the non-linear
# exponential model (NLS), while more parsimonious with only 3 
# parameters compared to the polynomial's 6, shows significantly 
# lower precision with an RSE of 6.60.

# Question 4

x_new_data <- data.frame(x = 1:10)

conf_int <- predict(chosen_polynomial, newdata = x_new_data, interval = "confidence")

pred_int <- predict(chosen_polynomial, newdata = x_new_data, interval = "prediction")

print("Confidence Intervals (Average):")
head(conf_int)

print("Prediction Intervals (Individual):")
head(pred_int)

#todo: graphe

# Exercice 5

# Question 1)

library(readr)
sp <- read_csv("sp500_history.csv")
head(sp)
sp["daily_return"] <- sp["Close"]-sp["Open"]
head(sp)

# Question 2)

# First we fit the normal model by using the MLE estimators of the 
# mean and std.

returns <- sp$daily_return
mu_hat <- mean(returns)
sigma_hat <- sd(returns)
cat("Estimated Mean (mu):", mu_hat, "\n")
cat("Estimated Std Dev (sigma):", sigma_hat, "\n")

# Let us first plot the distribution of the normal distribution with these estimators
# along with the samples. The fit looks pretty good at first glance. 

hist(returns, probability = TRUE, main = "Normal Fit to Daily Returns", col = "lightblue")
x_seq <- seq(min(returns), max(returns), length.out = 100)
y_norm <- dnorm(x_seq, mean = mu_hat, sd = sigma_hat)
lines(x_seq, y_norm, col = "red", lwd = 2)

# We can go further by using a QQ plot:

qqnorm(returns)
qqline(returns, col = "red")

# We see that the normal distribution is not a good fit. Indeed, we see that
# the empirical distribution has much more extreme outliers than a normal distribution
# typically would.

# Question 3

dcomponents_hand <- function(theta, x) {
  p <- length(theta$pi)
  n <- length(x)
  probs <- matrix(0, nrow = n, ncol = p)
  for (j in 1:p) {
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
  
  list(parameters = theta[[t + 1]], loglik = -deviance[t + 1] / 2, deviance = deviance[t+1])
}

#########################################################
# 2. DENSITY VISUALIZATION FUNCTION
#########################################################

calc_mixture_density <- function(x_grid, params) {
  y_dens <- rep(0, length(x_grid))
  for (j in 1:length(params$pi)) {
    y_dens <- y_dens + params$pi[j] * dnorm(x_grid, mean = params$mu[j], sd = params$sigma[j])
  }
  return(y_dens)
}

#########################################################
# 3. GENERATING THE 5 PLOTS (p = 2 to 6)
#########################################################

# Ensure 'returns' variable exists
# returns <- data_3$y  # Uncomment if needed

# Setup grid: 2 rows, 3 columns
par(mfrow = c(2, 3))
x_grid <- seq(min(returns), max(returns), length.out = 1000)

for (p in 2:6) {
  # Fit model
  fit <- mixture_gaussian_hand(returns, p = p)
  
  # Plot Histogram
  hist(returns, breaks = 50, probability = TRUE, 
       main = paste("GMM: p =", p), 
       col = "lightgray", border = "white",
       xlab = "Returns", ylab = "Density")
  
  # Overlay Density
  y_values <- calc_mixture_density(x_grid, fit$parameters)
  lines(x_grid, y_values, col = "darkgreen", lwd = 2)
  
  # Add Log-Likelihood info
  legend("topright", legend = paste("LL:", round(fit$loglik, 2)), bty = "n", cex = 0.8)
}

# Reset layout
par(mfrow = c(1, 1))

#Visually, the Mixture of Normals is a vastly superior fit compared to the single Normal model. While the single Normal ($p=1$) fails to capture the high peak at the center and the extreme 
#outliers in the tails (as seen in the Q-Q plot), the Mixture models effectively allocate different components to handle different 'volatility regimes.' Specifically, a $p=2$ or $p=3$ model provides a much better description of the empirical density without the excessive complexity of higher $p$ values.


# Question 4.a)

################################################################################
# MOTIVATION FOR THE LOCATION-SCALE STUDENT-T MODEL
################################################################################
# The shift from a Gaussian (Normal) framework to a Student-t distribution is 
# primarily driven by the empirical failure of the Normal distribution to 
# account for extreme events (outliers) in financial data.
#
# 1. HANDLING "FAT TAILS" (LEPTOKURTOSIS):
#    - The key advantage is the 'nu' (degrees of freedom) parameter.
#    - MODELING EXTREMES: Unlike the Normal distribution, the Student-t has 
#      "heavier" tails. This allows it to assign a realistic probability to 
#      large returns (like the +/- 150 outliers seen in the QQ-plot), which 
#      the Normal model would consider impossible.
#    - FLEXIBILITY: As nu -> infinity, the Student-t converges to a Normal 
#      distribution. This allows the data to "choose" its own tail thickness.
#
# 2. PARSIMONY (SIMPLICITY):
#    - While a Mixture of Normals (GMM) can capture tails, it is complex. 
#      A GMM with p=3 requires 8 parameters.
#    - The Student-t model captures heavy tails using only 3 parameters:
#        * m (Location): The center of the returns.
#        * a (Scale): The volatility/spread.
#        * nu (Degrees of Freedom): Controls the thickness of the tails.
#    - This makes the model more robust and less prone to overfitting 
#      compared to high-component mixtures or high-degree polynomials.
################################################################################

# Question 4.b)

log_lik_student <- function(theta, x) {
  nu <- theta[1]; m <- theta[2]; a <- theta[3]
  if (nu <= 0 || a <= 0) return(-Inf) # Constraints
  
  # Formula provided in your exercise
  # We use log-densities for numerical stability
  n <- length(x)
  log_c <- lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(pi * nu * a^2)
  log_pdf <- log_c - ((nu + 1) / 2) * log(1 + (x - m)^2 / (nu * a^2))
  
  return(sum(log_pdf))
}

# --- 2. Optimization (Direct Maximization) ---
# We minimize the negative log-likelihood
fit_t <- optim(par = c(5, mean(returns), sd(returns)), 
               fn = function(theta) -log_lik_student(theta, returns),
               method = "Nelder-Mead")

theta_mle <- fit_t$par
cat("MLE Estimates: nu =", theta_mle[1], "m =", theta_mle[2], "a =", theta_mle[3], "\n")

# --- 3. Model Comparison (using BIC as per course criteria) ---
n <- length(returns)
logL_t <- -fit_t$value
bic_t <- -2 * logL_t + 3 * log(n)
cat("BIC for Student-t Model:", bic_t, "\n")

dt_ls <- function(x, nu, m, a) {
  term1 <- gamma((nu + 1) / 2) / (gamma(nu / 2) * sqrt(pi * nu * a^2))
  term2 <- (1 + (x - m)^2 / (nu * a^2))^(-(nu + 1) / 2)
  return(term1 * term2)
}

# 2. Plot the Empirical Distribution (Histogram)
hist(returns, breaks = 60, probability = TRUE, 
     col = "lightgray", border = "white",
     main = "Student-t Fit vs. Empirical Distribution",
     xlab = "Daily Returns", ylab = "Density")

# 3. Overlay the Fitted Student-t Curve
x_grid <- seq(min(returns), max(returns), length.out = 1000)
y_t <- dt_ls(x_grid, theta_mle[1], theta_mle[2], theta_mle[3])

lines(x_grid, y_t, col = "purple", lwd = 3)

# 4. Add a legend with the parameter estimates
legend("topright", 
       legend = c(paste("nu (df) =", round(theta_mle[1], 2)),
                  paste("m (loc) =", round(theta_mle[2], 2)),
                  paste("a (scale) =", round(theta_mle[3], 2))),
       col = "purple", lwd = 2, bty = "n")


y_sorted <- sort(returns)
n <- length(returns)

# 2. Calculate theoretical probabilities (using the mid-point rule)
probs <- (1:n - 0.5) / n

# 3. Calculate theoretical quantiles for the Location-Scale Student-t
# We use the parameters (nu, m, a) from your optim() fit
theoretical_quantiles <- theta_mle[2] + theta_mle[3] * qt(probs, df = theta_mle[1])

# 4. Create the Q-Q Plot
plot(theoretical_quantiles, y_sorted, 
     main = "Q-Q Plot: Location-Scale Student-t",
     xlab = "Theoretical Student-t Quantiles",
     ylab = "Empirical Quantiles (Returns)",
     pch = 20, col = "darkorchid")

# 5. Add a 45-degree reference line
abline(0, 1, col = "red", lwd = 2)

# Optional: Add a grid for better readability
grid()

# The Student fit still struggles wit outliers although it seems better
# than the normal ddistribution fit.

# Question 5)

n_obs <- length(returns)

# 1. BIC for Normal Model (Question 2)
logL_norm <- sum(dnorm(returns, mean = mean(returns), sd = sd(returns), log = TRUE))
bic_norm <- -2 * logL_norm + 2 * log(n_obs)

# 2. BIC for GMM (p = 5) (Question 3)

fit_gmm5 <- mixture_gaussian_hand(returns, p = 5)
bic_gmm5 <- -2 * fit_gmm5$loglik + (3 * 5 - 1) * log(n_obs)

# 3. BIC for Student-t Model (Question 4)
bic_t <- -2 * logL_t + 3 * log(n_obs)

bic_comparison <- data.frame(
  Model = c("Normal (p=1)", "GMM (p=5)", "Student-t"),
  Parameters = c(2, 14, 3),
  LogLikelihood = c(logL_norm, fit_gmm5$loglik, logL_t),
  BIC = c(bic_norm, bic_gmm5, bic_t)
)

print(bic_comparison)

winner <- bic_comparison$Model[which.min(bic_comparison$BIC)]
cat("\nAccording to the BIC criterion, the best model is:", winner, "\n")



# --- Theoretical Quantiles Calculation ---

# 1. Normal (Q2)
theoretical_q_norm <- mu_hat + sigma_hat * qnorm(probs)

# 2. GMM p=5 (Q3)
# We need an additional function for GMM quantiles using interpolation
get_gmm_quantiles <- function(probs, params) {
  # Create a fine grid for interpolation
  x_grid_q <- seq(min(returns)*1.5, max(returns)*1.5, length.out = 10000)
  cdf_gmm <- rep(0, length(x_grid_q))
  for (j in 1:length(params$pi)) {
    cdf_gmm <- cdf_gmm + params$pi[j] * pnorm(x_grid_q, mean = params$mu[j], sd = params$sigma[j])
  }
  # Interpolate to find quantiles corresponding to probs
  return(approx(x = cdf_gmm, y = x_grid_q, xout = probs)$y)
}
theoretical_q_gmm5 <- get_gmm_quantiles(probs, fit_gmm5$parameters)

# 3. Student-t (Q4b)
theoretical_q_t <- theta_mle[2] + theta_mle[3] * qt(probs, df = theta_mle[1])


# --- Side-by-Side Plotting ---

par(mfrow = c(1, 3), mar = c(5, 4, 4, 1)) # 3 plots per row, adjust margins

# Plot 1: Normal Model
plot(theoretical_q_norm, y_sorted, 
     main = "Q-Q Plot: Normal",
     xlab = "Theoretical Normal Quantiles",
     ylab = "Empirical Quantiles",
     pch = 20, col = "red")
abline(0, 1, col = "black", lwd = 2)
grid()

# Plot 2: GMM p=5 Model
plot(theoretical_q_gmm5, y_sorted, 
     main = "Q-Q Plot: GMM p=5",
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

# The QQ plot and the BIC for the GMM is the best. We use the GMM.


# Question 6)

# Due to the large amount of data (4301), the central limit theorem allows us to use a t-test.

length(y_sorted)

test_mu <- t.test(sp$daily_return,mu=0)
print(test_mu)
print(paste("We cannot conclude as the p-value equals", test_mu$p.value))

