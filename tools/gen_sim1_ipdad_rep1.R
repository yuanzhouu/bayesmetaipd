# Generate Simulation Study 1, replicate 1 (matches Data/Generate_Simulation_1.R)
n_rep <- 1
L <- 40
n_sample <- 200
true_mu <- c(0.5, 1, 0.5, 0.25)
true_Sigma_theta <- rep(1, 4)
true_kappa <- 1
p_theta <- length(true_mu)
p_beta <- p_theta
i_rep <- 1
set.seed(i_rep + 33)

theta_l_mat <- array(0, c(L, p_theta))
dimnames(theta_l_mat)[[2]] <- c("theta0", "theta1", "theta2", "theta3")
theta_hat_l_mat <- array(0, c(L, p_theta))
dimnames(theta_hat_l_mat)[[2]] <- c("theta0", "theta1", "theta2", "theta3")
X_cube <- array(0, c(L, n_sample, p_theta))
dimnames(X_cube)[[3]] <- c("Int", "X1", "X2", "X1X2")
Y_mat <- array(0, c(L, n_sample))
beta_mat <- array(0, c(L, p_beta))
dimnames(beta_mat)[[2]] <- c("beta0", "beta1", "beta2", "beta3")
V_beta_cube <- array(0, c(L, p_beta, p_beta))
dimnames(V_beta_cube)[[2]] <- dimnames(V_beta_cube)[[3]] <- c("beta0", "beta1", "beta2", "beta3")
type_vec <- c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10))

for (i_study in 1:L) {
  theta_l_mat[i_study, ] <- rnorm(n = p_theta, mean = true_mu, sd = sqrt(true_Sigma_theta))
  X_cube[i_study, , 1] <- rep(1, n_sample)
  upsilon_l <- rnorm(1, mean = 0, sd = sqrt(0.1))
  X_cube[i_study, , 2] <- rnorm(n_sample, mean = upsilon_l, sd = sqrt(0.5))
  X_cube[i_study, , 3] <- ifelse(runif(n_sample) < 0.5, 1, 0)
  X_cube[i_study, , 4] <- X_cube[i_study, , 2] * X_cube[i_study, , 3]
  x.theta <- X_cube[i_study, , ] %*% theta_l_mat[i_study, ]
  Y_mat[i_study, ] <- x.theta + rnorm(n_sample, mean = 0, sd = sqrt(true_kappa))
  tempData <- data.frame(
    Y = Y_mat[i_study, ],
    X1 = X_cube[i_study, , 2],
    X2 = X_cube[i_study, , 3],
    X1X2 = X_cube[i_study, , 4]
  )
  theta_lm <- lm(Y ~ X1 + X2 + X1X2, data = tempData)
  theta_hat_l_mat[i_study, ] <- summary(theta_lm)$coefficients[, 1]
  if (type_vec[i_study] == 1) {
    beta_lm <- lm(Y ~ X1 + X2, data = tempData)
    beta_mat[i_study, 2:3] <- summary(beta_lm)$coefficients[2:3, 1]
    V_beta_cube[i_study, 2:3, 2:3] <- vcov(summary(beta_lm))[2:3, 2:3]
  }
  if (type_vec[i_study] == 2) {
    phi <- rep(0, n_sample)
    for (i in 1:n_sample) {
      if (tempData$X1[i] > 0 && tempData$X2[i] == 0) {
        phi[i] <- 1
      } else if (tempData$X1[i] > 0 && tempData$X2[i] == 1) {
        phi[i] <- 2
      } else if (tempData$X1[i] <= 0 && tempData$X2[i] == 0) {
        phi[i] <- 3
      } else {
        phi[i] <- 4
      }
    }
    ind.1 <- 1 * (phi == 1)
    ind.2 <- 1 * (phi == 2)
    ind.3 <- 1 * (phi == 3)
    ind.4 <- 1 * (phi == 4)
    beta_lm <- lm(tempData$Y ~ ind.1 + ind.2 + ind.3 + ind.4 - 1)
    beta_mat[i_study, ] <- summary(beta_lm)$coefficients[, 1]
    V_beta_cube[i_study, , ] <- vcov(summary(beta_lm))
  }
  if (type_vec[i_study] == 3) {
    beta_lm <- lm(Y ~ X1 + X2 + X1X2, data = tempData)
    beta_mat[i_study, 3:4] <- summary(beta_lm)$coefficients[3:4, 1]
    V_beta_cube[i_study, 3:4, 3:4] <- vcov(summary(beta_lm))[3:4, 3:4]
  }
}

SimulData <- list(list(
  theta_l_mat = theta_l_mat,
  theta_hat_l_mat = theta_hat_l_mat,
  X_cube = X_cube,
  Y_mat = Y_mat,
  type_vec = type_vec,
  beta_mat = beta_mat,
  V_beta_cube = V_beta_cube,
  true_mu = true_mu,
  true_kappa = true_kappa
))

root <- "C:/Users/Yuan Zhou/OneDrive - University of Cincinnati/UCsemester/survey_sampling/paper_repo"
outdir <- file.path(root, "Bayesian-Meta/Data")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
save.image(file = file.path(outdir, "SimulationData_1.RData"))

d1 <- SimulData[[1]]
sim1_ipdad_rep1 <- list(
  X_cube = d1$X_cube,
  Y_mat = d1$Y_mat,
  theta_l_mat = d1$theta_l_mat,
  beta_mat = d1$beta_mat,
  V_beta_cube = d1$V_beta_cube,
  type_vec = d1$type_vec,
  true_mu = true_mu,
  true_kappa = true_kappa,
  random_seed = .Random.seed
)
pkgdata <- file.path(root, "bayesmetaipd/data")
if (!dir.exists(pkgdata)) dir.create(pkgdata, recursive = TRUE)
save(sim1_ipdad_rep1, file = file.path(pkgdata, "sim1_ipdad_rep1.rda"), compress = "xz")
cat("saved sim1_ipdad_rep1\n")
