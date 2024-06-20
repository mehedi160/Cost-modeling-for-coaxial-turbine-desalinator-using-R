# Load necessary libraries
library(ggplot2)
library(pracma)

# Set plot theme
theme_set(theme_minimal(base_size = 15, base_family = "serif") + theme(text = element_text(size = 15)))

# Define the data
data <- matrix(c(1, 0.0292895, 0.0074324, 2.92895e-02,
                 2, 0.113752, 0.030239, 5.68760e-02,
                 3, 0.255394, 0.074549, 8.51314e-02,
                 4, 0.391834, 0.129176, 9.79584e-02,
                 5, 0.464199, 0.16741, 9.28397e-02,
                 6, 0.463196, 0.16681, 7.71993e-02,
                 7, 0.441931, 0.154578, 6.31330e-02,
                 8, 0.41121, 0.13852, 5.14012e-02,
                 9, 0.370715, 0.11956, 4.11905e-02,
                 10, 0.320892, 0.098771, 3.20892e-02), 
               byrow = TRUE, ncol = 4)

U_ <- c(0.7, 1.0, 1.3)

# Create an empty data frame to store the results
results <- data.frame()

# Loop through each U value and compute the corresponding SEC
for (U in U_) {
  rho_w <- 1024
  eta_gb <- 0.90
  rho_ss <- 8000
  sigma_y <- 250e6
  g <- 9.81
  mu_k <- 0.0015
  A_p <- 1.5 / 3.6e11
  b <- 1.0
  tsr_u <- 3
  tsr_d <- 4
  r_f <- 9
  l_ <- r_f
  r_0 <- r_f * 0.1
  r_in <- r_0 * 0.1
  eta <- 0.9
  d_gap <- r_0 / eta - r_0
  mu_a <- 0.00103
  K_ <- 0.801
  C_f <- 35
  Pi_0 <- 2.8e6
  
  p_r <- 0.5 * data[tsr_u, 2] * rho_w * pi * r_f^2 * U^3 + 0.5 * data[tsr_d, 2] * rho_w * pi * (b * r_f)^2 * ((1 - 2 * data[tsr_u, 3]) * U)^3
  T_u <- 0.5 * data[tsr_u, 4] * rho_w * pi * r_f^3 * U^2
  T_d <- 0.5 * data[tsr_d, 4] * rho_w * pi * (b * r_f)^3 * ((1 - 2 * data[tsr_u, 3]) * U)^2
  w_u <- tsr_u * U / r_f
  w_d <- tsr_d * (1 - 2 * data[tsr_u, 3]) * U / (b * r_f)
  
  a <- pi * rho_ss * rho_w * r_0^2 * mu_k * r_in * (r_0^2 - r_in^2) * (l_ + r_0) * g / sigma_y
  b1 <- 2 * pi * mu_a * r_0^3 * l_ / d_gap
  c <- mu_k * g * r_in * rho_w * pi * r_0^2 * l_
  d <- -eta_gb * (T_u * w_u + T_d * w_d)
  
  # Solving for w_c using uniroot
  w_c_func <- function(x) a * x^3 + b1 * x^2 + c * x + d
  w_c <- uniroot(w_c_func, c(0, 1000))$root
  
  del_p <- 0.5 * rho_w * w_c^2 * (r_0^2 - r_in^2) + rho_w * g * 50
  
  R <- seq(0.1, 1 - Pi_0 / del_p, length.out = 10)
  phi <- Pi_0 - 2 * del_p
  A_m <- 2 * pi * r_0 * l_
  Q_p <- sapply(R, function(R_) {
    phi^2 * A_p * A_m / (2 * Pi_0 * log((Pi_0 + phi) / (Pi_0 + phi * (1 - R_))) / R_ - 2 * phi)
  })
  
  SEC <- (T_u * w_u + T_d * w_d) / Q_p * 2.77778e-7
  
  # Store the results
  results <- rbind(results, data.frame(U = U, R = R, SEC = SEC))
}

# Plotting
p <- ggplot(results, aes(x = R, y = SEC, color = as.factor(U))) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  labs(x = expression(paste("Recovery factor, ", italic(R), " [%]")),
       y = expression(paste("Optimal SEC, ", gamma["CTD"]^"*", " [kWh/m"^3, "]")),
       color = expression(paste(U, "[m/s]"))) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_line(size = 0.5),
    panel.grid.minor = element_line(size = 0.25),
    axis.line = element_line(size = 1.5),
    axis.ticks = element_line(size = 1.5),
    legend.position = "best"
  )

# Saving the plot
ggsave("optimal_sec_R.png", plot = p, dpi = 300)
