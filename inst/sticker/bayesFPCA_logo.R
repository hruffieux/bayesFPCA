# Set seed for reproducibility
set.seed(1)

# Generate sample data for a smooth functional curve with credible intervals
q <- 100
time <- seq(1, q, length.out = q)
mean_curve <- sin(time / 10)+cos(time/8)
lower_band <- mean_curve - 0.5
upper_band <- mean_curve + 0.5

curve_data <- data.frame(
  time = time,
  mean = mean_curve,
  lower = lower_band,
  upper = upper_band
)

# Generate fewer points for added visual interest
nb_points <- 6
points_indices <- sample(1:q, nb_points)
points_data <- data.frame(
  time = time[points_indices],
  value = mean_curve[points_indices] + rnorm(nb_points, mean = 0, sd = 0.5)
)

# Ensure most points fall within the credible band
points_data$value[points_data$value < lower_band[points_indices]] <- lower_band[points_indices]
points_data$value[points_data$value > upper_band[points_indices]] <- upper_band[points_indices]

# Plot the functional curves with credible intervals and points
plot_curves <- ggplot() +
  geom_ribbon(data = curve_data, aes(x = time, ymin = lower, ymax = upper), fill = "gray35", alpha = 0.5) +
  geom_line(data = curve_data, aes(x = time, y = mean), color = "firebrick3", size = 1.5) + # darkblue
  geom_point(data = points_data, aes(x = time, y = value), color = "white", size = 1.5) +
  theme_void() +
  theme(legend.position = "none")


# Create the hex sticker
sticker(
  subplot = plot_curves,
  package = "bayesFPCA",
  p_size = 64,  # Larger text size
  s_x = 1,
  s_y = 0.89,
  s_width = 1.6,
  s_height = 1.5,
  p_x = 0.89,
  p_y = 1.39,
  u_color = "white",
  u_size = 1,
  h_fill = "gray75",  # lighter grey background
  h_color = "gray75",
  filename = "man/figures/bayesFPCA_logo.png",
  dpi = 1200
)
