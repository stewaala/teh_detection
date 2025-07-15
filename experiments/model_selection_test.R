library(data.table)
library(ggplot2)

is_Arithmetic <- FALSE
set.seed(42)

#––– Parameters ––––––––––––––––––––––––––––––––––––––––––––––––––
sample_sizes <- c(100, 500, 1000)
K_values  <- if (is_Arithmetic) c(0, 1.5, 3) else c(3, 5, 7)
n_rep        <- 100

#––– Compute BIC accuracy under the arithmetic DGP –––––––––––––––––––––––
results <- rbindlist(
  lapply(sample_sizes, function(n) {
    rbindlist(
      lapply(K_values, function(K) {
        n_correct <- sum(
          replicate(n_rep, {
            if (is_Arithmetic) {
              df <- simulate_arithmetic(n = n, K = K)
            }
            else {
              df <- simulate_geometric(n = n, K = K)
            }

            m_id  <- glm(Y ~ X1 + X2 + Z + A + Z:A,
                         family = gaussian(link = "identity"),
                         data   = df)
            m_log <- glm(Y ~ X1 + X2 + Z + A + Z:A,
                         family = Gamma(link = "log"),
                         data   = df)
            bic_id  <- BIC(m_id)
            bic_log <- BIC(m_log)
            if (is_Arithmetic) bic_id < bic_log else bic_id > bic_log
          })
        )
        data.table(n = n, K = K, pct = (n_correct / n_rep) * 100)
      })
    )
  })
)

#––– Plot ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
ggplot(results, aes(x = factor(n), y = pct, color = factor(K), group = factor(K))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = if (is_Arithmetic) "BIC Accuracy on Arithmetic DGP" else "BIC Accuracy on Geometric DGP",
    x     = "Sample size (n)",
    y     = "BIC % correct (identity model)",
    color = expression(K)
  ) +
  theme_minimal(base_size = 14)

