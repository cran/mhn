## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4,
  dpi = 72
)
set.seed(1)

## ----library, include = FALSE-------------------------------------------------
library(mhn)

## ----density-shape, echo = -c(2, 8), fig.cap = "MHN density as alpha crosses the log-concavity threshold; beta = 1, gamma = 1 fixed."----
x  <- seq(0.001, 4, length.out = 400)
op <- par(mar = c(4, 4, 2, 1))
plot(x, dmhn(x, alpha = 2,   beta = 1, gamma = 1), type = "l", lwd = 2,
     ylim = c(0, 1.4), ylab = "f(x)",
     main = "MHN density at beta = 1, gamma = 1")
lines(x, dmhn(x, alpha = 5,   beta = 1, gamma = 1), lwd = 2, col = "steelblue")
lines(x, dmhn(x, alpha = 1,   beta = 1, gamma = 1), lwd = 2, col = "tomato")
lines(x, dmhn(x, alpha = 0.5, beta = 1, gamma = 1), lwd = 2, col = "darkorange")
legend("topright", bty = "n", lwd = 2,
       col = c("darkorange", "tomato", "black", "steelblue"),
       legend = c("alpha = 0.5  (no interior mode)",
                  "alpha = 1    (truncated normal)",
                  "alpha = 2    (log-concave)",
                  "alpha = 5    (log-concave, shifted right)"))
par(op)

## ----rtdr-regions-------------------------------------------------------------
cases <- data.frame(
  region = c("(a)", "(b)", "(c)", "(d)"),
  alpha  = c(2.0,    0.7,    0.3,    0.3),
  gamma  = c(0.5,    1.0,    0.5,    5.0)
)
set.seed(42)
empirical <- mapply(
  function(a, g) mean(rmhn(5000, alpha = a, beta = 1, gamma = g, method = "rtdr")),
  cases$alpha, cases$gamma
)
theoretical <- mapply(mhn_mean, cases$alpha, beta = 1, cases$gamma)
data.frame(cases,
           empirical_mean   = round(empirical,   4),
           theoretical_mean = round(theoretical, 4))

## ----auto-dispatch------------------------------------------------------------
# Branch: alpha < 1 and gamma > 0  -->  expects RTDR
set.seed(7); auto1 <- rmhn(20, alpha = 0.7, beta = 1, gamma = 2)
set.seed(7); rtdr1 <- rmhn(20, alpha = 0.7, beta = 1, gamma = 2, method = "rtdr")
identical(auto1, rtdr1)

# Branch: alpha > 1 and gamma > 0  -->  expects Sun Algorithm 1
set.seed(7); auto2 <- rmhn(20, alpha = 3,   beta = 1, gamma = 2)
set.seed(7); sun2  <- rmhn(20, alpha = 3,   beta = 1, gamma = 2, method = "sun")
identical(auto2, sun2)

# Branch: gamma < 0 with samples_per_setup < 25  -->  expects Sun Algorithm 3
set.seed(7); auto3 <- rmhn(5,  alpha = 2,   beta = 1, gamma = -1)
set.seed(7); sun3  <- rmhn(5,  alpha = 2,   beta = 1, gamma = -1, method = "sun")
identical(auto3, sun3)

