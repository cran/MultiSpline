make_data <- function(G = 40, n = 20, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(G), each = n),
                  sch = rep(seq_len(G / 4), each = 4 * n))
  d$x   <- stats::rexp(nrow(d), 1 / 3)
  d$grp <- sample(c("a", "b", "c"), nrow(d), TRUE)
  u <- stats::rnorm(G)
  d$y  <- 2 * sin(d$x / 2) + u[d$id] + stats::rnorm(nrow(d))
  d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(-0.5 + sin(d$x / 2) + u[d$id]))
  d
}
