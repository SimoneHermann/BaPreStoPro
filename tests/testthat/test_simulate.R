


model <- set.to.class("Diffusion", parameter = list(phi = 0.5, gamma2 = 0.01))
t <- seq(0, 1, by = 0.01)
data1 <- simulate(model, seed = 123, t = t, y0 = 0.5, plot.series = FALSE)
data2 <- simulate(model, seed = 123, t = t, y0 = 0.5, plot.series = FALSE)

expect_equal(data1, data2)

