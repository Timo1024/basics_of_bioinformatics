y_axis <- function(x, t) {
    x * exp(1) ** (-x * (1 - t))
}

y_axis_1 <- function(x) {y_axis(x, 0.2)}
y_axis_2 <- function(x) {y_axis(x, 0.3)}
y_axis_3 <- function(x) {y_axis(x, 0.4)}
y_axis_4 <- function(x) {y_axis(x, 0.5)}
y_axis_5 <- function(x) {y_axis(x, 0.6)}

smoothness <- 1000

pdf("plot.pdf")

curve(y_axis_1, from = 0, to = 20, col = 2, n = smoothness,
    ylim = c(0, 1),
    ylab = " expected number of contigs scaled in units G/L",
    xlab = "coverage c")  # Draw Base R plot
curve(y_axis_2, from = 0, to = 20, col = 3, add = TRUE, n = smoothness)
curve(y_axis_3, from = 0, to = 20, col = 4, add = TRUE, n = smoothness)
curve(y_axis_4, from = 0, to = 20, col = 5, add = TRUE, n = smoothness)
curve(y_axis_5, from = 0, to = 20, col = 6, add = TRUE, n = smoothness)
# plot(x_axis, y_axis(c(0.2, 0.4, 0.6, 0.8)), type = "l")

dev.off()