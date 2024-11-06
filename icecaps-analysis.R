install.packages("imager")
library(imager)

icecap_1979 <- load.image("project/data/sea-ice-concentration-map-1979-aug.jpeg")
dim(icecap_1979)

icecap_1979_gray <- grayscale(icecap_1979)
plot(icecap_1979_gray)

# SVD analysis
icecap_1979_svd <- svd(icecap_1979_gray)

# Plot singular values
variance_percent_D <- 100 * (icecap_1979_svd$d^2) / sum(icecap_1979_svd$d^2)
cum_percent_D <- cumsum(variance_percent_D)
modeK <- 1:length(icecap_1979_svd$d)
K <- 30

plot(modeK[1:K], variance_percent_D[1:K],
    type = "o", col = "blue",
    xlab = "Mode number", pch = 16,
    ylab = "Percentage of mode variance",
    main = "Scree Plot of first 30 eigenvalues of B/W Icecap Image 1979"
)

dev.off()


# Scree plot of first 30 modes of A, incl. cumulative variance
par(mar = c(4, 4, 2, 4), mgp = c(2.2, 0.7, 0))
plot(1:K,
    variance_percent_D[1:K],
    ylim = c(0, 100),
    type = "o",
    ylab = "Percentage of Variance [%]",
    xlab = "EOF Mode Number",
    cex.lab = 1.2, cex.axis = 1.1, lwd = 2,
    main = "Scree Plot of first 30 eigenvalues of B/W Sunset Photo"
)
legend(3, 30,
    col = c("black"), lty = 1, lwd = 2.0,
    legend = c("Percentange Variance"), bty = "n",
    text.font = 2, cex = 1.0, text.col = "black"
)

par(new = TRUE)
plot(1:K, cum_percent_D[1:K],
    ylim = c(90, 100), type = "o",
    col = "blue", lwd = 2, axes = FALSE,
    xlab = "", ylab = ""
)
legend(3, 94.5,
    col = c("blue"), lty = 1, lwd = 2.0,
    legend = c("Cumulative Percentage Variance"), bty = "n",
    text.font = 2, cex = 1.0, text.col = "blue"
)
axis(4, col = "blue", col.axis = "blue", mgp = c(3, 0.7, 0))
mtext("Cumulative Variance [%]",
    col = "blue",
    cex = 1.2, side = 4, line = 2
)
