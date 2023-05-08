varcomp = read.csv("varcompggm.csv")
samples <- as.numeric(seq(from = 100, to = 500, by = 50))
matplot(samples, varcomp,
main = "", xlab = "Samples", ylab = "Average Variance",
type = "b", pch = 20, col = rainbow(2))
legend("topright", inset=0.01, legend = colnames(varcomp), pch = 20, col = rainbow(2))