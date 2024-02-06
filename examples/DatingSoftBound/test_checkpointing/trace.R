library(coda)

trace <- read.delim("complete_sample.txt", header=TRUE)

trace <- trace[, -1]

trace.mcmc <- coda::as.mcmc(trace)
