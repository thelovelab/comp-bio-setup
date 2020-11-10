cmd_args=commandArgs(TRUE)

dist <- cmd_args[1] # dist
n <- as.numeric(cmd_args[2]) # sample size
out <- cmd_args[3] # output filename

if (dist == "norm") {
  x <- rnorm(n)
} else if (dist == "unif") {
  x <- runif(n)
} else if (dist == "exp") {
  x <- rexp(n)
}
dat <- data.frame(dist=dist,x=x)
write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
