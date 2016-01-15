library('ggplot2')
library('parallel')
library('data.table')
library(foreach)
library(doMC)
registerDoMC(cores = 4)

sizes = c(1000, 5000, 10000, 50000, 1e5, 5e5, 1e6)
testFun <- function(x, size){
  set.seed(1)
  sum(rnorm(size))
}

results = data.frame(size=c(), type=c(), time=c())

for (size in sizes)
{
  beg = Sys.time()
  # z = apply(matrix(1:1000,1000,1), 1 , testFun, size = size)
  foreach(i=1:sizes) %dopar% sum(runif(100000))
  td = as.numeric(Sys.time() - beg, "secs")
  results = rbind(results, data.frame(size=size, type="quad_core", time=td))
}

# results = data.frame(size=c(), type=c(), time=c())

# cl <- makeCluster(4)

# for(size in sizes){
#   ## parallel computing
#   beg = Sys.time()
#   z=clusterApply(cl, 1:1000, testFun, size=size)
#   td = as.numeric(Sys.time() - beg, "secs")
#   results = rbind(results, data.frame(size=size, type="quad_core", time=td))
  
#   ## single threaded computing (to compare times and code)
#   beg = Sys.time()
#   z=lapply(1:1000, testFun, size=size)
#   td = as.numeric(Sys.time() - beg, "secs")
#   results = rbind(results, data.frame(size=size, type="single_core", time=td))
# }

# stopCluster(cl)

dt=data.table(results)

print(dt)
