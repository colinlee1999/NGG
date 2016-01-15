library('ggplot2')
library('parallel')
library('data.table')

sizes = c(1000, 5000, 10000, 50000, 1e5, 5e5, 1e6, 1e7)
testFun <- function(x, size){
  set.seed(1)
  sum(rnorm(size))
}
results = data.frame(size=c(), type=c(), time=c())

cl <- makeCluster(4)

for(size in sizes){
  ## parallel computing
  beg = Sys.time()
  z=clusterApply(cl, 1:1000, testFun, size=size)
  td = as.numeric(Sys.time() - beg, "secs")
  results = rbind(results, data.frame(size=size, type="quad_core", time=td))
  
  ## single threaded computing (to compare times and code)
  beg = Sys.time()
  z=lapply(1:1000, testFun, size=size)
  td = as.numeric(Sys.time() - beg, "secs")
  results = rbind(results, data.frame(size=size, type="single_core", time=td))
}

stopCluster(cl)

postscript('temp.ps')

ggplot(results, aes(x=size, y=time, group=type, color=type))+
  geom_point()+
  scale_x_log10()
dt=data.table(results)

dev.off()