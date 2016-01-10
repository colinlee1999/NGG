gen_data_by_NGG_lp <- function(G, n, psi, t_pi)
{
  data = matrix(, nrow = G, ncol = n)
  cluster_info = matrix(rep(0,G*3),G,3)

  colnames(cluster_info) = c("true_cluster","est_cluster","flag")

  for (row in 1:G)
  {
    # determine which cluster, using a uniform distribution to compare with t_pi
    temp = runif(1)
    if (temp<t_pi[1]) cluster = 1
    else if(temp<t_pi[1]+t_pi[2]) cluster = 2
    else cluster = 3

    cluster_info[row,1] = cluster

    alpha = psi[cluster * 4 - 3]
    beta = psi[cluster * 4 - 2]
    if (cluster<3)
    {
      xi = psi[cluster * 4 - 1]
      eta = psi[cluster * 4]  
    }

    tau_g = rgamma(1, alpha, beta)

    if (cluster == 1)
    {
      mu_g = rgamma(1, shape = xi, rate = eta)
    }
    else if (cluster == 2)
    {
      mu_g = - rgamma(1, shape = xi, rate = eta)
    }
    else if (cluster == 3)
    {
      mu_g = 0
    }

    data[row,] = rnorm(n,mean = mu_g, sd = sqrt(1/tau_g))
  }

  return (ExpressionSet(assayData = data, featureData = new("AnnotatedDataFrame", data = data.frame(cluster_info)))) 
}