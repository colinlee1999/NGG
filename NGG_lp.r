library(iCheck)

NGG_lp <- function(
  E_Set,
  b = c(2,2,2), 
  t_pi_prior = c(0.05, 0.05, 0.90),
  is_sim = 0, 
  verbose = 0,
  infinity = 1e100,
  converge_threshold = 1e-6,
  param_limit_min = c(-5,-5,-5,-5,-5,-5,-5,-5,-5,-5),
  param_limit_max = c(5,5,5,5,5,5,5,5,5,5),
  # param_limit_min = c(-3,-3,-3,-3,-3,-3,-3,-3,-3,-3),
  # param_limit_max = c(3,3,3,3,3,3,3,3,3,3),
  # param_limit_min = c(-2,-2,-2,-2,-2,-2,-2,-2,-2,-2),
  # param_limit_max = c(2,2,2,2,2,2,2,2,2,2),
  # param_limit_min = c(-6,-6,-6,-6,-6,-6,-6,-6,-6,-6),
  # param_limit_max = c(6,6,6,6,6,6,6,6,6,6),
  # param_limit_min = c(-10,-10,-10,-10,-10,-10,-10,-10,-10,-10),
  # param_limit_max = c(10,10,10,10,10,10,10,10,10,10),
  max_iteration_num_in_optim = 100,
  max_repeated_times = 500,
  M = 1000,
  limma_prior = 1
  )
{
  
  get_A_B <- function(
    beta, 
    sum_dgl_by_l_i,
    sum_dgl_square_by_l_i,
    n,
    mu_g)
  {
    result = beta+(sum_dgl_square_by_l_i - 2*mu_g*sum_dgl_by_l_i + n * mu_g^2)/2

    return(result)
  }

  get_params <- function(
    xi,
    eta,
    alpha,
    beta,
    n,
    lower,
    upper,
    M,
    sum_dgl_by_l_i,
    sum_dgl_square_by_l_i)
  {
    result = list(
      xi,
      eta,
      alpha,
      beta,
      n,
      lower,
      upper,
      M,
      0,
      0)

    names(result) = c(
      "xi",
      "eta",
      "alpha",
      "beta",
      "n",
      "lower",
      "upper",
      "M",
      "sum_dgl_by_l_i",
      "sum_dgl_square_by_l_i")

    return(result)
  }

  get_sample_gen_params <- function(
    xi,
    eta)
  {
    result = list(
      xi,
      eta)

    names(result) = c(
      "shape",
      "rate")

    return(result)
  }

  func_over_expressed <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    A = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = 1 / (A)^(n/2+alpha)
    return(result)
  }

  func_over_expressed_d_alpha <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    A = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = -log(A) / (A)^(n/2+alpha)
    return(result)
  }

  func_over_expressed_d_beta <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    A = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = -(n/2+alpha) / (A)^(n/2+alpha+1)
    return(result)
  }

  func_over_expressed_d_xi <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    A = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = (log(eta) + log(mu_g) - digamma(xi)) / (A)^(n/2+alpha)
    return(result)
  }

  func_over_expressed_d_eta <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    A = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = (xi/eta - mu_g) / (A)^(n/2+alpha)
    return(result)
  }

  func_under_expressed <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    B = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = 1 / (B)^(n/2+alpha)
    return(result)
  }

  func_under_expressed_d_alpha <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    B = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = -log(B) / (B)^(n/2+alpha)
    return(result)
  }

  func_under_expressed_d_beta <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    B = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = -(n/2+alpha) / (B)^(n/2+alpha+1)
    return(result)
  }

  func_under_expressed_d_xi <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    B = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = (log(eta) + log(-mu_g) - digamma(xi)) / (B)^(n/2+alpha)
    return(result)
  }

  func_under_expressed_d_eta <- function(
    t,
    params)
  {
    xi = params$"xi"
    eta = params$'eta'
    alpha = params$'alpha'
    beta = params$'beta'
    n = params$'n'
    sum_dgl_by_l_i = params$'sum_dgl_by_l_i'
    sum_dgl_square_by_l_i = params$'sum_dgl_square_by_l_i'

    mu_g = t / (1 - t^2)
    B = get_A_B(beta, sum_dgl_by_l_i, sum_dgl_square_by_l_i, n, mu_g)

    result = (xi/eta + mu_g) / (B)^(n/2+alpha)
    return(result)
  }

  func_gamma_cluster_1 <- function(
    t,
    params)
  {
    alpha = params$'shape'
    beta = params$'rate'
    result = beta^alpha/gamma(alpha) * 
              (t/(1-t^2))^(alpha-1) * 
              exp(-beta * (t/(1-t^2))) * 
              (t^2 + 1) / (t^2 - 1)^2
    return(result)
  }

  func_gamma_cluster_2 <- function(
    t,
    params)
  {
    alpha = params$'shape'
    beta = params$'rate'
    result = beta^alpha/gamma(alpha) * 
              (-t/(1-t^2))^(alpha-1) * 
              exp(-beta * (-t/(1-t^2))) * 
              (t^2 + 1) / (t^2 - 1)^2
    return(result)
  }

  func_synthesize <- function(
    t,
    func,
    params,
    sample_gen_func,
    sample_gen_params)
  {
    result = func(t, params) * sample_gen_func(t, sample_gen_params)
    return(result)
  }

  run_integral <- function(
    dat,
    func,
    params,
    sample_gen_func,
    sample_gen_params)
  {
    params$sum_dgl_by_l_i = dat[1]
    params$sum_dgl_square_by_l_i = dat[2]
    result = integrate(
              func_synthesize, 
              lower = params$'lower', 
              upper = params$'upper',
              func,
              params,
              sample_gen_func,
              sample_gen_params,
              subdivisions = params$'M'
              )
    return(result$value)
  }

  lf123 <- function(
    psi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n,
    M)
  {
    excerpt <- function(func_passed_in, func_passed_in_2, cluster_1)
    {
      xi = exp(delta)
      eta = exp(theta)
      alpha = exp(lambda)
      beta = exp(nu)
      if (cluster_1)
      {
        lower = 0
        upper = 1
      }
      else
      {
        lower = -1
        upper = 0
      }

      params = get_params(
        xi,
        eta,
        alpha,
        beta,
        n,
        lower,
        upper,
        M,
        0,
        0)

      sample_gen_params = get_sample_gen_params(
        xi,
        eta)

      integral_result = apply(
        cbind(
          sum_dgl_by_l,
          sum_dgl_square_by_l
          ),
        1,
        run_integral,
        func_passed_in,
        params,
        func_passed_in_2,
        sample_gen_params
        )

      result = alpha * log(beta) + lgamma(n/2 + alpha) - lgamma(alpha) - n/2 * log(2*pi) + log(integral_result)
      return(result)
    }

    # cluster 1
    lambda = psi[1]
    nu = psi[2]
    delta = psi[3]
    theta = psi[4]

    logf1 = excerpt(func_over_expressed, func_gamma_cluster_1, 1)

    #cluster 2
    lambda = psi[5]
    nu = psi[6]
    delta = psi[7]
    theta = psi[8]

    logf2 = excerpt(func_under_expressed, func_gamma_cluster_2, 0)

    #cluster 3
    lambda = psi[9]
    nu = psi[10]

    alpha = exp(lambda)
    beta = exp(nu)

    logf3 = alpha * log(beta) + lgamma(n/2 + alpha) - lgamma(alpha) - n/2 * log(2*pi) - (n/2 + alpha) * log (beta + sum_dgl_square_by_l/2)

    result = cbind(
      logf1,
      logf2,
      logf3)
    colnames(result) = c("logf1", "logf2", "logf3")

    return (result)
  }

  get_tilde_z <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n,
    M)
  {
    logf = lf123(
      psi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      M)
    max_logf = apply(logf, 1, max, na.rm = TRUE)
    t1 = t_pi[1] * exp(logf[,1] - max_logf)
    t2 = t_pi[2] * exp(logf[,2] - max_logf)
    t3 = t_pi[3] * exp(logf[,3] - max_logf)
    total = t1 + t2 + t3

    result = cbind(t1, t2, t3)/total
    return(result)
  }

  l_c <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b, 
    converge_threshold, 
    infinity,
    M)
  {
    # if (t_pi[1]<converge_threshold && t_pi[2]<converge_threshold) return (-infinity)
    # if (t_pi[1]<converge_threshold && t_pi[3]<converge_threshold) return (-infinity)
    # if (t_pi[2]<converge_threshold && t_pi[3]<converge_threshold) return (-infinity)

    logf = lf123(
      psi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n,
      M)
    if (verbose)
    {
      cat("max logf>>",max(logf),"\n")
      cat("min logf>>",min(logf),"\n")
      cat("psi>>\n",psi,"\n")  
    }    

    result = 0
    result = result + sum(tilde_z[,1] * logf[,1], na.rm = TRUE)
    result = result + sum(tilde_z[,2] * logf[,2], na.rm = TRUE)
    result = result + sum(tilde_z[,3] * logf[,3], na.rm = TRUE)
    result = result + sum(tilde_z[,1] * log(t_pi[1]), na.rm = TRUE)
    result = result + sum(tilde_z[,2] * log(t_pi[2]), na.rm = TRUE)
    result = result + sum(tilde_z[,3] * log(t_pi[3]), na.rm = TRUE)

    #with dirichlet distribution
    result = result + lgamma(b[1]+b[2]+b[3]) - lgamma(b[1]) - lgamma(b[2]) - lgamma(b[3])
    result = result + sum((b-1) * log(t_pi), na.rm = TRUE)

    return (result)
  }

  negative_l_c <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b, 
    converge_threshold, 
    infinity,
    M)
  {
    return (-l_c(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z, 
      b, 
      converge_threshold, 
      infinity,
      M))
  }

  gradient_l_c <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b, 
    converge_threshold, 
    infinity,
    M)
  {
    result = psi
    
    # cluster 1
    lambda = psi[1]
    nu = psi[2]
    delta = psi[3]
    theta = psi[4]

    cluster = 1

    alpha = exp(lambda)
    beta = exp(nu)
    xi = exp(delta)
    eta = exp(theta)

    lower = 0
    upper = 1

    params = get_params(
      xi,
      eta,
      alpha,
      beta,
      n,
      lower,
      upper,
      M,
      0,
      0)
    sample_gen_params = get_sample_gen_params(
      xi,
      eta)

    shared_denominator = apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_over_expressed,
            params,
            func_gamma_cluster_1,
            sample_gen_params
            )

    d_lambda = alpha * sum(
      tilde_z[,cluster] *
      (log(beta) + digamma(n/2+alpha) - digamma(alpha) +
       (
          apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_over_expressed_d_alpha,
            params,
            func_gamma_cluster_1,
            sample_gen_params
            )
        )/shared_denominator
      ))

    d_nu = beta * sum(
      tilde_z[,cluster] * 
      (alpha/beta +
       (
          apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_over_expressed_d_beta,
            params,
            func_gamma_cluster_1,
            sample_gen_params
            )
        )/shared_denominator
      ))

    d_delta = xi * sum(
      tilde_z[,cluster] *
      (
       (
          apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_over_expressed_d_xi,
            params,
            func_gamma_cluster_1,
            sample_gen_params
            )
        )/shared_denominator
      ))

    d_theta = eta * sum(
      tilde_z[,cluster] *
      (
       (
          apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_over_expressed_d_eta,
            params,
            func_gamma_cluster_1,
            sample_gen_params
            )
        )/shared_denominator
      ))

    result[1] = d_lambda
    result[2] = d_nu
    result[3] = d_delta
    result[4] = d_theta

    # cluster 2
    lambda = psi[5]
    nu = psi[6]
    delta = psi[7]
    theta = psi[8]
    
    cluster = 2

    alpha = exp(lambda)
    beta = exp(nu)
    xi = exp(delta)
    eta = exp(theta)

    lower = -1
    upper = 0

    params = get_params(
      xi,
      eta,
      alpha,
      beta,
      n,
      lower,
      upper,
      M,
      0,
      0)    

    sample_gen_params = get_sample_gen_params(
      xi,
      eta)

    shared_denominator = apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_under_expressed,
            params,
            func_gamma_cluster_2,
            sample_gen_params
            )

    d_lambda = alpha * sum(
      tilde_z[,cluster] *
      (log(beta) + digamma(n/2+alpha) - digamma(alpha) +
       (
          apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_under_expressed_d_alpha,
            params,
            func_gamma_cluster_2,
            sample_gen_params
            )
        )/shared_denominator
      ))

    d_nu = beta * sum(
      tilde_z[,cluster] *
      (alpha/beta +
       (
          apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_under_expressed_d_beta,
            params,
            func_gamma_cluster_2,
            sample_gen_params
            )
        )/shared_denominator
      ))

    d_delta = xi * sum(
      tilde_z[,cluster] *
      (
       (
          apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_under_expressed_d_xi,
            params,
            func_gamma_cluster_2,
            sample_gen_params
            )
        )/shared_denominator
      ))

    d_theta = eta * sum(
      tilde_z[,cluster] *
      (
       (
          apply(
            cbind(
              sum_dgl_by_l,
              sum_dgl_square_by_l
              ),
            1,
            run_integral,
            func_under_expressed_d_eta,
            params,
            func_gamma_cluster_2,
            sample_gen_params
            )
        )/shared_denominator
      ))

    result[5] = d_lambda
    result[6] = d_nu
    result[7] = d_delta
    result[8] = d_theta

    # cluster 3
    lambda = psi[9]
    nu = psi[10]
    
    cluster = 3

    alpha = exp(lambda)
    beta = exp(nu)

    d_lambda = alpha * sum(
      tilde_z[,cluster] * (
        log(beta) + digamma(n/2+alpha) - digamma(alpha)
        - log(beta+sum_dgl_square_by_l/2)
        )
      )

    d_nu = beta * sum(
      tilde_z[,cluster] * (
        alpha/beta - (n/2 + alpha) / (beta + sum_dgl_square_by_l/2)
        )
      )

    result[9] = d_lambda
    result[10] = d_nu

    if (verbose)
    {
      cat("gradient_l_c>>\n",result,"\n")  
    }    

    return (result)
  }

  gradient_negative_l_c <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b, 
    converge_threshold, 
    infinity,
    M)
  {
    return (-gradient_l_c(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z, 
      b, 
      converge_threshold, 
      infinity,
      M))
  }

  get_t_pi <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b,
    converge_threshold,
    infinity)
  {
    #with Dirichlet distribution
    denominator = length(sum_dgl_by_l) + sum(b, na.rm=TRUE) - 3

    t1 = (sum(tilde_z[,1], na.rm=TRUE) + b[1] - 1) / denominator
    t2 = (sum(tilde_z[,2], na.rm=TRUE) + b[2] - 1) / denominator
    t3 = (sum(tilde_z[,3], na.rm=TRUE) + b[3] - 1) / denominator

    return (c(t1,t2,t3))
  }

  # function body

  G = nrow(E_Set)
  n = ncol(E_Set)

  data_matrix_of_E_Set = exprs(E_Set)

  # 'sum_dgl_by_l' is an G * 1 matrix, the summation result of every row of 'data_matrix_of_E_Set'
  sum_dgl_by_l = apply(data_matrix_of_E_Set,1,sum, na.rm=TRUE)

  # 'sum_dgl_square_by_l' is an G * 1 matrix, the summation of every squared elements of 'data_matrix_of_E_Set' by row
  sum_dgl_square_by_l = apply(data_matrix_of_E_Set^2,1,sum, na.rm = TRUE)

  column_names = colnames(fData(E_Set))

  if (limma_prior)
  {
    result_limma = lmFitPaired(
      E_Set,
      probeID.var = column_names[1],
      gene.var = column_names[2],
      chr.var = column_names[3],
      verbose = FALSE)

    frame_unsorted = result_limma$frame.unsorted
    over_expressed_sub_script = frame_unsorted$pos[which(frame_unsorted$stats > 0 & frame_unsorted$p.adj < 0.05)]
    under_expressed_sub_script = frame_unsorted$pos[which(frame_unsorted$stats < 0 & frame_unsorted$p.adj < 0.05)]
    non_diff_sub_script = frame_unsorted$pos[which(frame_unsorted$p.adj >= 0.05)]

    t_pi_prior = c(
      length(under_expressed_sub_script)/G,
      length(over_expressed_sub_script)/G,
      0)
    t_pi_prior[3] = 1 - t_pi_prior[1] - t_pi_prior[2]

    over_expressed_E_Set = E_Set[over_expressed_sub_script]
    under_expressed_E_Set = E_Set[under_expressed_sub_script]
    non_diff_E_Set = E_Set[non_diff_sub_script]

    #########################################
    # cluster 1
    data_matrix_of_E_Set = exprs(over_expressed_E_Set)

    median_dgl_by_l = apply(data_matrix_of_E_Set, 1, median, na.rm=TRUE)
    delta_1 = log((median(median_dgl_by_l, na.rm=TRUE)^2)/(mad(median_dgl_by_l, na.rm=TRUE)^2))
    theta_1 = log(median(median_dgl_by_l, na.rm=TRUE)/(mad(median_dgl_by_l, na.rm=TRUE)^2))  

    temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)
    lambda_1 = log((median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2))
    nu_1 = log(median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2))
    # end of cluster 1
    #########################################

    #########################################
    # cluster 2
    data_matrix_of_E_Set = exprs(under_expressed_E_Set)

    median_dgl_by_l = -apply(data_matrix_of_E_Set, 1, median, na.rm=TRUE)
    delta_2 = log((median(median_dgl_by_l, na.rm=TRUE)^2)/(mad(median_dgl_by_l, na.rm=TRUE)^2))
    theta_2 = log(median(median_dgl_by_l, na.rm=TRUE)/(mad(median_dgl_by_l, na.rm=TRUE)^2))

    temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)
    lambda_2 = log((median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2))
    nu_2 = log(median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2))
    # end of cluster 2
    #########################################

    #########################################
    # cluster 3
    data_matrix_of_E_Set = exprs(non_diff_E_Set)

    temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)
    lambda_3 = log((median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2))
    nu_3 = log(median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2))
    # end of cluster 3
    #########################################
  }
  else
  {
    median_dgl_by_l = apply(data_matrix_of_E_Set, 1, median, na.rm=TRUE)
    sorted_median_dgl_by_l = sort(median_dgl_by_l)

    median_dgl_by_l = sorted_median_dgl_by_l[ceiling(G*(1 - t_pi_prior[1])):G]

    delta_1 = log((median(median_dgl_by_l, na.rm=TRUE)^2)/(mad(median_dgl_by_l, na.rm=TRUE)^2))
    temp = median(median_dgl_by_l, na.rm=TRUE)/(mad(median_dgl_by_l, na.rm=TRUE)^2)
    if (temp<=0)
    {
      theta_1 = - infinity
    }
    else
    {
      theta_1 = log(temp)
    }

    median_dgl_by_l = -sorted_median_dgl_by_l[1:floor(G*t_pi_prior[2])]
    delta_2 = log((median(median_dgl_by_l, na.rm=TRUE)^2)/(mad(median_dgl_by_l, na.rm=TRUE)^2))
    temp = median(median_dgl_by_l, na.rm=TRUE)/(mad(median_dgl_by_l, na.rm=TRUE)^2)
    if (temp<=0)
    {
      theta_2 = - infinity
    }
    else
    {
      theta_2 = log(temp)
    }

    # median_dgl_by_l = sorted_median_dgl_by_l[(floor(G*t_pi_prior[2])+1):(ceiling(G*(1 - t_pi_prior[1]))-1)]

    temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)
    lambda_1 = log((median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2))
    nu_1 = log(median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2))

    lambda_2 = lambda_1
    nu_2 = nu_1
    lambda_3 = lambda_1
    nu_3 = nu_1
  }    

  data_matrix_of_E_Set = exprs(E_Set)

  psi = c(lambda_1, nu_1, delta_1, theta_1,
          lambda_2, nu_2, delta_2, theta_2,
          lambda_3, nu_3)

  t_pi = t_pi_prior

  tilde_z = get_tilde_z(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l,
    n,
    M)

  #############################################################
  # for test
  if (0)
  {
    focus = 10
    precision = 0.001
    psi[focus] = psi[focus] - precision
    f1 = l_c(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z, 
      b, 
      converge_threshold, 
      infinity,
      M)
    psi[focus] = psi[focus] + precision
    f2 = l_c(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z, 
      b, 
      converge_threshold, 
      infinity,
      M)
    psi[focus] = psi[focus] + precision
    f3 = l_c(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z, 
      b, 
      converge_threshold, 
      infinity,
      M)

    d = (f3 - f1)/(2 * precision)
    print(c(f1,f2,f3))
    cat("numerical>>",d,'\n')

    psi[focus] = psi[focus] - precision
    true_d = gradient_l_c(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z, 
      b, 
      converge_threshold, 
      infinity,
      M)
    
    cat("true>>",true_d[focus],'\n')
    cat("diff>>",d - true_d[focus],'\n')

    stop()
  }
  # end of for test
  #############################################################

  if (verbose)
  {
    cat("starting EM>>\n")
  }

  mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c, 
    t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, sum_dgl_square_by_l = sum_dgl_square_by_l, 
    n = n, tilde_z = tilde_z, method = 'L-BFGS-B', 
    b = b, converge_threshold = converge_threshold, infinity = infinity, M = M,
    lower=param_limit_min, 
    upper=param_limit_max,
    control = list(maxit = max_iteration_num_in_optim))

  t_pi = get_t_pi(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z,
    b,
    converge_threshold,
    infinity)

  repeated_times = 0  

  while (repeated_times<max_repeated_times)
  {
    repeated_times = repeated_times + 1

    psi = mleinfo$par

    if (verbose)
    {
      print(c("repeated times:", repeated_times))
      cat("psi>>\n",psi,'\n')
      cat("t_pi>>\n",t_pi,'\n')
    }
    
    tilde_z = get_tilde_z(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n,
      M)

    last_mleinfo = mleinfo
    mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c, 
      t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, sum_dgl_square_by_l = sum_dgl_square_by_l, 
      n = n, tilde_z = tilde_z, method = 'L-BFGS-B', 
      b = b, converge_threshold = converge_threshold, infinity = infinity, M = M,
      lower=param_limit_min, 
      upper=param_limit_max,
      control = list(maxit = max_iteration_num_in_optim))

    t_pi = get_t_pi(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z,
      b,
      converge_threshold,
      infinity)

    #if (abs(last_mleinfo$value - mleinfo$value)<converge_threshold) break
    if (sum(abs(last_mleinfo$par - mleinfo$par))<converge_threshold) break
  }

  tilde_z = get_tilde_z(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n,
    M)

  fData(E_Set)$est_cluster = apply(tilde_z,1,which.max)
  if (is_sim)
  {
    fData(E_Set)$flag = (fData(E_Set)$est_cluster == fData(E_Set)$true_cluster)
  }

  result = list(
    E_Set = E_Set, 
    mleinfo = mleinfo, 
    t_pi = t_pi,
    repeated_times = repeated_times)
  invisible(result)
}
