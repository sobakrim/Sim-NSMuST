library(MASS)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(invgamma)
library(patchwork)
library(parallel)


############################ Main function #########################

SimulateParsimNS = function(L, ns, nt, p, params, spatial_coordinates, SpectralDensity, 
                            parallelize = FALSE, n_cores = 2, batch_size = 100) 
{
  # Simulates a nonstationary Multivariate space-time Gaussian field with a user-specified SpectralDensity function.
  #
  # Arguments:
  #   L: Integer number of Monte Carlo or spectral repetitions.
  #   ns: Integer number of spatial locations.
  #   nt: Integer number of time points.
  #   p: Integer number of variables.
  #   params: List containing model parameters and functions (length-scale, anisotropy, etc.).
  #   spatial_coordinates: Matrix or data frame with spatial coordinates.
  #   SpectralDensity: A function defining the spectral density to be used (e.g., Matern).
  #   parallelize: Logical flag indicating whether to parallelize the simulation.
  #   n_cores: Integer specifying the number of cores to use if parallelize = TRUE.
  #   batch_size: Integer specifying how many iterations per parallel batch.
  #
  # Returns:
  #   A 3D array of dimensions (nt, ns, p) containing the simulated field.
  
  Sig1 = SimulateGauIncrements(params)  
  aniso_precomp_list = lapply(
    params$anisotropies,
    precompute_aniso,         
    t_vals = params$t, 
    spatial_coordinates = spatial_coordinates
  )
  
  # Create aniso_5D: shape (p, nt, ns, 2, 2)
  aniso_5D = array(0, dim = c(p, nt, ns, 2, 2))
  det_array_all = array(0, dim = c(p, nt, ns))
  
  for (i_ in seq_len(p)) {
    aniso_list_t = aniso_precomp_list[[i_]]
    for (t_idx in seq_len(nt)) {
      for (s_idx in seq_len(ns)) {
        A = aniso_list_t[[t_idx]][[s_idx]]      # (2 x 2) matrix
        aniso_5D[i_, t_idx, s_idx, , ] = A
        det_array_all[i_, t_idx, s_idx] = det(A)
      }
    }
  }
  
  rm(aniso_precomp_list)
  gc()
  #### Matern param arrays: shape (p, nt, ns)
  matern_precomp_list = lapply(
    params$SpectralPars,
    precompute_matern_pars,
    t_vals = params$t,
    spatial_coordinates = spatial_coordinates
  )
  nu_3D = array(0, dim = c(p, nt, ns))
  r_3D  = array(0, dim = c(p, nt, ns))
  for(i_ in seq_len(p)) {
    nu_3D[i_, , ] = matern_precomp_list[[i_]]$nu
    r_3D[i_,  , ] = matern_precomp_list[[i_]]$r
  }
  rm(matern_precomp_list)
  gc()
  #### SigmaFun -> 4D array (nt, ns, 2,2) for p=2 dimension
  SigmaFun_precomp = precompute_sigma_fun(params$Sigma_fun, params$t, spatial_coordinates)
  
  SigmaFun_4D = array(0, dim = c(nt, ns, 2, 2))
  for(t_idx in seq_len(nt)) {
    for(s_idx in seq_len(ns)) {
      SigmaFun_4D[t_idx, s_idx, , ] = SigmaFun_precomp[[t_idx]][[s_idx]]
    }
  }
  rm(SigmaFun_precomp)
  gc()
  
  d = dim(spatial_coordinates)[2]
  n = nt * ns
  
  A_flat_all = array(NA_real_, dim = c(p, n * d, d))
  
  for (i_ in seq_len(p)) {
    A_mat = array(aniso_5D[i_, , , , ], dim = c(n, d, d))             
    A_flat_all[i_, , ] = matrix(aperm(A_mat, c(2, 1, 3)), nrow = n * d, ncol = d)
  }
  
  cst = (2*sqrt(pi))^(-2)
  #### single_iteration function
  single_iteration = function(iter) {
    
    Wl  = matrix(Sig1 %*% rnorm(nt * p), nrow = nt, ncol = p, byrow = F)
    
    Xil   = SimulateXi(1)         
    Vl    = rnorm(d)   
    phil  = runif(1, 0, 2*pi)      
    Al    = rnorm(p)    
    
    t_Vl_Vl   = sqrt(sum(Vl^2))
    temp_Wl   = (1 / sqrt(2)) * t_Vl_Vl * Wl 
    omegal    = sqrt(2 * Xil) * Vl   
    
    pdfXi_Xil = pdfXi(Xil)
    x = sqrt(2) * Vl
    pdfi      = GaussianSpectralDensity(x, Sigma =  diag(p)) 
    
    localZ = array(0, dim = c(nt, ns, p))
    temp   = spatial_coordinates %*% omegal  
    
    
    for(i_ in seq_len(p)) {
      
      spec_raw = SpectralDensity(Xil, nu_3D[i_, , ], r_3D[i_, , ])
      localSpec_mat = spec_raw / pdfXi_Xil
      
      x1 = x[1]; x2 =x[2]
      quad_form = x1^2*aniso_5D[i_, , , , ] [,,1,1] +
        2*x1*x2*aniso_5D[i_, , , , ] [,,1,2] +
        x2^2*aniso_5D[i_, , , , ] [,,2,2]
      
      det_array = matrix(det_array_all[i_, , ], nrow = nt, ncol = ns)
      
      gauss_mat = cst * sqrt(det_array) * exp(-quad_form/4) / pdfi[1]
      gauss_local_mat = sqrt(gauss_mat)
      
      scale_local_mat = matrix(crossprod(Al, matrix(aperm(SigmaFun_4D[, , i_, ]  , c(3, 1, 2)), nrow = p))  , nrow = nt, ncol = ns)
      scale_local_mat = scale_local_mat * sqrt(localSpec_mat)
      
      time_vec = temp_Wl[, i_]  
      cosarg_mat = time_vec + matrix(temp[,1], nrow = nt, ncol = ns, byrow = TRUE) + phil
      cos_mat = cos(cosarg_mat) 
      
      inc_mat = scale_local_mat * gauss_local_mat * cos_mat 
      localZ[,, i_] = inc_mat
    }
    
    localZ
  } 
  
  simulation_indices = seq_len(L)
  if(!parallelize) {
    Z_accum = array(0, dim = c(nt, ns, p))
    for (ell in seq_len(L)) {
      iter = single_iteration(ell)
      if(all(!is.finite(iter))){
        iter = 0
        L = L - 1
      }
      Z_accum = Z_accum + iter
    }
    Z = sqrt(2 / L) * Z_accum
  } else {
    library(parallel)
    num_batches = ceiling(L / batch_size)
    Z_accum = array(0, dim = c(nt, ns, p))
    valid_count = 0  # Track the number of valid simulations
    
    for (batch in seq_len(num_batches)) {
      batch_indices = seq((batch - 1) * batch_size + 1, min(batch * batch_size, L))
      Z_list = mclapply(batch_indices, single_iteration, mc.cores = n_cores)
      
      # Filter invalid simulations
      is_valid = sapply(Z_list, function(x) all(is.finite(x)))
      num_valid = sum(is_valid)
      num_invalid = length(batch_indices) - num_valid
      
      if (num_invalid > 0) {
        warning(sprintf("Removed %d invalid simulations in batch %d.", num_invalid, batch))
      }
      
      Z_list_valid = Z_list[is_valid]
      
      # Accumulate only valid results
      for (ell in seq_len(num_valid)) {
        Z_accum = Z_accum + Z_list_valid[[ell]]
      }
      
      valid_count = valid_count + num_valid
    }
    
    Z = sqrt(2 / valid_count) * Z_accum
  }
  return(Z)
}


############################ Helper functions #########################


C_M = function(x, kappa, nu) {
  # Computes the Matern covariance values for a vector of distances x, using kappa and nu.
  #
  # Arguments:
  #   x: Numeric vector of distances between points.
  #   kappa: Scalar range parameter (often noted as 'r' in Matern contexts).
  #   nu: Scalar smoothness parameter of the Matern function.
  #
  # Returns:
  #   Numeric vector of Matern covariance values for the input distances.
  
  matern_cov = ifelse(x == 0, 1, 
                      (2^(1 - nu) / gamma(nu)) * 
                        (kappa * abs(x))^nu * besselK(kappa * abs(x), nu))
  
  matern_cov[is.nan(matern_cov)] = 1
  
  return(matern_cov)
}


Gneiting_Matern_correlation = function(i, j, s1, s2, t1, t2, params, sigma_ij_x1x2, spatial_coordinates) {
  # Computes the Gneiting-Matern correlation between two fields i and j at locations s1, s2 and times t1, t2.
  #
  # Arguments:
  #   i: Integer specifying the index of the first field (variable).
  #   j: Integer specifying the index of the second field (variable).
  #   s1: Numeric vector of spatial coordinates (e.g., c(x, y)) for the first location.
  #   s2: Numeric vector of spatial coordinates for the second location.
  #   t1: Numeric scalar for the first time point.
  #   t2: Numeric scalar for the second time point.
  #   params: List containing relevant model parameters and functions (anisotropy, SpectralPars, etc.).
  #   sigma_ij_x1x2: Not used in the body here but presumably a placeholder for cross-covariance terms.
  #   spatial_coordinates: Complete matrix of all spatial coordinates (not fully used here).
  #
  # Returns:
  #   A numeric scalar representing the Gneiting-Matern correlation between fields i and j at locations s1, s2 and times t1, t2.
  
  Sigma_ii_x1 = solve(params$anisotropies[[i]](t1,s1))
  Sigma_jj_x2 = solve(params$anisotropies[[j]](t2,s2))
  Sigma_ij_x1x2 = ((Sigma_ii_x1 %*% t(Sigma_ii_x1)) + (Sigma_jj_x2 %*% t(Sigma_jj_x2))) / 2
  
  det_Sigma_ii_x1 = det(Sigma_ii_x1)^2
  det_Sigma_jj_x2  = det(Sigma_jj_x2)^2
  
  nu_ii_x1 = params$SpectralPars[[i]](t1,s1)$nu
  nu_jj_x2 = params$SpectralPars[[j]](t2,s2)$nu
  nu_ij_x1x2 = (nu_ii_x1 + nu_jj_x2) / 2
  
  kappa_ii_x1 = params$SpectralPars[[i]](t1,s1)$r
  kappa_jj_x2 = params$SpectralPars[[j]](t2,s2)$r
  kappa_ij_x1x2 = sqrt((kappa_ii_x1^2 + kappa_jj_x2^2)/2)
  
  a_param = params$a
  alpha_param = params$alpha
  gamma0 = function(h, a_param, alpha_param) {
    return(a_param * abs(h)^alpha_param)
  }
  
  gamma0_t1_t2 = gamma0(t1-t2, a_param, alpha_param)
  
  r_ii_x1 = params$l_funcs[[i]](t1)
  r_jj_x2 = params$l_funcs[[j]](t2)
  
  R_ij =  r_jj_x2 * r_ii_x1 * exp(- (t1-t2)^2 /2)
  Rii = (r_jj_x2^2 + r_ii_x1^2)/2
  gamma_ij = (gamma0_t1_t2 - R_ij + Rii)
  
  Sigma_ij_x1x2 = Sigma_ij_x1x2 + diag(rep(gamma_ij, ncol(Sigma_ij_x1x2)))
  
  spatial_diff = as.numeric(s1 - s2)
  
  distance = sqrt(t(spatial_diff) %*% solve(Sigma_ij_x1x2) %*% spatial_diff)
  
  CM_value = C_M(distance, kappa_ij_x1x2, nu_ij_x1x2)
  
  C_ij = (det_Sigma_ii_x1^(1/4) * det_Sigma_jj_x2^(1/4)) / det(Sigma_ij_x1x2)^(1/2) *
    (gamma(nu_ij_x1x2) / sqrt(gamma(nu_ii_x1) * gamma(nu_jj_x2))) * sigma_ij_x1x2 *
    (kappa_ii_x1^nu_ii_x1 * kappa_jj_x2^nu_jj_x2) / (kappa_ij_x1x2^(2 * nu_ij_x1x2)) *
    CM_value
  
  return(C_ij)
}


GaussianSpectralDensity = function(x, Sigma){
  # Computes the Gaussian spectral density for input vector x and covariance matrix Sigma.
  #
  # Arguments:
  #   x: Numeric vector at which to evaluate the density.
  #   Sigma: Covariance matrix for the Gaussian spectral density computation.
  #
  # Returns:
  #   Numeric scalar representing the Gaussian spectral density at x.
  
  d = ncol(Sigma)
  return((2*sqrt(pi))^(-d) * det(Sigma)^(1/2) * exp(- t(x) %*% Sigma %*% x / 4))
}


create_anisotropy_matrix = function(r1, r2, theta) {
  # Creates a 2x2 anisotropy matrix, scaled by (r1, r2) and rotated by theta.
  #
  # Arguments:
  #   r1: Numeric scalar controlling anisotropy in one direction.
  #   r2: Numeric scalar controlling anisotropy in the orthogonal direction.
  #   theta: Angle (in radians) to rotate the anisotropy matrix.
  #
  # Returns:
  #   A 2x2 matrix encoding the specified anisotropy and rotation.
  
  R_matrix = matrix(c(r1, 0, 0, r2), nrow = 2, byrow = TRUE)
  T_matrix = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, byrow = TRUE)
  A_matrix = R_matrix %*% T_matrix
  
  return(A_matrix)
}


MaternSpectralDensity = function(xi, params){
  # Evaluates the Matern spectral density at xi given parameters nu and r.
  #
  # Arguments:
  #   xi: Numeric scalar or vector where the Matern spectral density is evaluated.
  #   params: List containing 'nu' (smoothness) and 'r' (range) parameters for the Matern function.
  #
  # Returns:
  #   Numeric value(s) of the Matern spectral density evaluated at xi.
  
  nu = params$nu
  r = params$r
  result =  (1 / gamma(nu)) * (r / 2)^(2 * nu) *
    xi^(-nu - 1) * exp(-r^2 / (4 * xi))
  return(result)
}


CauchySpectralDensity = function(xi, params){
  # Evaluates the Cauchy spectral density at xi given parameters nu and r.
  #
  # Arguments:
  #   xi: Numeric scalar or vector where the Cauchy spectral density is evaluated.
  #   params: List containing 'nu' (related to shape) and 'r' (scale/range).
  #
  # Returns:
  #   Numeric value(s) of the Cauchy spectral density evaluated at xi.
  
  nu = params$nu
  r = params$r
  result = r^(-nu) * gamma(nu)^(-1) * xi^(nu-1) * exp(-xi/r)
  return(result)
}


SimulateGauIncrements = function(params) {
  # Constructs the Cholesky factor of the covariance matrix for Gaussian increments.
  #
  # Arguments:
  #   params: A list containing the fields configuration (p), time points (n, t_seq),
  #           length-scale functions (l_funcs), variogram parameters (a, alpha, b), etc.
  #
  # Returns:
  #   The Cholesky factor of the covariance matrix for Gaussian increments.
  
  p = params$p                # Number of fields
  n = params$nt                # Number of time points
  t_seq = params$t                # Vector of time points
  l_funcs = params$l_funcs    # List of length scale functions l_i(t)
  a = params$a          # Variogram parameter a
  alpha = params$alpha  # Variogram parameter 'alpha' in [0,2]
  b = params$b
  # Define gamma0(h)
  gamma0 = function(h, a_param, alpha_param) {
    return(a_param * abs(h)^alpha_param)
  }
  
  # Initialize Correlation matrix
  total_points = p * n
  Sigma = matrix(0, nrow = total_points, ncol = total_points)
  
  # Precompute gamma0(t) for all time points
  gamma0_t = gamma0(t_seq, a, alpha)
  
  # Compute gamma0(t1 - t2) for all pairs of time points
  delta_t = outer(t_seq, t_seq, "-")
  delta_t2 = delta_t^2
  gamma0_t1_t2 = gamma0(delta_t, a, alpha)
  
  # Precompute l_i(t) and r_ii(t, t)
  l_t_list = vector("list", p)
  r_ii_t_t_list = vector("list", p)
  for (i in 1:p) {
    l_t_list[[i]] = l_funcs[[i]](t_seq)
    r_ii_t_t_list[[i]] = l_t_list[[i]]
  }
  
  # Compute R_ii(t1, 0) and R_jj(0, t2)
  R_ii_t1_0_list = vector("list", p)
  R_jj_0_t2_list = vector("list", p)
  for (i in 1:p) {
    # Evaluate length scale functions at time t and at time 0
    l_i_t1 = l_t_list[[i]]           # Length scales at all time points for field i
    l_1_0 = l_funcs[[i]](0)       # Length scale at time 0 for field i
    
    r_ii_t1_t1 = l_i_t1            # r_ii(t1, t1)
    r_11_0_0 = l_1_0              # r_ii(0, 0)
    
    # Compute delta_t1_0 = t - 0 
    delta_t1_0 = t_seq
    
    # Compute R_ii(t1, 0)
    R_ii_t1_0 =  r_ii_t1_t1 * r_11_0_0 * exp(- delta_t1_0^2 /2)
    Rii = (r_11_0_0^2 + r_ii_t1_t1^2)/2
    R_ii_t1_0_list[[i]] = gamma0_t -   R_ii_t1_0 + Rii
  }
  for (i in 1:p) {
    # Evaluate length scale functions at time t and at time 0
    l_i_t1 = l_t_list[[i]]           # Length scales at all time points for field i
    l_1_0 = l_funcs[[i]](0)       # Length scale at time 0 for field i
    
    r_ii_t1_t1 = l_i_t1            # r_ii(t1, t1)
    r_11_0_0 = l_1_0              # r_ii(0, 0)
    
    # Compute delta_t1_0 = t - 0 (which is just t)
    delta_t1_0 = t_seq
    
    # Compute R_ii(t1, 0)
    R_ii_t1_0 =  r_ii_t1_t1 * r_11_0_0 * exp(- delta_t1_0^2 /2)
    Rii = (r_11_0_0^2 + r_ii_t1_t1^2)/2
    
    R_jj_0_t2_list[[i]] = gamma0_t -   R_ii_t1_0 + Rii
  }
  for (i in 1:p) {
    for (j in 1:p) {
      # Indices for placing values in Sigma
      rows = ((i - 1) * n + 1):(i * n)
      cols = ((j - 1) * n + 1):(j * n)
      
      l_i_t = l_t_list[[i]]
      l_j_t = l_t_list[[j]]
      
      # Compute r_ii(t1, t1) and r_jj(t2, t2)
      r_ii_t1_t1 = r_ii_t_t_list[[i]]  
      r_jj_t2_t2 = r_ii_t_t_list[[j]]  
      
      # Compute R_ij(t1, t2)
      R_ij =  r_jj_t2_t2 * r_ii_t1_t1 * exp(- delta_t2 /2)
      
      
      gamma0_t1 = matrix(R_ii_t1_0_list[[i]], nrow = n, ncol = n, byrow = F)
      gamma0_t2 = matrix(R_jj_0_t2_list[[j]], nrow = n, ncol = n, byrow = T)
      
      Rii = (r_jj_t2_t2^2 + r_ii_t1_t1^2)/2
      
      Cov = gamma0_t1 + gamma0_t2 - (gamma0_t1_t2 - R_ij + Rii)
      
      Sigma[rows, cols] = Cov
    }
  }
  chS = try(chol(Sigma), silent = T)
  if(is.character(chS)) chS = chol(Sigma + diag(1e-6, nrow = ncol(Sigma)))
  
  return(t(chS))
}


pdfXi = function(x){
  # Computes the inverse-gamma PDF (shape=1) for a numeric scalar/vector x.
  #
  # Arguments:
  #   x: Numeric vector (or scalar) at which the inverse-gamma density is evaluated.
  #
  # Returns:
  #   The value(s) of the inverse-gamma PDF at x (shape=1).
  
  return(dinvgamma(x, shape = 1))
}


SimulateXi = function(n){
  # Samples from the inverse-gamma distribution (shape=1).
  #
  # Arguments:
  #   n: Integer specifying how many samples to draw from the inverse-gamma distribution.
  #
  # Returns:
  #   A numeric vector of length n sampled from the inverse-gamma distribution (shape=1).
  
  return(rinvgamma(n, shape = 1))
}


precompute_aniso = function(aniso_fun, t_vals, spatial_coordinates) {
  # Precomputes (A %*% t(A)) where A = solve(aniso_fun(t, s)) for all t, s.
  #
  # Arguments:
  #   aniso_fun: A function that gives the anisotropy matrix as a function of time and space.
  #   t_vals: Numeric vector of time points.
  #   spatial_coordinates: Matrix or data frame of spatial coordinates.
  #
  # Returns:
  #   A list of length nt, where each element is another list (length ns) of A_matrices.
  
  nt = length(t_vals)
  ns = nrow(spatial_coordinates)
  A_big_list = vector("list", nt)
  for(t_idx in seq_len(nt)) {
    t_val = t_vals[t_idx]
    A_list_t = vector("list", ns)
    for(s_idx in seq_len(ns)) {
      s_coord = spatial_coordinates[s_idx, ]
      ans     = solve(aniso_fun(t_val, s_coord))  # p x p
      A_mat = ans %*% t(ans)
      A_list_t[[s_idx]] = A_mat
    }
    A_big_list[[t_idx]] = A_list_t
  }
  return(A_big_list)
}


precompute_matern_pars = function(matern_fun, t_vals, spatial_coordinates) {
  # Precomputes Matern parameters (nu and r) for each time and space location.
  #
  # Arguments:
  #   matern_fun: A function that returns a list of Matern parameters (e.g., nu, r) for each time/space.
  #   t_vals: Numeric vector of time points.
  #   spatial_coordinates: Matrix or data frame of spatial coordinates.
  #
  # Returns:
  #   A list with two matrices, nu_mat and r_mat, each of size (nt x ns).
  
  nt = length(t_vals)
  ns = nrow(spatial_coordinates)
  nu_mat = matrix(0, nt, ns)
  r_mat  = matrix(0, nt, ns)
  for(t_idx in seq_len(nt)) {
    for(s_idx in seq_len(ns)) {
      pars = matern_fun(t_vals[t_idx], spatial_coordinates[s_idx, ])
      nu_mat[t_idx, s_idx] = pars$nu
      r_mat[t_idx, s_idx]  = pars$r
    }
  }
  return(list(nu=nu_mat, r=r_mat))
}


precompute_sigma_fun = function(sigma_fun, t_vals, spatial_coordinates) {
  # Precomputes the output of sigma_fun for each time and space location.
  #
  # Arguments:
  #   sigma_fun: A function returning a covariance (or similar) matrix given time and location.
  #   t_vals: Numeric vector of time points.
  #   spatial_coordinates: Matrix or data frame of spatial coordinates.
  #
  # Returns:
  #   A list of length nt, where each element is a list of length ns, containing the sigma_fun output.
  
  nt = length(t_vals)
  ns = nrow(spatial_coordinates)
  S_big_list = vector("list", nt)
  for(t_idx in seq_len(nt)) {
    t_val = t_vals[t_idx]
    S_list_t = vector("list", ns)
    for(s_idx in seq_len(ns)) {
      s_coord = spatial_coordinates[s_idx, ]
      S_list_t[[s_idx]] = sigma_fun(t_val, s_coord)
    }
    S_big_list[[t_idx]] = S_list_t
  }
  return(S_big_list)
}


MaternSpectralDensity3D = function(Xi, nu_3D, r_3D) {
  # Computes the 3D Matern spectral density for Xi using arrays nu_3D and r_3D.
  #
  # Arguments:
  #   Xi: Numeric scalar or vector where the 3D Matern spectral density is evaluated.
  #   nu_3D: Numeric (or matrix/array) specifying the smoothness parameter in 3D.
  #   r_3D: Numeric (or matrix/array) specifying the range parameter in 3D.
  #
  # Returns:
  #   Numeric value(s) for the Matern spectral density in 3D evaluated at Xi.
  
  term1 = 1 / gamma(nu_3D)
  term2 = (r_3D / 2)^(2 * nu_3D)
  term3 = Xi^(-nu_3D - 1)
  term4 = exp(-(r_3D^2)/(4 * Xi))
  
  localSpec_mat = term1 * term2 * term3 * term4
  return(localSpec_mat)
}


MaternSpectralDensity_vec = function(xi, nu_vec, r_vec) {
  # Computes the Matern spectral density in a vectorized manner for multiple locations.
  #
  # Arguments:
  #   xi: Numeric scalar for which the Matern spectral density is evaluated.
  #   nu_vec: Numeric vector of smoothness parameters, each index corresponding to a location.
  #   r_vec: Numeric vector of range parameters, each index corresponding to a location.
  #
  # Returns:
  #   A numeric vector of the Matern spectral density evaluated at xi for each location.
  
  out = numeric(length(nu_vec))
  for(i in seq_along(nu_vec)){
    nu_val = nu_vec[i]
    r_val  = r_vec[i]
    pf = (1 / gamma(nu_val)) * ((r_val/2)^(2*nu_val))
    out[i] = pf * xi^(-nu_val - 1) * exp(- (r_val^2)/(4*xi))
  }
  out
}

############################ Plot functions #########################

Plot2VarsOverTime = function(Z_sim, 
                             time_steps, 
                             var1 = 1, 
                             var2 = 2, 
                             spatial_coords,
                             title_prefix = "Variable",
                             color_scale_fun = function(...) scale_fill_viridis_c(option = "H", ...),
                             theme_fun = theme_bw) 
{
  # Plots of two variables from a space-time simulation.
  #
  # Arguments:
  #   Z_sim: 3D array of simulated data, where dimensions typically represent time, space, and variables.
  #   time_steps: Numeric vector specifying which time indices to plot.
  #   var1: Integer index specifying the first variable to plot.
  #   var2: Integer index specifying the second variable to plot.
  #   spatial_coords: Matrix or data frame of spatial coordinates (e.g., x and y).
  #   title_prefix: Character string used as a title prefix for the plots.
  #   color_scale_fun: A function to define the color scale in the plots.
  #   theme_fun: A function to define the ggplot theme.
  #
  # Returns:
  #   A patchwork/ggplot object combining plots for the specified variables over the given time steps.
  
  x_coords = spatial_coords[, 1]
  y_coords = spatial_coords[, 2]
  
  df_list = list()
  for (ts_val in time_steps) {
    time_label = paste("Time", ts_val)
    
    df_v1 = data.frame(
      Z        = Z_sim[ts_val, , var1],
      X        = x_coords,
      Y        = y_coords,
      Variable = paste(title_prefix, var1),
      Time     = time_label
    )
    
    df_v2 = data.frame(
      Z        = Z_sim[ts_val, , var2],
      X        = x_coords,
      Y        = y_coords,
      Variable = paste(title_prefix, var2),
      Time     = time_label
    )
    df_list[[length(df_list) + 1]] = df_v1
    df_list[[length(df_list) + 1]] = df_v2
  }
  df = bind_rows(df_list)
  
  df = df %>%
    mutate(Time = factor(Time, levels = paste("Time", time_steps)))
  
  df_var1 = df %>% filter(Variable == paste(title_prefix, var1))
  df_var2 = df %>% filter(Variable == paste(title_prefix, var2))
  
  p_var1 = ggplot(df_var1, aes(x = X, y = Y, fill = Z)) +
    geom_raster() +
    color_scale_fun() +
    facet_wrap(~ Time, nrow = 1) +
    theme_fun() +
    labs(
      title = paste(title_prefix, var1),
      x = "",
      y = "y coordinate",
      fill = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      strip.text = element_text(size = 10)
    )
  
  p_var2 = ggplot(df_var2, aes(x = X, y = Y, fill = Z)) +
    geom_raster() +
    color_scale_fun() +
    facet_wrap(~ Time, nrow = 1) +
    theme_fun() +
    labs(
      title = paste(title_prefix, var2),
      x = "x coordinate",
      y = "y coordinate",
      fill = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      strip.text = element_text(size = 10)
    )
  
  p_var1 / p_var2
}
