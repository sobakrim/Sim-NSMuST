
source("Sim-NSMuST.R") ## load functions

# Define spatial and temporal resolutions
spatial_resolution = 0.05
temporal_resolution = 1

# Generate spatial grid for simulation
s1_seq = seq(0, 10, by = spatial_resolution)
s2_seq = seq(0, 10, by = spatial_resolution)
grid_spatial = expand.grid(s1 = s1_seq, s2 = s2_seq)
spatial_coordinates = as.matrix(grid_spatial) 
ns = nrow(spatial_coordinates)   

# Define temporal sequence
t_seq = seq(1, 4, by = temporal_resolution)
nt = length(t_seq)  

p = 2 # Number of variables
L = 50000 # Number of waves

# Define anisotropy functions
aniso_fun1 = function(t_val, s_coord) {
  diag(2)
}

aniso_fun2 = function(t_val, s_coord) {
  s1 = s_coord[1]
  s2 = s_coord[2]
  sum_s = s1 + s2
  ratio = 0.9 - 0.5 * (sum_s / 20) 
  angle = 2 * pi * (t_val - 1) / 5
  
  D = matrix(c(1,     0,
               0, ratio),
             nrow = 2, byrow = TRUE)
  
  R = matrix(c(cos(angle),  sin(angle),
               -sin(angle), cos(angle)),
             nrow = 2, byrow = TRUE)
  
  D %*% R
}

anisotropies = list(aniso_fun1, aniso_fun2)

# Define Matern covariance function parameters
matern_fun1 = function(t_val, s_coord) {
  list(nu=1, r=1)
}
matern_fun2 = function(t_val, s_coord) {
  s1 = s_coord[1]
  s2 = s_coord[2]
  sum_s = s1 + s2
  
  nu_val =  1 + 0.3 * (t_val - 1) * (sum_s / 10)
  
  r_val = 1
  
  list(nu = nu_val, r = r_val)
}
SpectralPars = list(matern_fun1, matern_fun2)

# Define covariance structure (between variables) function
Sigma_fun = function(t_val, s_coord) {
  s1 = s_coord[1]
  cval =9*s1/100
  out =  matrix(c(1, cval, cval, 1), 2,2, byrow=TRUE)
  return(t(chol(out)))
}

# Define simulation parameters
params = list(
  p = p,
  nt = nt,
  t = t_seq,
  l_funcs = list(
    function(t) { 0*t + 0.5 },
    function(t) { 0*t + 0.2 }
  ),
  a = 1/2,
  alpha = 1,
  delta = 0.4,
  b = 0.5,
  anisotropies = anisotropies,
  SpectralPars   = SpectralPars,
  Sigma_fun    = Sigma_fun
)

# Run simulation
t1 = Sys.time()
set.seed(14347)
Z = SimulateParsimNS(
  L                  = L,
  ns                 = ns,
  nt                 = nt,
  p                  = p,
  params             = params,
  spatial_coordinates = spatial_coordinates,
  SpectralDensity = MaternSpectralDensity3D,
  parallelize        = T, 
  batch_size = 100,
  n_cores = 5
)
t2 = Sys.time()
t2-t1

# Plot two variables over time
time_steps_vec = c(1,2,3,4)
pl = Plot2VarsOverTime(
  Z_sim = Z,
  time_steps = time_steps_vec,
  var1 = 1,
  var2 = 2,
  spatial_coords = spatial_coordinates
)
pl

ggplot2::ggsave(plot =pl,"~/Desktop/Simulation GM/plots/samples.pdf", width = 28, height = 18, units = "cm")
# Comparison between empirical and theoretical covariance 

# Resample spatial and temporal coordinates
set.seed(7145)
s1_seq = seq(0, 10, by = spatial_resolution)
s2_seq = seq(0, 10, by = spatial_resolution)
s1_seq = sample(s1_seq, 8)
s2_seq = sample(s2_seq, 8)
grid_spatial = expand.grid(s1 = s1_seq, s2 = s2_seq)
spatial_coordinates = as.matrix(grid_spatial) 
ns = nrow(spatial_coordinates) 

# Define new parameters for resampled data
params$nt = nt
params$t = t_seq


# Generate multiple realizations
num_realizations = 1000
Z_list = mclapply(1:num_realizations, mc.cores = parallel::detectCores(), function(i) SimulateParsimNS(
  L                  = L,
  ns                 = ns,
  nt                 = nt,
  p                  = p,
  params             = params,
  spatial_coordinates = spatial_coordinates,
  parallelize        = F, batch_size = L
))
Z_array = array(unlist(Z_list), dim = c(nt, ns, p, num_realizations))

# Compute empirical and theoretical correlations per pairs of space, time, and variables 

vars_df   = expand.grid(i = 1:p, j = 1:p) %>% subset(i <= j)
spaces_df = expand.grid(s1 = 1:ns, s2 = 1:ns) %>% subset(s1 <= s2)
times_df  = expand.grid(t1 = 1:nt, t2 = 1:nt) %>% subset(t1 <= t2)

selected_vars = lapply(seq_len(nrow(vars_df)), function(r) {
  c(vars_df$i[r], vars_df$j[r])
})

selected_spaces = lapply(seq_len(nrow(spaces_df)), function(r) {
  c(spaces_df$s1[r], spaces_df$s2[r])
})

selected_times = lapply(seq_len(nrow(times_df)), function(r) {
  c(times_df$t1[r], times_df$t2[r])
})


results = data.frame(
  i = integer(),
  j = integer(),
  s1 = integer(),
  s2 = integer(),
  t1 = integer(),
  t2 = integer(),
  empirical_corr = numeric(),
  theoretical_corr = numeric(),
  stringsAsFactors = FALSE
)

for (var_pair in selected_vars) {
  i = var_pair[1]
  j = var_pair[2]
  
  for (space_pair in selected_spaces) {
    s1 = spatial_coordinates[space_pair[1],]
    s2 = spatial_coordinates[space_pair[2],]
    
    for (time_pair in selected_times) {
      t1 = t_seq[time_pair[1]]
      t2 = t_seq[time_pair[2]]
      
      # sigma_ij_x1,x2
      sigma = params$Sigma_fun(t1, s1)  %*% t(params$Sigma_fun(t2, s2))
      sigma_ij_x1x2 = sigma[i,j]
      
      # Space, time, variable pair 
      Z_i_t1_s1 = Z_array[time_pair[1], space_pair[1], i, ]  
      Z_j_t2_s2 = Z_array[time_pair[2], space_pair[2], j, ]  
      
      # Empirical covariance 
      empirical_corr = cor(Z_i_t1_s1, Z_j_t2_s2)
      
      # Theoretical covariance 
    
      theoretical_corr = Gneiting_Matern_correlation(i, j, s1, s2, t1, t2, params,sigma_ij_x1x2, spatial_coordinates)
      results = rbind(results, data.frame(
        i = i,
        j = j,
        s1 = space_pair[1],
        s2 = space_pair[2],
        t1 = t1,
        t2 = t2,
        empirical_corr = empirical_corr,
        theoretical_corr = theoretical_corr
      ))
    }
  }
 
}

# Scatter plot comparing empirical and theoretical correlations (by pairs of variables)
pt = ggplot(results, aes(x = empirical_corr, y = theoretical_corr)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  coord_fixed() +
  facet_wrap(~ paste("i=", i, ", j=", j)) +  
  theme_bw() +
  labs(
    title = "",
    x = "Empirical correlation",
    y = "Theoretical correlation"
  ) + theme(plot.title = element_text(hjust = 0.5))

pt
ggplot2::ggsave(plot =pt,"~/Desktop/Simulation GM/plots/correlations.pdf", width = 28, height = 10, units = "cm")
