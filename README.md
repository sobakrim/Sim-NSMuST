<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>Simulation Algorithm for Nonstationary Multivariate Space-Time Gaussian Random Fields</title>
</head>
<body style="font-family: Arial, sans-serif;">

<hr />

<h2 id="overview">Overview</h2>
<p>
This repository provides R functions to simulate nonstationary Gaussian fields (multivariate and spatio-temporal) through Gaussian mixtures. 
Nonstationary Gaussian fields are a critical area of research in spatial analysis. Traditional simulation methods often rely on Cholesky decomposition, 
which can quickly become computationally prohibitive for moderately sized domains. By contrast, this simulation approach:
</p>
<ul>
  <li>Generates multivariate space-time Gaussian fields using Gneiting–Matérn or Gneiting–Cauchy models.</li>
  <li>Allows the Matérn (or Cauchy) parameters to vary with space and time.</li>
  <li>Incorporates anisotropy matrices that evolve over space and time.</li>
  <li>Uses parallel computing to accelerate the simulation process.</li>
</ul>

<hr />

<h2 id="collaborators">Collaborators</h2>
<ul>
  <li>Denis Allard &ndash; 
  <li>Lionel Benoit &ndash; 
</ul>

<hr />

<h2 id="example-usage">Example Usage</h2>
<p>The following is a more detailed example demonstrating how to set up a simulation environment and perform a basic run:</p>

<pre><code># Load required libraries
library(MASS)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(invgamma)
library(patchwork)

# Source the functions
source("Sim-NSMuST.R")  

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

# Number of variables (e.g., 2 variables/fields)
p = 2 

# Number of waves for spectral simulation
L = 50000

# Define anisotropy functions
aniso_fun1 = function(t_val, s_coord) {
  # Simple identity-based anisotropy (no real change)
  diag(2)
}

aniso_fun2 = function(t_val, s_coord) {
  s1 = s_coord[1]
  s2 = s_coord[2]
  sum_s = s1 + s2
  ratio = 0.9 - 0.5 * (sum_s / 20)  # Variation based on location
  angle = 2 * pi * (t_val - 1) / 5 # Variation based on time
  
  D = matrix(c(1, 0,
               0, ratio),
             nrow = 2, byrow = TRUE)
  
  R = matrix(c(cos(angle),  sin(angle),
               -sin(angle), cos(angle)),
             nrow = 2, byrow = TRUE)
  
  D %*% R
}

# Combine them in a list, one anisotropy function per variable
anisotropies = list(aniso_fun1, aniso_fun2)

# Define Matern covariance function parameters for each variable
matern_fun1 = function(t_val, s_coord) {
  list(nu = 1, r = 1)  # Constant parameters
}
matern_fun2 = function(t_val, s_coord) {
  s1 = s_coord[1]
  s2 = s_coord[2]
  sum_s = s1 + s2
  
  # Slightly more complex function of space & time
  nu_val = 1 + 0.3 * (t_val - 1) * (sum_s / 10)
  r_val = 1
  
  list(nu = nu_val, r = r_val)
}
SpectralPars = list(matern_fun1, matern_fun2)

# Define the cross-covariance (between variables) function
Sigma_fun = function(t_val, s_coord) {
  s1 = s_coord[1]
  # Just a simple function that depends on s1
  cval = 9*s1/100
  out = matrix(c(1, cval,
                 cval, 1), 
               2, 2, byrow = TRUE)
  # Return the Cholesky factor
  return(t(chol(out)))
}

# Define additional simulation parameters
params = list(
  p = p,
  nt = nt,
  t = t_seq,
  # Time-varying length-scale functions for each variable
  l_funcs = list(
    function(t) { 0*t + 0.5 },
    function(t) { 0*t + 0.2 }
  ),
  a = 1/2,
  alpha = 1,
  delta = 0.4,
  b = 0.5,
  
  # Include anisotropy and spectral parameter definitions
  anisotropies = anisotropies,
  SpectralPars = SpectralPars,
  
  # Covariance function
  Sigma_fun = Sigma_fun
)

# Run the simulation
set.seed(14347)  
Z = SimulateParsimNS(
  L                  = L,
  ns                 = ns,
  nt                 = nt,
  p                  = p,
  params             = params,
  spatial_coordinates = spatial_coordinates,
  SpectralDensity     = MaternSpectralDensity3D,
  parallelize         = TRUE, 
  batch_size          = 100,
  n_cores             = 5
)

# Plot two variables over time at selected time steps
time_steps_vec = c(1,2,3,4)
pl = Plot2VarsOverTime(
  Z_sim         = Z,
  time_steps    = time_steps_vec,
  var1          = 1,
  var2          = 2,
  spatial_coords = spatial_coordinates
)

pl
</code></pre>
</body>
</html>
