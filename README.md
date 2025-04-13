<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>Simulation algorithm for non-stationsanry multivariate space-time Gaussian Random Fields - README</title>
</head>
<body style="font-family: Arial, sans-serif;">

<p>
This repository contains R functions for simulating nonstationary Gaussian fields in space and time, using Gaussian mixtures. 
</p>

<hr />

<h2 id="overview">Overview</h2>
<p>
Simulating <strong>nonstationary</strong> (spatio-temporal) Gaussian fields is a key topic in geostatistics and spatial
analysis. Traditional stationary approaches may not capture local anisotropies or varying smoothness parameters. This
repository demonstrates how to:
</p>
<ul>
  <li>Combine time-varying length-scale (or range) functions.</li>
  <li>Introduce anisotropy matrices that change over time and space.</li>
  <li>Use different spectral densities (e.g., Matern, Cauchy) for flexible modeling.</li>
  <li>Efficiently generate Gaussian increments consistent with these nonstationary assumptions.</li>
</ul>

<hr />

<h2 id="example-usage">Example Usage</h2>
<p>Below is a toy example showing how to set up a minimal environment and run a simulation:</p>

<pre><code># Load required libraries
library(MASS)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(invgamma)
library(patchwork)

# Source the nonstationary Gaussian field functions
source("Sim-NSMuST.R")

# Example length-scale function for a single field
example_length_scale &lt;- function(t) {
  # Suppose it grows linearly over time
  return(0.5 + 0.1 * t)
}

# Example anisotropy function (creates a simple scale matrix)
example_aniso_func &lt;- function(t, s) {
  r1 &lt;- 1 + 0.01 * t
  r2 &lt;- 1
  theta &lt;- 0.2 * t
  create_anisotropy_matrix(r1, r2, theta)
}

# Example Matern parameters
example_matern_func &lt;- function(t, s) {
  list(nu = 1.0 + 0.05 * t, r = 0.3 + 0.01 * t)
}

# Spatial coordinates (fake grid)
nx &lt;- 10
ny &lt;- 10
x_coords &lt;- seq(0, 1, length.out = nx)
y_coords &lt;- seq(0, 1, length.out = ny)
grid &lt;- expand.grid(x_coords, y_coords)

# Setup params
params &lt;- list(
  p = 1,
  nt = 10,
  t = seq(0, 9, length.out = 10),
  l_funcs = list(example_length_scale),
  a = 0.2,
  alpha = 1.0,
  b = 1,
  anisotropies = list(example_aniso_func),
  MaternPars = list(example_matern_func),
  Sigma_fun = function(time, coord) {
    # Return a simple 2x2 covariance matrix
    matrix(c(1, 0, 0, 1), ncol = 2)
  }
)

# Run simulation
result &lt;- SimulateParsimNS(
  L = 50,                  # number of Monte Carlo iterations
  ns = nrow(grid),
  nt = 10,
  p = 1,
  params = params,
  spatial_coordinates = grid
)

# Plotting the result for time steps 1 to 4
Plot2VarsOverTime(
  Z_sim = result,
  time_steps = 1:4,
  var
