import dadi
import numpy

# A two population model (M2) consisting of FOUR parameters: 
# 3 population sizes (ancestral and 2 recent) and 
# a divergence time
def M2_fs((nuA, nu1, nu2, T_ms), (n1,n2), pts):
# Define the grid we'll use
xx = dadi.Numerics.default_grid(pts)

# Convert from units used in ms (scaling by 4 N0) to units used in dadi 
# (scaling by 2 N0)
T = T_ms*2

# Here theta0 is set to nuA because we define the ancestral population size to also be 
# the reference population size.
phi = dadi.PhiManip.phi_1D(xx, theta0 = nuA)

# The divergence
phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

# Steady state during time T
phi = dadi.Integration.two_pops(phi, xx, T, nu1 = nu1, nu2 = nu2, theta0 = nuA)

fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx , xx))
return fs
