import dadi
import numpy

# A two population model (M6) consisting of SEVEN parameters: 
# 3 population sizes (ancestral and 2 recent), 
# the time since divergence to the admixture event (T1), 
# the time since the latter to the present (T2) and 
# 2 parameters representing bi-directional admixture rates
def M6_fs((nuA, nu1, nu2, T1_ms, f21, f12, T2_ms), (n1,n2), pts):
# Define the grid we'll use
xx = dadi.Numerics.default_grid(pts)

# Convert from units used in ms (scaling by 4 N0) to units used in dadi 
# (scaling by 2 N0)
T1 = T1_ms*2
T2 = T2_ms*2

# Here theta0 is set to nuA because we define the ancestral population size to also be
# the reference population size.
phi = dadi.PhiManip.phi_1D(xx, theta0 = nuA)

# The divergence
phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

# Steady state during time T1
phi = dadi.Integration.two_pops(phi, xx, T1, nu1 = nu1, nu2 = nu2, theta0 = nuA)

# The "pulse" admixture event
phi = dadi.PhiManip.phi_2D_admix_2_into_1(phi, f21, xx, xx)
phi = dadi.PhiManip.phi_2D_admix_1_into_2(phi, f12, xx, xx)

# Steady state during time T2
phi = dadi.Integration.two_pops(phi, xx, T2, nu1 = nu1, nu2 = nu2, theta0 = nuA)

fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx , xx))
return fs
