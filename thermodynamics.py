import numpy as np
from scipy.optimize import newton

unit_atm = 101300.0
latent_heat = 2.25e6
gas_const = 8.3147
R_d   = gas_const / 28.9e-3
R_vap = gas_const / 18.0e-3
C_v = 2.5 * R_d
C_p = R_d + C_v
THETA_vap = latent_heat * R_d / C_p / R_vap

p_ref = unit_atm
p_ground = 101300.0
kappa = R_d / C_p
zeroK = 273.15

# This reference is rather important since Claus-Claperon approximation
# is not that good when it comes to water
vapor_pressure_ref   = 2338.8
vapor_pressure_T_ref = zeroK + 20.0

def saturated_vapor_pressure(T):
	global R_vap, vapor_pressure_ref, vapor_pressure_T_ref
	return vapor_pressure_ref * np.exp(- latent_heat / R_vap * (1.0/T - 1.0 / vapor_pressure_T_ref) )

def inv_saturated_vapor_pressure(p_vap):
	global R_vap, latent_heat, vapor_pressure_ref, vapor_pressure_T_ref
	return 1.0 / ( 1.0 / vapor_pressure_T_ref  - np.log(p_vap / vapor_pressure_ref) * R_vap / latent_heat) if p_vap != 0.0 else np.nan
	
def saturated_vapor_mass(T, p):  # kg / kg
	global R_vap, R_d
	es = saturated_vapor_pressure(T)
	#return ( R_d / R_vap ) * (es / p)
	# analytic sol
	return ( R_d / R_vap ) / (p / es - 1.0)   # analytic

def inv_saturated_vapor_mass(w, p):
	global R_vap, R_d
	#return inv_saturated_vapor_pressure(w * p * R_vap / R_d)

	# analytic sol
	return inv_saturated_vapor_pressure(p / (1.0 + R_d / R_vap / w))

def cal_dew(T, p, RH):
	return inv_saturated_vapor_pressure(RH * saturated_vapor_pressure(T)) 
	



def cal_theta(T, p):
	return T / EXENER(p)

def cal_theta_e(T, p):
	global latent_heat, C_p
	return cal_theta(T, p) + latent_heat * saturated_vapor_mass(T, p) / C_p

def EXENER(p):
	global kappa, p_ref
	return (p / p_ref) ** kappa

def theta_e_helper(theta_e, T, p):
	return theta_e - cal_theta_e(T, p)

#def theta_e_fprime_helper(T, p):
#	global latent_heat, R_vap, C_p
#	return - 1.0 / EXENER(p) - (THETA_vap * latent_heat / R_vap / T**2.0) * saturated_vapor_pressure(T) / p

def theta_e_fprime_helper(T, p):
	global latent_heat, R_vap, C_p
	x = p / saturated_vapor_pressure(T)
	return - 1.0 / EXENER(p) - (THETA_vap * latent_heat * x / (x - 1.0)**2.0 / R_vap / T**2.0)


def solve_T_given_theta_e_and_p(theta_e, p):
	# initial guess is important
	return newton((lambda T: theta_e_helper(theta_e, T, p)), zeroK, fprime=(lambda T: theta_e_fprime_helper(T, p)), tol=1e-2, maxiter=50)
	#return newton((lambda T: theta_e_helper(theta_e, T, p)), zeroK, fprime=None, tol=1e-2, maxiter=500)

