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

def saturated_vapor_pressure(T):
	global R_vap
	return unit_atm * np.exp(- latent_heat / R_vap * (1.0/T - 1.0 / (zeroK + 100.0)) )

def saturated_water_mass(T, p):
	global R_vap, R_d
	es = saturated_vapor_pressure(T)
	return ( R_d / R_vap ) * (es / p)

def cal_theta(T, p):
	return T / EXENER(p)

def cal_theta_e(T, p):
	global latent_heat, C_p
	return cal_theta(T, p) + latent_heat * saturated_water_mass(T, p) / C_p

def EXENER(p):
	global kappa, p_ref
	return (p / p_ref) ** kappa

def theta_e_helper(theta_e, T, p):
	return theta_e - cal_theta_e(T, p)

def theta_e_fprime_helper(T, p):
	global latent_heat, R_vap, C_p
	return - 1.0 / EXENER(p) - (THETA_vap * latent_heat / R_vap / T**2.0) * saturated_vapor_pressure(T) / p

def solve_T_given_theta_e_and_p(theta_e, p):
	return newton((lambda T: theta_e_helper(theta_e, T, p)), theta_e, fprime=(lambda T: theta_e_fprime_helper(T, p)), tol=1e-2, maxiter=50)
