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

p_ref = 100000.0
kappa = R_d / C_p
zeroK = 273.15

# This reference is rather important since Claus-Claperon approximation
# is not that good when it comes to water
vapor_pressure_ref   = 2338.8
vapor_pressure_T_ref = zeroK + 20.0

def vapor_pressure_to_mixing_ratio(p, p_vap):
	global R_d, R_vap
	return R_d / R_vap * (p / p_vap - 1.0)**(-1.0)

def saturated_vapor_pressure(T):
	global R_vap, vapor_pressure_ref, vapor_pressure_T_ref
	return vapor_pressure_ref * np.exp(- latent_heat / R_vap * (1.0/T - 1.0 / vapor_pressure_T_ref) )

def inv_saturated_vapor_pressure(p_vap):
	global R_vap, latent_heat, vapor_pressure_ref, vapor_pressure_T_ref
	return 1.0 / ( 1.0 / vapor_pressure_T_ref  - np.log(p_vap / vapor_pressure_ref) * R_vap / latent_heat) if p_vap != 0.0 else np.nan
	
def saturated_vapor_mass(T, p):  # kg / kg
	es = saturated_vapor_pressure(T)
	#return ( R_d / R_vap ) * (es / p)
	# analytic sol
#	return ( R_d / R_vap ) / (p / es - 1.0)   # analytic
	return vapor_pressure_to_mixing_ratio(p, es)   # analytic

def inv_saturated_vapor_mass(w, p):
	global R_vap, R_d
	#return inv_saturated_vapor_pressure(w * p * R_vap / R_d)

	# analytic sol
	return inv_saturated_vapor_pressure(p / (1.0 + R_d / R_vap / w))

def cal_dew(T, p, RH):
	return inv_saturated_vapor_pressure(RH * saturated_vapor_pressure(T)) 
	

def cal_LCL_helper(p, theta, q):
	"""

	p* = e_vs(T*) (1 + R_d / (q * R_v) )
	T* = EXNER(p*) theta

	residue = p* - e_vs(EXNER(p*) theta) (1 + R_d / (q * R_v) )

	"""
	global R_d, R_vap
	return p - saturated_vapor_pressure(EXNER(p) * theta) * (1.0 + R_d / q / R_vap)

def cal_LCL(T, p, RH, mode='RH'):
	"""
	[T], [p] and [RH] are informations of air parcel at the initial point of adiabatic lifting.
	
	[RH] will be treated as mixing ratio if [mode] is set to 'q'
	"""
	global R_d, R_vap
	if mode == 'RH':
		q = vapor_pressure_to_mixing_ratio(p, saturated_vapor_pressure(T) * RH)
	elif mode == 'q': # input is already q
		q = RH
	else:
		raise ValueError

	theta = T_to_theta(T, p)
	
	print("theta: %.2f, T: %.2f, p: %.2f, RH: %.2f, q: %.5f " % (theta, T, p, RH, q))
	
	return newton( (lambda pp: cal_LCL_helper(pp, theta, q)), p, fprime=None, tol=1e-3, maxiter=50 )
	

def T_to_theta(T, p):
	return T / EXNER(p)

def theta_to_T(theta, p):
	return theta * EXNER(p)


def cal_theta_es(T, p):
	return cal_theta_e(T, p, saturated_vapor_mass(T, p))

def cal_theta_e(T, p, q):
	global latent_heat, C_p
	return T_to_theta(T, p) + latent_heat * q / C_p


def EXNER(p):
	global kappa, p_ref
	return (p / p_ref) ** kappa

def theta_es_helper(theta_e, T, p):
	return theta_e - cal_theta_es(T, p)

def theta_es_fprime_helper(T, p):
	global latent_heat, R_vap, C_p
	x = p / saturated_vapor_pressure(T)
	return - 1.0 / EXNER(p) - (THETA_vap * latent_heat * x / (x - 1.0)**2.0 / R_vap / T**2.0)


def solve_T_given_theta_es_and_p(theta_e, p):
	# initial guess is important
	return newton((lambda T: theta_es_helper(theta_e, T, p)), zeroK, fprime=(lambda T: theta_es_fprime_helper(T, p)), tol=1e-3, maxiter=50)

