import numpy as np
import csv
from thermodynamics import *

class Sounding:
	def __init__(self, fname, truncate=True):
		self.fname = fname
		self.data = {}
		self.default_header = ['TIME', 'HGT', 'PRE', 'TEMP', 'RH', 'WS', 'WD', 'ASC']
		self.load()
		self.process(truncate)

	def load(self):
		for key in self.default_header:
			self.data[key] = []

		with open(self.fname, 'r') as csvfile:
			rdr = csv.reader(csvfile)
			for row in rdr:
				if len(row) != 8:
					continue
				for i, key in enumerate(self.default_header):
					self.data[key].append(float(row[i]))

		for key in self.default_header:
			self.data[key] = np.array(self.data[key])

	def process(self, truncate):
		trunc_i = len(self.data['PRE']-1)
		if truncate:
			for i in range(len(self.data['PRE']) - 1):
				if self.data['PRE'][i] < 500.0 and self.data['PRE'][i] < self.data['PRE'][i+1]:
					trunc_i = i
					break

			for key in self.default_header:
				self.data[key] = self.data[key][0:trunc_i + 1]

		self.data_length = len(self.data['PRE'])
		# Transform to SI unit
		self.data['TEMP'] += zeroK
		self.data['PRE']  *= 100.0
		self.data['RH']   /= 100.0

		# Dew point
		dew = np.zeros((len(self.data['RH']),))
		for i in range(len(dew)):
			dew[i] = cal_dew(self.data['TEMP'][i], self.data['PRE'][i], self.data['RH'][i])

		self.data['DEW'] =  dew

		ref_q = RH_to_q(self.data['TEMP'][0], self.data['PRE'][0], self.data['RH'][0])
		self.data['LCL'], self.data['PARCEL_T'] = self.genPARCEL_T(self.data['TEMP'][0], self.data['PRE'][0], ref_q)
		self.data['CAPE_int'] = self.genCAPE_int(self.data['PARCEL_T'])

		self.data['LCL_p10'], self.data['PARCEL_T_p10'] = self.genPARCEL_T(self.data['TEMP'][0] + 10, self.data['PRE'][0], ref_q)
		self.data['CAPE_int_p10'] = self.genCAPE_int(self.data['PARCEL_T_p10'])

	def genPARCEL_T(self, ref_T, ref_p, ref_q):
		# LCL
		print(ref_q)
		ref_theta = T_to_theta(ref_T, ref_p)
		lcl = cal_LCL(ref_T, ref_p, ref_q, mode='q')

		# parcel line
		parcel = np.zeros((self.data_length,))
		LCL_theta_e = cal_theta_e(ref_T, ref_p, saturated_vapor_mass(theta_to_T(ref_theta, lcl), lcl))
		for i in range(len(parcel)):
			if self.data['PRE'][i] > lcl: # below LCL, dry line
				parcel[i] = theta_to_T(ref_theta, self.data['PRE'][i])
			else: # above LCL, wet line
				parcel[i] = solve_T_given_theta_es_and_p(LCL_theta_e, self.data['PRE'][i])
		
		return lcl, parcel
	

	def genCAPE_int(self, parcel_T):
		
		CAPE_int = np.zeros((self.data_length,))

		CAPE_int[0] = 0.0
		for i in range(self.data_length - 1):
			dH = self.data['HGT'][i + 1] - self.data['HGT'][i]
			CAPE_int[i+1] = CAPE_int[i] + ((parcel_T[i] - self.data['TEMP'][i])/self.data['TEMP'][i] + (parcel_T[i+1] - self.data['TEMP'][i+1]) / self.data['TEMP'][i]) * dH / 2.0
		
		return CAPE_int * 9.8
