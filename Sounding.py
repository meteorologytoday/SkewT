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
		if truncate:
			trunc_i = 0
			for i in range(len(self.data['PRE']) - 1):
				if self.data['PRE'][i] < 500.0 and self.data['PRE'][i] < self.data['PRE'][i+1]:
					trunc_i = i
					break

			for key in self.default_header:
				self.data[key] = self.data[key][0:trunc_i + 1]

		# Transform to SI unit
		self.data['TEMP'] += zeroK
		self.data['PRE']  *= 100.0
		self.data['RH']   /= 100.0

		# Dew point
		dew = np.zeros((len(self.data['RH']),))
		for i in range(len(dew)):
			dew[i] = cal_dew(self.data['TEMP'][i], self.data['PRE'][i], self.data['RH'][i])
		self.data['DEW'] =  dew
