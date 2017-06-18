import numpy as np
import csv
from thermodynamics import *

class Sounding:
	def __init__(self, fname):
		self.fname = fname
		self.data = {}
		self.header = ['TIME', 'HGT', 'PRE', 'TEMP', 'RH', 'WS', 'WD', 'ASC']
		self.load()
		self.process()

	def load(self):
		for key in self.header:
			self.data[key] = []

		with open(self.fname, 'r') as csvfile:
			rdr = csv.reader(csvfile)
			for row in rdr:
				if len(row) != 8:
					continue
				
				for i, key in enumerate(self.header):
					self.data[key].append(float(row[i]))

		for key in self.header:
			self.data[key] = np.array(self.data[key])

	def process(self):
		# Transform to SI unit
		self.data['TEMP'] += zeroK
		print(self.data['TEMP'])
		self.data['PRE']  *= 100.0

		self.data['RH']   /= 100.0

		# Dew point
		dew = np.zeros((len(self.data['RH']),))
		for i in range(len(dew)):
			dew[i] = cal_dew(self.data['TEMP'][i], self.data['PRE'][i], self.data['RH'][i])
		self.data['DEW'] =  dew
