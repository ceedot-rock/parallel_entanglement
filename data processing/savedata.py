####################################################################
# Class: savedata.py

# Description: This class converts .mat files created by the SaveData Matlab
# class into a convenient data structure for working in Python. This can be
# useful for performing analysis in Python on any data generated using the
# Bassett Lab Matlab framework
####################################################################

import numpy as np #type: ignore
import scipy.io as spio #type: ignore
import os.path

class SaveData():
	"""
	Creates data object from .mat file created by Matlab SaveData class. Requires numpy and scipy

	Input parameters
	-----------
		filename: String containing path to file to import. Works with our without including .mat ending
		convert: Boolean deciding whether to convert data into cleaner ndarrays. Optional, default True

	Example usage
	----------
		#Import data. Make sure /bassett-lab/analysis-scripts/ is added to your path or use sys.path.append()
		from python_analysis.save_data.savedata import SaveData

		#Load in data
		sd = SaveData('S:/Projects/Materials Exploration/EuGaN/Data/MLSB1A/1acryostat_test_062024')
		
		#Prints overview of all relevant parameters and their data types. Usage should be intuitive from here.
		print(sd)
	"""

	def __init__(self, filename, convert = True):

		#Check that file exists
		exists = os.path.isfile(filename) or os.path.isfile(filename + '.mat')
		if not exists:
			raise Exception('File not found. Try inputting the full file path.')

		self.converted = False
		self.verbose = False

		#Loads file
		file = spio.loadmat(filename, simplify_cells = True)

		#Pull data from dict into class parameters
		self.matheader = file['__header__']
		self.directory = file['directory']
		self.name = file['name']
		self.nSweeps = file['nSweeps']
		self.columnNames = file['columnNames'].tolist() if isinstance(file['columnNames'], np.ndarray) else [file['columnNames']]
		self.sDataHeader = file['sDataHeader'].tolist()
		self.sData = file['sData']
		self.data = file['data']

		#Adds extraInfo if it exists
		if 'extraInfo' in file.keys():
			if isinstance(file['extraInfo'],dict):
			#if isinstance(file['extraInfo'], dict) or file['extraInfo'].size != 0:
				self.extraInfo = file['extraInfo']
			elif isinstance(file['extraInfo'],list):
				self.extraInfo = file['extraInfo']

		if convert:
			try:
				self.convert_data()
			except AssertionError as msg:
				print(msg)
				print("Creating object without converting data")

		#Converts sData into dict for convenience
		self.sDataDict = {}
		for index, key in enumerate(self.sDataHeader):
			self.sDataDict[key] = np.transpose(self.sData)[index]


	def convert_data(self):
		'''Converts data field into clean [nSweeps x nPoints] numpy matrices.
		Handles various cases differently. Called automatically by class
		constructor unless convert = False'''
		
		#Checks if data already converted and exits
		if self.converted:
			print('Data already converted!')
			return

		#Conversion if obtained a dict
		if isinstance(self.data, list):
			newdata = {}
			datakeys = self.data[0].keys()
			for key in datakeys:
				#Get amount of points per sweep
				singleton = not isinstance(self.data[0][key], np.ndarray)

				#Get number of points and pre-allocate array
				if singleton:
					nPoints = 1
					values = np.empty((self.nSweeps, 1))
				else:
					nPoints = self.data[0][key].size
					values = np.empty((self.nSweeps,) + self.data[0][key].shape)

				#Iterate over each sweep
				for iSweep in range(self.nSweeps):
					points = self.data[iSweep][key]

					#Confirm that size remains consistent
					points_size = self.data[iSweep][key].size if isinstance(self.data[0][key], np.ndarray) else 1
					assert points_size == nPoints, "Conversion failed! Number of points was not constant for each sweep."

					#Add to matrix
					values[iSweep, :] = points

				#Add to dictionary
				newdata[key] = values

			self.data = newdata
			self.converted = True
			return self.data

		#Conversion if obtained a numpy array
		elif isinstance(self.data, np.ndarray):

			#Handles case where you have [nSweeps,] ndarray containing separate ndarrays
			if self.data.size == self.nSweeps and not isinstance(self.data[0], str):

				#the number of points collected each sweep might not be constant, so get largest amount of them
				nPoints_list = []
				for iSweep in range(self.nSweeps):
					iSweep_shape = np.shape(self.data[iSweep])
					if len(iSweep_shape) == 1:
						nPoints_list.append(self.data[iSweep].size)
					elif len(iSweep_shape) == 2:
						nPoints_list.append(iSweep_shape[0])
					else:
						raise ValueError("Unexpected shape per sweep. Shape should be (N counts collected,) for a single counter, or (N counts collected,N counters)")
				nPoints = np.max(nPoints_list)
				nPoints_where = np.where(nPoints_list == nPoints)[0][0]

				#Pre-allocate array
				newdata = np.full((self.nSweeps,) + np.shape(self.data[nPoints_where]),fill_value=np.nan)
				for iSweep in range(self.nSweeps):
					points = self.data[iSweep]
					
					#Add to matrix
					newdata[iSweep, 0:len(points)] = points

				self.data = newdata
				self.converted = True
				return self.data
			
			else:
				if self.verbose:
					print("Conversion may not be necessary. Keeping data array as is.")
				return self.data

		elif isinstance(self.data, dict):
			if self.verbose:
				print("Conversion may not be necessary. Keeping data array as is.")
			self.converted = True
			return self.data

		else:
			raise UserWarning("Received unexpected data type for conversion!")


	def print_extra_info(self):
		string = ''
		assert hasattr(self, 'extraInfo'), "There is no extra info!"
		assert isinstance(self.extraInfo, dict), "Expected extraInfo to be a dictionary! This case not implemented"
		string += "sd.extraInfo: Dictionary with the following fields:\n"
		for key in self.extraInfo.keys():
			if isinstance(self.extraInfo[key], dict):
				string += "\tsd.extraInfo['{}']: Dictionary with the following fields:\n".format(key)
				for more_keys in self.extraInfo[key].keys():
					if isinstance(self.extraInfo[key][more_keys], np.ndarray):
						string += "\t\tsdextraInfo['{}']['{}']: Numpy array of shape {}\n".format(key, more_keys, self.extraInfo[key][more_keys].shape)
					else:
						string += "\t\tsd.extraInfo['{}']['{}']: {}\n".format(key, more_keys, self.extraInfo[key][more_keys])
			elif isinstance(self.extraInfo[key], np.ndarray):
				string += "\tsd.extraInfo['{}']: Numpy array of size {}\n".format(key, self.extraInfo[key].shape)
			else:
				string += "\tsd.extraInfo['{}']: '{}'\n".format(key, self.extraInfo[key])
		print(string)

	def __str__(self):
		'''Defines output printed to command line when using print() command with this object'''
		string = "SaveData object with the following properties, using 'sd' as stand-in for your variable name:\n"
		string += "\tsd.directory: Dictionary with the following fields:\n"
		for key in self.directory.keys():
			string += "\t\tsd.directory['{}']: '{}'\n".format(key, self.directory[key])
		string += "\tsd.name: '{}'\n".format(self.name)
		string += "\tsd.nSweeps: {}\n".format(self.nSweeps)
		string += "\tsd.columnNames: List of length {}\n".format(len(self.columnNames))
		string += "\tsd.sDataDict: Dictionary with the following fields:\n"
		for key in self.sDataDict.keys():
			string += "\t\tsd.sDataDict['{}']: {}\n".format(key, self.sDataDict[key])

		if hasattr(self, 'extraInfo'):
			string += "\tsd.extraInfo: Additional information. Use sd.print_extra_info() for details.\n"

		if isinstance(self.data, np.ndarray):
			string += "\tsd.data: Numpy array of shape {}".format(self.data.shape)
			return string
		if self.converted:
			string +="\tsd.data: Dictionary with the following fields:\n"
			for key in self.data.keys():
				if isinstance(self.data[key], np.ndarray):
					string += "\t\tsd.data['{}']: Numpy array of shape {}\n".format(key, self.data[key].shape)
				else:
					string += "\t\tsd.data['{}'']: {}\n".format(key, self.data[key])
		else:
			string +="\tsd.data: List of length {} containing a dictionary of values for each sweep: \n".format(len(self.data))
			for i in range(len(self.data)):
				string +="\t\tSweep {}:\n".format(i)
				for key in self.data[i].keys():
					if isinstance(self.data[i][key], np.ndarray):
						string += "\t\t\tsd.data[{}]['{}']: Numpy array of shape {}\n".format(i, key, self.data[i][key].shape)
					else:
						string += "\t\t\tsd.data[{}]['{}']: {}\n".format(i, key, self.data[i][key])
		return string 