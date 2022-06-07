# This script analyzes different transformation techniques
# Copyright (C) 2020  Urban Marhl
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


class EvokedMaps:
	def __init__(self):
		import numpy as np
		self.data = np.array(()) 
		self.times = np.array(())
		self.xyz = np.array(()) 
		self.times = np.array(()) 
		self.names = np.array(())	
	
	def add_names(self, names):
		self.names = names
	
	def add_channels(self, xyz):
		self.xyz = xyz
	
	def add_times(self, times):
		self.times = times
	
	def add_data(self, data):
		self.data = data

	def add_mne_Evoked(self, Evoked):
		import numpy as np
		if self.data.size == 0:
			self.data = Evoked.data
			self.times = Evoked.times
		else:
			self.data = np.hstack((self.data, Evoked.data))
			self.times = np.hstack((self.times, Evoked.times))	
		channels = []
		for i in range(len(Evoked.info['ch_names'])):
			channels.append(Evoked.info["chs"][i]['loc'])
		channels = np.array(channels)
		self.xyz = channels
		self.names = Evoked.ch_names

	def copy(self):
		import copy
		copied_evoked = copy.deepcopy(self)
		return copied_evoked


class covariance_matrix:
	def __init__(self):
		import numpy as np
		self.data = np.array(()) 
		self.names = np.array(())
		self.std = np.array(())

	def channel_switching(self, ch_indx1, ch_indx2):
		self.data[[ch_indx1, ch_indx2]] = self.data[[ch_indx2, ch_indx1]]
		self.data[:, [ch_indx1, ch_indx2]] = self.data[:, [ch_indx2, ch_indx1]]
		self.names[[ch_indx1, ch_indx2]] = self.names[[ch_indx2, ch_indx1]]
		self.std[[ch_indx1, ch_indx2]] = self.std[[ch_indx2, ch_indx1]]
	
	def update_std(self):
		import numpy as np
		for k, j in enumerate(self.std):
			self.std[k] = np.sqrt(self.data[k, k])

	def copy(self):
		import copy
		copied = copy.deepcopy(self)
		return copied


class corrcoef_matrix:
	def __init__(self):
		import numpy as np
		self.data = np.array(()) 
		self.names = np.array(())
		self.std = np.array(())

	def channel_switching(self, ch_indx1, ch_indx2):
		self.data[[ch_indx1, ch_indx2]] = self.data[[ch_indx2, ch_indx1]]
		self.data[:,[ch_indx1, ch_indx2]] = self.data[:,[ch_indx2, ch_indx1]]
		self.names[[ch_indx1, ch_indx2]] = self.names[[ch_indx2, ch_indx1]]
		self.std[[ch_indx1, ch_indx2]] = self.std[[ch_indx2, ch_indx1]]
	
	def calculate_from_covmatrix(self, covmatrix):
		import numpy as np
		self.data = covmatrix.data / np.array((np.dot(np.matrix(covmatrix.std).T, np.matrix(covmatrix.std))))
		self.names = covmatrix.names
		self.std = covmatrix.std

	def copy(self):
		import copy
		copied = copy.deepcopy(self)
		return copied


class information_index:
	def __init__(self):
		import numpy as np
		self.data = np.array(()) 

	def calculate_index(self, cov_mat, start_indx=0):
		import numpy as np
		cal_ind = np.zeros((len(cov_mat.data)))
		for i in range(start_indx, len(cov_mat.data)):
			cal_ind[i] = np.sum((cov_mat.data[start_indx:, i]**2)/(cov_mat.std[i]**2))
		self.data = cal_ind

	def calculate_index_2(self, cov_mat, start_indx=0):
		import numpy as np
		cal_ind = np.zeros((len(cov_mat.data)))
		for i in range(start_indx, len(cov_mat.data)):
			cal_ind[i] = np.trace((1.0/(cov_mat.std[i]**2))*np.array(np.dot(np.matrix(cov_mat.data[i, start_indx:]).T, np.matrix(cov_mat.data[i, start_indx:]))))
		self.data = cal_ind

	def calculate_index_3(self, cov_mat, start_indx=0):
		import numpy as np
		cor_mat = corrcoef_matrix()
		cor_mat.calculate_from_covmatrix(cov_mat)
		cal_ind = np.zeros((len(cor_mat.data)))
		for i in range(start_indx, len(cov_mat.data)):
			a = (cor_mat.data[start_indx:, i])	
			b = (cov_mat.std[start_indx:])
			cal_ind[i] = np.sum((cor_mat.data[start_indx:, i]**2)*(cor_mat.std[start_indx:]**2))
		self.data = cal_ind

class leadsel_matrix:
	def __init__(self):
		import numpy as np
		self.data = np.array(()) 
		self.unchosen = np.array(())
		self.chosen = np.array(())
		self.error = np.array(())
	
	def calculate_lsm(self, cov_mat, no_chosen_ch):
		import numpy as np
		KMM = cov_mat.data[:no_chosen_ch,:no_chosen_ch]
		KMM_inv = np.linalg.inv(KMM)
		KUM = cov_mat.data[no_chosen_ch:,:no_chosen_ch]
		self.data = np.dot(KUM, KMM_inv)
		self.chosen = cov_mat.names[:no_chosen_ch]
		self.unchosen = cov_mat.names[no_chosen_ch:]


def get_leadsel_matrix(cov_matrix, no_best_ch, opm_sensors=False):
	import numpy as np
	work_cov = cov_matrix.copy()
	work_cov_orig = cov_matrix.copy()
	inf_ind = information_index()
	lsm = leadsel_matrix()
	RMS_error = []
	j=0
	RMS_error.append(np.sum(work_cov.std[:])/np.square(len(work_cov.std)-j))

	i=0
	while i < no_best_ch:
	# for i in range(no_best_ch):
		inf_ind.calculate_index_3(work_cov, i)
		best = np.argsort(inf_ind.data)[::-1]

		work_cov.channel_switching(best[0], i)
		work_cov_orig.channel_switching(best[0], i)

		j=i+1

		# PRVI PRISTOP
		K11I = np.linalg.inv(np.mat(work_cov.data[i, i]))
		K22 = np.mat(work_cov.data[j:, j:])
		K12 = np.mat(work_cov.data[j:, i])
		K12T = K12.T
		work_cov.data[j:, j:] = K22 - np.dot(np.dot(K12T, K11I), K12)
		work_cov.update_std()

		RMS_error.append(np.sum(work_cov.std[j:] / np.square(len(work_cov.std) - j)))
		i = i + 1

		if opm_sensors == True:
			if (best[0] % 2) == 0:
				best_pair = best[0] + 1
			else:
				best_pair = best[0] - 1

			work_cov.channel_switching(best_pair, i)
			work_cov_orig.channel_switching(best_pair, i)
			j = i + 1

			K11I = np.linalg.inv(np.mat(work_cov.data[i, i]))
			K22 = np.mat(work_cov.data[j:, j:])
			K12 = np.mat(work_cov.data[j:, i])
			K12T = K12.T
			work_cov.data[j:, j:] = K22 - np.dot(np.dot(K12T, K11I), K12)
			work_cov.update_std()

			RMS_error.append(np.sum(work_cov.std[j:] / np.square(len(work_cov.std) - j)))
			i = i + 1

		# DRUGI PRISTOP
#		work_cov = work_cov_orig.copy()
#		K11I = np.linalg.inv(np.mat(work_cov_orig.data[:j, :j]))
#		K22 = np.mat(work_cov_orig.data[j:, j:])
#		K12 = np.mat(work_cov_orig.data[j:, :j])
#		K12T = K12.T
#		work_cov.data[j:, j:] = K22 - np.dot(np.dot(K12, K11I), K12T)
#		work_cov.update_std()

	lsm.error = RMS_error
	print(np.min(RMS_error))
	lsm.calculate_lsm(work_cov_orig, no_best_ch)

	return lsm 


def get_mne_EvokedMaps(Evoked):
	evoked_maps = EvokedMaps()
	evoked_maps.add_mne_Evoked(Evoked)
	return evoked_maps


def get_covariance_matrix(EvokedMaps, time_inter=None):
	import numpy as np
	#QUESTION SHOULD STD BE ON WHOLE INTERVAL OR NOT?
	
	cov_matrix = covariance_matrix()
	if len(EvokedMaps.data) == len(EvokedMaps.names):
		cov_matrix.names = np.array(EvokedMaps.names)
	else:
		cov_matrix.names = ["CH" + item for item in np.arange(len(EvokedMaps.data)).astype(str)]

	if time_inter == None:
		cov_matrix.data = np.cov(EvokedMaps.data, bias = True)
		cov_matrix.std = np.std(EvokedMaps.data, axis=1)
	elif isinstance(time_inter, float):
		time_ind = np.argmin(abs(EvokedMaps.times - time_inter))
		evoked_copy = EvokedMaps.copy()
		evoked_copy.data = EvokedMaps.data[:,time_ind]	
		cov_matrix.data = np.cov(evoked_copy.data, bias = True)
	elif isinstance(time_inter, list):  # or type(time_inter).__module__ == np.__name__:
		if len(time_inter) == 1:
			time_ind = np.argmin(abs(EvokedMaps.times - time_inter))
			evoked_copy = EvokedMaps.copy()
			evoked_copy.data = EvokedMaps.data[:,time_ind]	
			cov_matrix.data = np.cov(evoked_copy.data, bias = True)
			cov_matrix.std = np.std(evoked_copy.data, axis=1)
		elif len(time_inter) == 2:
			time_idx =np.where((EvokedMaps.times>=time_inter[0])&(EvokedMaps.times<=time_inter[1]))[0].tolist()
			cov_matrix.data = np.cov(EvokedMaps.data[:, time_idx]
, bias = True)
			cov_matrix.std = np.std(EvokedMaps.data[:, time_idx], axis=1)
		else:
			raise ValueError('Variable time_inter should be a list with len(time_inter)=2.')
	return cov_matrix	


def get_corrcoef_matrix(EvokedMaps, time_inter=None):
	import numpy as np
	#QUESTION SHOULD STD BE ON WHOLE INTERVAL OR NOT?

	corr_matrix = corrcoef_matrix()
	if len(EvokedMaps.data) == len(EvokedMaps.names):
		corr_matrix.names = np.array(EvokedMaps.names)
	else:
		corr_matrix.names = [ "CH" + item for item in np.arange(len(EvokedMaps.data)).astype(str)]

	if time_inter == None:
		corr_matrix.data = np.corrcoef(EvokedMaps.data, bias = True)
		corr_matrix.std = np.std(EvokedMaps.data, axis=1)
	elif isinstance(time_inter, float):
		time_ind = np.argmin(abs(EvokedMaps.times - time_inter))
		evoked_copy = EvokedMaps.copy()
		evoked_copy.data = EvokedMaps.data[:,time_ind]	
		corr_matrix.data = np.corrcoef(evoked_copy.data, bias = True)
		corr_matrix.std = np.std(evoked_copy.data, axis=1)
	elif isinstance(time_inter, list):
		if len(time_inter) == 1:
			time_ind = np.argmin(abs(EvokedMaps.times - time_inter))
			evoked_copy = EvokedMaps.copy()
			evoked_copy.data = EvokedMaps.data[:,time_ind]	
			corr_matrix.data = np.corrcoef(evoked_copy.data, bias = True)
			corr_matrix.std = np.std(evoked_copy.data, axis=1)
		elif len(time_inter) == 2:
			time_ind_first = np.argmin(abs(EvokedMaps.times - time_inter[0]))
			time_ind_last = np.argmin(abs(EvokedMaps.times - time_inter[1]))
			evoked_copy = EvokedMaps.copy()
			evoked_copy.data = EvokedMaps.data[:,time_ind_first:time_ind_last+1]	
			corr_matrix.data = np.corrcoef(evoked_copy.data, bias = True)
			corr_matrix.std = np.std(evoked_copy.data, axis=1)
		else:
			raise ValueError('Variable time_inter should be a list with len(time_inter)=2.')
	return corr_matrix	


def create_evokedmaps_lsm(evokedmaps, LSM, time_inter=None):
	import numpy as np
	lsm_evoked = evokedmaps.copy()

	chosen_indx = np.zeros(len(LSM.chosen), dtype=np.int) 
	for i in range(len(LSM.chosen)):
		chosen_indx[i] = np.where(np.array(evokedmaps.names)==LSM.chosen[i])[0]

	unchosen_indx = np.zeros(len(LSM.unchosen), dtype=np.int) 
	for i in range(len(LSM.unchosen)):
		unchosen_indx[i] = np.where(np.array(evokedmaps.names)==LSM.unchosen[i])[0]
	
	unchosen_mag = np.dot(LSM.data, lsm_evoked.data[chosen_indx, :])

	for j, i in enumerate(unchosen_indx):
		lsm_evoked.data[i,:] = unchosen_mag[j,:]
	return lsm_evoked 

