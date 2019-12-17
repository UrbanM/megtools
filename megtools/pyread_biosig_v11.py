def write_squid_hdrflt(data_path, evoked, block_name, template):
	from shutil import copy
	import numpy as np

	hdr_file = data_path + "/processed/"+block_name+".flt.hdr"
	flt_file = data_path + "/processed/"+block_name+".flt"
	copy(template, hdr_file)

	number_of_samples = evoked.data.shape[1]
	number_of_gradiometers = evoked.data.shape[0]

	# with is like your try .. finally block in this case
	with open(hdr_file, 'r') as file:
		# read a list of lines into data
		data = file.readlines()

	print(data)
	search1 = 'name_of_data_file='
	search2 = 'number_of_samples='

	for i in range(len(data)):
		if search1 in data[i]:
			block_name1 = block_name.replace('/', '')
			data[i] = search1 + block_name1 + ".flt\n"

	for i in range(len(data)):
		if search2 in data[i]:
			data[i] = search2 + str(number_of_samples) + "\n"

	data.append('\nmatrix_of_rotation_translation=4\n')
	data.append('*-------------------------------------------------------------------------\n')
	rot_matrix = evoked.info['dev_head_t']['trans']
	data.append("%.4f" % rot_matrix[0, 0] + "  " + "%.4f" % rot_matrix[0, 1]+ "  " +"%.4f" % rot_matrix[0, 2] + "  " + "%.4f" % rot_matrix[0, 3] +"\n")
	data.append("%.4f" % rot_matrix[1, 0] + "  " + "%.4f" % rot_matrix[1, 1]+ "  " +"%.4f" % rot_matrix[1, 2] + "  " + "%.4f" % rot_matrix[1, 3] +"\n")
	data.append("%.4f" % rot_matrix[2, 0] + "  " + "%.4f" % rot_matrix[2, 1]+ "  " +"%.4f" % rot_matrix[2, 2] + "  " + "%.4f" % rot_matrix[2, 3] +"\n")
	data.append("%.4f" % rot_matrix[3, 0] + "  " + "%.4f" % rot_matrix[3, 1]+ "  " +"%.4f" % rot_matrix[3, 2] + "  " + "%.4f" % rot_matrix[3, 3] +"\n")
	data.append("}")

	with open(hdr_file, 'w') as file:
		file.writelines(data)

	print(template[:-4])

	no_sa, no_ch, no_se, no_gr, no_mo = imp_hdr_param(hdr_file)
	flt_data = np.fromfile(template[:-4], dtype=np.float32)
	nl = float(len(flt_data)) / float(no_ch)
	nl = int(nl)
	flt_data = flt_data.reshape(nl, no_ch)

	flt_data = flt_data[:number_of_samples, :]
	for i in range(number_of_gradiometers):
		flt_data[:, i] = evoked.data[i]
	flt_data = flt_data.ravel()
	print(flt_data)

	flt_data.tofile(flt_file)

	return


def imp_samp_freq(fn):
	import re
	import numpy
	import string

	num_lines = sum(1 for line in open(fn))

	search = 'sampling_exponent='
	i = 0
	F = open(fn, 'r')
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1
	sa_ex = float(f.partition("=")[2])

	search = 'sampling_step='
	i = 0
	F = open(fn, 'r')
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1
	sa_st = float(f.partition("=")[2])

	return sa_st * 10 ** (sa_ex)


def imp_hdr_param(fn):
	import re
	import numpy
	import string

	num_lines = sum(1 for line in open(fn))

	search = 'number_of_samples='
	i = 0
	F = open(fn, 'r')
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1
	no_sa = int(''.join(list(filter(str.isdigit, f))))

	search = 'number_of_channels='
	i = 0
	F = open(fn, 'r')
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1
	no_ch = int(''.join(list(filter(str.isdigit, f))))

	search = 'number_of_sensors='
	i = 0
	F = open(fn, 'r')
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1
	no_se = int(''.join(list(filter(str.isdigit, f))))

	search = 'number_of_groups='
	i = 0
	F = open(fn, 'r')
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1
	no_gr = int(''.join(list(filter(str.isdigit, f))))

	search = 'number_of_modules='
	i = 0
	F = open(fn, 'r')
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1
	no_mo = int(''.join(list(filter(str.isdigit, f))))

	return no_sa, no_ch, no_se, no_gr, no_mo


def imp_param_channels(fn):
	# STRUCTUR OF channels
	# seq id u name calib grd grd_name grp n_sensors
	#
	import re
	import numpy

	search = 'parameter_of_channels={'
	end = '}'
	F = open(fn, 'r')
	channels = []
	num_lines = sum(1 for line in open(fn))

	i = 0
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1

	while i < num_lines:
		f = F.readline()
		if end in f:
			break
		if '\t' in f:
			None
		else:
			ff = str_list = list(filter(None, re.split('  |m\r\n| ', f)))
			if len(ff) > 6:
				channels.append(ff)
		i = i + 1
	channels = numpy.array(channels)
	return channels


def imp_param_sensors(fn):
	# STRUCTUR OF sensors
	# id  name	 type mod	x	 y	 z	 a	 b	 c	area
	#
	import re
	import numpy

	search = 'parameter_of_sensors={'
	end = '}'
	F = open(fn, 'r')
	sensors = []
	num_lines = sum(1 for line in open(fn))
	i = 0
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1

	while i < num_lines:
		f = F.readline()
		if end in f:
			break
		ff = list(filter(None, re.split('  |m\r\n| |\r\n|\n', f)))
		sensors.append(ff)
		i = i + 1
	sensors = numpy.array(sensors)
	return sensors


def imp_param_groups(fn):
	# STRUCTUR OF groups
	# id  u name	 unit   exp  calib
	#
	import re
	import numpy

	search = 'parameter_of_groups={'
	end = '}'
	F = open(fn, 'r')
	groups = []
	num_lines = sum(1 for line in open(fn))

	i = 0
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1

	while i < num_lines:
		f = F.readline()
		if end in f:
			break
		ff = filter(None, re.split('  |m\r\n| ', f))
		groups.append(ff)
		i = i + 1
	groups = numpy.array(groups)
	groups[:, 0] = groups[:, 0].astype(int)
	return groups


def imp_param_modules(fn):
	# STRUCTUR OF MODULES
	# id  name	  x	  y	  z	  a	  b	  c	  unit exp name
	#
	import re
	import numpy

	search = 'parameter_of_modules={'
	end = '}'
	F = open(fn, 'r')
	modules = []
	num_lines = sum(1 for line in open(fn))

	i = 0
	while i < num_lines:
		f = F.readline()
		if search in f:
			break
		i = i + 1

	while i < num_lines:
		f = F.readline()
		if end in f:
			break
		ff = filter(None, re.split('  |m\r\n| ', f))
		modules.append(ff)
		i = i + 1
	modules = numpy.array(modules)
	modules[:, 0] = modules[:, 0].astype(int)
	return modules


def imp_bin_data(fn1, fn2):
	import numpy
	no_sa, no_ch, no_se, no_gr, no_mo = imp_hdr_param(fn2)
	data = numpy.fromfile(fn1, dtype=numpy.float32)
	nl = float(len(data)) / float(no_ch)
	nl = int(nl)
	data = data.reshape(nl, no_ch)
	return data


def find_bad_channels(data, min_ch, max_ch):
	import numpy as np
	import vector_functions_v12 as vfun

	odklon = np.zeros(max_ch - min_ch)
	bad_ch = []

	j = 0
	for i in range(min_ch, max_ch):
		odklon[j] = vfun.root_mean_square(data[i, :])
		j += 1
	jmax = j
	mean_odklon = np.mean(odklon)

	j = 0
	for i in range(min_ch, max_ch):
		if odklon[j] < (mean_odklon / 10.0):
			bad_ch.append(i)
		j += 1

	return bad_ch


def import_squid(header_name, value_name):
	import numpy as np
	import meg_filtering_v11 as mfil
	min_chan = 0
	max_chan = 125
	min_sens = 0
	max_sens = 250

	ch_info = imp_param_channels(header_name)  # seq id u name calib grd grd_name grp n_sensors
	data = np.transpose(imp_bin_data(value_name, header_name))
	locations = imp_param_sensors(header_name)
	sfreq = (imp_samp_freq(header_name))

	# order = 5
	# cutoff = 100.0  # desired cutoff frequency of the filter, Hz
	# fs = 1 / sfreq  # sample rate, Hz
	# for i in range(min_chan, max_chan+1, 1):
	#     data[i, :] = mfil.butter_filtering(data[i, :], order, fs, cutoff)

	names = ch_info[min_chan:max_chan, 6]

	bad_channels = find_bad_channels(data, min_chan, max_chan)
	# bad_channels = []

	xyz1 = np.array(locations[min_sens:max_sens:2, 4:10], dtype=np.float32)
	xyz2 = np.array(locations[min_sens + 1:max_sens:2, 4:10], dtype=np.float32)

	mag = np.array(data[min_chan:max_chan, :], dtype=np.float32)
	mag = mag * 10.0 ** (-15.0)

	ii = np.argmax(xyz1[:, 0])
	a = (xyz1[:, 0] - xyz1[ii, 0]) ** 2
	b = (xyz1[:, 1] - xyz1[ii, 1]) ** 2
	c = (xyz1[:, 2] - xyz1[ii, 2]) ** 2

	# lenghts = np.sqrt(a + b + c)
	# nearest = np.argsort(lenghts)
	#
	# # nearest = np.delete(nearest, bad_channels, 0)
	# nearest = nearest[0:] #int(len(lenghts))]
	# nearest = np.sort(nearest)

	return xyz1[:], xyz2[:], mag[:], ch_info[:, 3], sfreq


def import_opm(header_name, value_name, bad_opm):
	import numpy as np
	import meg_filtering_v11 as mfil

	min_chan = 128
	max_chan = 158
	min_sens = 253
	max_sens = 283

	ch_info = imp_param_channels(header_name)  # seq id u name calib grd grd_name grp n_sensors
	data = np.transpose(imp_bin_data(value_name, header_name))
	locations = imp_param_sensors(header_name)
	sfreq = (imp_samp_freq(header_name))

	order = 2
	cutoff_low = 3.0  # desired cutoff frequency of the filter, Hz
	cutoff_high = 40.0
	fs = 1 / sfreq  # sample rate, Hz
	for i in range(min_chan, max_chan, 1):
		data[i, :] = mfil.butter_bandpass_filter(data[i, :], cutoff_low, cutoff_high, fs, order)

	names = ch_info[:, 6]
	# bad_channels = find_bad_channels(data, min_chan, max_chan)

	xyz1 = np.array(locations[min_sens:max_sens, 4:10], dtype=np.float32)
	mag = np.array(data[min_chan:max_chan, :], dtype=np.float32)
	mag = mag * 10.0 ** (-15.0)

	ch_info = np.array(ch_info[min_chan:max_chan, 3])

	bads = []
	if bad_opm is not None:
		for i in bad_opm:
			for j in range(len(ch_info)):
				if i == ch_info[j]:
					bads.append(j)

	# ii = np.argmax(xyz1[:, 0])
	# a = (xyz1[:, 0] - xyz1[ii, 0]) ** 2
	# b = (xyz1[:, 1] - xyz1[ii, 1]) ** 2
	# c = (xyz1[:, 2] - xyz1[ii, 2]) ** 2
	#
	# lenghts = np.sqrt(a + b + c)
	# nearest = np.argsort(lenghts)

	# nearest = np.delete(nearest, bad_channels, 0)
	# nearest = nearest[0:40]

	return xyz1, mag, ch_info, sfreq, bads


def import_triger(header_name, value_name):
	import numpy as np
	trig_channel = 158
	data = np.transpose(imp_bin_data(value_name, header_name))
	trig = np.array(data[trig_channel, :], dtype=np.float32)
	return trig


def import_time(header_name):
	import numpy as np
	sfreq = (imp_samp_freq(header_name))
	sam_num = imp_hdr_param(header_name)[0]
	time = np.arange(sam_num) * sfreq
	return time


def import_sensors(header_name, system):
	import numpy as np

	if system == "squid":
		min_chan = 0
		max_chan = 125
		min_sens = 0
		max_sens = 250

		ch_info = imp_param_channels(header_name)  # seq id u name calib grd grd_name grp n_sensors
		locations = imp_param_sensors(header_name)
		sfreq = (imp_samp_freq(header_name))
		xyz1 = np.array(locations[min_sens:max_sens:2, 4:10], dtype=np.float32)
		xyz2 = np.array(locations[min_sens + 1:max_sens:2, 4:10], dtype=np.float32)

	elif system == "opm":
		min_chan = 128
		max_chan = 158
		min_sens = 253
		max_sens = 283

		ch_info = imp_param_channels(header_name)  # seq id u name calib grd grd_name grp n_sensors
		locations = imp_param_sensors(header_name)
		sfreq = (imp_samp_freq(header_name))
		xyz1 = np.array(locations[min_sens:max_sens, 4:10], dtype=np.float32)
		xyz2 = np.array(locations[min_sens:max_sens, 4:10], dtype=None)

	else:
		print("wrong system")

	return xyz1[:], xyz2[:], ch_info[:, 3], sfreq


def create_topomap(data_path ,system="squid"):
	import mne
	import numpy as np
	import matplotlib.pyplot as plt
	default_squid_path = data_path + "ET160_template.0100.flt.hdr"
	xyz1, xyz2, ch_info, sfreq= import_sensors(default_squid_path, system)

	# print(ch_info)

	def map_et_coord(xx, yy, zz):
		# projected helmet coordinates onto a plane

		cc = np.where((xx.any() == 0.0) and (yy.any() == 0.0))

		x2 = xx * xx
		z2 = zz * zz
		y2s = (yy - 0.13526) ** (2.0)

		fact=np.sqrt((x2+z2+y2s)/(x2+z2))

		xx = -xx * fact
		yy = zz * fact

		if cc[0] != -1:
			xx[cc] = 0.0
			yy[cc] = 0.0

		return xx, yy

	xx, yy = map_et_coord(xyz1[:, 0], xyz1[:, 1], xyz1[:, 2])

	pos = np.zeros((len(ch_info), 4))

	for i in range(0,len(xx),1):
		pos[i, 0] = xx[i]
		pos[i, 1] = yy[i]
		pos[i, 2] = 0.04
		pos[i, 3] = 0.03

	ids = np.arange(len(ch_info))

	lout = mne.channels.Layout((min(xx), max(xx), min(yy), max(yy)), pos=pos, names=ch_info, ids=ids ,kind="Vectorview-all")
	# plt.show()
	lout.save(data_path+"ET160.lout")
	return


def create_topomap_uniform(xyz, names, data_path):
	import mne
	import numpy as np
	import matplotlib.pyplot as plt
	import pymeg_visualize_v13 as pvis

	xyz[:, 0] = xyz[:, 0] - np.mean(xyz[:, 0])
	xyz[:, 1] = xyz[:, 1] - np.mean(xyz[:, 1])
	xyz[:, 2] = xyz[:, 2] - np.mean(xyz[:, 2])

	xy = pvis.squid_xy_reconstruciton2(xyz[:, 0:3])

	pos = np.zeros((len(xy), 4))

	for i in range(0, len(xy), 1):
		pos[i, 0] = xy[i, 0]
		pos[i, 1] = xy[i, 1]
		pos[i, 2] = 0.04
		pos[i, 3] = 0.03

	ids = np.arange(len(xy))

	lout = mne.channels.Layout((min(xy[:, 0]), max(xy[:, 0]), min(xy[:, 1]), max(xy[:, 1])), pos=pos, names=names, ids=ids ,kind="Vectorview-all")
	plt.show()
	lout.save(data_path+"uniform_76_OPM.lout")
	return


def imp_sensor_holders(fn):
	# STRUCTUR OF SENSOR HOLDERS
	# ECHO:
	#
	import re
	import numpy as np

	F = open(fn, 'r')
	sensors = []
	orientations = []
	num_lines = sum(1 for line in open(fn))

	i = 0
	while i < num_lines:
		f = F.readline()
		f = f.replace('"','')
		f = f.replace(':', ',')
		f = f.replace('\n', '')
		ff = list(filter(None, re.split(", ", f)))
		if ff[1] == "POS":
			sensors.append([ff[3], ff[4], ff[5]])
			sensors.append([ff[3], ff[4], ff[5]])
		elif ff[1] == "DIRz":
			orientations.append([ff[3], ff[4], ff[5]])
		elif ff[1] == "DIRy":
			orientations.append([ff[3], ff[4], ff[5]])
		i = i + 1

	holders = [a + b for a,b in zip(sensors,orientations)]

	return np.array(holders).astype(np.float)


def imp_sensor_occupied(fn):
	# STRUCTUR OF SENSOR HOLDERS
	# ECHO:
	#
	import re
	import numpy as np

	F = open(fn, 'r')
	sensors = []
	orientations = []
	num_lines = sum(1 for line in open(fn))

	i = 0
	j = 0
	while i < num_lines:
		f = F.readline()
		f = f.replace('"','')
		f = f.replace(':', ',')
		f = f.replace('\n', '')
		ff = list(filter(None, re.split(", ", f)))
		if j == 0:
			sensors.append([ff[1], ff[2], ff[3]])
			sensors.append([ff[1], ff[2], ff[3]])
			j = 1
		elif j == 1:
			orientations.append([ff[1], ff[2], ff[3]])
			j = 2
		elif j == 2:
			orientations.append([ff[1], ff[2], ff[3]])
			j = 0
		i = i + 1

	holders = [a + b for a,b in zip(sensors,orientations)]

	return np.array(holders).astype(np.float)


def import_coregistration(landmark_trans, zebris_file):
	import re
	import numpy as np
	import vector_functions_v12 as vfun

	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D

	elp = np.zeros((3,3))
	elp_zebris = np.zeros((8,3))
	
	if landmark_trans is None:
		return None	
	else:
		if zebris_file is None:
			return None
		else:
			F = open(landmark_trans, 'r')
			num_lines = sum(1 for line in open(landmark_trans))
			i = 0
			while i < num_lines:
				f = F.readline()
				f = f.replace('"','')
				f = f.replace(':', ',')
				f = f.replace('\n', '')
				ff = list(filter(None, re.split(" ", f)))
				if(ff[0] == 'Na'): elp[0] = ff[1:4]
				if(ff[0] == 'LP'): elp[1] = ff[1:4]
				if(ff[0] == 'RP'): elp[2] = ff[1:4]
				i += 1

			zeb=[]
			F = open(zebris_file, 'r')
			num_lines = sum(1 for line in open(zebris_file))
			i = 0
			while i < num_lines:
				f = F.readline()
				f = f.replace('"','')
				f = f.replace(':', ',')
				f = f.replace('\n', '')
				f = f.replace('\t\t', ' ')
				f = f.replace('\t', '')
				ff = list(filter(None, re.split(" ", f)))
				ffd = np.array(ff[1:4], dtype=float) 
				if(ff[0] == 'Coil1'): elp_zebris[3] += ffd
				elif(ff[0] == 'Coil2'): elp_zebris[4] += ffd
				elif(ff[0] == 'Coil3'): elp_zebris[5] += ffd
				elif(ff[0] == 'Coil4'): elp_zebris[6] += ffd
				elif(ff[0] == 'Coil5'): elp_zebris[7] += ffd
				elif(ff[0] == 'Na'): elp_zebris[0] += ffd
				elif(ff[0] == 'LP'): elp_zebris[1] += ffd
				elif(ff[0] == 'RP'): elp_zebris[2] += ffd
				elif(i>2): zeb.append(ff[1:4])
				i += 1

			zeb = np.array(zeb, dtype=float)

			elp_zebris[:, :] /= 2
			elp = np.array(elp, dtype=float)

			centroid_elp = np.mean(elp, axis=0)
			centroid_elp_zebris = np.mean(elp_zebris[:3,:], axis=0)

			elp_trans = np.copy(elp) - centroid_elp
			elp_zebris_trans = np.copy(elp_zebris) - centroid_elp_zebris
			zeb_trans = np.copy(zeb) - centroid_elp_zebris

			# R, t = vfun.rigid_transform_3D(elp_zebris_trans[0:3,:].T, elp_trans[0:3].T)

			# elp_zebris_rot = np.copy(elp_zebris_trans)
			# for i in range(len(elp_zebris_trans)):
			# 	elp_zebris_rot[i] = np.dot(R, elp_zebris_trans[i])
			#
			# zeb_rot = np.copy(zeb_trans)
			# for i in range(len(zeb)):
			# zeb_rot[i] = np.dot(R, zeb_trans[i])
			# elp_zebris_rot = elp_zebris_rot + centroid_elp
			# zeb_rot = zeb_rot + centroid_elp

			R, t = vfun.rigid_transform_3D(elp_trans[0:3].T, elp_zebris_trans[0:3, :].T)
			elp_rot = np.copy(elp_trans)
			for i in range(len(elp_trans)):
				elp_rot[i] = np.dot(R, elp_trans[i])
			elp_rot = elp_rot + centroid_elp_zebris

			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')
			# ax.scatter(elp[:, 0], elp[:, 1], elp[:, 2], c="red")
			# ax.scatter(elp_trans[:, 0], elp_trans[:, 1], elp_trans[:, 2], c="red")
			ax.scatter(elp_rot[:, 0], elp_rot[:, 1], elp_rot[:, 2], c="red")
			ax.scatter(elp_zebris[:3, 0], elp_zebris[:3, 1], elp_zebris[:3, 2], c="green")
			# ax.scatter(elp_zebris_trans[:3, 0], elp_zebris_trans[:3, 1], elp_zebris_trans[:3, 2], c="green")
			# ax.scatter(elp_zebris_rot[:, 0], elp_zebris_rot[:, 1], elp_zebris_rot[:, 2], c="blue")
			# ax.scatter(zeb_rot[:, 0], zeb_rot[:, 1], zeb_rot[:, 2], c="black")
			plt.show()

			# elp_zebris[3:8] = vfun.rotateandtranslate(elp_zebris[3:8], R, t)

			return elp_zebris*10**(-3), zeb*10**(-3)


def import_coregistration2(zebris_file):
	import re
	import numpy as np
	import vector_functions_v12 as vfun

	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D

	elp_zebris = np.zeros((8, 3))

	zeb = []
	F = open(zebris_file, 'r')
	num_lines = sum(1 for line in open(zebris_file))
	i = 0
	while i < num_lines:
		f = F.readline()
		f = f.replace('"', '')
		f = f.replace(':', ',')
		f = f.replace('\n', '')
		f = f.replace('\t\t', ' ')
		f = f.replace('\t', '')
		ff = list(filter(None, re.split(" ", f)))
		ffd = np.array(ff[1:4], dtype=float)
		if (ff[0] == 'Coil1'):
			elp_zebris[3] += ffd
		elif (ff[0] == 'Coil2'):
			elp_zebris[4] += ffd
		elif (ff[0] == 'Coil3'):
			elp_zebris[5] += ffd
		elif (ff[0] == 'Coil4'):
			elp_zebris[6] += ffd
		elif (ff[0] == 'Coil5'):
			elp_zebris[7] += ffd
		elif (ff[0] == 'Na'):
			elp_zebris[0] += ffd
		elif (ff[0] == 'LP'):
			elp_zebris[1] += ffd
		elif (ff[0] == 'RP'):
			elp_zebris[2] += ffd
		elif (i > 2):
			zeb.append(ff[1:4])
		i += 1

	zeb = np.array(zeb, dtype=float)
	elp_zebris[:, :] /= 2

	return elp_zebris * 10 ** (-3), zeb * 10 ** (-3)


def import_opm_trans(fn,name):
	import re
	import numpy as np
	F = open(fn, 'r')
	num_lines = sum(1 for line in open(fn))
	i = 0
	while i < num_lines:
		f = F.readline()
		f = f.replace('"', '')
		f = f.replace(':', ',')
		f = f.replace('\n', '')
		f = f.replace('\t\t', ' ')
		f = f.replace('\t', '')
		if f == name:
			a=F.readline()
			aa = list(filter(None, re.split(" ", a)))
			ad = np.array(aa, dtype=float)
			b=F.readline()
			bb= list(filter(None, re.split(" ", b)))
			bd = np.array(bb, dtype=float)
			i = num_lines
		i+=1
	print(ad[0], bd)
	return ad, bd

