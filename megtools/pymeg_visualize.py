def squid_butterflyplot(mag, names, freq):
	import matplotlib.pyplot as plt
	import numpy as np
	import matplotlib.ticker as ticker

	mag = mag * 10 ** (15)

	i, j = mag.shape
	x = np.arange(0, j) * freq * 10 ** (3)
	x = x - 200

	# def fmt(x, pos):
	# 	a, b = '{:.1e}'.format(x).split('e')
	# 	b = int(b)
	# 	return r'${} \times 10^{{{}}}$'.format(a, b)
	# majorformatter = ticker.FuncFormatter(fmt)

	fig = plt.figure(figsize=(10, 9))
	# fig.suptitle("Magnetic field t=%1.4f [s]" % time)

	ax = fig.add_subplot(1, 1, 1)
	for ii in range(0, i, 1):
		# if names[ii] == 'gr1s':
		ax.plot(x, mag[ii, :], lw=1.0)

	# ax.yaxis.set_major_formatter(majorformatter)
	ax.set_xlim(min(x), 800)
	ax.set_xlabel('$t$[ms]', fontsize=20)
	ax.set_ylabel('$B$[fT]', fontsize=20)
	plt.xticks(size=20)
	plt.yticks(size=20)

	fig.tight_layout()
	# fn = "PICS/butterflyplot.png"
	# plt.savefig(fn, dpi=600)  # save the figure to file
	plt.show()
	# plt.close()
	# plt.draw()
	# plt.pause(0.001)
	return 0


def simple_butterflyplot(mag, times):
	import matplotlib.pyplot as plt
	import numpy as np
	import matplotlib.ticker as ticker

	mag = mag * 10 ** (15)

	i, j = mag.shape
	# x = np.arange(0, j) * freq * 10 ** (3)
	# x = x - 200

	# def fmt(x, pos):
	# 	a, b = '{:.1e}'.format(x).split('e')
	# 	b = int(b)
	# 	return r'${} \times 10^{{{}}}$'.format(a, b)
	# majorformatter = ticker.FuncFormatter(fmt)

	fig = plt.figure(figsize=(10, 9))
	# fig.suptitle("Magnetic field t=%1.4f [s]" % time)

	ax = fig.add_subplot(1, 1, 1)
	for ii in range(0, i, 1):
		# if names[ii] == 'gr1s':
		ax.plot(times, mag[ii, :], lw=1.0)

	# ax.yaxis.set_major_formatter(majorformatter)
	# ax.set_xlim(min(x), 800)
	ax.set_xlabel('$t$[s]', fontsize=20)
	ax.set_ylabel('$B$[fT]', fontsize=20)
	plt.xticks(size=20)
	plt.yticks(size=20)

	fig.tight_layout()
	# fn = "PICS/butterflyplot.png"
	# plt.savefig(fn, dpi=600)  # save the figure to file
	plt.show()
	# plt.close()
	# plt.draw()
	# plt.pause(0.001)
	return 0


def opm_butterflyplot(mag, names, freq, bads):
	import matplotlib.pyplot as plt
	import numpy as np

	fig = plt.figure(figsize=(10, 9))

	mag = mag * 10 ** (15)

	i, j = mag.shape
	x = np.arange(0, j) * freq * 10 ** (3)
	x = x - 200

	ax1 = plt.subplot(211)
	for ii in range(0, i, 1):
		if ii in bads:
			None
		else:
			if names[ii].find('y') != -1:
				ax1.plot(x, mag[ii, :], lw=1.0)
	ax1.set_xlim(-100, 800)
	ax1.set_xlabel('$t$[ms]', fontsize=20)
	ax1.set_ylabel('$B$[fT]', fontsize=20)
	ax1.set_title("Tangential component", fontsize=16)
	plt.xticks(size=20)
	plt.yticks(size=20)

	ax2 = plt.subplot(212)
	for ii in range(0, i, 1):
		if ii in bads:
			None
		else:
			if names[ii].find('z') != -1:
				ax2.plot(x, mag[ii, :], lw=1.0)
	ax2.set_xlim(-100, 800)
	ax2.set_xlabel('$t$[ms]', fontsize=20)
	ax2.set_ylabel('$B$[fT]', fontsize=20)
	ax2.set_title("Radial component", fontsize=16)
	plt.xticks(size=20)
	plt.yticks(size=20)
	fig.tight_layout()
	fn = "PICS/butterflyplot_opm.png"
	plt.savefig(fn, dpi=600)  # save the figure to file
	plt.show()
	plt.close()
	return 0


def squid_contourplot(mag, xyz, names_ch, names_sen, ti):
	import matplotlib.pyplot as plt
	import matplotlib.mlab as ml
	from mpl_toolkits.mplot3d import Axes3D
	from scipy.interpolate import griddata
	import numpy as np

	k = len(names_sen)
	t = 0.002 * ti

	xyz = np.array(xyz, dtype=np.float32)
	mag = np.array(mag, dtype=np.float32)

	list_s = []
	for ii in range(0, k, 1):
		if names_sen[ii].find('S') == -1:
			list_s.append(ii)

	xyz_s = xyz[list_s, :]

	xyz_s = squid_xy_reconstruciton(xyz_s)

	x_min = min(xyz_s[:, 0])
	x_max = max(xyz_s[:, 0])
	y_min = min(xyz_s[:, 1])
	y_max = max(xyz_s[:, 1])
	nx, ny = 200, 200

	xi = np.linspace(x_min, x_max, nx)
	yi = np.linspace(y_min, y_max, ny)
	xi, yi = np.meshgrid(xi, yi)

	zi = griddata((xyz_s[:, 0], xyz_s[:, 1]), mag[:, ti], (xi, yi), method='cubic')  # ne pozabi pri mag na ti

	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'))
	plt.clim(np.min(mag), np.max(mag))
	plt.colorbar()
	plt.scatter(xyz_s[:, 0], xyz_s[:, 1], marker='o', c='b', s=5, zorder=10)
	ii = 0
	for i in names_ch:
		plt.text(xyz_s[ii, 0], xyz_s[ii, 1], i, fontsize=6)
		ii += 1
	plt.xlim(x_min, x_max)
	plt.ylim(y_min, y_max)
	plt.title("SourceDetection DATA t=%1.4f [s]" % t)
	fn = "PICS/pic_%04d.png" % ti
	plt.savefig(fn, dpi=300)  # save the figure to file
	# plt.show()
	plt.close()

	return 0


def mag_contourplot(mag, xyz, names_sen, time):
	import matplotlib.pyplot as plt
	from scipy.interpolate import griddata
	import vector_functions_v12 as vfun
	import numpy as np

	k = len(names_sen)
	xyz = np.array(xyz, dtype=np.float32)
	mag = np.array(mag, dtype=np.float32)

	list_s = []
	for ii in range(0, k, 1):
		if names_sen[ii].find('S') == -1:
			list_s.append(ii)

	xyz_s = xyz[:, :]

	xyz_s[:,0] = xyz_s[:,0] - np.mean(xyz_s[:,0])
	xyz_s[:, 1] = xyz_s[:, 1] - np.mean(xyz_s[:, 1])
	xyz_s[:, 2] = xyz_s[:, 2] - np.mean(xyz_s[:, 2])
	xyz_s = vfun.map_et_coord(xyz_s[:,0], xyz_s[:,1], xyz_s[:,2])
	# xyz_s = squid_xy_reconstruciton(xyz_s)

	x_min = min(xyz_s[:, 0])
	x_max = max(xyz_s[:, 0])
	y_min = min(xyz_s[:, 1])
	y_max = max(xyz_s[:, 1])
	nx, ny = 200, 200

	xi = np.linspace(x_min, x_max, nx)
	yi = np.linspace(y_min, y_max, ny)
	xi, yi = np.meshgrid(xi, yi)

	zi = griddata((xyz_s[:, 0], xyz_s[:, 1]), mag[:], (xi, yi), method='linear')  # ne pozabi pri mag na ti

	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'))
	plt.clim(np.min(mag), np.max(mag))
	plt.colorbar()
	plt.scatter(xyz_s[:, 0], xyz_s[:, 1], marker='o', c='b', s=5, zorder=10)
	ii = 0
	# for i in names_ch:
	# 	plt.text(xyz_s[ii, 0], xyz_s[ii, 1], i, fontsize=6)
	# 	ii += 1
	plt.xlim(x_min, x_max)
	plt.ylim(y_min, y_max)
	plt.title("SourceDetection DATA t=%1.4f [s]" % time)
	# fn = "PICS/pic_%04d.png" % ti
	# plt.savefig(fn, dpi=600)  # save the figure to file
	plt.show()
	plt.close()

	return 0


def opm_contourplot(mag, xyz, names, time):
	import matplotlib.pyplot as plt
	from scipy.interpolate import griddata
	import numpy as np

	i = len(mag)
	t = time

	xyz = np.array(xyz, dtype=np.float32)
	mag = np.array(mag, dtype=np.float32)

	list_z = []
	for ii in range(0, i, 1):
		if names[ii].find('z') != -1:
			list_z.append(ii)

	list_y = []
	for ii in range(0, i, 1):
		if names[ii].find('y') != -1:
			list_y.append(ii)

	xyz_z = xyz[list_z, :]
	xyz_y = xyz[list_y, :]

	mag_z = mag[list_z]
	mag_y = mag[list_y]

	names_z = []
	names_y = []
	for i in list_z:
		names_z.append(names[i])

	for i in list_y:
		names_y.append(names[i])

	# x_min = min(xyz_z[:, 0])
	# x_max = max(xyz_z[:, 0])
	y_min = min(xyz_z[:, 1])
	y_max = max(xyz_z[:, 1])
	z_min = min(xyz_z[:, 2])
	z_max = max(xyz_z[:, 2])
	ny, nz = 200, 200

	xi = np.linspace(y_min, y_max, ny)
	yi = np.linspace(z_min, z_max, nz)
	xi, yi = np.meshgrid(xi, yi)
	zi = griddata((xyz_z[:, 1], xyz_z[:, 2]), mag_z[:], (xi, yi), method='cubic')

	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'))
	plt.clim(np.min(mag_z), np.max(mag_z))
	plt.colorbar()
	plt.scatter(xyz_z[:, 1], xyz_z[:, 2], marker='o', c='b', s=5, zorder=10)
	ii = 0
	for i in names_z:
		plt.annotate(i, (xyz_z[ii, 1], xyz_z[ii, 2]))
		ii += 1
	plt.xlim(y_min, y_max)
	plt.ylim(z_min, z_max)
	plt.title("Radialna komponenta t=%1.4f [s]" % t)
	# fn = "PICS/Radial/radial_%04d.png" % ti
	# plt.savefig(fn, dpi=600)  # save the figure to file
	# plt.close()

	xi = np.linspace(y_min, y_max, ny)
	yi = np.linspace(z_min, z_max, nz)
	xi, yi = np.meshgrid(xi, yi)
	zi = griddata((xyz_y[:, 1], xyz_y[:, 2]), mag_y[:], (xi, yi), method='cubic')

	plt.subplots(nrows=1, ncols=1)  # create figure & 1 axis
	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'))
	plt.clim(np.min(mag_y), np.max(mag_y))
	plt.colorbar()
	plt.scatter(xyz_y[:, 1], xyz_y[:, 2], marker='o', c='b', s=5, zorder=10)
	ii = 0
	for i in names_y:
		plt.annotate(i, (xyz_y[ii, 1], xyz_y[ii, 2]))
		ii += 1
	plt.xlim(y_min, y_max)
	plt.ylim(z_min, z_max)
	plt.title("Tangecialna komponenta t=%1.4f [s]" % t)
	# fn = "PICS/Tangential/tangential_%04d.png" % ti
	# plt.savefig(fn, dpi=600)  # save the figure to file
	# plt.close()

	plt.show()

	return 0


def squid_xy_reconstruciton(xyz):
	import numpy as np
	import math

	xy = np.zeros((xyz.shape[0], 2))
	ind = np.argmax(xyz[:, 2])

	for i in range(0, len(xyz), 1):
		if i != ind:
			dist = np.sqrt(
				(xyz[i, 0] - xyz[ind, 0]) ** 2 + (xyz[i, 1] - xyz[ind, 1]) ** 2 + (xyz[i, 2] - xyz[ind, 2]) ** 2)
			dx = xyz[i, 0] - xyz[ind, 0]
			dy = xyz[i, 1] - xyz[ind, 1]
			theta = math.atan2(dy, dx)
			xy[i, 0] = dist * np.cos(theta) + xyz[ind, 0]
			xy[i, 1] = dist * np.sin(theta) + xyz[ind, 1]
		else:
			xy[i, 0] = xyz[ind, 0]
			xy[i, 1] = xyz[ind, 1]

	return xy


def dipole_represent_nifti(out, num, nifti_path, rotation, translation):
	import head_reconstruct_v12 as hrec
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	import numpy as np

	def rotate_via_numpy(xr, yr, radians):
		"""Use numpy to build a rotation matrix and take the dot product."""
		c, s = np.cos(radians), np.sin(radians)
		j = np.matrix([[c, s], [-s, c]])
		mm = np.dot(j, [xr, yr])
		return mm.item(0), mm.item(1)

	img, values = hrec.import_nifti(nifti_path)

	M = img.affine[:3, :3]
	abc = img.affine[:3, 3]
	shape = img.shape

	values1 = values
	values1 = np.flip(values1, 0)

	def f(i, j, k):
		""" Return X, Y, Z coordinates for i, j, k """
		return M.dot([i, j, k]) + abc

	ii = 0

	fig = plt.figure(figsize=(3*num, 8))

	for i in range(num):

		xs = 0
		ys = 0
		zs = 0
		dist = 100.0

		xx = out[0 + i * 6]
		yy = out[1 + i * 6]
		zz = out[2 + i * 6]

		yy, zz = rotate_via_numpy(yy, zz, np.radians(rotation[0]))
		xx, yy = rotate_via_numpy(xx, yy, np.radians(rotation[2]))

		xx = xx + translation[0]
		yy = yy + translation[1]
		zz = zz + translation[2]

		for l in range(shape[0]):
			j = 0
			while (j < shape[1]):
				k = 0
				while (k < shape[2]):
					res1 = f(l, j, k)[0] - xx
					res2 = f(l, j, k)[1] - yy
					res3 = f(l, j, k)[2] - zz
					if abs(res1) > 10.0:
						j = shape[1]
						k = shape[2]
						l += 1
					elif abs(res2) > 10.0:
						k = shape[2]
						j += 1
					else:
						temp_dist = np.sqrt(res1 * res1 + res2 * res2 + res3 * res3)
						if temp_dist < dist:
							xs = l
							ys = j
							zs = k
							dist = temp_dist
					k += 1
				j += 1

		if i == 0:
			pos = 1
			title = "first dipole"
		if i == 1:
			pos = 2
			title = "second dipole"

		xsi = shape[0] - xs

		ax = fig.add_subplot(3, num, pos, facecolor=(.18, .31, .31))
		circle = plt.Circle((xsi, ys), 4, color='red')
		# arr = plt.Arrow(0.5, 0.5, 10, 10, width=5.0, color='red')
		# ax.add_artist(arr)
		# plt.title(title, fontsize=20)
		ax.pcolormesh(np.transpose(values1[:, :, zs]), cmap=cm.gray)
		ax.add_artist(circle)
		# ax.set_aspect('equal')
		plt.xticks([])
		plt.yticks([])
		pos += num

		# axial plot
		ax = fig.add_subplot(3, num, pos)
		circle2 = plt.Circle((ys, zs), 4, color='red')
		# plt.title('Axial')
		ax.pcolormesh(np.transpose(values1[xs, :, :]), cmap=cm.gray)
		ax.add_artist(circle2)
		# ax.set_aspect('equal')
		plt.xticks([])
		plt.yticks([])
		pos += num

		# coronal plot
		ax = fig.add_subplot(3, num, pos)
		circle3 = plt.Circle((xsi, zs), 4, color='red')
		# plt.title('Coronal')
		ax.pcolormesh(np.transpose(values1[:, ys, :]), cmap=cm.gray)
		ax.add_artist(circle3)
		# ax.set_aspect('equal')
		plt.xticks([])
		plt.yticks([])
		pos += num
	#
	fig.tight_layout()
	fn = "PICS/location.png"
	plt.savefig(fn, dpi=600)  # save the figure to file
	plt.show()
	# plt.close()

	return 0


def dipole_represent_dicom(out, num, dicom_path, time, tint):
	import matplotlib.pyplot as plt
	import pydicom as dicom
	import numpy as np
	import matplotlib.cm as cm
	import os

	def rotate_via_numpy(xr, yr, radians):
		"""Use numpy to build a rotation matrix and take the dot product."""
		c, s = np.cos(radians), np.sin(radians)
		j = np.matrix([[c, s], [-s, c]])
		mm = np.dot(j, [xr, yr])
		return mm.item(0), mm.item(1)

	PathDicom = dicom_path
	lstFilesDCM = []  # create an empty list
	for dirName, subdirList, fileList in os.walk(PathDicom):
		for filename in fileList:
			if ".dcm" in filename.lower():  # check whether the file's DICOM
				lstFilesDCM.append(os.path.join(dirName, filename))
				lstFilesDCM.sort()

	RefDs = dicom.read_file(lstFilesDCM[1])
	ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))
	ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))

	x = np.zeros(ConstPixelDims[0])
	y = np.zeros(ConstPixelDims[1])
	z = np.zeros(len(lstFilesDCM))

	x[:] = -(RefDs.ImagePositionPatient[1] + np.arange(ConstPixelDims[0]) * float(RefDs.PixelSpacing[0]))
	y[:] = RefDs.ImagePositionPatient[2] - np.arange(ConstPixelDims[1]) * float(RefDs.PixelSpacing[1])

	for ii in range(0, len(lstFilesDCM), 1):
		dataset = dicom.read_file(lstFilesDCM[ii])
		z[ii] = dataset.ImagePositionPatient[0]

	# The array is sized based on 'ConstPixelDims'
	ArrayDicom = np.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

	# loop through all the DICOM files
	for filenameDCM in lstFilesDCM:
		# read the file
		ds = dicom.read_file(filenameDCM)
		# store the raw image data
		ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = ds.pixel_array

	fig = plt.figure(figsize=(4 * num, 13))
	# fig.suptitle("Source location t=%1.4f [s]" % time)

	for i in range(num):
		xx = out[1 + i * 6]
		yy = out[2 + i * 6]
		zz = out[0 + i * 6]

		if i == 0:
			pos = 1
			title = "first dipole"
		if i == 1:
			pos = 2
			title = "second dipole"

		xx, yy = rotate_via_numpy(xx, yy, np.radians(-12.0))

		ax = fig.add_subplot(3, num, pos, aspect='equal')
		lengths_z = (z - zz) * (z - zz)
		zi = np.argsort(lengths_z)[0]
		circle = plt.Circle((xx, yy), 4, color='red')
		# arr = plt.Arrow(0.5, 0.5, 10, 10, width=5.0, color='red')
		# ax.add_artist(arr)
		plt.title(title, fontsize=20)
		ax.pcolormesh(x, y, ArrayDicom[:, :, zi], cmap=cm.gray)
		ax.add_artist(circle)
		ax.set_xlabel('$x$[mm]', fontsize=20)
		ax.set_ylabel('$z$[mm]', fontsize=20)
		plt.xticks(size=20)
		plt.yticks(size=20)
		pos += num

		# axial plot
		ax = fig.add_subplot(3, num, pos, aspect='equal')
		lengths_y = (y - yy) * (y - yy)
		yi = np.argsort(lengths_y)[0]
		circle2 = plt.Circle((zz, xx), 4, color='red')
		# plt.title('Axial')
		ax.pcolormesh(z, x, ArrayDicom[yi, :, :], cmap=cm.gray)
		ax.add_artist(circle2)
		ax.set_xlabel('$y$[mm]', fontsize=20)
		ax.set_ylabel('$x$[mm]', fontsize=20)
		plt.xticks(size=20)
		plt.yticks(size=20)
		pos += num

		# coronal plot
		ax = fig.add_subplot(3, num, pos, aspect='equal')

		lengths_x = (x - xx) * (x - xx)
		xi = np.argsort(lengths_x)[0]
		circle3 = plt.Circle((zz, yy), 4, color='red')
		# plt.title('Coronal')
		ax.pcolormesh(z, y, ArrayDicom[:, xi, :], cmap=cm.gray)
		ax.add_artist(circle3)
		ax.set_xlabel('$y$[mm]', fontsize=20)
		ax.set_ylabel('$z$[mm]', fontsize=20)
		plt.xticks(size=20)
		plt.yticks(size=20)
		pos += num

	fig.tight_layout()
	fn = "PICS/location_voja_squid_%04d.png" % tint
	plt.savefig(fn, dpi=600)  # save the figure to file
	plt.show()
	# plt.close()

	return 0


def dipole_represent_dicom_opm(out, num, dicom_path, time, tint):
	import matplotlib.pyplot as plt
	import pydicom as dicom
	import numpy as np
	import matplotlib.cm as cm
	import os

	def rotate_via_numpy(xr, yr, radians):
		"""Use numpy to build a rotation matrix and take the dot product."""
		c, s = np.cos(radians), np.sin(radians)
		j = np.matrix([[c, s], [-s, c]])
		mm = np.dot(j, [xr, yr])
		return mm.item(0), mm.item(1)

	PathDicom = dicom_path
	lstFilesDCM = []  # create an empty list
	for dirName, subdirList, fileList in os.walk(PathDicom):
		for filename in fileList:
			if ".dcm" in filename.lower():  # check whether the file's DICOM
				lstFilesDCM.append(os.path.join(dirName, filename))
				lstFilesDCM.sort()

	RefDs = dicom.read_file(lstFilesDCM[1])
	ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))
	ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))

	x = np.zeros(ConstPixelDims[0])
	y = np.zeros(ConstPixelDims[1])
	z = np.zeros(len(lstFilesDCM))

	x[:] = -(RefDs.ImagePositionPatient[1] + np.arange(ConstPixelDims[0]) * float(RefDs.PixelSpacing[0]))
	y[:] = RefDs.ImagePositionPatient[2] - np.arange(ConstPixelDims[1]) * float(RefDs.PixelSpacing[1])

	for ii in range(0, len(lstFilesDCM), 1):
		dataset = dicom.read_file(lstFilesDCM[ii])
		z[ii] = dataset.ImagePositionPatient[0]

	# The array is sized based on 'ConstPixelDims'
	ArrayDicom = np.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

	# loop through all the DICOM files
	for filenameDCM in lstFilesDCM:
		# read the file
		ds = dicom.read_file(filenameDCM)
		# store the raw image data
		ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = ds.pixel_array

	fig = plt.figure(figsize=(10, 10))
	# fig.suptitle("Source location t=%1.4f [s]" % time)

	for i in range(num):
		xx = out[1 + i * 6]
		yy = out[2 + i * 6]
		zz = out[0 + i * 6]

		# if i == 0:
		# 	pos = 1
		# 	title = "first dipole"
		# if i == 1:
		# 	pos = 2
		# 	title = "second dipole"
		pos = 1

		# xx, yy = rotate_via_numpy(xx, yy, np.radians(-12.0))
		yy, zz = rotate_via_numpy(yy, zz, np.radians(3.0))

		ax = fig.add_subplot(2, 2, pos, aspect='equal')
		lengths_x = (x - xx) * (x - xx)
		xi = np.argsort(lengths_x)[0]
		circle3 = plt.Circle((zz, yy), 4, color='red')
		# plt.title('Coronal')
		ax.pcolormesh(z, y, ArrayDicom[:, xi, :], cmap=cm.gray)
		ax.add_artist(circle3)
		ax.set_xlabel('$y$[mm]', fontsize=20)
		ax.set_ylabel('$z$[mm]', fontsize=20)
		plt.xticks(size=20)
		plt.yticks(size=20)
		pos += num

		# axial plot
		ax = fig.add_subplot(2, 2, pos, aspect='equal')
		lengths_y = (y - yy) * (y - yy)
		yi = np.argsort(lengths_y)[0]
		circle2 = plt.Circle((zz, xx), 4, color='red')
		# plt.title('Axial')
		ax.pcolormesh(z, x, ArrayDicom[yi, :, :], cmap=cm.gray)
		ax.add_artist(circle2)
		ax.set_xlabel('$y$[mm]', fontsize=20)
		ax.set_ylabel('$x$[mm]', fontsize=20)
		plt.xticks(size=20)
		plt.yticks(size=20)
		pos += num

		# coronal plot
		ax = fig.add_subplot(2, 1, 2, aspect='equal')
		lengths_z = (z - zz) * (z - zz)
		zi = np.argsort(lengths_z)[0]
		circle = plt.Circle((xx, yy), 4, color='red')
		# arr = plt.Arrow(0.5, 0.5, 10, 10, width=5.0, color='red')
		# ax.add_artist(arr)
		# plt.title(title, fontsize=20)
		ax.pcolormesh(x, y, ArrayDicom[:, :, zi], cmap=cm.gray)
		ax.add_artist(circle)
		ax.set_xlabel('$x$[mm]', fontsize=20)
		ax.set_ylabel('$z$[mm]', fontsize=20)
		plt.xticks(size=20)
		plt.yticks(size=20)
		pos += num

	# fig.tight_layout()
	# fn = "PICS/location_voja_squid_%04d.png" % tint
	# plt.savefig(fn, dpi=600)  # save the figure to file
	plt.show()
	plt.close()

	return 0


def dipole_represent_dicom_opm_horizontal(out, num, dicom_path, time, tint):
	import matplotlib.pyplot as plt
	import pydicom as dicom
	import numpy as np
	import matplotlib.cm as cm
	import os

	def rotate_via_numpy(xr, yr, radians):
		"""Use numpy to build a rotation matrix and take the dot product."""
		c, s = np.cos(radians), np.sin(radians)
		j = np.matrix([[c, s], [-s, c]])
		mm = np.dot(j, [xr, yr])
		return mm.item(0), mm.item(1)

	PathDicom = dicom_path
	lstFilesDCM = []  # create an empty list
	for dirName, subdirList, fileList in os.walk(PathDicom):
		for filename in fileList:
			if ".dcm" in filename.lower():  # check whether the file's DICOM
				lstFilesDCM.append(os.path.join(dirName, filename))
				lstFilesDCM.sort()

	RefDs = dicom.read_file(lstFilesDCM[1])
	ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))
	ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))

	x = np.zeros(ConstPixelDims[0])
	y = np.zeros(ConstPixelDims[1])
	z = np.zeros(len(lstFilesDCM))

	x[:] = -(RefDs.ImagePositionPatient[1] + np.arange(ConstPixelDims[0]) * float(RefDs.PixelSpacing[0]))
	y[:] = RefDs.ImagePositionPatient[2] - np.arange(ConstPixelDims[1]) * float(RefDs.PixelSpacing[1])

	for ii in range(0, len(lstFilesDCM), 1):
		dataset = dicom.read_file(lstFilesDCM[ii])
		z[ii] = dataset.ImagePositionPatient[0]

	# The array is sized based on 'ConstPixelDims'
	ArrayDicom = np.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

	# loop through all the DICOM files
	for filenameDCM in lstFilesDCM:
		# read the file
		ds = dicom.read_file(filenameDCM)
		# store the raw image data
		ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = ds.pixel_array

	fig = plt.figure(figsize=(15, 4))
	# fig.suptitle("Source location t=%1.4f [s]" % time)

	for i in range(num):
		xx = out[1 + i * 6]
		yy = out[2 + i * 6]
		zz = out[0 + i * 6]

		# if i == 0:
		# 	pos = 1
		# 	title = "first dipole"
		# if i == 1:
		# 	pos = 2
		# 	title = "second dipole"
		pos = 1

		xx, yy = rotate_via_numpy(xx, yy, np.radians(-12.0))

		ax = fig.add_subplot(1, 3, pos, aspect='equal')
		lengths_x = (x - xx) * (x - xx)
		xi = np.argsort(lengths_x)[0]
		circle3 = plt.Circle((zz, yy), 4, color='red')
		# plt.title('Coronal')
		ax.pcolormesh(z, y, ArrayDicom[:, xi, :], cmap=cm.gray)
		ax.add_artist(circle3)
		ax.set_xlabel('$y$[mm]', fontsize=20)
		ax.set_ylabel('$z$[mm]', fontsize=20)
		plt.xticks(size=20)
		plt.yticks(size=20)
		pos += num

		# axial plot
		ax = fig.add_subplot(1, 3, pos, aspect='equal')
		lengths_y = (y - yy) * (y - yy)
		yi = np.argsort(lengths_y)[0]
		circle2 = plt.Circle((zz, xx), 4, color='red')
		# plt.title('Axial')
		ax.pcolormesh(z, x, ArrayDicom[yi, :, :], cmap=cm.gray)
		ax.add_artist(circle2)
		ax.set_xlabel('$y$[mm]', fontsize=20)
		ax.set_ylabel('$x$[mm]', fontsize=20)
		plt.xticks(size=20)
		plt.yticks(size=20)
		pos += num

		# coronal plot
		ax = fig.add_subplot(1, 3, pos, aspect='equal')
		lengths_z = (z - zz) * (z - zz)
		zi = np.argsort(lengths_z)[0]
		circle = plt.Circle((xx, yy), 4, color='red')
		# arr = plt.Arrow(0.5, 0.5, 10, 10, width=5.0, color='red')
		# ax.add_artist(arr)
		# plt.title(title, fontsize=20)
		ax.pcolormesh(x, y, ArrayDicom[:, :, zi], cmap=cm.gray)
		ax.add_artist(circle)
		ax.set_xlabel('$x$[mm]', fontsize=20)
		ax.set_ylabel('$z$[mm]', fontsize=20)
		plt.xticks(size=20)
		plt.yticks(size=20)
		pos += num

	fig.tight_layout()
	fn = "PICS/fig3_" + str(round(time, 3)).replace(".", "") + ".png"
	plt.savefig(fn)  # save the figure to file
	# plt.show()
	plt.close()

	return 0


def twin_squid_contourplot(mag, mag2, xyz, names_sen, time):
	import matplotlib.pyplot as plt
	import matplotlib.ticker as ticker
	from scipy.interpolate import griddata
	import numpy as np

	k = len(names_sen)
	xyz = np.array(xyz, dtype=np.float32)
	mag = np.array(mag, dtype=np.float32)
	mag2 = np.array(mag2, dtype=np.float32)

	mag = mag * 10 ** 15
	mag2 = mag2 * 10 ** 15

	list_s = []
	for ii in range(0, k, 1):
		if names_sen[ii].find('S') == -1:
			list_s.append(ii)

	xyz_s = xyz[list_s, :]
	xyz_s = squid_xy_reconstruciton(xyz_s)

	x_min = min(xyz_s[:, 0])
	x_max = max(xyz_s[:, 0])
	y_min = min(xyz_s[:, 1])
	y_max = max(xyz_s[:, 1])
	nx, ny = 200, 200

	xi = np.linspace(x_min, x_max, nx)
	yi = np.linspace(y_min, y_max, ny)
	xi, yi = np.meshgrid(xi, yi)

	zi = griddata((xyz_s[:, 0], xyz_s[:, 1]), mag[:], (xi, yi), method='cubic')  # ne pozabi pri mag na ti
	zi2 = griddata((xyz_s[:, 0], xyz_s[:, 1]), mag2[:], (xi, yi), method='cubic')  # ne pozabi pri mag na ti

	fig = plt.figure(figsize=(10, 3.5))
	# fig.suptitle("Magnetic field t=%1.4f [s]" % time)

	ax = fig.add_subplot(1, 2, 1, aspect='equal')
	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'))
	plt.clim(np.min(mag), np.max(mag))

	# def fmt(x, pos):
	# 	a, b = '{:.1e}'.format(x).split('e')
	# 	b = int(b)
	# 	return r'${} \times 10^{{{}}}$'.format(a, b)

	# cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
	cbar = plt.colorbar()
	ax.scatter(xyz_s[:, 0], xyz_s[:, 1], marker='o', c='b', s=5, zorder=10)
	ii = 0
	# for i in names_ch:
	# 	ax.text(xyz_s[ii, 0], xyz_s[ii, 1], i, fontsize=6)
	# 	ii += 1
	plt.xlim(x_min, x_max)
	plt.ylim(y_min, y_max)
	tick_spacing = 100
	ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
	ax.set_xlabel('$\Delta x$[mm]', fontsize=13)
	ax.set_ylabel('$\Delta y$[mm]', fontsize=13)
	plt.title("Measured", fontsize=13)

	cbar.ax.set_xlabel('$B$[fT]', fontsize=13)
	# cbar.ax.xaxis.set_label_position('top')
	cbar.ax.tick_params(labelsize=13)
	plt.xticks(size=13)
	plt.yticks(size=13)

	ax = fig.add_subplot(1, 2, 2, aspect='equal')
	plt.contour(xi, yi, zi2)
	plt.pcolormesh(xi, yi, zi2, cmap=plt.get_cmap('rainbow'))
	plt.clim(np.min(mag), np.max(mag))
	# cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
	cbar = plt.colorbar()
	ax.scatter(xyz_s[:, 0], xyz_s[:, 1], marker='o', c='b', s=5, zorder=10)
	ii = 0
	# for i in names_ch:
	# 	ax.text(xyz_s[ii, 0], xyz_s[ii, 1], i, fontsize=6)
	# 	ii += 1
	plt.xlim(x_min, x_max)
	plt.ylim(y_min, y_max)
	tick_spacing = 100
	ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
	ax.set_xlabel('$\Delta x$[mm]', fontsize=13)
	ax.set_ylabel('$\Delta y$[mm]', fontsize=13)
	plt.title("Calculated", fontsize=13)

	cbar.ax.set_xlabel('$B$[fT]', fontsize=13)
	# cbar.ax.xaxis.set_label_position('top')
	cbar.ax.tick_params(labelsize=13)
	plt.xticks(size=13)
	plt.yticks(size=13)

	fig.tight_layout()
	fn = "PICS/fig5_" + str(round(time, 3)).replace(".", "") + ".png"
	plt.savefig(fn)  # save the figure to file

	plt.show()
	plt.close()

	return 0


def quadro_opm_contourplot(mag, mag2, xyz, names, time, tint):
	import matplotlib.pyplot as plt
	from scipy.interpolate import griddata
	import matplotlib.ticker as ticker
	import numpy as np

	i = len(mag)
	t = time

	size = 20

	xyz = np.array(xyz, dtype=np.float32)
	mag = np.array(mag, dtype=np.float32)
	mag2 = np.array(mag2, dtype=np.float32)

	mag = mag * 10.0 ** (15.0)
	mag2 = mag2 * 10.0 ** (15.0)

	list_z = []
	for ii in range(0, i, 1):
		if names[ii].find('z') != -1:
			list_z.append(ii)

	list_y = []
	for ii in range(0, i, 1):
		if names[ii].find('y') != -1:
			list_y.append(ii)

	xyz_z = xyz[list_z, :]
	xyz_y = xyz[list_y, :]

	mag_z = mag[list_z]
	mag_y = mag[list_y]

	mag2_z = mag2[list_z]
	mag2_y = mag2[list_y]

	vmin_y = np.min([np.min(mag_y), np.min(mag2_y)])
	vmax_y = np.max([np.max(mag_y), np.max(mag2_y)])

	vmin_z = np.min([np.min(mag_z), np.min(mag2_z)])
	vmax_z = np.max([np.max(mag_z), np.max(mag2_z)])

	names_z = []
	names_y = []

	for i in list_z:
		names_z.append(names[i])

	for i in list_y:
		names_y.append(names[i])

	x_min = min(xyz_z[:, 0])
	x_max = max(xyz_z[:, 0])
	y_min = min(xyz_z[:, 1])
	y_max = max(xyz_z[:, 1])
	z_min = min(xyz_z[:, 2])
	z_max = max(xyz_z[:, 2])

	ny, nz = 200, 200

	xi = np.linspace(y_min, y_max, ny)
	yi = np.linspace(z_min, z_max, nz)
	xi, yi = np.meshgrid(xi, yi)

	# def fmt(x, pos):
	# 	a, b = '{:.1e}'.format(x).split('e')
	# 	b = int(b)
	# 	return r'${} \times 10^{{{}}}$'.format(a, b)

	fig = plt.figure(figsize=(12, 10))
	# fig.suptitle("Magnetic field t=%1.4f [s]" % time)

	ax = fig.add_subplot(2, 1, 1, aspect='auto')
	ax.set_title("Measured", fontsize=25)
	ax.tick_params(labelcolor=(1., 1., 1., 0.0), top=False, bottom=False, left=False, right=False)
	ax._frameon = False

	ax = fig.add_subplot(2, 1, 2, aspect='auto')
	ax.set_title("Fitted", fontsize=25)
	ax.tick_params(labelcolor=(1., 1., 1., 0.0), top=False, bottom=False, left=False, right=False)
	ax._frameon = False

	fig.tight_layout()

	zi = griddata((xyz_z[:, 1], xyz_z[:, 2]), mag_z[:], (xi, yi), method='cubic')
	ax = fig.add_subplot(2, 2, 1, aspect='auto')
	ax.set_title("Radial", fontsize=20)
	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_z, vmax=vmax_z)
	# plt.clim(np.min(mag_z), np.max(mag_z))
	# cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
	cbar = plt.colorbar()
	cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
	cbar.ax.xaxis.set_label_position('top')
	cbar.ax.tick_params(labelsize=size)
	plt.scatter(xyz_z[:, 1], xyz_z[:, 2], marker='o', c='b', s=5, zorder=10)
	ax.set_xlabel('$y$[mm]', fontsize=size)
	ax.set_ylabel('$z$[mm]', fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	plt.xticks(np.arange(-125.0, 26.0, 50.0))
	ii = 0
	for i in names_z:
		j = str(i).replace("O", "").replace("z", "").replace("y", "")
		plt.annotate(j, (xyz_z[ii, 1], xyz_z[ii, 2]))
		ii += 1

	zi = griddata((xyz_z[:, 1], xyz_z[:, 2]), mag_y[:], (xi, yi), method='cubic')
	ax = fig.add_subplot(2, 2, 2, aspect='auto')
	ax.set_title("Tangential", fontsize=20)
	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_y, vmax=vmax_y)
	# plt.clim(np.min(mag_z), np.max(mag_z))
	# cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
	cbar = plt.colorbar()
	cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
	cbar.ax.xaxis.set_label_position('top')
	cbar.ax.tick_params(labelsize=size)
	plt.scatter(xyz_z[:, 1], xyz_z[:, 2], marker='o', c='b', s=5, zorder=10)
	ax.set_xlabel('$y$[mm]', fontsize=size)
	ax.set_ylabel('$z$[mm]', fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	plt.xticks(np.arange(-125.0, 26.0, 50.0))
	ii = 0
	for i in names_y:
		j = str(i).replace("O", "").replace("z", "").replace("y", "")
		plt.annotate(j, (xyz_z[ii, 1], xyz_z[ii, 2]))
		ii += 1

	zi = griddata((xyz_z[:, 1], xyz_z[:, 2]), mag2_z[:], (xi, yi), method='cubic')
	ax = fig.add_subplot(2, 2, 3, aspect='auto')
	ax.set_title("Radial", fontsize=20)
	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_z, vmax=vmax_z)
	# plt.clim(np.min(mag_z), np.max(mag_z))
	# cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
	cbar = plt.colorbar()
	cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
	cbar.ax.xaxis.set_label_position('top')
	cbar.ax.tick_params(labelsize=size)
	plt.scatter(xyz_z[:, 1], xyz_z[:, 2], marker='o', c='b', s=5, zorder=10)
	ax.set_xlabel('$y$[mm]', fontsize=size)
	ax.set_ylabel('$z$[mm]', fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	plt.xticks(np.arange(-125.0, 26.0, 50.0))
	ii = 0
	for i in names_z:
		j = str(i).replace("O", "").replace("z", "").replace("y", "")
		plt.annotate(j, (xyz_z[ii, 1], xyz_z[ii, 2]))
		ii += 1

	zi = griddata((xyz_z[:, 1], xyz_z[:, 2]), mag2_y[:], (xi, yi), method='cubic')
	ax = fig.add_subplot(2, 2, 4, aspect='auto')
	ax.set_title("Tangential", fontsize=20)
	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_y, vmax=vmax_y)
	# plt.clim(np.min(mag_z), np.max(mag_z))
	# cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
	cbar = plt.colorbar()
	cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
	cbar.ax.xaxis.set_label_position('top')
	cbar.ax.tick_params(labelsize=size)
	plt.scatter(xyz_z[:, 1], xyz_z[:, 2], marker='o', c='b', s=5, zorder=10)
	ax.set_xlabel('$y$[mm]', fontsize=size)
	ax.set_ylabel('$z$[mm]', fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	plt.xticks(np.arange(-125.0, 26.0, 50.0))
	ii = 0
	for i in names_y:
		j = str(i).replace("O", "").replace("z", "").replace("y", "")
		plt.annotate(j, (xyz_z[ii, 1], xyz_z[ii, 2]))
		ii += 1

	fig.tight_layout()
	fn = "PICS/voja_opm_%04d.png" % tint
	plt.savefig(fn, dpi=600)  # save the figure to file

	plt.show()
	plt.close()

	# plt.xlim(y_min, y_max)
	# plt.ylim(z_min, z_max)
	# plt.title("Radialna komponenta t=%1.4f [s]" % t)
	# # fn = "PICS/Radial/radial_%04d.png" % ti
	# # plt.savefig(fn, dpi=600)  # save the figure to file
	# # plt.close()
	#
	# xi = np.linspace(y_min, y_max, ny)
	# yi = np.linspace(z_min, z_max, nz)
	# xi, yi = np.meshgrid(xi, yi)
	# zi = griddata((xyz_y[:, 1], xyz_y[:, 2]), mag_y[:], (xi, yi), method='cubic')
	#
	# plt.subplots(nrows=1, ncols=1)  # create figure & 1 axis
	# plt.contour(xi, yi, zi)
	# plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'))
	# plt.clim(np.min(mag_y), np.max(mag_y))
	# plt.colorbar()
	# plt.scatter(xyz_y[:, 1], xyz_y[:, 2], marker='o', c='b', s=5, zorder=10)
	# ii = 0
	# for i in names_y:
	# 	plt.annotate(i, (xyz_y[ii, 1], xyz_y[ii, 2]))
	# 	ii += 1
	# plt.xlim(y_min, y_max)
	# plt.ylim(z_min, z_max)
	# plt.title("Tangecialna komponenta t=%1.4f [s]" % t)
	# # fn = "PICS/Tangential/tangential_%04d.png" % ti
	# # plt.savefig(fn, dpi=600)  # save the figure to file
	# # plt.close()
	#
	# plt.show()

	return 0


def quadro_opm_contourplot_horizontal(mag, mag2, xyz, names, time, tint):
	import matplotlib.pyplot as plt
	from scipy.interpolate import griddata
	import matplotlib.ticker as ticker
	import numpy as np

	i = len(mag)
	t = time

	size = 16

	xyz = np.array(xyz, dtype=np.float32)
	mag = np.array(mag, dtype=np.float32)
	mag2 = np.array(mag2, dtype=np.float32)

	mag = mag * 10.0 ** (15.0)
	mag2 = mag2 * 10.0 ** (15.0)

	list_z = []
	for ii in range(0, i, 1):
		if names[ii].find('z') != -1:
			list_z.append(ii)

	list_y = []
	for ii in range(0, i, 1):
		if names[ii].find('y') != -1:
			list_y.append(ii)

	xyz_z = xyz[list_z, :]
	xyz_y = xyz[list_y, :]

	xyz = xyz_z

	mag_z = mag[list_z]
	mag_y = mag[list_y]

	mag2_z = mag2[list_z]
	mag2_y = mag2[list_y]

	# vmin_y = np.min([np.min(mag_y), np.min(mag2_y)])
	# vmax_y = np.max([np.max(mag_y), np.max(mag2_y)])
	#
	# vmin_z = np.min([np.min(mag_z), np.min(mag2_z)])
	# vmax_z = np.max([np.max(mag_z), np.max(mag2_z)])

	vmin_y = np.min(mag_y)
	vmax_y = np.max(mag_y)

	vmin_z = np.min(mag_z)
	vmax_z = np.max(mag_z)


	names_z = []
	names_y = []

	for i in list_z:
		names_z.append(names[i])

	for i in list_y:
		names_y.append(names[i])

	x_min = min(xyz_z[:, 0])
	x_max = max(xyz_z[:, 0])
	y_min = min(xyz_z[:, 1])
	y_max = max(xyz_z[:, 1])
	z_min = min(xyz_z[:, 2])
	z_max = max(xyz_z[:, 2])

	ny, nz = 200, 200

	xi = np.linspace(y_min, y_max, ny)
	yi = np.linspace(z_min, z_max, nz)
	xi, yi = np.meshgrid(xi, yi)

	aspect_ratio = (z_max - z_min) / (y_max - y_min)

	fig = plt.figure(figsize=(2*10, 10 * aspect_ratio * 0.9 * 2))
	fig.tight_layout()

	ax = fig.add_subplot(2, 2, 1)
	# ax.set_aspect(aspect_ratio)
	zi = griddata((xyz[:, 1], xyz[:, 2]), mag_z[:], (xi, yi), method='cubic')
	plt.contour(xi, yi, zi, aspect=aspect_ratio)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_z, vmax=vmax_z)
	cbar = plt.colorbar(ticks=[vmin_z, 0, vmax_z])
	cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
	# cbar.ax.xaxis.set_label_position('right')
	cbar.ax.tick_params(labelsize=size)
	plt.scatter(xyz[:, 1], xyz[:, 2], marker='o', c='b', s=5, zorder=10)
	ax.set_xlabel('$y$[mm]', fontsize=size)
	ax.set_ylabel('$z$[mm]', fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	plt.xlim(y_min - 1, y_max + 1)
	plt.ylim(z_min - 1, z_max + 1)
	# ax.set_title("Measured radial component", fontsize=16)
	ax.set_title("RADIAL COMPONENT", fontsize=20)

	ax = fig.add_subplot(2, 2, 3)
	zi = griddata((xyz[:, 1], xyz[:, 2]), mag_y[:], (xi, yi), method='cubic')
	plt.contour(xi, yi, zi, aspect=aspect_ratio)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_y, vmax=vmax_y)
	cbar = plt.colorbar(ticks=[vmin_y, 0, vmax_y])
	cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
	cbar.ax.tick_params(labelsize=size)
	plt.scatter(xyz[:, 1], xyz[:, 2], marker='o', c='b', s=5, zorder=10)
	ax.set_xlabel('$y$[mm]', fontsize=size)
	ax.set_ylabel('$z$[mm]', fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	plt.xlim(y_min - 1, y_max + 1)
	plt.ylim(z_min - 1, z_max + 1)
	# ax.set_title("Measured tangential component", fontsize=16)
	ax.set_title("TANGENTIAL COMPONENT", fontsize=20)

	ax = fig.add_subplot(2, 2, 2)
	zi = griddata((xyz[:, 1], xyz[:, 2]), mag2_z[:], (xi, yi), method='cubic')
	plt.contour(xi, yi, zi, aspect=aspect_ratio)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_z, vmax=vmax_z)
	cbar = plt.colorbar(ticks=[vmin_z, 0, vmax_z])
	cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
	cbar.ax.tick_params(labelsize=size)
	plt.scatter(xyz[:, 1], xyz[:, 2], marker='o', c='b', s=5, zorder=10)
	ax.set_xlabel('$y$[mm]', fontsize=size)
	ax.set_ylabel('$z$[mm]', fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	plt.xlim(y_min - 1, y_max + 1)
	plt.ylim(z_min - 1, z_max + 1)
	# ax.set_title("Fitted radial component", fontsize=16)
	ax.set_title("RADIAL COMPONENT", fontsize=20)

	ax = fig.add_subplot(2, 2, 4)
	zi = griddata((xyz[:, 1], xyz[:, 2]), mag2_y[:], (xi, yi), method='cubic')
	plt.contour(xi, yi, zi, aspect=aspect_ratio)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_y, vmax=vmax_y)
	cbar = plt.colorbar(ticks=[vmin_y, 0, vmax_y])
	cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
	cbar.ax.tick_params(labelsize=size)
	plt.scatter(xyz[:, 1], xyz[:, 2], marker='o', c='b', s=5, zorder=10)
	ax.set_xlabel('$y$[mm]', fontsize=size)
	ax.set_ylabel('$z$[mm]', fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	plt.xlim(y_min - 1, y_max + 1)
	plt.ylim(z_min - 1, z_max + 1)
	# ax.set_title("Fitted tangential component", fontsize=16)
	ax.set_title("TANGENTIAL COMPONENT", fontsize=20)

	fig.tight_layout()
	fn = "PICS/voja_opm_%04d.png" % tint
	plt.savefig(fn)  # save the figure to file

	plt.show()
	plt.close()

	return 0


def plot_opm_data(mag_z, mag_y, xyz, time):
	import matplotlib.pyplot as plt
	from scipy.interpolate import griddata
	import numpy as np

	size = 16

	no_arrays = 0
	if isinstance(mag_z, (list, tuple, np.ndarray)):
		no_arrays += 1
	if isinstance(mag_y, (list, tuple, np.ndarray)):
		no_arrays += 1

	xyz = np.array(xyz, dtype=np.float32)
	x_min = min(xyz[:, 0])
	x_max = max(xyz[:, 0])
	y_min = min(xyz[:, 1])
	y_max = max(xyz[:, 1])
	z_min = min(xyz[:, 2])
	z_max = max(xyz[:, 2])
	ny, nz = 200, 200

	xi = np.linspace(y_min, y_max, ny)
	yi = np.linspace(z_min, z_max, nz)
	xi, yi = np.meshgrid(xi, yi)

	aspect_ratio = (z_max - z_min) / (y_max - y_min)

	fig = plt.figure(figsize=(10 * no_arrays, 10 * aspect_ratio * 0.9))
	fig.tight_layout()

	if isinstance(mag_z, (list, tuple, np.ndarray)):
		mag_z = np.array(mag_z, dtype=np.float32)
		mag_z = mag_z * 10.0 ** (15.0)

		vmin_z = np.min(mag_z)
		vmax_z = np.max(mag_z)

		ax = fig.add_subplot(1, no_arrays, 1)

		zi = griddata((xyz[:, 1], xyz[:, 2]), mag_z[:], (xi, yi), method='cubic')
		plt.contour(xi, yi, zi, aspect=aspect_ratio)
		plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_z, vmax=vmax_z)
		cbar = plt.colorbar(ticks=[vmin_z, 0, vmax_z])
		cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
		cbar.ax.tick_params(labelsize=size)
		plt.scatter(xyz[:, 1], xyz[:, 2], marker='o', c='b', s=5, zorder=10)
		ax.set_xlabel('$y$[mm]', fontsize=size)
		ax.set_ylabel('$z$[mm]', fontsize=size)
		plt.xticks(size=size)
		plt.yticks(size=size)
		plt.xlim(y_min - 1, y_max + 1)
		plt.ylim(z_min - 1, z_max + 1)
		ax.set_title("Radial component", fontsize=16)

	if isinstance(mag_y, (list, tuple, np.ndarray)):
		mag_y = np.array(mag_y, dtype=np.float32)
		mag_y = mag_y * 10.0 ** (15.0)

		vmin_y = np.min(mag_y)
		vmax_y = np.max(mag_y)

		ax = fig.add_subplot(1, no_arrays, no_arrays)

		zi = griddata((xyz[:, 1], xyz[:, 2]), mag_y[:], (xi, yi), method='cubic')
		plt.contour(xi, yi, zi, aspect=aspect_ratio)
		plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'), vmin=vmin_y, vmax=vmax_y)
		cbar = plt.colorbar(ticks=[vmin_y, 0, vmax_y])
		cbar.ax.set_xlabel('$B$[fT]', fontsize=size)
		cbar.ax.tick_params(labelsize=size)
		plt.scatter(xyz[:, 1], xyz[:, 2], marker='o', c='b', s=5, zorder=10)
		ax.set_xlabel('$y$[mm]', fontsize=size)
		ax.set_ylabel('$z$[mm]', fontsize=size)
		plt.xticks(size=size)
		plt.yticks(size=size)
		plt.xlim(y_min - 1, y_max + 1)
		plt.ylim(z_min - 1, z_max + 1)
		ax.set_title("Tangential component", fontsize=16)

	fig.tight_layout()
	fn = "PICS/fig2_" + str(round(time, 3)).replace(".", "") + ".png"
	plt.savefig(fn)  # save the figure to file


def mag_contourplot_simulation(mag, xyz):
	import matplotlib.pyplot as plt
	from scipy.interpolate import griddata
	import numpy as np

	xyz = np.array(xyz, dtype=np.float32)
	mag = np.array(mag, dtype=np.float32)

	xyz_s = squid_xy_reconstruciton(xyz)

	x_min = min(xyz_s[:, 0])
	x_max = max(xyz_s[:, 0])
	y_min = min(xyz_s[:, 1])
	y_max = max(xyz_s[:, 1])
	nx, ny = 200, 200

	xi = np.linspace(x_min, x_max, nx)
	yi = np.linspace(y_min, y_max, ny)
	xi, yi = np.meshgrid(xi, yi)

	zi = griddata((xyz_s[:, 0], xyz_s[:, 1]), mag[:], (xi, yi), method='cubic')  # ne pozabi pri mag na ti

	plt.contour(xi, yi, zi)
	plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('rainbow'))
	plt.clim(np.min(mag), np.max(mag))
	plt.colorbar()
	plt.scatter(xyz_s[:, 0], xyz_s[:, 1], marker='o', c='b', s=5, zorder=10)
	ii = 0
	plt.xlim(x_min, x_max)
	plt.ylim(y_min, y_max)

	plt.show()
	plt.close()
	return 0


def plot_vfield(x, y, z, dx, dy, dz):
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.quiver(x, y, z, dx, dy, dz)
	# cont = ax.contour(x, y, p, cmap='gist_earth')
	# ax.clabel(cont)
	ax.set(aspect=1, title='Streamplot with contours')
	plt.show()
	return


def plot_sens_orient(xyz):
	from matplotlib import pyplot
	from mpl_toolkits.mplot3d import Axes3D
	fig = pyplot.figure()
	ax = Axes3D(fig)
	ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], alpha=0.6, c='red', edgecolors='none', s=30)
	ax.quiver(xyz[:, 0], xyz[:, 1], xyz[:, 2], 10 * xyz[:, 3], 10 * xyz[:, 4], 10 * xyz[:, 5])
	ax.grid(False)
	ax.axis('off')
	pyplot.show()
	return


def simple_plot(x, y, name = "default", yaxis = "y-os", xaxis="x-os"):
	import matplotlib.pyplot as plt
	size = 16
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.xlim(min(x), max(x))
	plt.ylim(min(y), max(y))
	# ax.scatter(x, y)
	ax.plot(x, y, '-')
	ax.set_ylabel(yaxis, fontsize=size)
	ax.set_xlabel(xaxis, fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	fig.tight_layout()
	# fn = "PICS/" + name + ".png"  # "PICS/fig3.png"
	# plt.savefig(fn)  # save the figure to file
	plt.show(block=False)

	return 0


def plot_channel(x, y, ch_name, yaxis, xaxis, x_range):
	import matplotlib.pyplot as plt
	size = 16
	fig = plt.figure()
	ax = fig.add_subplot(111)

	if isinstance(x_range, list):  # type() == list():
		plt.xlim(x_range)
	else:
		plt.xlim(min(x), max(x))
	# plt.ylim(min(y), max(y))
	# ax.scatter(x, y)
	ax.plot(x, y, '-')
	ax.set_ylabel(yaxis, fontsize=size)
	ax.set_xlabel(xaxis, fontsize=size)
	plt.xticks(size=size)
	plt.yticks(size=size)
	fig.tight_layout()
	# fn = "PICS/" + name + ".png"  # "PICS/fig3.png"
	# plt.savefig(fn)  # save the figure to file
	plt.show(block=False)

	return 0


def beamformers_represent_nifti(neu_act, sourcespace, nifti_path, rotation, translation, threshold):
	import head_reconstruct_v11 as hrec
	import vector_functions_v12 as vfun
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	import matplotlib
	import numpy as np

	img, values = hrec.import_nifti(nifti_path)

	M = img.affine[:3, :3]
	abc = img.affine[:3, 3]
	shape = img.shape

	# values1 = values
	# values1 = np.flip(values1, 0)

	def f(i, j, k):
		""" Return X, Y, Z coordinates for i, j, k """
		return M.dot([i, j, k]) + abc

	ii = 0

	max_col = np.max(neu_act)
	min_col = np.min(neu_act)

	scale = min_col + np.linspace(0, 20, 21)*((max_col-min_col)/21.0)
	# print(scale)

	ii = np.argmax(neu_act)
	xx = sourcespace[ii, 0]
	yy = sourcespace[ii, 1]
	zz = sourcespace[ii, 2]

	yy, zz = vfun.rotate_via_numpy(yy, zz, -np.radians(rotation[0]))
	xx, yy = vfun.rotate_via_numpy(xx, yy, -np.radians(rotation[2]))

	xx = xx - translation[0]
	yy = yy - translation[1]
	zz = zz - translation[2]

	dist = pow(10.0, 6.0)
	for l in range(shape[0]):
		j = 0
		while (j < shape[1]):
			k = 0
			while (k < shape[2]):
				res1 = f(l, j, k)[0] - xx
				res2 = f(l, j, k)[1] - yy
				res3 = f(l, j, k)[2] - zz
				if abs(res1) > 10.0:
					j = shape[1]
					k = shape[2]
					l += 1
				elif abs(res2) > 10.0:
					k = shape[2]
					j += 1
				else:
					temp_dist = np.sqrt(res1 * res1 + res2 * res2 + res3 * res3)
					if temp_dist < dist:
						xs = l
						ys = j
						zs = k
						dist = temp_dist
				k += 1
			j += 1

	XX = np.zeros(shape[0])
	YY = np.zeros(shape[1])
	ZZ = np.zeros(shape[2])
	for l in range(shape[0]):
		XX[l] = f(l, ys, zs)[0]
	for k in range(shape[1]):
		YY[k] = f(xs, k, zs)[1]
	for j in range(shape[2]):
		ZZ[j] = f(xs, ys, j)[2]


	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	ax.pcolormesh(YY, XX, values[:, :, zs], cmap=cm.gray)

	xx = sourcespace[:, 0]
	yy = sourcespace[:, 1]
	zz = sourcespace[:, 2]
	yy, zz = vfun.rotate_via_numpy(yy, zz, -np.radians(rotation[0]))
	xx, yy = vfun.rotate_via_numpy(xx, yy, -np.radians(rotation[2]))
	xx = xx - translation[0]
	yy = yy - translation[1]
	zz = zz - translation[2]

	for ii in range(len(neu_act)):
		if abs(zz[ii] - f(xs, ys, zs)[2]) < 5.0 and neu_act[ii] > threshold:
			rectangle = matplotlib.patches.Rectangle((yy[ii], xx[ii]), 10, 10, color='red')
			ax.add_artist(rectangle)
	plt.xticks([])
	plt.yticks([])

	fig.tight_layout()

	plt.show()

	return 0


def sensors_visualize(wolke_trans, head_name, cortex_name, sensors_ocupied = None): #, sensorholders, header_name, value_name, vertices_name):
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	import numpy as np
	import vector_functions_v12 as vfun
	import pyread_biosig_v11 as pbio

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	# vertices = np.loadtxt(vertices_name, dtype=np.float)
	#

	xyzw = np.loadtxt(wolke_trans, dtype=np.float)
	xyzh = np.loadtxt(head_name, dtype=np.float)
	xyzc = np.loadtxt(cortex_name, dtype=np.float)
	# sensors = pbio.imp_sensor_occupied(sensors_ocupied)

	# xyzw[:, 0] = xyzw[:, 0]

	# xyzh[:,1], xyzh[:,2] = vfun.rotate_via_numpy(xyzh[:,1], xyzh[:,2], np.radians(rotation[0]))
	# xyzh[:,0], xyzh[:,1] = vfun.rotate_via_numpy(xyzh[:,0], xyzh[:,1], np.radians(rotation[2]))
	#
	# xyzh[:,0] = xyzh[:,0] + translation[0]
	# xyzh[:,1] = xyzh[:,1] + translation[1]
	# xyzh[:,2] = xyzh[:,2] + translation[2]


	# xyz1, xyz2, mag, names_sen, mes_step = pbio.import_squid(header_name, value_name)
	# xyz1[:, 0:3] = xyz1[:, 0:3] 
	# xyz2[:, 0:3] = xyz2[:, 0:3]
	# 
	# ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c='red')
	# ax.scatter(xyz1[:, 0], xyz1[:, 1], xyz1[:, 2], c='gray')
	# ax.scatter(sensorholders[:, 0], sensorholders[:, 1], sensorholders[:, 2], c='green')
	# ax.quiver(sensorholders[:, 0], sensorholders[:, 1], sensorholders[:, 2], 10 * sensorholders[:, 3], 10 * sensorholders[:, 4], 10 * sensorholders[:, 5], color='green')

	ax.scatter(xyzh[::60, 0], xyzh[::60, 1], xyzh[::60, 2], c='blue')
	ax.scatter(xyzw[::2, 0], xyzw[::2, 1], xyzw[::2, 2], c='gray')
	# ax.scatter(sensors[:, 0], sensors[:, 1], sensors[:, 2], c='orange')
	ax.scatter(xyzc[:, 0], xyzc[:, 1], xyzc[:, 2], c='red')
	plt.show()
	return
