def trans_cart_spher(dip1):
	import math
	r = math.sqrt(dip1[0]**2 + dip1[1]**2 + dip1[2]**2)
	theta = math.acos(dip1[2]/r)
	phi = math.atan(dip1[1]/dip1[0])
	n_dip = [None]*3
	n_dip[0] = r
	n_dip[1] = theta
	n_dip[2] = phi

	return n_dip


def trans_spher_cart(dip1):
	import math
	x = dip1[0] * math.sin(dip1[1]) * math.cos(dip1[2])
	y = dip1[0] * math.sin(dip1[1]) * math.sin(dip1[2])
	z = dip1[0] * math.cos(dip1[1])
	n_dip = [None]*3
	n_dip[0] = x
	n_dip[1] = y
	n_dip[2] = z
	return n_dip


def rotation_x(alpha):
	from math import cos, sin
	matrix = [[1, 0, 0], [0, cos(alpha), -sin(alpha)], [0, sin(alpha), cos(alpha)]]
	return matrix


def rotation_y(alpha):
	from math import cos, sin
	matrix = [[cos(alpha), 0, sin(alpha)], [0, 1, 0], [-sin(alpha), 0, cos(alpha)]]
	return matrix


def rotation_z(alpha):
	from math import cos, sin
	matrix = [[cos(alpha), -sin(alpha), 0], [sin(alpha), cos(alpha), 0], [0, 0, 1]]
	return matrix


def rm_radial_component(dip):
	import numpy as np
	import random
	from math import cos, sin, pi
	from numpy import linalg as LA

	radius, theta, phi = trans_cart_spher(dip[0:3])
	theta = -theta
	phi = -phi

	dipole_dir = dip[3:6]

	# transformacija nazaj v prvotni koordinatni sistem !!ROTACIJSKE MATRIKE!!
	dipole_dir = np.dot(rotation_z(phi), dipole_dir)
	dipole_dir = np.dot(rotation_y(theta), dipole_dir)

	dipole_dir[2] = 0.0
	#
	dipole_dir = np.dot(rotation_y(-theta), dipole_dir)
	dipole_dir = np.dot(rotation_z(-phi), dipole_dir)

	dipole = np.concatenate((dip[0:3], dipole_dir), axis=0)

	return dipole


def create_net_sphere_hex(circles, direction, lenght, radius):
	import math
	import numpy as np

	# INITIALIZING
	positions = []
	spher_positions = []
	orientations = []
	directions = []
	dirr = [None]*3
	positions_2d = []

	for b in range(0, 3, 1):
		if direction[b] > 0:
			for c in range(0, circles, 1):
				if c > 0:
					d_phi = 2 * math.pi / (c * 6)
					d_theta = math.pi / ((circles) * 2)   # math.atan(lenght / radius)
					dmax = (c * 6)
				else:
					d_phi = 0.0
					d_theta = 0.0
					dmax = 1

				theta = d_theta * c
				for d in range(0, dmax, 1):
					phi = d * d_phi
					positions.append([radius * math.sin(theta) * math.cos(phi),
									  radius * math.sin(phi) * math.sin(theta),
									  radius * math.cos(theta)])
					spher_positions.append([radius, theta, phi])
					positions_2d.append([radius * c * math.cos(phi), radius * c * math.sin(phi)])
					dirr[0] = [math.cos(theta) * math.cos(phi), math.cos(theta) * math.sin(phi), -math.sin(theta)]
					dirr[1] = [-math.sin(phi), math.cos(phi), 0]
					dirr[2] = [math.sin(theta) * math.cos(phi), math.sin(phi) * math.sin(theta), math.cos(theta)]
					orientations.append(dirr[b])
					directions.append(b)

	return np.array(positions), np.array(orientations), np.array(spher_positions), np.array(directions), np.array(positions_2d)


def create_2d_square_sp(n1, n2, d, z):
	import numpy
	source_space = numpy.zeros((n1*n2, 6))
	ii = 0
	for i in range(0, n1, 1):
		for j in range(0, n2, 1):
			source_space[ii, 0] = i * d
			source_space[ii, 1] = j * d
			source_space[ii, 2] = z
			ii += 1

	return source_space


def create_2d_square_me(n1, n2, d, z):
	import numpy
	source_space = numpy.zeros((n1*n2, 6))
	ii = 0
	for i in range(0, n1, 1):
		for j in range(0, n2, 1):
			source_space[ii, 0] = i * d
			source_space[ii, 1] = j * d
			source_space[ii, 2] = z
			source_space[ii, 3] = 0.0
			source_space[ii, 4] = 0.0
			source_space[ii, 5] = 1.0
			ii += 1

	return source_space


def root_mean_square(field):
	import numpy as np
	rms = np.sqrt(np.mean(np.square(field)))
	return rms


def R_squared(izr, izm):
	import numpy as np
	# FOR LINEAR FIT ONLY
	SE_L = np.sum((izm - izr)*(izm - izr))
	mean_izm = np.average(izm)
	mean_izr = np.average(izr)
	SE_Y = np.sum((izm - mean_izr) * (izm - mean_izr))
	RR = 1 - SE_L / SE_Y
	if RR < 0: RR = 0.0
	return RR


def CHI_squared(izr, izm):
	# Pearson's chi-squared test
	import numpy as np
	CHICHI = np.sum((izm - izr)**2/izr)
	return CHICHI


def RMSerr(izr, izm):
	# defined as vojkos
	# rmserr @ root mean square error
	import numpy as np
	rms = np.sum((izm - izr)**2)
	rms = rms / len(izm)
	return rms


def rel_err(izr, izm):
	# defined as vojkos
	import numpy as np
	rms1 = np.sqrt(np.sum(izr * izr) / len(izr))
	rms2 = np.sqrt(np.sum(izm * izm) / len(izm))

	rel = np.sum((izm - izr)**2)
	rel = rel / len(izm)

	# rmserr @ root mean square error
	relerr = np.sqrt(rel)/rms2
	return relerr


def rel_err_vojko(izm, izr):
	import numpy as np
	rms1 = np.sqrt(np.mean(np.square(izm - izr)))
	rms2 = np.sqrt(np.mean(np.square(izm)))

	relerr = rms1/rms2
	return relerr


def corr_coeff_vojko(izm, izr):
	import numpy as np
	covariance = np.sum((izm - np.mean(izm))*(izr - np.mean(izr)))
	std1 = np.sqrt(np.sum((izm - np.mean(izm))**2))
	std2 = np.sqrt(np.sum((izr - np.mean(izr))**2))
	return covariance/(std1*std2)


def create_rand_dipole_quart(rad, theta, phi):
	import random
	import numpy

	dipcheck = 0
	while dipcheck == 0:
		dip = numpy.zeros(6)
		dip[0] = rad * random.uniform(-1.0, 1.0)   # rad * random.random()
		dip[1] = rad * random.uniform(-1.0, 1.0)   # theta * random.random()
		dip[2] = rad * random.uniform(-1.0, 1.0)  # phi * random.random()

		dip[0:3] = trans_cart_spher(dip[0:3])
		if dip[0] <= rad and dip[1] <= theta and dip[2] <= phi and dip[2] > 0:
			dip[0:3] = trans_spher_cart(dip[0:3])
			dipcheck = 1

	dip[3] = random.uniform(-1.0, 1.0)
	dip[4] = random.uniform(-1.0, 1.0)
	dip[5] = random.uniform(-1.0, 1.0)

	norma = numpy.sqrt(dip[3]**2 + dip[4]**2 + dip[5]**2)

	dip[3] /= norma
	dip[4] /= norma
	dip[5] /= norma

	return dip


def create_rand_dipole_sphere(x0, y0, z0, rad_min, rad_max):
	import random
	import numpy

	dipcheck = 0
	while dipcheck == 0:
		dip = numpy.zeros(6)
		dip[0] = rad_max * random.uniform(-1.0, 1.0)   # rad * random.random()
		dip[1] = rad_max * random.uniform(-1.0, 1.0)   # theta * random.random()
		dip[2] = rad_max * random.uniform(-1.0, 1.0)  # phi * random.random()

		if rad_min <= numpy.sqrt(dip[0]**2 + dip[1]**2 + dip[2]**2) <= rad_max:
			dipcheck = 1
			premik = [x0, y0, z0]
			dip[0:3] = dip[0:3] + premik

	dip[3] = random.uniform(-1.0, 1.0)
	dip[4] = random.uniform(-1.0, 1.0)
	dip[5] = random.uniform(-1.0, 1.0)
	norma = numpy.sqrt(dip[3]**2 + dip[4]**2 + dip[5]**2)
	dip[3] /= norma
	dip[4] /= norma
	dip[5] /= norma

	return dip


def dist_two_points(point1, point2):
	import numpy as np
	dist = np.sqrt((point1[0]-point2[0])**2+(point1[1]-point2[1])**2+(point1[2]-point2[2])**2)
	return dist


def arrange_leftright_dip(dipole):
	import numpy as np
	if len(dipole) > 7:
		if dipole[0] < dipole[6]:
			leftdipole = dipole[0:6]
			rightdipole = dipole[6:12]
			dipole = np.hstack((leftdipole, rightdipole))
		else:
			leftdipole = dipole[0:6]
			rightdipole = dipole[6:12]
			dipole = np.hstack((leftdipole, rightdipole))

		return dipole, leftdipole, rightdipole
	else:
		return dipole


def dotproduct(v1, v2):
	import math
	return sum((a*b) for a, b in zip(v1, v2))


def length(v):
	import math
	return math.sqrt(dotproduct(v, v))


def angle(v1, v2):
	import math
	return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


def rotate_via_numpy(xx, yy, radians):
	import numpy as np
	"""Use numpy to build a rotation matrix and take the dot product."""
	c, s = np.cos(radians), np.sin(radians)
	j = np.array([[c, s], [-s, c]])
	mm = np.dot(j, [xx, yy])
	return mm[0], mm[1]


def check_nearest_neighbours(lists, point):
	import numpy as np
	dx = lists[:, 0] - point[0]
	dy = lists[:, 1] - point[1]
	dz = lists[:, 2] - point[2]
	distances = np.sqrt(dx * dx + dy * dy + dz * dz)
	return distances


def create_custom_opm_holders(sensorholders, circles, directions, data_path, save = "nosave"):
	import math
	import numpy as np

	minz = min(sensorholders[:, 2])
	maxz = max(sensorholders[:, 2])
	dz = (maxz - minz)/(circles-1)

	x0 = np.mean(sensorholders[:, 0])
	y0 = np.mean(sensorholders[:, 1])

	a = (sensorholders[:, 0]) ** 2
	b = (sensorholders[:, 1]) ** 2
	c = (sensorholders[:, 2]) ** 2
	lenghts = np.sqrt(a + b + c)

	dr = (3.2/4.0)*max(lenghts)

	z = minz
	indexes = np.zeros(circles)
	for j in range(circles):
		for i in range(len(sensorholders)):
			aa = z - sensorholders[i, 2]
			if abs(aa) < dz/2:
				indexes[j] += 1
		z += dz
	indexes /= 2

	z = minz
	xyz=[]
	for j in range(circles):
		dtheta =  2 * math.pi / indexes[j]
		if (j % 2) == 0:
			theta = dtheta/2
		else:
			theta = 0

		for i in range(int(indexes[j])):
			dl = math.sqrt(dr*dr - 1.18*(z-minz)*(z-minz))
			xx = dl*math.cos(theta) + x0
			yy = dl*math.sin(theta) + y0
			zz = z
			s_xyz = np.array([x0, y0, minz])-np.array([xx, yy, zz])
			s_xyz = s_xyz/np.linalg.norm(s_xyz)

			s_xyz_p = np.array((-s_xyz[1], s_xyz[0], 0))
			s_xyz_p = s_xyz_p / np.linalg.norm(s_xyz_p)

			xyz.append([dl * math.cos(theta) + x0, dl * math.sin(theta) + y0, z, s_xyz[0], s_xyz[1], s_xyz[2]])
			xyz.append([dl * math.cos(theta) + x0, dl * math.sin(theta) + y0, z, s_xyz_p[0], s_xyz_p[1], s_xyz_p[2]])
			theta += dtheta
		z += dz

	if save == "save":
		np.savetxt(data_path + "uniform_opm_holders.txt", xyz)
	return np.array(xyz)


def create_custom_opm_holders_elipsoid(sensorholders, circles, directions, save_path, save = "nosave"):
	import math
	import numpy as np

	minz = min(sensorholders[:, 2])
	maxz = max(sensorholders[:, 2])
	dz = (maxz - minz)/(circles-1)

	x0 = np.mean(sensorholders[:, 0])
	y0 = np.mean(sensorholders[:, 1])

	z = minz
	indexes = np.zeros(circles)
	for j in range(circles):
		for i in range(len(sensorholders)):
			aa = z - sensorholders[i, 2]
			if abs(aa) < dz/2:
				indexes[j] += 1
		z += dz
	indexes /= 2

	a = (sensorholders[:, 0]) ** 2
	b = (sensorholders[:, 1]) ** 2
	c = (sensorholders[:, 2]) ** 2
	lenghts = np.sqrt(a + b + c)

	dr = (3.2/4.0)*max(lenghts)
	
	ratio = 1.15
	scalling = 95.0
	
	coefs = (scalling, scalling*ratio, scalling)  # Coefficients in a0/c x**2 + a1/c y**2 + a2/c z**2 = 1 

	z = minz
	xyz=[]
	for j in range(circles):
		dtheta =  2 * math.pi / indexes[j]
		if (j % 2) == 0:
			theta = dtheta/2
		else:
			theta = 0

		for i in range(int(indexes[j])):
			dl = math.sqrt(dr*dr - 1.18*(z-minz)*(z-minz))
			
			# (this is the equation of an ellipsoid):
			phi = np.arccos((z-minz)/coefs[2])
			print(phi)
			xx = coefs[0] * np.sin(phi) * np.cos(theta) + x0
			yy = coefs[1] * np.sin(phi) * np.sin(theta) + y0
			zz = z

			s_xyz = np.array([x0, y0, minz])-np.array([xx, yy, zz])
			s_xyz = s_xyz/np.linalg.norm(s_xyz)

			s_xyz_p = np.array((-s_xyz[1], s_xyz[0], 0))
			s_xyz_p = s_xyz_p / np.linalg.norm(s_xyz_p)

			xyz.append([xx, yy, zz, s_xyz[0], s_xyz[1], s_xyz[2]])
			xyz.append([xx, yy, zz, s_xyz_p[0], s_xyz_p[1], s_xyz_p[2]])
			theta += dtheta
		z += dz

	if save == "save":
		np.savetxt(save_path + "uniform_opm_holders_elipsoid.txt", xyz)
	return np.array(xyz)


def map_et_coord(xx, yy, zz):
	import  numpy as np
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

	return np.array([xx, yy])


def rigid_transform_3D(A, B):
	# Input: expects 3xN matrix of points
	# Returns R,t
	# R = 3x3 rotation matrix
	# t = 3x1 column vector

	import numpy as np
	from math import sqrt
	assert len(A) == len(B)

	num_rows, num_cols = A.shape;

	if num_rows != 3:
		raise Exception("matrix A is not 3xN, it is {}x{}".format(num_rows, num_cols))

	[num_rows, num_cols] = B.shape;
	if num_rows != 3:
		raise Exception("matrix B is not 3xN, it is {}x{}".format(num_rows, num_cols))

	# find mean column wise
	centroid_A = np.mean(A, axis=1)
	centroid_B = np.mean(B, axis=1)

	# subtract mean
	Am = np.copy(A)
	Bm = np.copy(B)
	for i in range(num_cols):
		Am[:,i] = A[:,i] - centroid_A
		Bm[:,i] = B[:,i] - centroid_B

	# dot is matrix multiplication for array
	H = np.dot(Am, np.transpose(Bm))

	# find rotation
	U, S, Vt = np.linalg.svd(H)
	R = np.dot(Vt.T, U.T)

	# special reflection case
	if np.linalg.det(R) < 0:
		print("det(R) < R, reflection detected!, correcting for it ...\n");
		Vt[2,:] *= -1
		R = np.dot(Vt.T, U.T)

	t = -np.dot(R,centroid_A) + centroid_B

	return R, t


def create_rot_matrix(v1, v2):
	import numpy as np

	a, b = (v1 / np.linalg.norm(v1)).reshape(3), (v2 / np.linalg.norm(v2)).reshape(3)

	I = np.identity(3)
	v = np.cross(a, b)
	s = np.linalg.norm(v)
	c = np.dot(a, b)
	vx = np.array([[0., -v[2], v[1]], [v[2], 0., -v[0]], [-v[1], v[0], 0.]])
	rot_matrix = I + vx + np.dot(vx, vx) * (1.0 - c) / (s * s)

	return rot_matrix
