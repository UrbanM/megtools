def source_detection(mag, xyz1, xyz2, estimates, num, bad, system):
	# ####### system = "grad" : gradiometer
	# ####### system = "mag" : magnetometer

	import numpy as np
	from scipy import optimize
	import megtools.vector_functions as vfun

	# type = grad, mag

	xyz1 = np.array(xyz1, dtype=np.float32)
	if system == "grad":
		xyz2 = np.array(xyz2, dtype=np.float32)

	mag = np.array(mag, dtype=np.float32)
	x0 = np.array(estimates, dtype=np.float32)

	x0[0:6] = vfun.rm_radial_component(x0[0:6])
	if num == 2:
		x0[6:12] = vfun.rm_radial_component(x0[6:12])

	ii = np.argmax(xyz1[:, 0])
	a = (xyz1[:, 0] - xyz1[ii, 0]) ** 2
	b = (xyz1[:, 1] - xyz1[ii, 1]) ** 2
	c = (xyz1[:, 2] - xyz1[ii, 2]) ** 2

	lenghts = np.sqrt(a + b + c)
	nearest = np.argsort(lenghts)[:]
	bad_channels = bad

	# nearest = np.delete(nearest, bad_channels, 0)
	# nearest = nearest[0:40]


	# if system == "grad":
	# 	argument = (xyz1[:, :], xyz2[:, :], mag[:], num, 1, 1, "grad", None)

	if system == "grad":
		argument = (xyz1[nearest, :], xyz2[nearest, :], mag[nearest], num, 1, 1, "grad")

	elif system == "mag":
		argument = (xyz1[nearest, :], None, mag[nearest], num, 1, 1, "mag")
	else:
		print("Wrong system")
		return

	result_1 = optimize.leastsq(
		magfield,
		x0,
		args=argument
	)

	result = result_1[0]

	result[0:6] = vfun.rm_radial_component(result[0:6])
	if num == 2:
		result[6:12] = vfun.rm_radial_component(result[6:12])

	return result


def magfield(x0, mer1, mer2, pol, num, res, model, system):
	# ####### model = 0 : infinite conductor (Maxwell)
	# ####### model = 1 : spherical model (Sarwas)

	# ####### system = "grad" : gradiometer
	# ####### system = "mag" : magnetometer

	# ####### res = 1 : calculate ressiduals
	# ####### res = 0 : won't calculate ressiduals

	def magfieldinfinite():
		import numpy as np

		if system == "grad":
			isys = 1
		else:
			isys = 0

		for ii in range(0, num, 1):
			i = 0
			while i <= isys:
				if i == 0:
					mer = mer1
				else:
					mer = mer2

				mu0 = 1.25663706 * 10.0 ** (-6.0)

				aa1 = mer[:, 0] - x0[0]
				aa2 = mer[:, 1] - x0[1]
				aa3 = mer[:, 2] - x0[2]

				norm_aa = aa1 * aa1 + aa2 * aa2 + aa3 * aa3
				norm_aa = np.sqrt(norm_aa)
				norm_aa = norm_aa * norm_aa * norm_aa

				kross_q_a1 = (aa3 * x0[4] - aa2 * x0[5]) / norm_aa
				kross_q_a2 = (aa3 * x0[3] - aa1 * x0[5]) / norm_aa
				kross_q_a3 = (aa2 * x0[3] - aa1 * x0[4]) / norm_aa

				if i == 0:
					skalarno_smermerilnika1 = mu0 * (kross_q_a1 * mer[:, 3] + kross_q_a2 * mer[:, 4] + kross_q_a3 * mer[:, 5])
					skalarno_smermerilnika2 = 0.0
				else:
					skalarno_smermerilnika2 = mu0 * (kross_q_a1 * mer[:, 3] + kross_q_a2 * mer[:, 4] + kross_q_a3 * mer[:, 5])
				i += 1

			if ii == 0:
				razlika1 = skalarno_smermerilnika1 - skalarno_smermerilnika2
				razlika2 = 0.0
			if ii == 1:
				razlika2 = skalarno_smermerilnika1 - skalarno_smermerilnika2

			vsota = razlika1 + razlika2
		return vsota

	def magfieldspherical():
		import numpy as np
		mu0 = 1.25663706 * 10.0 ** (-6.0)

		if system == "grad":
			isys = 1
		else:
			isys = 0

		for ii in range(0, num, 1):
			i = 0
			while i <= isys:
				if i == 0:
					mer = mer1
				else:
					mer = mer2

				aa1 = mer[:, 0] - x0[0 + (ii * 6)]
				aa2 = mer[:, 1] - x0[1 + (ii * 6)]
				aa3 = mer[:, 2] - x0[2 + (ii * 6)]

				norm_aa = aa1 * aa1 + aa2 * aa2 + aa3 * aa3
				norm_aa = np.sqrt(norm_aa)

				rr1 = mer[:, 0]
				rr2 = mer[:, 1]
				rr3 = mer[:, 2]
				norm_rr = rr1 * rr1 + rr2 * rr2 + rr3 * rr3
				norm_rr = np.sqrt(norm_rr)

				skalar_aa_rs = aa1 * rr1 + aa2 * rr2 + aa3 * rr3
				f = norm_aa * (norm_rr * norm_aa + skalar_aa_rs)

				gradF1 = (f / (norm_aa ** 2) + norm_aa + norm_rr) * aa1 + ((norm_aa ** 2) / norm_rr + norm_aa) * rr1
				gradF2 = (f / (norm_aa ** 2) + norm_aa + norm_rr) * aa2 + ((norm_aa ** 2) / norm_rr + norm_aa) * rr2
				gradF3 = (f / (norm_aa ** 2) + norm_aa + norm_rr) * aa3 + ((norm_aa ** 2) / norm_rr + norm_aa) * rr3

				KONST = (mu0 / (4 * np.pi)) / (f * f)

				kross_q_r0s1 = x0[4 + (ii * 6)] * x0[2 + (ii * 6)] - x0[1 + (ii * 6)] * x0[5 + (ii * 6)]
				kross_q_r0s2 = -x0[3 + (ii * 6)] * x0[2 + (ii * 6)] + x0[0 + (ii * 6)] * x0[5 + (ii * 6)]
				kross_q_r0s3 = x0[3 + (ii * 6)] * x0[1 + (ii * 6)] - x0[0 + (ii * 6)] * x0[4 + (ii * 6)]

				mixed = kross_q_r0s1 * rr1 + kross_q_r0s2 * rr2 + kross_q_r0s3 * rr3

				B1 = (f * kross_q_r0s1 - mixed * gradF1)
				B2 = (f * kross_q_r0s2 - mixed * gradF2)
				B3 = (f * kross_q_r0s3 - mixed * gradF3)

				if i == 0:
					skalarno_smermerilnika1 = KONST * (B1 * mer[:, 3] + B2 * mer[:, 4] + B3 * mer[:, 5])
					skalarno_smermerilnika2 = 0.0
				else:
					skalarno_smermerilnika2 = KONST * (B1 * mer[:, 3] + B2 * mer[:, 4] + B3 * mer[:, 5])
				i += 1
			if ii == 0:
				razlika1 = skalarno_smermerilnika1 - skalarno_smermerilnika2
				razlika2 = 0.0
			if ii == 1:
				razlika2 = skalarno_smermerilnika1 - skalarno_smermerilnika2

		vsota = razlika1 + razlika2
		return vsota

	if model == 0:
		result = magfieldinfinite()
	elif model == 1:
		result = magfieldspherical()

	if res == 1:
		result = result - pol

	return result


def calculate_leadfield(source_space, meters, meters1, model, system):
	import numpy

	# ####### model = 0 : infinite conductor (Maxwell)
	# ######## model = 1 : spherical model (Sarwas)

	def calculate_basic_field(x0, mer):
		mu0 = 1.25663706 * 10.0 ** (-6.0)

		aa1 = mer[0] - x0[0]
		aa2 = mer[1] - x0[1]
		aa3 = mer[2] - x0[2]

		norm_aa = aa1 * aa1 + aa2 * aa2 + aa3 * aa3
		norm_aa = numpy.sqrt(norm_aa)
		norm_aa = norm_aa * norm_aa * norm_aa

		kross_a_p1 = (aa2 * mer[5] - aa3 * mer[4]) / norm_aa
		kross_a_p2 = (-aa1 * mer[5] + aa3 * mer[3]) / norm_aa
		kross_a_p3 = (aa1 * mer[4] - aa2 * mer[3]) / norm_aa

		# vec_merilnika = mu0 * numpy.hstack((kross_a_p1, kross_a_p2, kross_a_p3))
		return mu0 * kross_a_p1, mu0 * kross_a_p2, mu0 * kross_a_p3

	def calculate_spherical_field(x0, mer):
		import numpy as np
		mu0 = 1.25663706 * 10.0 ** (-6.0)

		aa1 = mer[0] - x0[0]
		aa2 = mer[1] - x0[1]
		aa3 = mer[2] - x0[2]

		norm_aa = aa1 * aa1 + aa2 * aa2 + aa3 * aa3
		norm_aa = np.sqrt(norm_aa)

		rr1 = mer[0]
		rr2 = mer[1]
		rr3 = mer[2]

		mp1 = mer[3]
		mp2 = mer[4]
		mp3 = mer[5]

		rp1 = x0[0]
		rp2 = x0[1]
		rp3 = x0[2]

		norm_rr = rr1 * rr1 + rr2 * rr2 + rr3 * rr3
		norm_rr = np.sqrt(norm_rr)

		skalar_aa_rr = aa1 * rr1 + aa2 * rr2 + aa3 * rr3
		f = norm_aa * (norm_rr * norm_aa + skalar_aa_rr)

		gradF1 = ((norm_aa * norm_aa / norm_rr) + (skalar_aa_rr / norm_aa) + (2 * norm_aa) + (2 * norm_rr)) * rr1 - \
		         (norm_aa + (2 * norm_rr) + (skalar_aa_rr / norm_aa)) * x0[0]
		gradF2 = ((norm_aa * norm_aa / norm_rr) + (skalar_aa_rr / norm_aa) + (2 * norm_aa) + (2 * norm_rr)) * rr2 - \
		         (norm_aa + (2 * norm_rr) + (skalar_aa_rr / norm_aa)) * x0[1]
		gradF3 = ((norm_aa * norm_aa / norm_rr) + (skalar_aa_rr / norm_aa) + (2 * norm_aa) + (2 * norm_rr)) * rr3 - \
		         (norm_aa + (2 * norm_rr) + (skalar_aa_rr / norm_aa)) * x0[2]

		KONST = (mu0 / (4 * np.pi)) / (f * f)

		kross_rp_mp1 = rp2 * mp3 - rp3 * mp2
		kross_rp_mp2 = - rp1 * mp3 + rp3 * mp1
		kross_rp_mp3 = rp1 * mp2 - rp2 * mp1

		skalar_gradF_mp = gradF1 * mp1 + gradF2 * mp2 + gradF3 * mp3

		kross_rp_rr1 = rp2 * rr3 - rp3 * rr2
		kross_rp_rr2 = - rp1 * rr3 + rp3 * rr1
		kross_rp_rr3 = rp1 * rr2 - rp2 * rr1

		B1 = (f * kross_rp_mp1 - skalar_gradF_mp * kross_rp_rr1)
		B2 = (f * kross_rp_mp2 - skalar_gradF_mp * kross_rp_rr2)
		B3 = (f * kross_rp_mp3 - skalar_gradF_mp * kross_rp_rr3)

		return KONST * B1, KONST * B2, KONST * B3

	lead_field = numpy.zeros((meters.shape[0], source_space.shape[0], 3))
	for j in range(meters.shape[0]):
		for i in range(source_space.shape[0]):
			if system == "OPM":
				if model == 0:
					lead_field[j, i, :] = calculate_basic_field(source_space[i, :], meters[j, :])
				elif model == 1:
					lead_field[j, i, :] = calculate_spherical_field(source_space[i, :], meters[j, :])
				else:
					print("Nisi izbral pravi model, dodaj 0 ali 1")
			elif system == "SQUID":
				if model == 0:
					lead_field1 = calculate_basic_field(source_space[i, :], meters[j, :])
					lead_field2 = calculate_basic_field(source_space[i, :], meters1[j, :])
					lead_field[j, i, :] = numpy.subtract(lead_field1, lead_field2)
				elif model == 1:
					lead_field1 = calculate_spherical_field(source_space[i, :], meters[j, :])
					lead_field2 = calculate_spherical_field(source_space[i, :], meters1[j, :])
					lead_field[j, i, :] = numpy.subtract(lead_field1, lead_field2)
				else:
					print("Nisi izbral pravi model, dodaj 0 ali 1")
			else:
				print("Nisi izbral pravega sistema, dodaj OPM ali SQUID")

	return lead_field


def calculate_lead_matrix(source_space, lead_field, reg_value=1/40):
	import numpy
	
	lead_field_matrix = numpy.zeros((lead_field.shape[0], lead_field.shape[0]))
	for i in range(lead_field.shape[0]):
		for j in range(lead_field.shape[0]):
			for k in range(source_space.shape[0]):
				lead_field_matrix[i, j] += numpy.dot(lead_field[i, k, :], lead_field[j, k, :])

	U, D, VT = numpy.linalg.svd(lead_field_matrix[:, :])

	UT = numpy.transpose(U)
	V = numpy.transpose(VT)

	DI = numpy.linalg.inv(D * numpy.identity(len(D)))

	for i in range(len(D)):
		if D[i] < reg_value * max(D):
			DI[i, :] = 0.0

	VDIUT_matrix = numpy.dot(numpy.dot(V, DI), UT)

	return lead_field_matrix, VDIUT_matrix

def calculate_field_lead_matrix(VDIUT_matrix, lead_field_matrix, mag, lead_field):
	import numpy
	
	weight = numpy.dot(VDIUT_matrix, mag)

	p = numpy.zeros((lead_field.shape[1], 3))
	for j in range(lead_field.shape[1]):
		for k in range(mag.shape[0]):
			p[j, 0] += lead_field[k, j, 0] * weight[k]
			p[j, 1] += lead_field[k, j, 1] * weight[k]
			p[j, 2] += lead_field[k, j, 2] * weight[k]

	mag_field = numpy.dot(numpy.transpose(weight), lead_field_matrix)

	return weight, mag_field, p


def calculate_field(lead, sp_dipoles):
	import numpy
	a, b, c = numpy.shape(lead)

	mag_field = numpy.zeros(a)
	for i in range(0, a, 1):
		for j in range(0, b, 1):
			mag_field[i] += numpy.dot(lead[i, j, :], sp_dipoles[j, :])

	return mag_field


def import_sourcespace(vertices, faces, sp_anlges, sp_dangles, leftright):
	import numpy
	import megtools.vector_functions as vfun

	vertices_temp = numpy.copy(vertices)

	sourcespace = []
	chosen_faces = []

	x_povp = numpy.mean(vertices_temp[:, 0])
	y_povp = numpy.mean(vertices_temp[:, 1])
	z_povp = numpy.mean(vertices_temp[:, 2])

	vertices_temp[:, 0] = vertices_temp[:, 0] - x_povp
	vertices_temp[:, 1] = vertices_temp[:, 1] - y_povp
	vertices_temp[:, 2] = vertices_temp[:, 2] - z_povp

	min_th = sp_anlges[0] - sp_dangles[0]
	max_th = sp_anlges[0] + sp_dangles[0]
	min_phi = sp_anlges[1] - sp_dangles[1]
	max_phi = sp_anlges[1] + sp_dangles[1]

	for j, i in enumerate(faces):
		x1 = (vertices_temp[i[0], 0] + vertices_temp[i[1], 0] + vertices_temp[i[2], 0]) / 3.0
		y1 = (vertices_temp[i[0], 1] + vertices_temp[i[1], 1] + vertices_temp[i[2], 1]) / 3.0
		z1 = (vertices_temp[i[0], 2] + vertices_temp[i[1], 2] + vertices_temp[i[2], 2]) / 3.0
		if leftright == "left":
			vi = vfun.trans_cart_spher([x1, y1, z1])
			if min_th < vi[1] < max_th and min_phi < vi[2] < max_phi and x1 < 0:
				sourcespace.append([x1, y1, z1])
				chosen_faces.append(j)
		if leftright == "right":
			vi = vfun.trans_cart_spher([x1, y1, z1])
			if min_th < vi[1] < max_th and min_phi < vi[2] < max_phi and x1 >= 0:
				sourcespace.append([x1, y1, z1])
				chosen_faces.append(j)
		if leftright == "both":
			vi = vfun.trans_cart_spher([x1, y1, z1])
			if min_th < vi[1] < max_th and min_phi < vi[2] < max_phi:
				sourcespace.append([x1, y1, z1])
				chosen_faces.append(j)

	sourcespace = numpy.array(sourcespace)
	sourcespace[:, 0] = sourcespace[:, 0] + x_povp
	sourcespace[:, 1] = sourcespace[:, 1] + y_povp
	sourcespace[:, 2] = sourcespace[:, 2] + z_povp

	return sourcespace, numpy.array(chosen_faces)

def create_beamformers_sourcespace(freesurfer_path, name, rotation, translation):
	import numpy as np
	import vector_functions_v12 as vfun
	import mne
	import head_reconstruct_v12 as hrec

	img, values = hrec.import_nifti(freesurfer_path + name + "/mri/wmparc.mgz")

	M = img.affine[:3, :3]
	abc = img.affine[:3, 3]

	def f(i, j, k):
		""" Return X, Y, Z coordinates for i, j, k """
		return M.dot([i, j, k]) + abc

	img1 = values[:, :, :]
	a, b, c = img1.shape

	x_scatters = []
	y_scatters = []
	z_scatters = []
	no_point = []

	nn = 10

	for kk in range(0, c, nn):
		for jj in range(0, b, nn):
			for ii in range(0, a, nn):
				if img1[ii, jj, kk] > 0:
					x_scatters.append(f(ii, jj, kk)[0])
					y_scatters.append(f(ii, jj, kk)[1])
					z_scatters.append(f(ii, jj, kk)[2])

	xyz = [list(a) for a in zip(x_scatters, y_scatters, z_scatters)]
	xyz = np.array(xyz)

	cube = np.array((10.0, 10.0, 10.0))

	xyz_min = np.array((np.min(xyz[:, 0]), np.min(xyz[:, 1]), np.min(xyz[:, 2])))
	xyz_max = np.array((np.max(xyz[:, 0]), np.max(xyz[:, 1]), np.max(xyz[:, 2])))

	a = xyz_min[0]
	b = xyz_min[1]
	c = xyz_min[2]

	xyz_temp = np.array((a, b, c))

	source_space=[]

	while xyz_temp[0] <= xyz_max[0]:
		xyz_temp[1] = b
		while xyz_temp[1] <= xyz_max[1]:
			xyz_temp[2] = c
			while xyz_temp[2] <= xyz_max[2]:
				da_ne = 0
				lenx = len(xyz)
				i=0
				while(i<lenx):
					if(xyz_temp[0] - (cube[0] / 2.0)) < xyz[i,0] < (xyz_temp[0] + (cube[0] / 2.0)):
						if (xyz_temp[1] - (cube[1] / 2.0)) < xyz[i,1] < (xyz_temp[1] + (cube[1] / 2.0)):
							if (xyz_temp[2] - (cube[2] / 2.0)) < xyz[i,2] < (xyz_temp[2] + (cube[2] / 2.0)):
								xyz = np.delete(xyz,i,0)
								i-=1
								lenx-=1
								da_ne = 1
					i+=1
				if da_ne == 1:
					source_space.append([xyz_temp.item(0), xyz_temp.item(1), xyz_temp.item(2)])
				xyz_temp[2] += cube[2]
			xyz_temp[1] += cube[1]
		xyz_temp[0] += cube[0]

	source_space = np.array(source_space)

	source_space[:, 1], source_space[:, 2] = vfun.rotate_via_numpy(source_space[:, 1], source_space[:, 2], np.radians(rotation[0]))
	source_space[:, 0], source_space[:, 1] = vfun.rotate_via_numpy(source_space[:, 0], source_space[:, 1], np.radians(rotation[2]))

	source_space[:, 0] = source_space[:, 0] + translation[0]
	source_space[:, 1] = source_space[:, 1] + translation[1]
	source_space[:, 2] = source_space[:, 2] + translation[2]

	return source_space

def calculate_covariancematrix(mag):
	import numpy as np
	num = mag.shape[1]
	sumag = np.sum(mag, axis=1)/float(num)

	cov_matrix = np.matmul(np.matrix(mag[:, 0]-sumag).T, np.matrix(mag[:, 0]-sumag))
	for i in range(1,num):
		cov_matrix = cov_matrix + np.matmul(np.matrix(mag[:, i]-sumag).T, np.matrix(mag[:, i]-sumag))

	cov_matrix = (1.0/(float(num)-1.0))*cov_matrix
	return cov_matrix





