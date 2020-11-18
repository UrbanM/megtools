def plot_butterfly(evoked1, time=None, time_avg=None, realpicks=None):
	import mne
	import matplotlib.pyplot as plt
	from matplotlib import rc
	rc('text', usetex=True)

	if realpicks != None:
#		chosen = evoked.ch_names[realpicks]
		evoked = evoked1.copy()
		evoked.pick(realpicks)
	else:
		evoked = evoked1.copy()
		evoked.pick_types(meg=True, exclude='bads')

#	evoked = evoked.copy()
#	evoked.pick_types(meg=True, exclude='bads')
	
	if time != None:
		evoked_copy = evoked.copy().crop(tmin = time[0], tmax = time[1], include_tmax=True)

	mag = evoked_copy.data
	times = evoked_copy.times

	mag = mag * 10 ** (15)
	i, j = mag.shape

	fig = plt.figure(figsize=(10, 5))

	ax = fig.add_subplot(1, 1, 1)
	for ii in range(0, i, 1):
		ax.plot(times, mag[ii, :], lw=1.0, c='black')

	if time_avg != None and isinstance(time_avg, list):
		plt.axvspan(time_avg[0], time_avg[1], facecolor='r', alpha=0.5)

	ax.set_xlim(min(times), max(times))
	ax.set_xlabel('$t$[s]', fontsize=20)
	ax.set_ylabel('$B$[fT]', fontsize=20)
	plt.xticks(size=20)
	plt.yticks(size=20)

	fig.tight_layout()

	return fig


def map_et_coord(xx, yy, zz):
	import numpy as np
	# projected helmet coordinates onto a plane

	cc = np.where((xx.any() == 0.0) and (yy.any() == 0.0))

	x2 = xx * xx
	z2 = zz * zz
	y2s = (yy - 0.13526) ** (2.0)

	fact= 0.1 *np.sqrt((x2+z2+y2s)/(x2+z2))

	xx = -xx * fact
	yy = zz * fact

	if cc[0] != -1:
		xx[cc] = 0.0
		yy[cc] = 0.0

	return xx, yy


def plotxyz(evoked, name=None, freesurfer_path=None):
	import numpy as np
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	
	unit_v = np.array([0.0, 0.0, 1.0])
	xyz = np.zeros((len(evoked.ch_names), 12))
	vector = np.zeros((len(evoked.ch_names), 3))
	for i in range(len(evoked.ch_names)):
		xyz[i,:] = evoked.info['chs'][i]['loc']
		rot_matrix = np.array(([xyz[i,3], xyz[i,4], xyz[i,5]], [xyz[i,6], xyz[i,7], xyz[i,8]], [xyz[i,9], xyz[i,10], xyz[i,11]]))
		vector[i,:] = np.dot(unit_v, rot_matrix)

	if name != None and freesurfer_path != None:
		import mne
		xyzh = mne.read_surface(freesurfer_path +"/"+name+"/"+"bem/watershed/"+ name +"_outer_skull_surface")[0]/1000.0
		ax.scatter(xyzh[:, 0], xyzh[:, 1], xyzh[:, 2], s=4, c='green', alpha=0.1)

	ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c='gray', alpha=0.5)
	ax.quiver(xyz[:, 0], xyz[:, 1], xyz[:, 2], vector[:, 0], vector[:, 1], vector[:, 2], length=0.01)

	X = xyz[:, 0]
	Y = xyz[:, 1]
	Z = xyz[:, 2]

	max_range = np.array([X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]).max() / 2.0

	mid_x = (X.max() + X.min()) * 0.5
	mid_y = (Y.max() + Y.min()) * 0.5
	mid_z = (Z.max() + Z.min()) * 0.5
	ax.set_xlim(mid_x - max_range, mid_x + max_range)
	ax.set_ylim(mid_y - max_range, mid_y + max_range)
	ax.set_zlim(mid_z - max_range, mid_z + max_range)
	ax.axis('off')
	
	return

def plot_topo_v2(evoked, data_path, block_name, system, time, subject_dir, realpicks = None, halve=False):
	import mne
	import matplotlib.pyplot as plt
	import matplotlib.cbook as cbook
	import numpy as np
	import megtools.pyread_biosig as pbio
	import megtools.pymeg_visualize as pvis
	import megtools.vector_functions as vfun
	from scipy.interpolate import griddata
	import matplotlib.mlab as ml
	from mpl_toolkits.mplot3d import Axes3D
	from numpy import linalg as LA

	from matplotlib import rc
	# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
	## for Palatino and other serif fonts use:
	# rc('font',**{'family':'serif','serif':['Palatino']})
	rc('text', usetex=True)



	# import os
	# os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2015/bin/x86_64-darwin'
	# print(os.getenv("PATH"))

	if realpicks != None:
#		chosen = evoked.ch_names[realpicks]
		evoked = evoked.copy()
		evoked.pick(realpicks)
	else:
		evoked = evoked.copy()
		evoked.pick_types(meg=True, exclude='bads')

	if isinstance(time, list):
		max_i_time = np.argsort(abs(evoked.times - max(time)))[0]
		min_i_time = np.argsort(abs(evoked.times - min(time)))[0]
		i_time = np.arange(min_i_time, max_i_time+1)
	elif isinstance(time, float):
		i_time = np.argsort(abs(evoked.times - time))[0]

	if system == "OPM":
		sens_all_path = data_path + "/sensorholders/" + block_name[0:4] + "_sensor_pos_ori_all.txt"
		xyz_all = pbio.imp_sensor_holders(sens_all_path)
		xyz_all[:, 0:3] = xyz_all[:, 0:3] / 1000.0
		xyz_all[:, 3:6] = -xyz_all[:, 3:6]

		rotation, translation = pbio.import_opm_trans(data_path + "/MEG/opm_trans.txt", block_name[0:4] )
		translation = translation /1000.0

		unit_v = np.array([0.0, 0.0, 1.0])
		xyz_occupied = np.zeros((len(evoked.ch_names), 6))
		channelnames = []
		for i in range(len(evoked.ch_names)):
			channelnames.append(evoked.info['chs'][i]['ch_name'])
			xyz = evoked.info['chs'][i]['loc']
			rot_matrix = np.array(([xyz[3], xyz[4], xyz[5]], [xyz[6], xyz[7], xyz[8]], [xyz[9], xyz[10], xyz[11]]))
			vector = np.dot(unit_v, rot_matrix)
			xyz_occupied[i, 0:3] = xyz[0:3]
			xyz_occupied[i, 3:6] = vector
		channelnames = np.array(channelnames)

		surf = mne.read_surface(subject_dir + block_name[0:4]  + "/surf/" + "lh.white", read_metadata=True)
		xyz_occupied[:, 0:3] = xyz_occupied[:, 0:3] + (surf[2]['cras']/1000.0)

		xyz_occupied[:, 0] = xyz_occupied[:, 0] + translation[0]
		xyz_occupied[:, 1] = xyz_occupied[:, 1] + translation[1]
		xyz_occupied[:, 2] = xyz_occupied[:, 2] + translation[2]

		xyz_occupied[:, 1], xyz_occupied[:, 2] = vfun.rotate_via_numpy(xyz_occupied[:, 1], xyz_occupied[:, 2], np.radians(rotation[0]))
		xyz_occupied[:, 0], xyz_occupied[:, 2] = vfun.rotate_via_numpy(xyz_occupied[:, 0], xyz_occupied[:, 2], np.radians(rotation[1]))
		xyz_occupied[:, 0], xyz_occupied[:, 1] = vfun.rotate_via_numpy(xyz_occupied[:, 0], xyz_occupied[:, 1], np.radians(rotation[2]))
		xyz_occupied[:, 4], xyz_occupied[:, 5] = vfun.rotate_via_numpy(xyz_occupied[:, 4], xyz_occupied[:, 5], np.radians(rotation[0]))
		xyz_occupied[:, 3], xyz_occupied[:, 5] = vfun.rotate_via_numpy(xyz_occupied[:, 3], xyz_occupied[:, 5], np.radians(rotation[1]))
		xyz_occupied[:, 3], xyz_occupied[:, 4] = vfun.rotate_via_numpy(xyz_occupied[:, 3], xyz_occupied[:, 4], np.radians(rotation[2]))

#		fig = plt.figure()
#		ax = fig.gca(projection='3d')
#		ax.scatter(xyz_occupied[:, 0], xyz_occupied[:, 1], xyz_occupied[:, 2], c='gray', alpha=0.5)
#		ax.quiver(xyz_occupied[:, 0], xyz_occupied[:, 1], xyz_occupied[:, 2], xyz_occupied[:, 3], xyz_occupied[:, 4], xyz_occupied[:, 5], length=0.01)
#		ax.scatter(xyz_all[:, 0], xyz_all[:, 1], xyz_all[:, 2], c='gray', alpha=0.5)
#		ax.quiver(xyz_all[:, 0], xyz_all[:, 1], xyz_all[:, 2], xyz_all[:, 3], xyz_all[:, 4], xyz_all[:, 5], length=0.01, color="red")
#		plt.show()

		xx, yy = map_et_coord(-xyz_all[:,0], xyz_all[:,2], xyz_all[:,1])
		xx =(np.array(xx).reshape((-1, 1)) )
		yy =(np.array(yy).reshape((-1, 1)) )
		xy = np.concatenate((xx,yy), axis=1)

		xy = pvis.squid_xy_reconstruciton(xyz_all[:, 0:3])
#		xy_orig = xy

		center = (np.mean(xy[:,0]), np.mean(xy[:,1]))
		radius = 1.2*(np.max(xy[:,0])-np.min(xy[:,0]))/2.0

		picks = []
		orients = []  # tan or rad
		signs = []  # -1 or 1, if -1 it is the new generation

		# write prittier code
		for i, j in enumerate(xyz_occupied):
			dist = np.sqrt((xyz_all[:,0]-j[0])**2+(xyz_all[:,1]-j[1])**2+(xyz_all[:,2]-j[2])**2)
			pick1 = np.argsort(dist)[0]
			pick2 = np.argsort(dist)[1]

			norm1 = LA.norm(np.cross(j[0:3], j[3:6]))
			norm2 = LA.norm(np.cross(xyz_all[pick1, 0:3], xyz_all[pick1, 3:6]))
			norm3 = LA.norm(np.cross(xyz_all[pick2, 0:3], xyz_all[pick2, 3:6]))

			if abs(norm1-norm2) < abs(norm1-norm3):
				picks.append(pick1)
				dist1 = vfun.dist_two_points(j[3:6], xyz_all[pick1, 3:6])
				dist2 = vfun.dist_two_points(-j[3:6], xyz_all[pick1, 3:6])
				if dist2 < dist1:
					signs.append(-1)
				else:
					signs.append(1)
				
				if(norm2<norm3):
					orients.append("rad")
				else:
					orients.append("tan")
			else:
				picks.append(pick2)
				dist1 = vfun.dist_two_points(j[3:6], xyz_all[pick2, 3:6])
				dist2 = vfun.dist_two_points(-j[3:6], xyz_all[pick2, 3:6])
				if dist2 < dist1:
					signs.append(-1)
				else:
					signs.append(1)
				
				if(norm2<norm3):
					orients.append("tan")
				else:
					orients.append("rad")

		rad_picks = np.where(np.array(orients) == "rad")[0].tolist()
		tan_picks = np.where(np.array(orients) == "tan")[0].tolist()
	
		xy_rad = xy[np.array(picks)[rad_picks],:]
		xy_tan = xy[np.array(picks)[tan_picks],:]
		
		mag = evoked.data[:, i_time]
		mag = mag*(10.0**15.0)
		
		if isinstance(time, list):
			mag = evoked.data[:, i_time]
			mag = np.average(mag, axis=1)
			mag = mag*(10.0**15.0)
		
		mag = mag*signs
		#print(signs)
		mag_rad = mag[rad_picks]
		mag_tan = mag[tan_picks]

		plot_both=0
		if len(rad_picks) > 0 and len(tan_picks) > 0:
			plot_both = 1
			if halve!=False:
				fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(7, 5))
			else:
				fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12, 5))

		dx = 0.0
		dy = 0.0

		if len(rad_picks) > 0:
			if plot_both==0:
				fig, ax1 = plt.subplots()
			x_min = min(xy_rad[:, 0])
			x_max = max(xy_rad[:, 0])
			y_min = min(xy_rad[:, 1])
			y_max = max(xy_rad[:, 1])
			nx, ny = 200, 200

			xi = np.linspace(x_min, x_max, nx)
			yi = np.linspace(y_min, y_max, ny)
			xi, yi = np.meshgrid(xi, yi)
			xi_rad = xi
			yi_rad = yi
			zi_rad = griddata((xy_rad[:, 0], xy_rad[:, 1]), mag_rad, (xi, yi), method='cubic')

			ax1.set(adjustable='box', aspect='equal')
			patch, xy_circle = cut_circle_patch(center, radius, ax1, 350, halve=halve)
			ax1.plot(xy_circle[:, 0], xy_circle[:, 1], '-k', linewidth=2)
			im11 = ax1.pcolormesh(xi_rad - dx, yi_rad + dx, zi_rad, cmap=plt.get_cmap('hot'))
			im21 = ax1.contour(xi_rad - dx, yi_rad + dx, zi_rad, colors="black")
			im = ax1.scatter(xy_rad[:, 0] - dx, xy_rad[:, 1] + dy, s=2, c="black")
			clb1 = fig.colorbar(im11, shrink=0.8, extend='both', ax=ax1)
			clb1.ax.tick_params(labelsize=30)
			clb1.ax.set_title('$B [\mathrm{fT}]$', fontsize=30)
			ax1.set_title("$\mathrm{radial}$", fontsize = 30)
			ax1.axis('off')

		if len(tan_picks) > 0:
			if plot_both==0:
				fig, ax2 = plt.subplots()
			x_min = min(xy_tan[:, 0])
			x_max = max(xy_tan[:, 0])
			y_min = min(xy_tan[:, 1])
			y_max = max(xy_tan[:, 1])
			nx, ny = 200, 200

			xi = np.linspace(x_min, x_max, nx)
			yi = np.linspace(y_min, y_max, ny)
			xi, yi = np.meshgrid(xi, yi)
			xi_tan = xi
			yi_tan = yi
			zi_tan = griddata((xy_tan[:, 0], xy_tan[:, 1]), mag_tan, (xi, yi), method='cubic')

			ax2.set(adjustable='box', aspect='equal')
			patch, xy_circle = cut_circle_patch(center, radius, ax2, 350, halve=halve)
			ax2.plot(xy_circle[:,0], xy_circle[:,1],'-k', linewidth=2)
			im12 = ax2.pcolormesh(xi_tan - dx, yi_tan + dx, zi_tan, cmap=plt.get_cmap('hot'))
			im22 = ax2.contour(xi_tan-dx, yi_tan+dx, zi_tan, colors="black")
			im = ax2.scatter(xy_tan[:,0]-dx, xy_tan[:,1]+dy, s=2, c="black")
			clb2 = fig.colorbar(im12, shrink=0.8, extend='both', ax=ax2)
			clb2.ax.set_title('$B [\mathrm{fT}]$', fontsize=30)
			clb2.ax.tick_params(labelsize=30)
			ax2.set_title("$\mathrm{tangential}$", fontsize=30)
			ax2.axis('off')

		fig.tight_layout()
#		plt.rc('font', family='serif')

	return fig, picks, i_time, rad_picks, tan_picks



def plot_topo(evoked, data_path, block_name, system, time, position=None, multi=None, endalign=None):
	#implement to remove bad channels
	import mne
	import matplotlib.pyplot as plt
	import matplotlib.cbook as cbook
	import numpy as np
	import megtools.pyread_biosig as pbio
	import megtools.pymeg_visualize as pvis
	import megtools.vector_functions as vfun
	from scipy.interpolate import griddata
	import matplotlib.mlab as ml
	from matplotlib import rc
	rc('text', usetex=True)

	evoked.pick_types(meg=True, exclude='bads')
	
	if isinstance(time, list):
		max_i_time = np.argsort(abs(evoked.times - max(time)))[0]
		min_i_time = np.argsort(abs(evoked.times - min(time)))[0]
		i_time = np.arange(min_i_time, max_i_time+1)
		print(i_time)
	elif isinstance(time, float):
		i_time = np.argsort(abs(evoked.times - time))[0]

	if system == "SQUID":
		lout_dict = {}
		lout = mne.channels.read_layout("ET160.lout", data_path+"/MEG/")
		lout_names = lout.names
		lout_pos = lout.pos
		for i,j in enumerate(lout_names):
			lout_dict[j]=lout_pos[i]

		xy_all = []
		xy  = []
		chs_names_evoked = evoked.ch_names
		for i,j in enumerate(lout_names):
			if j in chs_names_evoked:
				xy.append(lout_dict[j])
			xy_all.append(lout_dict[j])
		xy = np.array(xy)
		xy_all = np.array(xy_all)
		center = (np.mean(xy_all[:,0]), np.mean(xy_all[:,1]))
		radius = 1.08*(np.max(xy_all[:,0])-np.min(xy_all[:,0]))/2.0
		no_ch_all = len(xy_all)
		no_ch = len(xy)
		
		mag_avg = evoked.data[:, i_time]
		if isinstance(time, list):
			mag_avg = np.average(mag_avg, axis=1)
		
		if position == "right":
			chosen = np.argsort(xy[:, 0])[::-1][0:int(0.6*no_ch)]
			xy = xy[chosen,:]
			mag = mag_avg[chosen]
		elif position == "left":
			chosen = np.argsort(xy[:, 0])[0:int(0.6*no_ch)]
			xy = xy[chosen,:]
			mag = mag_avg[chosen]
		else:
			chosen=[]
			mag = mag_avg[:]

		mag = mag*(10.0**15.0)

		x_min = min(xy[:, 0])
		x_max = max(xy[:, 0])
		y_min = min(xy[:, 1])
		y_max = max(xy[:, 1])
		nx, ny = 200, 200

		xi = np.linspace(x_min, x_max, nx)
		yi = np.linspace(y_min, y_max, ny)
		xi, yi = np.meshgrid(xi, yi)

		zi = griddata((xy[:, 0], xy[:, 1]), mag, (xi, yi), method='cubic')

		fig, ax = plt.subplots()
		im1 = plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap('hot'))
		im2 = plt.contour(xi, yi, zi, colors="black")
		im = plt.scatter(xy[:,0], xy[:,1], s=2, c="black")
		clb = fig.colorbar(im1, shrink=0.8, extend='both')
		clb.ax.set_title('$B [\mathrm{fT}]$', fontsize = 30)
		clb.ax.tick_params(labelsize=30)
		patch, xy_circle = cut_circle_patch(center, radius, ax, 320)
		im1.set_clip_path(patch)
		plt.plot(xy_circle[:,0], xy_circle[:,1],'-k', linewidth=2)
#		plt.tight_layout()
		plt.rc('font', family='serif')
		plt.axis('off')
	return fig, i_time, chosen

def cut_circle_patch(xy_s, radius, ax, max_angle, halve=False):
	import matplotlib.patches as patches
	import numpy as np
	no_angles = 360
	dphi =(2 * np.pi)/no_angles
	phi = 0
	xy = np.zeros((no_angles, 2))
	for i in range(no_angles):
		phi = dphi*i
		if	np.cos(phi) < np.cos(np.deg2rad(max_angle)):
			xy[i, 0] = np.sin(phi)*(radius)
			xy[i, 1] = np.cos(phi)*(radius)
		else:
			xy[i, 0] = np.sin(phi)*(radius)
			xy[i, 1] = np.cos(np.deg2rad(max_angle))*(radius)

	xy[:,0] += xy_s[0]
	xy[:,1] += xy_s[1]

	if halve!=False:
		idxs = np.where(xy[:,0]>xy_s[0])[0]
		xy=xy[idxs,:]
	elif halve=="left":
		idxs = np.where(xy[:,0]<xy_s[0])[0]
		xy=xy[idxs,:]

	patch = patches.Polygon(xy, transform=ax.transData)
	return patch, xy

def calculate_topomap(evoked, system):
	import mne
	import numpy as np
	import matplotlib.pyplot as plt
	default_squid_path = data_path + "ET160_template.0100.flt.hdr"
	xyz1, xyz2, ch_info, sfreq= import_sensors(default_squid_path, system)

	print(ch_info)

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
