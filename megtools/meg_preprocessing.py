def find_events_AEF(values, times, triger_range, plot=False):
	import megtools.meg_filtering as mfil
	import megtools.pymeg_visualize as pvis
	import matplotlib.pyplot as plt
	import numpy as np

	values = -values
	values1 = np.zeros(len(values))

	# n, bins, patches = plt.hist(values, 20)
	# pvis.plot_channel(times, values, "TrigChannel", "signal", "time", None)
	# plt.show()

	for j in range(0, len(values)):
		if values[j] > triger_range[0] and  values[j] < triger_range[1]:
			values1[j] = 1

	spikes = []
	spike_count = 0
	for j in range(0, len(values1)):
		if values1[j] > 0 and spike_count < 3:
			spike_count += 1
			if spike_count == 3: spikes.append([j-2, 0, 1])
		if values1[j] == 0 and spike_count != 0:
			spike_count = 0
	
	if plot==True:
		plt.plot(np.array(spikes), np.ones(len(spikes)), 'o')
		plt.plot(np.arange(0,len(values1)), values1, '-')
		pvis.plot_channel(times, values, "TrigChannel", "signal", "time", None)
		plt.show()
	return np.array(spikes)


def calculate_gfp(evoked=None, data=None, times=None):
	import numpy as np
	if evoked==None and data==None:
		print("To arguments add evoked or data+times!")
		return
	if evoked!=None and data!=None:
		print("Cant have both, use only evoked or data+times!")
		return
	
	if evoked!=None:
		import mne
		evoked1 = evoked.copy()
		evoked1.pick(picks="all",exclude="bads")
		GFP = np.std(evoked1.data, axis=0)
		return GFP, evoked1.times

	if data!=None:
		print("Not yet implemented for only data!")
		return


def find_M100_peak(evoked, prefered_time=0.10, time_range=0.02, show=False, savefig=False):
	import matplotlib.pyplot as plt
	import scipy.signal as ssig
	import numpy as np
	import megtools.pymeg_visualize as pvis

#	evoked.plot(gfp=True, show=False)

	GFP, times = calculate_gfp(evoked=evoked)

	peaks = ssig.find_peaks(GFP)[0]
	avg_of_peaks = np.average(GFP[peaks]) 

	high_peaks = peaks
	max_peak = 0
#	for i in peaks:
#		if GFP[i] > 1.0*avg_of_peaks:
#			high_peaks.append(i)
	
	first = 1
	for i in high_peaks:
		if first == 1:
			if prefered_time==None:
				max_peak = i
				first = 0
			elif abs(times[i]-prefered_time)<=time_range:
				max_peak = i
				first = 0
		else:
			if prefered_time==None:
				max_peak = i
				first = 0
			elif abs(times[i]-prefered_time)<=time_range:
				if GFP[i] > GFP[max_peak]:
					max_peak = i

	if show==True or savefig!=False:
		from matplotlib import rc

		plt.rcParams.update({
			"font.family": "serif"
			#    "font.serif": [],                    # use latex default serif font
			#    "font.sans-serif": ["DejaVu Sans"],  # use a specific sans-serif font
		})
		# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
		rc('text', usetex=True)
		GFP = GFP * 10**(15)
		fig, (ax2, ax) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6,4))
		
		ax.plot(times, GFP, c="black")
		ax.scatter(times[peaks], GFP[peaks], c="r", s=40)
		ax.set_ylabel("GFP$\,\mathrm{[fT]}$", fontsize=25) #, labelpad=20)
		ax.set_xlabel("$t \,\mathrm{[s]}$", fontsize=25)
		ax.tick_params(axis='both', which='major', labelsize=25)
		if max_peak != 0:
			ax.scatter(times[max_peak], GFP[max_peak], c="g", s=50)
		ax.yaxis.set_label_coords(-0.20,0.5)

		
		ax2.set_ylabel("$B\,\mathrm{[fT]}$", fontsize=25)#, labelpad=5)
		evoked1=evoked.copy()
		evoked1.pick("all", exclude="bads")
		ax2.plot(evoked1.times, evoked1.data.T*10**(15), c="black")
		ax2.tick_params(axis='both', which='major', labelsize=25)
		ax2.yaxis.set_label_coords(-0.20,0.5)
		plt.xlim((0.0,0.4))
		plt.tight_layout()
#		plt.savefig("test.png")

	if savefig!=False:
		plt.savefig(savefig, dpi=600)
	if show==True:
		plt.show()
	if show==True or savefig!=False:
		plt.close()

		# plt1 = pvis.simple_plot(times, GFP, yaxis="$G \,\mathrm{[fT]}$", xaxis="$t \,\mathrm{[s]}$", usetex=True, ratio=0.3, size=20, c="black", xrange=[0.0,0.4])
	return times[max_peak]

