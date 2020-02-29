def find_events_AEF(values, times, triger_range):
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

#	plt.plot(np.array(spikes), np.ones(len(spikes)), 'o')
#	plt.plot(np.arange(0,len(values1)), values1, '-')
#	pvis.plot_channel(times, values, "TrigChannel", "signal", "time", None)
#	plt.show()
	return np.array(spikes)
