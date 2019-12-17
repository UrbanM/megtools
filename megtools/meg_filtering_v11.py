def butter_lowpass(cutoff, fs, order=5):
    from scipy.signal import butter
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    from scipy.signal import lfilter

    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def butter_bandpass(lowcut, highcut, fs, order=5):
    import pymeg_visualize_v14 as pvis
    from scipy.signal import butter, freqz
    import matplotlib.pyplot as plt
    import numpy as np

    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')

    w, h = freqz(b, a, worN=2000)
    # plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)
    # plt.xlim(0,70)
    # pvis.simple_plot(b, a, "filter", "signal", "time")

    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    from scipy.signal import lfilter
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def calculate_FFT(field, T):
    from scipy.fftpack import fft
    import numpy as np
    import pymeg_visualize_v13 as pvis

    N = np.shape(field)[0]
    yf = fft(field)
    xf = np.linspace(0.0, 1.0 / (2.0 * T), N/2 )

    pvis.simple_plot(xf, 2.0/N * np.abs(yf[:int(N/2)]), "FFT", "FFT", "freq" )

    return yf, xf

def channel_linearization(trig, mag, time):
    import numpy as np
    from scipy import interpolate

    points = []
    no_neigh = 20

    for i in range(no_neigh, len(trig)-no_neigh, 1):
        if trig[i] > -0.1 and trig[i] < 0.1:
            j = 1
            ok = 1
            while(j<no_neigh):
                if trig[i-j] < -0.1 or trig[i] > 0.1:
                    i = i + no_neigh - j
                    ok = 0
                    break
                j+=1
            if ok == 1:
                points.append(i)

    tck = interpolate.splrep(time[points], mag[0,points], s=0)

    mag_new = interpolate.splev(time, tck, der=0)

    return mag_new, points

def channel_linearization(trig, mag, time):
    import numpy as np
    from scipy import interpolate

    points = []
    no_neigh = 5
    points.append(0)

    for i in range(no_neigh, len(trig) -1, 1):
        if trig[i] > -0.1 and trig[i] < 0.1:
            j = 1
            ok = 1
            while(j<no_neigh):
                if trig[i-j] < -0.1 or trig[i-j] > 0.1:
                    i = i + no_neigh - j
                    ok = 0
                    break
                j+=1
            if ok == 1:
                points.append(i)
    points.append(np.shape(mag)[1]-1)


    mag_new = mag
    mag_new[0,0] = np.mean(mag[0,0:no_neigh])
    mag_new[0,points[-1]] = np.mean(mag[0, -no_neigh:])

    f = interpolate.interp1d(time[points], mag_new[0,points])
    mag_new1 = f(time)

    tck = interpolate.splrep(time[points], mag_new[0,points], s=0)

    mag_new = interpolate.splev(time, tck, der=0)

    # mag[0] = mag[0] - mag_new1
    return mag_new1, points


def average_epochs(trig, epochs):
    search = 0
    end = 0
    no_neigh = 30
    start_end = []
    for i in range(len(trig)):
        if search == 0:
            if trig[i] > 0.1 or trig[i] < -0.1:
                i_start = i
                print(i_start)
                search = 1
        else:
            if trig[i] > -0.1 and trig[i] < 0.1:
                j = 1
                ok = 1
                while (j < 10):
                    if trig[i - j] < -0.1 or trig[i - j] > 0.1:
                        ok = 0
                        break
                    j += 1

                if ok == 1:
                    if trig[i - j] < -0.1 or trig[i - j] > 0.1:
                        i_end = i - j
                        print(i_end)
                        ok = 0
                        search = 0
                        start_end.append((i_start, i_end))
    print(start_end)
    return


def Implement_Notch_Filter(time, band, freq, ripple, order, filter_type, data):
    from scipy.signal import iirfilter, lfilter
    fs   = 1/time
    nyq  = fs/2.0
    low  = freq - band/2.0
    high = freq + band/2.0
    low  = low/nyq
    high = high/nyq
    b, a = iirfilter(order, [low, high], rp=ripple, btype='bandstop',
                     analog=False, ftype=filter_type)
    filtered_data = lfilter(b, a, data)
    return filtered_data

def nearest_neighbour_filter(data,k_nearest):
    import numpy as np
    nei_avg = np.array(data)
    max_ind = len(data)
    print(data)

    for i in range(max_ind):
        i_max = i + k_nearest
        i_min = i - k_nearest
        if(i_max > max_ind): i_max = max_ind
        if(i_min < 0): i_min = 0
        nei_avg[i] = np.mean(data[i_min: i_max])
    return nei_avg
