def recc_corrected_transformed_new(corrected, transformed, time, time2=None, scaling=None):
    import megtools.vector_functions as vfun
    import numpy as np

    corrected = corrected.copy()
    transformed = transformed.copy()

    # corrected = corrected.resample(500.0)
    # corrected = corrected.crop(tmin=time_crop[0], tmax=time_crop[1], include_tmax=True)
    bads = corrected.info['bads']

    # transformed = transformed.resample(500.0)
    # transformed = transformed.crop(tmin=time_crop[0], tmax=time_crop[1], include_tmax=True)
    transformed.info['bads'] = bads

    corr_names = corrected.info["ch_names"]
    tran_names = transformed.info["ch_names"]

    channels_in_both = set(corr_names).intersection(tran_names)

    corrected.pick_channels(channels_in_both)
    transformed.pick_channels(channels_in_both)

    corrected.pick_types(meg=True, exclude="bads")
    transformed.pick_types(meg=True, exclude="bads")

    mag_corrected = corrected.data
    mag_transformed = transformed.data

    if isinstance(time, list):
        max_i_time = np.argsort(abs(corrected.times - max(time)))[0]
        min_i_time = np.argsort(abs(corrected.times - min(time)))[0]
        i_time = np.arange(min_i_time, max_i_time + 1)

        if time2 is not None:
            max_i_time = np.argsort(abs(transformed.times - max(time2)))[0]
            min_i_time = np.argsort(abs(transformed.times - min(time2)))[0]
            i_time2 = np.arange(min_i_time, max_i_time)
        else:
            i_time2 = i_time

        mag_transformed = mag_transformed[:, i_time2]
        mag_transformed = np.average(mag_transformed, axis=1)
        mag_corrected = mag_corrected[:, i_time]
        mag_corrected = np.average(mag_corrected, axis=1)

    elif isinstance(time, float):
        i_time = np.argsort(abs(corrected.times - time))[0]
        if time2 is not None:
            i_time2 = np.argsort(abs(transformed.times - time2))[0]
        else:
            i_time2 = i_time
        mag_transformed = mag_transformed[:, i_time2]
        mag_corrected = mag_corrected[:, i_time]

    if scaling is not None:
        mag_transformed = mag_transformed*scaling

    # def rel_test(izr, izm):
    # 	rel = np.sum((izm - izr) ** 2)
    # 	rel = rel / np.sum(izm ** 2)
    # 	return rel

    # print(rel_test(mag_corrected, mag_transformed))
    print(vfun.rel_err_vojko(mag_corrected, mag_transformed))
    print(vfun.corr_coeff_vojko(mag_corrected, mag_transformed))
    return vfun.rel_err_vojko(mag_corrected, mag_transformed), vfun.corr_coeff_vojko(mag_corrected, mag_transformed)
