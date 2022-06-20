def evokedmaps_zeros(evokedmaps, LSM):
    import numpy as np
    lsm_evoked = evokedmaps.copy()

    unchosen_indx = np.zeros(len(LSM.unchosen), dtype=int)
    for i in range(len(LSM.unchosen)):
        unchosen_indx[i] = np.where(np.array(evokedmaps.names) == LSM.unchosen[i])[0]

    for j, i in enumerate(unchosen_indx):
        lsm_evoked.data[i, :] = np.zeros(np.shape(lsm_evoked.data[i, :]))
    return lsm_evoked

