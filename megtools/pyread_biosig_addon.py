# def create_opm_object(sensorholder_path, name, subject_dir, return_holders=0, chosen_ch=[], sfreq=500.0):
#     import megtools.vector_functions as vfun
#     import mne
#     import megtools.pyread_biosig as pbio
#     import numpy as np
# #	import edit_multikit as emul
# #	import compare_measurements as cmea
#
#     sensorholder_file = sensorholder_path + name + '_sensor_pos_ori_all.txt'
#
#     rotation, translation = pbio.import_opm_trans(sensorholder_path + "opm_trans_v2.txt", name)
#     translation = translation /1000.0
#
#     holders = pbio.imp_sensor_holders(sensorholder_file)
#     holders[:, 0:3] = holders[:, 0:3] / 1000.0
#
#     #!!!!! be very carefull what comes first.
# #	holders[:, 0] = holders[:, 0] - translation[0]
# #	holders[:, 1] = holders[:, 1] - translation[1]
# #	holders[:, 2] = holders[:, 2] - translation[2]
#
#     holders[:, 1], holders[:, 2] = vfun.rotate_via_numpy(holders[:, 1], holders[:, 2], np.radians(-rotation[0]))
#     holders[:, 0], holders[:, 2] = vfun.rotate_via_numpy(holders[:, 0], holders[:, 2], np.radians(-rotation[1]))
#     holders[:, 0], holders[:, 1] = vfun.rotate_via_numpy(holders[:, 0], holders[:, 1], np.radians(-rotation[2]))
#     holders[:, 4], holders[:, 5] = vfun.rotate_via_numpy(holders[:, 4], holders[:, 5], np.radians(-rotation[0]))
#     holders[:, 3], holders[:, 5] = vfun.rotate_via_numpy(holders[:, 3], holders[:, 5], np.radians(-rotation[1]))
#     holders[:, 3], holders[:, 4] = vfun.rotate_via_numpy(holders[:, 3], holders[:, 4], np.radians(-rotation[2]))
#
#     holders[:, 0] = holders[:, 0] - translation[0]
#     holders[:, 1] = holders[:, 1] - translation[1]
#     holders[:, 2] = holders[:, 2] - translation[2]
#
#     surf = mne.read_surface(subject_dir + name + "/surf/" + "lh.white", read_metadata=True)
#     holders[:, 0:3] = holders[:, 0:3] - (surf[2]['cras']/1000.0)
#
#     print(holders.shape)
#
#     ch_types = len(holders) * ["mag"]
#     ch_names = [str(a) + str(b) for a, b in zip(ch_types, [str(x + 1) for x in range(len(holders) + 1)])]
#     info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types=ch_types)
#     evoked = mne.EvokedArray(np.zeros((len(ch_names), 2)), info)
#
#     unit_v = np.array([0.0, 0.0, 1.0])
#     for i in range(len(ch_names)):
#         rot_mat = vfun.create_rot_matrix(holders[i, 3:6], unit_v)
#         evoked.info['chs'][i]['coil_type'] = 9999
#         evoked.info['chs'][i]['scanno'] = i + 1
#         evoked.info['chs'][i]['logno'] = i + 1
#         evoked.info['chs'][i]['kind'] = 1
#         evoked.info['chs'][i]['range'] = 1.0
#         evoked.info['chs'][i]['cal'] = 3.7000000285836165e-10
#         evoked.info['chs'][i]['unit'] = 112
#         evoked.info['chs'][i]['loc'] = np.array(
#             [holders[i, 0], holders[i, 1], holders[i, 2], rot_mat[0, 0], rot_mat[0, 1], rot_mat[0, 2],
#             rot_mat[1, 0], rot_mat[1, 1], rot_mat[1, 2], rot_mat[2, 0], rot_mat[2, 1], rot_mat[2, 2]])
#
#     if len(chosen_ch) > 0:
#         print(chosen_ch)
#         evoked.pick_channels(chosen_ch, ordered=True)
#
# #	evoked = emul.reorient_evoked(evoked)
#
#     # if name == "tisa":
#         # evoked = cmea.correct_tisa_sign(evoked)
#     if return_holders == 1:
#         holders2 = np.zeros((len(ch_names), 6))
#         for i in range(len(ch_names)):
#             a = evoked.info['chs'][i]['loc'][3:12]
#             b = np.reshape(a, (3, 3))
#             holders2[i, 0:3] = evoked.info['chs'][i]['loc'][0:3]
#             holders2[i, 3:6] = np.dot(unit_v, b)
#         return evoked.info, holders2
#     else:
#         return evoked.info


def create_opm_info_all(sensorholder_path, name, subject_dir, gen12, components=None):
    import numpy as np

    if components is None:
        components = ["rad", "tan"]

    xyz1 = import_sensor_pos_ori_all(sensorholder_path, name, subject_dir, gen12=gen12)
    xyz2 = []
    for i in range(0, len(xyz1), 2):
        xyz2.append(xyz1[i])
        xyz2.append(xyz1[i + 1])
        xyz = xyz1[i][0:3]
        dir = np.cross(xyz1[i][3:6], xyz1[i + 1][3:6])
        xyz2.append(np.hstack((xyz, dir)))
    xyz2 = np.array(xyz2).astype(float)
    xyz = xyz2.copy()
    xyz[0::3, 3:6] = -xyz[0::3, 3:6]
    xyz[1::3, 3:6] = -xyz[1::3, 3:6]

    info = create_opms(xyz, components=components)
    return info


def create_opms(holders, components=None):
    import mne
    import numpy as np
    import megtools.vector_functions as vfun

    if components is None:
        components = ["rad", "tan", "ver"]

    ch_types = int(len(holders) * (len(components) / 3)) * ["mag"]
    ch_numbers = []

    ch_numbers1 = []
    ch_numbers2 = []
    ch_numbers3 = []
    if "rad" in components:
        ch_numbers1 = range(0, len(holders), 3)
        ch_names1 = ["rad" + '{:03}'.format(a) for a in np.arange(int(len(holders) / 3))]
    if "tan" in components:
        ch_numbers2 = range(1, len(holders), 3)
        ch_names2 = ["tan" + '{:03}'.format(a) for a in np.arange(int(len(holders) / 3))]
    if "ver" in components:
        ch_numbers3 = range(2, len(holders), 3)
        ch_names3 = ["ver" + '{:03}'.format(a) for a in np.arange(int(len(holders) / 3))]

    # for i in range(len(list(ch_numbers1))):
    ch_numbers = list(ch_numbers1) + list(ch_numbers2) + list(ch_numbers3)
    ch_numbers = sorted([i for i in ch_numbers])

    ch_names = []
    ii_temp = int(len(ch_numbers) / len(components))
    print(ii_temp)
    for i in range(ii_temp):
        if "rad" in components:
            ch_names.append(ch_names1[i])
        if "tan" in components:
            ch_names.append(ch_names2[i])
        if "ver" in components:
            ch_names.append(ch_names3[i])

    # ch_names = [str(a) + str(b) for a, b in zip(ch_types, [str(x + 1) for x in range(len(holders) + 1)])]

    info = mne.create_info(ch_names=ch_names, sfreq=1, ch_types=ch_types)

    unit_v = np.array([0.0, 0.0, 1.0])
    for j, i in enumerate(ch_numbers):
        rot_mat = vfun.create_rot_matrix(holders[i, 3:6], unit_v)
        info['chs'][j]['coil_type'] = 9999
        info['chs'][j]['scanno'] = j + 1
        info['chs'][j]['logno'] = j + 1
        info['chs'][j]['kind'] = 1
        info['chs'][j]['range'] = 1.0
        info['chs'][j]['cal'] = 3.7000000285836165e-10
        info['chs'][j]['unit'] = 112
        info['chs'][j]['loc'] = np.array(
            [holders[i, 0], holders[i, 1], holders[i, 2], rot_mat[0, 0], rot_mat[0, 1], rot_mat[0, 2],
             rot_mat[1, 0], rot_mat[1, 1], rot_mat[1, 2], rot_mat[2, 0], rot_mat[2, 1], rot_mat[2, 2]])
    return info


def import_sensor_pos_ori_all(sensorholder_path, name, subject_dir, gen12=1):
    import megtools.pyread_biosig as pbio
    import megtools.vector_functions as vfun
    import mne
    import numpy as np

    if gen12 == 2:
        sensorholder_file = sensorholder_path + name + '_sensor_pos_ori_all_v2.txt'
        opm_trans_path = sensorholder_path + "opm_trans_v2.txt"
    else:
        sensorholder_file = sensorholder_path + name + '_sensor_pos_ori_all.txt'
        opm_trans_path = sensorholder_path + "opm_trans.txt"

    rotation, translation = pbio.import_opm_trans(opm_trans_path, name)
    translation = translation / 1000.0

    holders = pbio.imp_sensor_holders(sensorholder_file)
    holders[:, 0:3] = holders[:, 0:3] / 1000.0

    if gen12 == 2:
        # !!!!! be very carefull what comes first.
        holders[:, 0] = holders[:, 0] - translation[0]
        holders[:, 1] = holders[:, 1] - translation[1]
        holders[:, 2] = holders[:, 2] - translation[2]

    holders[:, 1], holders[:, 2] = vfun.rotate_via_numpy(holders[:, 1], holders[:, 2], np.radians(-rotation[0]))
    holders[:, 0], holders[:, 2] = vfun.rotate_via_numpy(holders[:, 0], holders[:, 2], np.radians(-rotation[1]))
    holders[:, 0], holders[:, 1] = vfun.rotate_via_numpy(holders[:, 0], holders[:, 1], np.radians(-rotation[2]))
    holders[:, 4], holders[:, 5] = vfun.rotate_via_numpy(holders[:, 4], holders[:, 5], np.radians(-rotation[0]))
    holders[:, 3], holders[:, 5] = vfun.rotate_via_numpy(holders[:, 3], holders[:, 5], np.radians(-rotation[1]))
    holders[:, 3], holders[:, 4] = vfun.rotate_via_numpy(holders[:, 3], holders[:, 4], np.radians(-rotation[2]))

    if gen12 != 2:
        holders[:, 0] = holders[:, 0] - translation[0]
        holders[:, 1] = holders[:, 1] - translation[1]
        holders[:, 2] = holders[:, 2] - translation[2]

        surf = mne.read_surface(subject_dir + name + "/surf/" + "lh.white", read_metadata=True)
        holders[:, 0:3] = holders[:, 0:3] - (surf[2]['cras'] / 1000.0)

    return holders


def geometry_initialization(block_name, subject_dir, workdir_path):
    # Steps for calculating BEM, sourcepaces etc..
    import mne
    name = block_name[0:4]
    geometry_files = dict()
    # geometry_files["fname_trans"] = workdir_path + block_name + '-trans.fif'
    geometry_files["src_path"] = workdir_path + name + '-oct6-src.fif'
    geometry_files["bem_surface"] = workdir_path + name + '-5120-5120-5120-bem.fif'
    geometry_files["fname_bem"] = workdir_path + name + '-5120-5120-5120-bem-sol.fif'

    # fname_trans = geometry_files["fname_trans"]
    src_path = geometry_files["src_path"]
    bem_surface = geometry_files["bem_surface"]
    fname_bem = geometry_files["fname_bem"]

    # try:
    #     with open(fname_trans) as f:
    #         None
    # except IOError:
    #     mne.gui.coregistration(subjects_dir=subject_dir, subject=name)
    # # mne.gui.coregistration(subjects_dir=subject_dir, subject=name)

    try:
        with open(src_path) as f:
            src = src_path
    except IOError:
        src = mne.setup_source_space(name, subjects_dir=subject_dir, spacing='oct6')
        mne.write_source_spaces(src_path, src)

    try:
        with open(bem_surface) as f:
            None
    except IOError:
        model = mne.make_bem_model(name, subjects_dir=subject_dir)
        mne.write_bem_surfaces(bem_surface, model)

    try:
        with open(fname_bem) as f:
            None
    except IOError:
        bem_sol = mne.make_bem_solution(model)
        mne.write_bem_solution(fname_bem, bem_sol)

    return geometry_files
