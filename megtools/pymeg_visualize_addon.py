def plot_magnetometers3(subject_dir, name, xyz1, rot_mat, magnetometer_number, coil_def, filename="rand_name"):
    # plot sensors
    import os
    # import megtools.my3Dplot as m3p
    import numpy as np

    accuracy = 2

    surface1 = os.path.join(subject_dir, name, 'bem', 'watershed', name + '_inner_skull_surface')
    surface2 = os.path.join(subject_dir, name, 'bem', 'watershed', name + '_outer_skull_surface')
    surface3 = os.path.join(subject_dir, name, 'bem', 'watershed', name + '_outer_skin_surface')

    surfaces = [surface1, surface2, surface3]

    with open(coil_def) as f:
        content = f.readlines()

    idx_range = []

    for j, i in enumerate(content):
        if i.split()[0].isdigit():
            if magnetometer_number == int(i.split()[1]) and accuracy == int(i.split()[2]):
                magnetometer_gradiometer = int(i.split()[0])
                print(i, j)
                idx_range = [j+1, j+int(i.split()[3])+1]

    arr = []
    for i in range(idx_range[0], idx_range[1]):
        arr.append(content[i].split())

    arr = np.array(arr, dtype=float)

    uniques, counts = np.unique(arr[:, 0], return_counts=True)

    locations = xyz1[:, 0:3]
    locations2 = []
    directions2 = []
    unit_v = np.array([0., 0., 1.])
    directions = np.zeros(np.shape(locations))

    elements = []
    temp_list_ind = 0
    for i in range(len(locations)):
        for jj, ii in enumerate(uniques):
            temp_list_el = []
            for j in range(len(arr)):
                if arr[j, 0] == ii:
                    arr_temp = np.dot(arr[j, 1:4], rot_mat[i])
                    locations2.append(locations[i]+arr_temp)
                    directions2.append(np.dot(arr[j, 4:7], rot_mat[i]))
                    temp_list_el.append(temp_list_ind)
                    temp_list_ind += 1
            elements.append(temp_list_el)

            # print(ii)
        # for j in range(len(arr)):
        #     arr_temp = np.dot(arr[j, 1:4], rot_mat[i])
        #     locations2.append(locations[i]+arr_temp)
        #     directions2.append(np.dot(arr[j, 4:7], rot_mat[i]))

        # directions[i] = np.dot(unit_v, rot_mat[i])

    sensors = np.hstack((np.array(locations2), np.array(directions2)))

    p1 = plot_sensors_pyvista1(surfaces, sensors=sensors, elements=elements, arrow_color="black", grad=magnetometer_gradiometer)

    # p1.show(screenshot=name + 'opm_rad.png')
    # p1 = m3p.plot_sensors_pyvista(surfaces, sensors=[])
    p1.show(screenshot=filename+'.png')
    # p1.show()
    p1.close()

    return


def plot_sensors_pyvista1(surfaces, sensors, sensors2=[], elements=[], arrow_color="black", grad=0):
    import pyvista as pv
    import mne
    import numpy as np
    import random

    pv.set_plot_theme("document")
    p = pv.Plotter()

    # for i in sensors:
    #     sphere = pv.Sphere(center=i[0:3] * 10 ** 3, radius=5)
    #     p.add_mesh(sphere, color="black")
    #     arrow = pv.Arrow(start=i[0:3] * 10 ** 3, direction=i[3:6], scale=15)
    #     p.add_mesh(arrow, color=arrow_color)

    for cunt, element in enumerate(elements):
            pts = np.zeros((len(element),3))
            if grad == 3 or grad == 2:
                if cunt % 2:
                    color = (random.uniform(0., 1.), random.uniform(0., 1.), random.uniform(0., 1.))
            else:
                color = "gray"
            for i, j in enumerate(element):
                pts[i,0] = sensors[j, 0]*10**3
                pts[i,1] = sensors[j, 1]*10**3
                pts[i,2] = sensors[j, 2]*10**3
            if len(element) == 1:
                sphere = pv.Sphere(center=pts, radius=1)
                p.add_mesh(sphere, color="black")
                arrow = pv.Arrow(start=pts[i, 0:3], direction=sensors[j, 3:6], scale=20)
                p.add_mesh(arrow, color=arrow_color)
            if len(element) == 4:
                faces = np.array([4, 0, 1, 2, 3])
                mesh = pv.PolyData(pts, faces)
                p.add_mesh(mesh, color=color)
            if len(element) == 6:
                faces = np.array([6, 0, 3, 5, 1, 4, 2])
                mesh = pv.PolyData(pts, faces)
                p.add_mesh(mesh, color=color)
            if len(element) == 8:
                faces = np.array([[4, 0, 1, 3, 2], [4, 4, 5, 7, 6], [4, 2, 3, 7, 6], [4, 0, 1, 5, 4], [4, 1, 3, 7, 5],
                                  [4, 0, 2, 6, 4]])
                mesh = pv.PolyData(pts, faces)
                p.add_mesh(mesh, color="red")
                arrow = pv.Arrow(start=np.mean(pts, axis=0), direction=sensors[j, 3:6], scale=20)
                p.add_mesh(arrow, color=arrow_color)

    step = 1.0 / (len(surfaces) + 1)
    opacities = np.linspace(1 - step, 0, num=len(surfaces), endpoint=False)
    for i, surface in enumerate(surfaces):
        rr_mm, tris = mne.read_surface(surface)
        tres = np.ones((len(tris), 1), dtype=int) * 3
        tris = np.hstack([tres, tris])
        gray1 = (0.5, 0.5, 0.5)
        polygon = pv.PolyData(rr_mm, tris)
        p.add_mesh(polygon, color=gray1, opacity=opacities[i] - 0.3)

    camera = pv.Camera()
    # p.camera.zoom(1.6)
    p.camera.zoom(3.0)

    return p


def plot_magnetometers31(subject_dir, name, evoked, magnetometer_number, coil_def, filename="rand_name"):
    # plot sensors
    import os
    import megtools.my3Dplot as m3p
    import numpy as np

    xyz1 = []
    rot_mat = []
    for j, i in enumerate(evoked.info.ch_names):
        xyz1.append(evoked.info['chs'][j]['loc'][0:3])
        rot_mat.append(evoked.info['chs'][j]['loc'][3:12].reshape((3, 3)))
    xyz1 = np.array(xyz1)
    rot_mat = np.array(rot_mat)

    accuracy = 2

    surface1 = os.path.join(subject_dir, name, 'bem', 'watershed', name + '_inner_skull_surface')
    surface2 = os.path.join(subject_dir, name, 'bem', 'watershed', name + '_outer_skull_surface')
    surface3 = os.path.join(subject_dir, name, 'bem', 'watershed', name + '_outer_skin_surface')

    surfaces = [surface1, surface2, surface3]

    with open(coil_def) as f:
        content = f.readlines()

    for j, i in enumerate(content):
        if (i.split()[0].isdigit()):
            if magnetometer_number == int(i.split()[1]) and accuracy == int(i.split()[2]):
                magnetometer_gradiometer = int(i.split()[0])
                print(i, j)
                idx_range = [j+1, j+int(i.split()[3])+1]


    arr = []
    for i in range(idx_range[0], idx_range[1]):
        arr.append(content[i].split())

    arr = np.array(arr, dtype=float)

    uniques, counts = np.unique(arr[:, 0], return_counts=True)

    locations = xyz1[:, 0:3]
    locations2 = []
    directions2 = []
    unit_v = np.array([0., 0., 1.])
    directions = np.zeros(np.shape(locations))

    elements = []
    temp_list_ind = 0
    for i in range(len(locations)):
        for jj, ii in enumerate(uniques):
            temp_list_el = []
            for j in range(len(arr)):
                if arr[j, 0] == ii:
                    arr_temp = np.dot(arr[j, 1:4], rot_mat[i])
                    locations2.append(locations[i]+arr_temp)
                    directions2.append(np.dot(arr[j, 4:7], rot_mat[i]))
                    temp_list_el.append(temp_list_ind)
                    temp_list_ind += 1
            elements.append(temp_list_el)

            # print(ii)
        # for j in range(len(arr)):
        #     arr_temp = np.dot(arr[j, 1:4], rot_mat[i])
        #     locations2.append(locations[i]+arr_temp)
        #     directions2.append(np.dot(arr[j, 4:7], rot_mat[i]))

        # directions[i] = np.dot(unit_v, rot_mat[i])

    sensors = np.hstack((np.array(locations2), np.array(directions2)))

    p1 = plot_sensors_pyvista1(surfaces, sensors=sensors, elements=elements, arrow_color="black", grad=magnetometer_gradiometer)

    # p1.show(screenshot=name + 'opm_rad.png')
    # p1 = m3p.plot_sensors_pyvista(surfaces, sensors=[])
    p1.show(screenshot=filename+'.png')
    # p1.show()
    p1.close()

    return
