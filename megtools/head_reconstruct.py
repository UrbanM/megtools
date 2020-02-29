#!/usr/bin/env python
# -*- coding: utf-8 -*-


def convert_dicom_nifti(dicom_path, nifti_path):
    import dicom2nifti

    dicom2nifti.convert_directory(dicom_path, nifti_path)

    return 0


def sort_dicom():
    import pydicom
    import os

    PathDicom = "/home/marhl/Sync/IMFM/Programi/HeadReconstruction/DATA/MEG-Anjo/"
    secondDir = "/home/marhl/Sync/IMFM/Programi/HeadReconstruction/DATA/MEG-Anjo/SliceSet/"
    lstFilesDCM = []  # create an empty list
    lstFilenamesDCM = []  # create an empty list
    for dirName, subdirList, fileList in os.walk(PathDicom):
        for filename in fileList:
            # if ".dcm" in filename.lower():  # check whether the file's DICOM
            lstFilesDCM.append(os.path.join(dirName, filename))
            lstFilenamesDCM.append(filename)

    new_names = []
    for i in lstFilenamesDCM:
        i = secondDir + i
        new_names.append(i)

    lenght = len(lstFilesDCM)
    for i in range(0, lenght, 1):
        # read the file
        dataset = pydicom.read_file(lstFilesDCM[i])

        de = 'SeriesDescription'
        a = dataset.data_element(de).value
        if a == 't1_mpr_ns_sag_pat2_iso_asy_meg':
            os.rename(lstFilesDCM[i], new_names[i])
    return 0


def import_dicom(import_path):
    import pydicom as dicom
    import numpy as np
    import os

    lstFilesDCM = []  # create an empty list
    for dirName, subdirList, fileList in os.walk(import_path):
        for filename in fileList:
            if ".dcm" in filename.lower():  # check whether the file's DICOM
                lstFilesDCM.append(os.path.join(dirName, filename))
                lstFilesDCM.sort()

    RefDs = dicom.read_file(lstFilesDCM[1])
    # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
    ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))

    # Load spacing values (in mm)
    x = -(RefDs.ImagePositionPatient[1] + np.arange(ConstPixelDims[0]) * float(RefDs.PixelSpacing[0]))
    y = RefDs.ImagePositionPatient[2] - np.arange(ConstPixelDims[1]) * float(RefDs.PixelSpacing[1])
    z = np.zeros(len(lstFilesDCM))
    for ii in range(0, len(lstFilesDCM), 1):
        dataset = dicom.read_file(lstFilesDCM[ii])
        z[ii] = dataset.ImagePositionPatient[0]
    z = np.sort(z)

    # The array is sized based on 'ConstPixelDims'
    ArrayDicom = np.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

    # loop through all the DICOM files
    for filenameDCM in lstFilesDCM:
        # read the file
        ds = dicom.read_file(filenameDCM)
        # store the raw image data
        ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = ds.pixel_array

    return x, y, z, ArrayDicom


def import_nifti(import_path):
    import nibabel as nib
    import numpy as np
    import matplotlib.pyplot as plt

    img = nib.load(import_path)

    img_data = img.get_fdata()
    shape = img.shape
    # dimension = img.header.get_sform()

    M = img.affine[:3, :3]
    abc = img.affine[:3, 3]

    def f(i, j, k):
        """ Return X, Y, Z coordinates for i, j, k """
        return M.dot([i, j, k]) + abc

    x = []
    y = []
    z = []

    # location = np.array((shape[0]*shape[1]*shape[2], 3))
    # ii=0
    # for i in range(shape[0]):
    #     for j in range(shape[1]):
    #         for k in range(shape[2]):
    #             location[] = f(i,j,k)
    #     print(i)
    #
    # print(location)

    # img_data = np.moveaxis(img_data, 0, -1)
    # img_data = np.flip(img_data, 1)

    # y_n = z
    # z_n = y
    # z = z_n
    # y = y_n

    return img, img_data


def slice_mri_nifti(import_path, export_path, rotation, translation):
    import numpy as np
    import cv2
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    import vector_functions_v12 as vfun

    img, values = import_nifti(import_path)

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

    nn = 1
    # nn = int(input("raizrisi vsake koliko tock:"))
    # Slicanje po konst. x slicih
    for ii in range(0, a, nn):
        img = np.uint8(img1[ii, :, :])
        edges = cv2.Canny(img, 40, 700)

        for jj in range(0, c, nn):
            first = -1
            last = -1
            for kk in range(0, b, nn):
                if edges[kk, jj] > 0:
                    if first > 0:
                        last = kk
                    else:
                        first = kk
            if first >= 0:
                x_scatters.append(f(ii, first, jj)[0])
                y_scatters.append(f(ii, first, jj)[1])
                z_scatters.append(f(ii, first, jj)[2])
                stevec = ii * a * b + jj * b + first
                no_point.append(stevec)

            if last >= 0:
                x_scatters.append(f(ii, last, jj)[0])
                y_scatters.append(f(ii, last, jj)[1])
                z_scatters.append(f(ii, last, jj)[2])
                stevec = ii * a * b + jj * b + last
                no_point.append(stevec)

    # Slicanje po konst. z slicih
    for kk in range(0, c, nn):
        img = np.uint8(img1[:, :, kk])
        edges = cv2.Canny(img, 60, 600)

        for jj in range(0, b, nn):
            first = -1
            last = -1
            for ii in range(0, a, nn):
                if edges[ii, jj] > 0:
                    if first > 0:
                        last = ii
                    else:
                        first = ii
            if first >= 0:
                x_scatters.append(f(first, jj, kk)[0])
                y_scatters.append(f(first, jj, kk)[1])
                z_scatters.append(f(first, jj, kk)[2])
                stevec = first * a * b + kk * b + jj
                no_point.append(stevec)

            if last >= 0:
                x_scatters.append(f(last, jj, kk)[0])
                y_scatters.append(f(last, jj, kk)[1])
                z_scatters.append(f(last, jj, kk)[2])
                stevec = last * a * b + kk * b + jj
                no_point.append(stevec)

    x_scatters = np.array(x_scatters)
    y_scatters = np.array(y_scatters)
    z_scatters = np.array(z_scatters)
    no_point = np.array(no_point)

    scatters = np.stack((x_scatters, y_scatters, z_scatters, no_point), 1)

    # scatters = scatters[np.argsort(scatters[:, 3])]
    #
    # # Check if two of the same elements in the array and delete them
    # ii = 1
    # while ii < scatters.shape[0]:
    #     if scatters[ii, 3] == scatters[ii - 1, 3]:
    #         scatters = np.delete(scatters, ii, 0)
    #     ii += 1

    # scatters = np.delete(scatters, 3, 1)

    scatters[:,1], scatters[:,2] = vfun.rotate_via_numpy(scatters[:,1], scatters[:,2], np.radians(rotation[0]))
    scatters[:,0], scatters[:,1] = vfun.rotate_via_numpy(scatters[:,0], scatters[:,1], np.radians(rotation[2]))

    scatters[:,0] = scatters[:,0] + translation[0]
    scatters[:,1] = scatters[:,1] + translation[1]
    scatters[:,2] = scatters[:,2] + translation[2]

    # Write binnary file
    scatters.tofile(export_path)  # save the data
    np.savetxt(export_path + ".txt", scatters)

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(scatters[:, 0], scatters[:, 1], scatters[:, 2], c='r', marker='o', s=1)
    # ax.set_xlabel('x os rezine', fontsize=18)
    # ax.set_ylabel('posamezne rezine', fontsize=18)
    # ax.set_zlabel('y os rezine', fontsize=18)
    # plt.show()

    return 0


def simplify_cortex(xyz):
    from scipy.spatial import distance
    import vector_functions_v12 as vfun

    # 1 najdi izhodisce, izracunaj povprecje
    x_povp = np.mean(xyz[:, 0])
    y_povp = np.mean(xyz[:, 1])
    z_povp = np.mean(xyz[:, 2])

    min_dist = 0.5

    # dst = distance.euclidean(xyz[1, :], xyz[0, :])

    dx = xyz[:, 0] - x_povp
    dy = xyz[:, 1] - y_povp
    dz = xyz[:, 2] - z_povp
    distance = np.sqrt(dx*dx + dy*dy + dz*dz)
    max_ind = np.argmax(distance)

    radius = distance[max_ind]
    phi_int = 15
    rho_int = 15
    dphi = 2*np.pi/phi_int
    drho = np.pi/rho_int

    phi = [dphi * n for n in range(phi_int)]
    rho = [drho * n for n in range(rho_int)]

    rhophi = np.zeros((phi_int*rho_int, 2))
    ii = 0
    for j in rho:
        for i in phi:
            rhophi[ii, :] = (i, j)
            ii += 1
    #2 napravi krog-n kotnik

    xyz=[]
    rr = radius
    while len(rhophi) > 0 and rr > 0.0:
        for i in rhophi:
            xyz_one = vfun.trans_spher_cart([radius, i[0], i[1]])
            # xyz_one[0] = xyz_one[0] + x_povp
            # xyz_one[1] = xyz_one[1] + y_povp
            # xyz_one[2] = xyz_one[2] + z_povp

            print(xyz_one)

            # dx = xyz[:, 0] - xyz_one[0]
            # dy = xyz[:, 1] - xyz_one[1]
            # dz = xyz[:, 2] - xyz_one[2]
            distance = np.sqrt(dx * dx + dy * dy + dz * dz)

    #3 zacni manjsat ta n kotnik

    return


def simplify_cortex_imported(xyz, datapath):
    from scipy.spatial import distance
    import numpy as np
    import spherical_functions_v10 as sfun

    def check_nearest_neighbours(lists, point):
        dx = lists[:, 0] - point[0]
        dy = lists[:, 1] - point[1]
        dz = lists[:, 2] - point[2]
        distances = np.sqrt(dx * dx + dy * dy + dz * dz)
        return distances

    # 1 najdi izhodišče, izračunaj povprečje
    x_povp = np.mean(xyz[:, 0])
    y_povp = np.mean(xyz[:, 1])
    z_povp = np.mean(xyz[:, 2])

    min_dist = 5.0

    xyz[:, 0] = xyz[:, 0] - x_povp
    xyz[:, 1] = xyz[:, 1] - y_povp
    xyz[:, 2] = xyz[:, 2] - z_povp

    distance = check_nearest_neighbours(xyz, np.zeros(3))

    n = 1
    max_ind = np.argsort(distance)[::-1][:n]
    max_radius = distance[max_ind]

    filename = datapath + "/trikotizacija/dodeca_482.tri"
    vert, face = sfun.import_spherical(filename)
    face = face - 1
    vert1 = vert

    fixed = np.zeros(len(vert))
    rad = np.linspace(max_radius, 0.0, num=20)

    for j in rad:
        k = 0
        for i in vert:
            if fixed[k] == 0:
                i = i * j
                distance = check_nearest_neighbours(xyz, i)
                ll = 0
                for n in distance:
                    if n < min_dist:
                        ll += 1
                if ll > 0:
                    vert1[k] = i
                    fixed[k] = 1
            k += 1

    vert1[:, 0] = vert1[:, 0] + x_povp
    vert1[:, 1] = vert1[:, 1] + y_povp
    vert1[:, 2] = vert1[:, 2] + z_povp

    return vert1, face

def create_cortex_vertices(import_path, export_path, data_path, name, rotation, translation):
    import mne
    import numpy as np
    import my3Dplot_v11 as mplt
    import vector_functions_v12 as vfun

    vertices1, faces1, metadata1 = mne.read_surface(import_path + name + "/surf/lh.pial", read_metadata=True)
    vertices2, faces2, metadata2 = mne.read_surface(import_path + name + "/surf/rh.pial", read_metadata=True)
    # print(mne.surface.complete_surface_info(surf))
    vertices1 = vertices1 + metadata1["cras"]
    vertices2 = vertices2 + metadata2["cras"]

    vertices = np.vstack((vertices1, vertices2))
    faces = np.vstack((faces1, faces2))

    # np.savetxt('vertices.txt', vertices, delimiter=' ')
    # np.savetxt('faces.txt', faces, fmt='%i', delimiter=' ')

    # vertices = vertices[0::100]

    vert, face = simplify_cortex_imported(vertices, data_path)

    vert[:,1], vert[:,2] = vfun.rotate_via_numpy(vert[:,1], vert[:,2], np.radians(rotation[0]))
    vert[:,0], vert[:,1] = vfun.rotate_via_numpy(vert[:,0], vert[:,1], np.radians(rotation[2]))

    vert[:,0] = vert[:,0] + translation[0]
    vert[:,1] = vert[:,1] + translation[1]
    vert[:,2] = vert[:,2] + translation[2]

    # vertices_nth = vertices[0::100]
    # mplt.TriangularVisualizationVTK(vert, face)
    # mplt.scatter3D(vert)
    # mplt.simple_3D_plt(vert)
    # mplt.line3D(vert)

    np.savetxt(export_path + 'vertices_' + name +'.txt', vert, delimiter=' ')
    np.savetxt(export_path + 'faces_' + name +'.txt', face, fmt='%i', delimiter=' ')

    return vert, face



