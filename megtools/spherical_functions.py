def create_icosahedron():
    import numpy

    choice = [-1, 1]
    phi = (1.0/2.0)*(1+numpy.sqrt(5.0))
    phichoice = choice*phi

    vertices = numpy.zeros((12, 3))
    k=0
    for i in choice:
        for j in phichoice:
            k += 1
    # aces = numpy.zeros(12)

    return vertices


def import_spherical(filename):
    import os
    import re
    import numpy

    # path = os.path.dirname(os.path.realpath(__file__))
    # complete_path = path + filename
    fin = open(filename, 'r')

    # import vertice
    num_vertices = int(fin.readline())
    vertices = []
    ii = 1
    jj = 0
    while ii < num_vertices:
        line = fin.readline()
        line = list(filter(None, re.split('  |m\r\n| |\r\n|\n', line)))
        vertices.append(line[1:4])
        ii = int(line[0])
    vertices = numpy.array(vertices).astype(float)

    # import faces
    num_faces = int(fin.readline())
    faces = []
    ii = 0
    jj = 0
    while ii < num_faces:
        line = fin.readline()
        line = list(filter(None, re.split('  |m\r\n| |\r\n|\n', line)))
        faces.append(line[1:4])
        ii = int(line[0])
    faces = numpy.array(faces).astype(int)

    return vertices, faces
