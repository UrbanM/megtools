def line3D(xyz):
    # pylint: disable=no-member
    """ a simple visuals for 3D plot """

    import numpy as np
    from vispy import app, gloo, visuals, scene

    # Define a simple vertex shader. We use $template variables as placeholders for
    # code that will be inserted later on.
    vertex_shader = """
    void main()
    {
        vec4 visual_pos = vec4($position, 1);
        vec4 doc_pos = $visual_to_doc(visual_pos);
    
        gl_Position = $doc_to_render(doc_pos);
    }
    """

    fragment_shader = """
    void main() {
      gl_FragColor = $color;
    }
    """

    # now build our visuals
    class Plot3DVisual(visuals.Visual):
        """ template """

        def __init__(self, x, y, z):
            """ plot 3D """
            visuals.Visual.__init__(self, vertex_shader, fragment_shader)

            # build Vertices buffer
            data = np.c_[x, y, z]
            v = gloo.VertexBuffer(data.astype(np.float32))

            # bind data
            self.shared_program.vert['position'] = v
            self.shared_program.frag['color'] = (1.0, 0.0, 0.0, 1.0)

            # config
            self.set_gl_state('opaque', clear_color=(1, 1, 1, 1))
            self._draw_mode = 'line_strip'

        def _prepare_transforms(self, view):
            """ This method is called when the user or the scenegraph has assigned
            new transforms to this visual """
            # Note we use the "additive" GL blending settings so that we do not
            # have to sort the mesh triangles back-to-front before each draw.
            tr = view.transforms
            view_vert = view.view_program.vert
            view_vert['visual_to_doc'] = tr.get_transform('visual', 'document')
            view_vert['doc_to_render'] = tr.get_transform('document', 'render')

    # build your visuals, that's all
    Plot3D = scene.visuals.create_visual_node(Plot3DVisual)

    # The real-things : plot using scene
    # build canvas
    canvas = scene.SceneCanvas(keys='interactive', show=True)

    # Add a ViewBox to let the user zoom/rotate
    view = canvas.central_widget.add_view()
    view.camera = 'turntable'
    view.camera.fov = 50
    view.camera.distance = 5

    # plot ! note the parent parameter
    p1 = Plot3D(xyz[:, 0], xyz[:, 1], xyz[:, 2], parent=view.scene)

    # run
    app.run()


def scatter3D(xyz):
    import sys
    import numpy as np
    from vispy import scene
    from vispy.scene import visuals

    #
    # Make a canvas and add simple view
    #
    canvas = scene.SceneCanvas(keys='interactive', show=True, bgcolor='w')

    view = canvas.central_widget.add_view()
    view.camera = 'turntable'
    view.camera.fov = 60

    scatter = visuals.Markers(parent=view.scene)
    # scatter.antialias = 0
    scatter.set_data(pos=xyz, edge_color=(0.0, 0.0, 0.0, 1.0), face_color=(0.6, 0.5, 0.4, 1.0), size=30)
    scatter.set_gl_state(depth_test=True, blend=True) #blend_func=('src_alpha', 'one_minus_src_alpha'))

    # Add axes
    axis = visuals.XYZAxis(parent=view.scene)
    canvas.app.run()


def TriangularVisualizationVTK(vertices, faces, chosen_faces, sources, dipoles):
    import vtk
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import numpy as np

    abs_values = np.sqrt(dipoles[:, 0] ** 2 + dipoles[:, 1] ** 2 + dipoles[:, 2] ** 2)

    pal = cmx.hot
    cNorm = colors.Normalize(vmin=min(abs_values), vmax=max(abs_values))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pal)
    scalarMap.get_clim()

    colors = vtk.vtkNamedColors()

    triangles = []
    mapper = []
    actor = []
    points = vtk.vtkPoints()

    # Create triangles
    ii = 0
    j = 0
    for i in range(0, len(faces), 1):
        points.InsertNextPoint(vertices[faces[i, 0]])
        points.InsertNextPoint(vertices[faces[i, 1]])
        points.InsertNextPoint(vertices[faces[i, 2]])

        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, ii)
        triangle.GetPointIds().SetId(1, ii + 1)
        triangle.GetPointIds().SetId(2, ii + 2)

        triangles.append(vtk.vtkCellArray())
        triangles[i].InsertNextCell(triangle)

        # Create a polydata object
        trianglePolyData = vtk.vtkPolyData()

        # Add the geometry and topology to the polydata
        trianglePolyData.SetPoints(points)
        trianglePolyData.SetPolys(triangles[i])

        # Create mapper and actor
        mapper.append(vtk.vtkOpenGLPolyDataMapper())
        mapper[i].SetInputData(trianglePolyData)
        actor.append(vtk.vtkOpenGLActor())

        if i in chosen_faces:
            color = scalarMap.to_rgba(abs_values[j])
            color = [color[0], color[1], color[2]]
            j += 1
        else:
            color = colors.GetColor3d("silver")

        actor[i].GetProperty().SetColor(color)
        actor[i].GetProperty().SetOpacity(1.00)
        actor[i].SetMapper(mapper[i])
        ii += 3

    # # Create arrows
    arr_mapper = []
    arr_actor = []
    arrowSource = vtk.vtkArrowSource()
    arrowSource.SetShaftRadius(2.0)
    arrowSource.SetTipRadius(4.0)
    # arrowSource.SetTipLength(2.0)

    for i in range(0, len(dipoles), 1):
        startPoint = sources[i] - 200*dipoles[i]/2
        endPoint = sources[i] + 200*dipoles[i]/2

        # Compute a basis
        normalizedX = [0] * 3
        normalizedY = [0] * 3
        normalizedZ = [0] * 3

        # The X axis is a vector from start to end
        vtk.vtkMath.Subtract(endPoint, startPoint, normalizedX)
        length = 100 * vtk.vtkMath.Norm(normalizedX)
        vtk.vtkMath.Normalize(normalizedX)

        # The Z axis is an arbitrary vector cross X
        arbitrary = [20.0, 20.0, 20.0]
        vtk.vtkMath.Cross(normalizedX, arbitrary, normalizedZ)
        vtk.vtkMath.Normalize(normalizedZ)

        # The Y axis is Z cross X
        vtk.vtkMath.Cross(normalizedZ, normalizedX, normalizedY)

        matrix = vtk.vtkMatrix4x4()
        # Create the direction cosine matrix
        matrix.Identity()
        for j in range(0, 3):
            matrix.SetElement(j, 0, normalizedX[j])
            matrix.SetElement(j, 1, normalizedY[j])
            matrix.SetElement(j, 2, normalizedZ[j])

        # Apply the transforms
        transform = vtk.vtkTransform()
        transform.Translate(startPoint)
        transform.Concatenate(matrix)
        transform.Scale(length, 0.2, 0.2)

        # Transform the polydata
        transformPD = vtk.vtkTransformPolyDataFilter()
        transformPD.SetTransform(transform)
        transformPD.SetInputConnection(arrowSource.GetOutputPort())

        # Create a mapper and actor
        arr_mapper.append(vtk.vtkPolyDataMapper())
        arr_actor.append(vtk.vtkActor())

        arr_mapper[i].SetInputConnection(arrowSource.GetOutputPort())
        arr_actor[i].SetUserMatrix(transform.GetMatrix())
        # arr_mapper[i].SetInputConnection(transformPD.GetOutputPort())

        arr_actor[i].SetMapper(arr_mapper[i])
        arr_actor[i].GetProperty().SetColor(colors.GetColor3d("Red"))

    # vtk.vtkScalarBarActor()

    # Create a renderer, render window, and an interactor
    renderer = vtk.vtkOpenGLRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName("Triangle")
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # # Add the actors to the scene
    for ii in range(0, len(actor), 1):
        renderer.AddActor(actor[ii])

    for ii in range(0, len(dipoles), 1):
        renderer.AddActor(arr_actor[ii])

    renderer.SetBackground(colors.GetColor3d("White"))

    # Render and interact
    renderWindow.Render()
    renderWindowInteractor.Start()
    return


def SimpleTriangularVisualizationVTK(vertices, faces):
    import vtk

    colors = vtk.vtkNamedColors()
    triangles = []
    mapper = []
    actor = []
    points = vtk.vtkPoints()
    # Create a triangle
    ii = 0
    for i in range(0, len(faces), 1):
        points.InsertNextPoint(vertices[faces[i, 0]])
        points.InsertNextPoint(vertices[faces[i, 1]])
        points.InsertNextPoint(vertices[faces[i, 2]])

        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, ii)
        triangle.GetPointIds().SetId(1, ii + 1)
        triangle.GetPointIds().SetId(2, ii + 2)

        triangles.append(vtk.vtkCellArray())
        triangles[i].InsertNextCell(triangle)

        # Create a polydata object
        trianglePolyData = vtk.vtkPolyData()

        # Add the geometry and topology to the polydata
        trianglePolyData.SetPoints(points)
        trianglePolyData.SetPolys(triangles[i])

        # Create mapper and actor
        mapper.append(vtk.vtkPolyDataMapper())
        mapper[i].SetInputData(trianglePolyData)
        actor.append(vtk.vtkActor())
        actor[i].GetProperty().SetColor(colors.GetColor3d("Red"))
        actor[i].SetMapper(mapper[i])
        ii += 3

    # Create a renderer, render window, and an interactor
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName("Triangle")
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Add the actors to the scene
    for ii in range(0, len(actor), 1):
        renderer.AddActor(actor[ii])
    renderer.SetBackground(colors.GetColor3d("White"))

    # Render and interact
    renderWindow.Render()
    renderWindowInteractor.Start()
    return


def simple_3D_plt(xyz):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c='r', marker='o', s=1)
    plt.show()
    return


def cortex_connectivity_plt(source_space, activity, conectivity_matrix):
    import vtk
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import numpy as np
    import vector_functions_v12 as vfun

    USER_MATRIX = True
    colors = vtk.vtkNamedColors()

    sphere_radius = 4.0 * activity/(np.max(activity))
    sphere_radius = np.clip(sphere_radius, 0.5, 4.0)
    cylinder_radius = 0.5 * conectivity_matrix/(np.max(conectivity_matrix))

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d("White"))

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(renderer)

    # An interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renwin)

    highest = np.argsort(activity)[-50:]
    highest1 = np.sort(highest)

    # for i in range(len(source_space)):
    #     if activity[i] > 0.95*np.max(activity):
    for i in range(len(source_space)):
        if(activity[i]>0.02):
            source = vtk.vtkSphereSource()

            source.SetRadius(sphere_radius[i])
            source.SetCenter(source_space[i, 0], source_space[i, 1], source_space[i, 2])
            source.SetPhiResolution(11)
            source.SetThetaResolution(21)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            r = vtk.vtkMath.Random(.4, 1.0)
            g = vtk.vtkMath.Random(.4, 1.0)
            b = vtk.vtkMath.Random(.4, 1.0)
            actor.GetProperty().SetColor(colors.GetColor3d("Red"))

            renderer.AddActor(actor)

    for i in highest:
        for j in highest:
            if conectivity_matrix[i][j] >= 0.50*np.max(conectivity_matrix): # and activity[i] > 0.80*np.max(activity) and activity[j] > 0.80*np.max(activity):
                source = vtk.vtkCylinderSource()
                source.SetRadius(cylinder_radius[i][j])
                source.SetResolution(50)

                startPoint = source_space[i]
                endPoint = source_space[j]

                # Compute a basis
                normalizedX = [0] * 3
                normalizedY = [0] * 3
                normalizedZ = [0] * 3

                # The X axis is a vector from start to end
                vtk.vtkMath.Subtract(endPoint, startPoint, normalizedX)
                length = vtk.vtkMath.Norm(normalizedX) #vfun.dist_two_points(startPoint, endPoint)
                vtk.vtkMath.Normalize(normalizedX)

                # The Z axis is an arbitrary vector cross X
                rng = vtk.vtkMinimalStandardRandomSequence()
                rng.SetSeed(8775070)  # For testing.8775070

                arbitrary = [0]*3
                for l in range(0, 3):
                    rng.Next()
                    arbitrary[l] = rng.GetRangeValue(-10, 10)
                vtk.vtkMath.Cross(normalizedX, arbitrary, normalizedZ)
                vtk.vtkMath.Normalize(normalizedZ)

                # The Y axis is Z cross X
                vtk.vtkMath.Cross(normalizedZ, normalizedX, normalizedY)
                matrix = vtk.vtkMatrix4x4()

                # Create the direction cosine matrix
                matrix.Identity()
                for k in range(0, 3):
                    matrix.SetElement(k, 0, normalizedX[k])
                    matrix.SetElement(k, 1, normalizedY[k])
                    matrix.SetElement(k, 2, normalizedZ[k])

                # Apply the transforms
                transform = vtk.vtkTransform()
                transform.Translate(startPoint)
                transform.Concatenate(matrix)
                transform.RotateZ(-90.0)  # align cylinder to x axis
                transform.Scale(1.0, length, 1.0)  # scale along the height vector
                transform.Translate(0, .5, 0)  # translate to start of cylinder
                # transform.Scale(length, 0.5, 0.5)

                # Transform the polydata
                transformPD = vtk.vtkTransformPolyDataFilter()
                transformPD.SetTransform(transform)
                transformPD.SetInputConnection(source.GetOutputPort())

                # Create a mapper and actor
                arr_mapper = vtk.vtkPolyDataMapper()
                arr_actor = vtk.vtkActor()

                if USER_MATRIX:
                    arr_mapper.SetInputConnection(source.GetOutputPort())
                    arr_actor.SetUserMatrix(transform.GetMatrix())
                else:
                    arr_mapper.SetInputConnection(transformPD.GetOutputPort())
                arr_actor.SetMapper(arr_mapper)
                arr_actor.GetProperty().SetColor(colors.GetColor3d("Black"))

                # arr_mapper.SetInputConnection(source.GetOutputPort())
                # arr_actor.SetUserMatrix(transform.GetMatrix())
                #
                # arr_actor.SetMapper(arr_mapper)
                # arr_actor.GetProperty().SetColor(colors.GetColor3d("Black"))

                renderer.AddActor(arr_actor)
            # print(j)
        # print(i)

    # Start
    interactor.Initialize()
    interactor.Start()

    return
