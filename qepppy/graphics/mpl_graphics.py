import numpy as np

from .. import utils


def draw_cell(ax, cell, center=[0,0,0]):
    """
    Draw a cell centered on 'center'.
    Params:
     - ax: matplotlib 3D axis object
     - center: center for plotting the cell
    """
    V = cell
    for n1 in range(3):
        orig = np.array(center)
        v0 = V[n1]
        for n2 in range(4):
            v = np.array(np.vstack((orig, orig + v0)))
            if n2 == n1:
                orig = V[(n2+1)%3] + V[(n2+2)%3]
            else:
                orig = V[n2%3]
            ax.plot(v[:,0], v[:,1], v[:,2], color='black', linewidth=0.5)
        ax.plot(v[:,0], v[:,1], v[:,2], color='black', linewidth=0.5)

def draw_Wigner_Seitz(
    ax, recipr, draw_corners=True,
    edge_args={
        'color': 'black'
    },
    corner_args={
        'color': 'green'
    }
    ):
    """
    Draw the Wigner Seitz cell for the lattice.
    Params:
     - ax: matplotlib 3D axis object
    """
    res = []
    dim = len(recipr)
    from scipy.spatial import KDTree, Voronoi
    L = utils.generate_repetition_grid([[-1,0,1]]*dim, vect_matrix=recipr)

    vor = Voronoi(L)
    P = vor.vertices
    R = vor.ridge_vertices

    tree       = KDTree(L)
    d,i        = tree.query([0,]*dim)
    dist, ind  = tree.query(P, k=L.shape[0])
    w          = (np.abs(dist.T - dist.T[0]) > 1E-5).T
    closest    = ind.copy()
    closest[w] = -1

    cond    = np.where([not i in a for a in closest])[0]
    P[cond] = np.zeros(dim)


    for i1, e in enumerate(R):
        for i2, r in enumerate(e):
            if r in cond:
                R[i1][i2] = -1

    if draw_corners:
        app = ax.scatter(*P.T, **corner_args)
        res.append(app)

    for vert in R:
        vert.append(vert[0])
        v = np.asarray(vert)
        if np.all(v >= 0):
            app = ax.plot(*P[v].T, **edge_args)
            res.extend(app)

    return res

def draw_plane(ax, v1, v2, origin=[0,0,0], **kwargs):
    a, b = np.meshgrid([0,1], [0,1], indexing='ij')
    app = np.dot(np.array([a.flatten(), b.flatten()]).T, [v1,v2])
    xx,yy,zz = (app + origin).T.reshape(3,2,2)

    return ax.plot_surface(xx,yy,zz, **kwargs)

def draw_sphere(ax, radius=1, center=[0,0,0], long=10, lat=7, **kwargs):
    center = np.array(center).reshape(-1,3)

    u = np.linspace(0, 2 * np.pi, long)
    v = np.linspace(0, np.pi, lat)
    x0 = radius * np.outer(np.cos(u), np.sin(v))
    y0 = radius * np.outer(np.sin(u), np.sin(v))
    z0 = radius * np.outer(np.ones(np.size(u)), np.cos(v))

    res = []
    for (xoff,yoff,zoff) in center:
        app = ax.plot_surface(x0+xoff, y0+yoff, z0+zoff, **kwargs)
        res.append(app)

    return res


def draw_cylinder(ax, radius=1, axis=[0,0,1], start=[0,0,0], color='b'):
    x0, y0, z0 = start

    norm   = np.linalg.norm(axis)
    axis  /= norm
    c_teta = axis[2]
    s_teta = np.sqrt(1 - c_teta**2)
    c_phi  = axis[0]/s_teta if s_teta else 1
    s_phi  = axis[1]/s_teta if s_teta else 0

    u = np.linspace(0, 2*np.pi, 8)
    v = np.linspace(0, 1*norm, 2)

    x = radius * np.outer(np.cos(u)*s_phi + np.sin(u)*c_phi*c_teta, np.ones(np.size(v)))   + axis[0] * v + x0
    y = radius * np.outer(-np.cos(u)*c_phi + np.sin(u)*s_phi*c_teta, np.ones(np.size(v)) ) + axis[1] * v + y0
    z = radius * np.outer(-np.sin(u)*s_teta, np.ones(np.size(v)))                          + axis[2] * v + z0

    return ax.plot_surface(x, y, z, color=color)

def draw_bond(ax, start, end, color1, color2, graph_lvl=0):
    if graph_lvl == 0:
        v = np.vstack((start, end))
        ax.plot(v[:,0], v[:,1], v[:,2], color='black', linewidth=1.5)
    elif graph_lvl == 1 or graph_lvl == 2:
        mid = (start + end) / 2
        v = np.vstack((start, mid))
        ax.plot(v[:,0], v[:,1], v[:,2], color=color1, linewidth=3.5)
        v = np.vstack((mid, end))
        ax.plot(v[:,0], v[:,1], v[:,2], color=color2, linewidth=3.5)
    elif graph_lvl == 3:
        mid = (start + end) / 2
        if color1 != color2:
            draw_cylinder(ax, radius=0.15, axis=(end-start)/2, start=start,  color=color1)
            draw_cylinder(ax, radius=0.15, axis=(end-start)/2, start=mid, color=color2)
        else:
            draw_cylinder(ax, radius=0.15, axis=end-start, start=start,  color=color1)
    else:
        raise ValueError("arg 'graph_lvl' must be <= 3")

def draw_atom(ax, X,Y,Z, color='k', name='None', radius=1, graph_lvl=0):
    if graph_lvl == 0 or graph_lvl == 1:
        ax.scatter(
            X, Y, Z,
            s=80*radius,
            marker='o',
            depthshade=False,
            c=color,
            label=name
            )
    elif graph_lvl == 2 or graph_lvl == 3:
        ax.scatter(
            X, Y, Z,
            s=10,
            marker='o',
            depthshade=False,
            c=color,
            label=name
            )
        for x,y,z in zip(X,Y,Z):
            draw_sphere(ax, radius=radius*0.3, center=[x,y,z], color=color)
    else:
        raise ValueError("arg 'graph_lvl' must be <= 3")
