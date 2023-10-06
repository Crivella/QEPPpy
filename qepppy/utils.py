import numpy as np


def recipr_base(base):
    return np.linalg.inv(base).T * 2 * np.pi

def generate_repetition_grid(R, vect_matrix=None):
    from itertools import product
    res = np.array(list(product(*R)))
    if not vect_matrix is None:
        res = res.dot(vect_matrix)

    return res

def xyz_mesh(shape, base=None, rep=1, reverse=False):
    _, n1,n2,n3 = shape
    try:
        r1,r2,r3 = rep
    except (TypeError, ValueError):
        r1 = r2 = r3 = rep
    a = np.linspace(0, r1, n1*r1 + 1)[1:] #[:-1] #+ .5/n1
    b = np.linspace(0, r2, n2*r2 + 1)[1:] #[:-1] #+ .5/n2
    c = np.linspace(0, r3, n3*r3 + 1)[1:] #[:-1] #+ .5/n3

    if reverse:
        a, c = c, a

    # Specific order to obtain the array with shape (n1,n2,n3) as the data grid
    # The 'b,a,c' order is because for a 3d meshgrid the resulting shape is (1,2,3) --> (2,1,3)
    # The 'y,x,z' order is because of how the 3d meshgrid output behaves:
    #    x,y,z=np.meshgrid(1,2,3)
    #       will cause the x to change value along axis=1
    #                       y to change value along axis=0
    #                       z to change value along axis=2
    # Since the FFT grid has the axis=0,1,2 corresponding to x,y,z i need to do the proper remapping
    # y,x,z = np.meshgrid(b,a,c)

    # OR juse use 'ij' indexing :)
    x,y,z = np.meshgrid(a,b,c, indexing='ij')

    if reverse:
        x,z = z,x

    if not base is None:
        XYZ  = np.dot(
            base.T,
            [x.flatten(),y.flatten(),z.flatten()],
            )
    else:
        XYZ = [x,y,z]

    return np.array(XYZ).reshape(3,*x.shape)

def lowdin_ortho(base):
    """
    https://arxiv.org/abs/1105.3571v1
    O = (M . M^T)^{-1/2} . M
    Params:
     - base: np.array with states along axis=0. If the states are
             multidimensional (eg: FFT grids), flatten them before. After the
             orthonormalization, restore the shape
    """
    from scipy.linalg import inv, sqrtm

    shape   = base[0].shape
    base    = np.array([a.flatten() for a in base])
    overlap = np.dot(base.conj(), base.T)
    print('OVERLAP: \n', overlap)

    base = np.dot(
        inv(sqrtm(overlap)).conj(),
        base
        )

    overlap = np.dot(base.conj(), base.T)
    print('OVERLAP: \n', overlap)

    return base.reshape(-1,*shape)

def remap_plane(
    invT,
    Xlim, Ylim, Zlim,
    shape, rep, fixaxis=2
    ):
    lim = (Xlim, Ylim, Zlim)
    shape = np.array(shape)

    # xmin, xmax = Xlim
    # ymin, ymax = Ylim
    # zmin, zmax = Zlim

    m1 = np.linspace(*lim[(fixaxis+1)%3], shape[0] * rep[0])
    m2 = np.linspace(*lim[(fixaxis+2)%3], shape[1] * rep[1])
    a,b = np.meshgrid(m1,m2)

    fix = np.ones(a.shape) * (sum(lim[fixaxis]))/2
    rec = np.round(
            np.dot(
            invT,
            [a.flatten(),b.flatten(),fix.flatten()]
            ),
            decimals=8
        ) % 1

    rec = np.array(rec) * shape.reshape(3,1)
    rec = rec.astype(dtype='int')
    i,j,k = rec

    return a, b, tuple((i,j,k))

def remap_generic_plane(
        cell: np.ndarray,
        v1: np.ndarray, v2: np.ndarray,
        shape_plane: tuple,
        shape_out: tuple,
        offset=None
    ) -> tuple[np.ndarray, np.ndarray]:
    """Remap a generic plane defined by two vectors v1 and v2 in a grid defined by the cell vectors 'cell'
    such that cell = [a1,a2,a3]

    Args:
        cell (np.ndarray): matrix with the cell vectors as rows
        v1 (np.ndarray): 1st vector defining the plane
        v2 (np.ndarray): 2nd vector defining the plane
        shape_plane (tuple): shape of the plane grid (arbitrary)
        shape_out (tuple): shape of the output grid (!!! this is used to generate the return indexes hence it should
            match the shape of the grid on which you wish to remap the plane)
        offset (np.ndarray, optional): By default the plane is defined from the origin (0,0,0). Use this to offset the
            origin of the plane (cartesian coordinates). Defaults to None.

    Returns:
        points_cart (np.ndarray): Points of the plane in cartesian coordinates
        rec (np.ndarray): Array of indexes (int) if the plane is remapped on a grid with shape 'shape_out'
            e.g.:  grid[*rec] will give the values of grid at the points of the plane
    """
    shape_out = np.array(shape_out)
    cell_inv = np.linalg.inv(cell)
    Xpc, Ypc = np.meshgrid(
        np.linspace(0,1, shape_plane[0]),
        np.linspace(0,1, shape_plane[1]),
        indexing='ij'
    )
    cc = np.hstack((v1.reshape(-1,1), v2.reshape(-1,1)))

    points_crys = np.vstack((Xpc.flatten(), Ypc.flatten()))
    points_cart = np.dot(cc, points_crys)

    if not offset is None:
        points_cart += offset.reshape(3,1)

    rec = np.round(
        np.dot(cell_inv, points_cart),
        decimals=8
    )
    rec = rec % 1
    rec = np.array(rec) * shape_out.reshape(3,1)
    rec = rec.astype(dtype='int')

    return points_cart, rec
