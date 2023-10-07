from itertools import product
from typing import Annotated, Literal

import numpy as np
import numpy.typing as npt


def recipr_base(base: npt.ArrayLike) -> npt.ArrayLike:
    """Return the reciprocal base of the given lattice.
    Makes use of the relation that A^T . B = 2*pi * I  where I is the identity matrix

    Args:
        base (npt.ArrayLike): Matrix of basis vectors (as row or columns). If the array are given as row/column 
            the returned array will also be row/column.

    Returns:
        npt.ArrayLike: Matrix of reciprocal basis vectors as a row/column depending on the input.
    """
    return np.linalg.inv(base).T * 2 * np.pi

def generate_repetition_grid(
        R: list[list[int]], 
        vect_matrix: Annotated[npt.ArrayLike, Literal['N', 'N']] = None
    ) -> Annotated[npt.NDArray[np.float64], Literal['M', 'N']]:
    """Generate a repetition grid using the cartesian product of the given list of repetition vectors.
    Optionally, the repetition vectors can be multiplied by a matrix (eg: the reciprocal base) to obtain
    a grid in cartesian coordinates.
    Example: R = [[0,1,2], [0,1,2], [0,1,2]]

    Args:
        R (list[list[int]]): List of list of indexes where the i-th list contains the indexes of the i-th axis
        vect_matrix (npt.ArrayLike, optional): ArrayLike of shape N*N where (N = number of axes) where the basis
            vector are the rows of the matrix. Defaults to None.

    Returns:
        npt.NDArray[np.float64]: Matrix of output vectors. The shape is (len(R[0])*len(R[1])*...*len(R[N-1]), N)
            The vectors are the rows of the matrix and will be in cartesian coordinates if vect_matrix is not given.
    """
    res = np.array(list(product(*R)))
    if not vect_matrix is None:
        res = res.dot(vect_matrix)

    return res

def xyz_mesh(
        shape: tuple[int, int, int],
        base: Annotated[npt.NDArray[np.float64], Literal[3, 3]] = None,
        rep: int | tuple[int, int, int] = 1
        ) -> Annotated[npt.NDArray[np.float64], Literal[3, 'shape[0]', 'shape[1]', 'shape[2]']]:
    """Generate an omogeneous XYZ meshgrid of points, optionally in cartesian coordinates.

    Args:
        shape (tuple[int, int, int]): (n1, n2, n3) where n1,n2,n3 are the number of points along the 3 axis.
        base (npt.NDArray, optional): Matrix of basis vectors as rows (3*3). Defaults to None.
        rep (int | tuple[int], optional): Number of repetitions. The i-th coordinate of the points will be indexed from
           ` 0 -> shape[i]*rep[i]` if rep is a tuple else from `0 to shape[i]*rep` for all i if rep is an integer. 
           Defaults to 1.

    Returns:
        npt.NDArray: Array of shape (3, *shape) containing the coordinates of the points (in cartesian coordinates if
            base is given).
    """
    n1,n2,n3 = shape
    try:
        r1,r2,r3 = rep
    except (TypeError, ValueError):
        r1 = r2 = r3 = rep
    a = np.linspace(0, r1, n1*r1 + 1)[1:] #[:-1] #+ .5/n1
    b = np.linspace(0, r2, n2*r2 + 1)[1:] #[:-1] #+ .5/n2
    c = np.linspace(0, r3, n3*r3 + 1)[1:] #[:-1] #+ .5/n3

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

    if not base is None:
        XYZ  = np.dot(
            base.T,
            [x.flatten(),y.flatten(),z.flatten()],
            )
    else:
        XYZ = [x,y,z]

    return np.array(XYZ).reshape(3,*x.shape)

def lowdin_ortho(states: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
    """
    Run the Lowdin orthogonalization on the given base.
        https://arxiv.org/abs/1105.3571v1
        O = (M . M^T)^{-1/2} . M
    Args:
        base: np.ndarray with states along axis=0 (as rows). If the states are
             multidimensional (eg: FFT grids), flatten them before. After the
             orthonormalization, restore the shape
    """
    from scipy.linalg import inv, sqrtm

    shape   = states[0].shape
    # states    = np.array([a.flatten() for a in states])
    states = states.reshape(states.shape[0], -1)
    overlap = np.dot(states.conj(), states.T)
    # print('OVERLAP: \n', overlap)

    states = np.dot(
        inv(sqrtm(overlap)).conj(),
        states
        )

    # overlap = np.dot(states.conj(), states.T)
    # print('OVERLAP: \n', overlap)

    return states.reshape(-1,*shape)

def remap_plane(
    cell: Annotated[npt.ArrayLike, Literal[3, 3]],
    Xlim: tuple[int, int], Ylim: tuple[int, int], Zlim: tuple[int, int],
    shape: tuple[int, int, int], 
    rep: tuple[int, int] = (1,1), 
    fixaxis: int = 2
    ) -> tuple[
        npt.NDArray[np.float64], npt.NDArray[np.float64], 
        tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.int32]]
        ]:
    """Remap a plane onto grid defined by the cell vectors 'cell' such that cell = [a1,a2,a3].
    The plane is defined by 2 of the 3 basis vector of the grid when keeping the other constant.

    Args:
        cell (Annotated): ArrayLike of the cell vectors as rows [a1, a2, a3]
        Xlim (tuple[int, int]): tuple of min/max along the X axis for the plane (cartesian coordinates).
        Ylim (tuple[int, int]): tuple of min/max along the Y axis for the plane (cartesian coordinates).
        Zlim (tuple[int, int]): tuple of min/max along the Z axis for the plane (cartesian coordinates).
        shape (tuple[int, int, int]): shape of the grid. !!! This should match the shape of the grid you want to remap
            the plane on. For example you have a volume_data_grid as a result of `data = F(x,y,z)` and you want to collect
            data along a specific plane (e.g. to extract an isocontour)
        rep (tuple[int, int], optional): Number repetitions of the plane axis (e.g. you want the plane defined not on 
            `a1,a2` but `2*a1,3*a2`). Defaults to (1,1).
            Using rep will not change the density of the points on the plane.
        fixaxis (int, optional): Index of the axis to be kept constant for the plane. The value of this constant is not 
         necessarily 0, but will be given by np.average(?lim) where ? is the fixed axis. Defaults to 2 (Z).

    Returns:
        npt.NDArray : Array of shape (shape[0]*rep[0], shape[1]*rep[1]) containing the X cartesian coordinates of
            the points of the plane.
        npt.NDArray : Array of shape (shape[0]*rep[0], shape[1]*rep[1]) containing the Y cartesian coordinates of
        tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.int32]] : Tuple of 3 arrays containing the 
            indexes of the points of the plane in the grid defined by 'shape' and 'rep'. 
            The indexes are given in the order (i,j,k)
            E.G.: grid[*rec] will give the values of grid at the points of the plane
    """
    lim = (Xlim, Ylim, Zlim)
    shape = np.array(shape)

    m1 = np.linspace(*lim[(fixaxis+1)%3], shape[0] * rep[0])
    m2 = np.linspace(*lim[(fixaxis+2)%3], shape[1] * rep[1])
    a,b = np.meshgrid(m1,m2)

    fix = np.ones(a.shape) * (sum(lim[fixaxis]))/2
    rec = np.round(
            np.dot(
            np.linalg.inv(cell.T),
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
