import numpy as np
from tps import ThinPlateSpline # for the simulated IFF
from scipy.interpolate import griddata

from zernike_polynomials import generate_zernike_matrix as assemble_zern_mat
from hexagonal_geometry import semi_structured_point_cloud


def matmul(mat, vec):
    """
    Simple function to perform matrix multiplication
    for both block matrices (saved as an array of
    matrices) and regular matrices
    """
    
    if len(np.shape(mat)) == 3: # Array of matrices
        id_ctr = 0 
        n_hex, n_pix, N = np.shape(mat)
        res = np.zeros(n_pix*n_hex)
        for n, mat_n in enumerate(mat):
            ids = np.arange(id_ctr, id_ctr+N)
            res[n_pix*n:n_pix*(n+1)] = mat_n @ vec[ids]
            id_ctr += N
    else:
        res = mat @ vec
        
    return res

def define_capsens_matrix(mask, pix_scale, act_coords, r_in, r_out, capsens_coords = None):
    """
    Determine the pixel CapSens matrix to estimate the gap measured by
    the capacitive sensors for a give shape as meas_gap = CSMAT @ shape

    Parameters
    ----------
    mask : ndarray(bool)
        The mask representing the mirror on which to operate.
    pix_scale : float
        The number of pixels per meter.
    act_coords : ndarray(float) [Nacts,2]
        The local actuator coordinates (in meters).
    r_in : float
        The CapSens inner radius (in meters).
    r_out : float
        The CapSens outer radius (in meters).
    capsens_coords : ndarray(float) [Nacts,2], optional
        The local actuator coordinates (in meters).
        Useful if there is an eccentricity between 
        CapSens and actuators. The default is act_coords.

    Returns
    -------
    CSMAT : ndarray(float) [Nacts,Npixels]
        The obtained CapSens pixel matrix, normalized by the
        pixel area (i.e. number of pixels per CapSens)

    """
    
    if capsens_coords is None:
        capsens_coords = act_coords
        
    X,Y = np.shape(mask)
    d = lambda i,j: np.sqrt(i**2+j**2)

    # Rounding to int is less accurate but gives 
    # more consistent results in terms of pixel area
    pix_act_coords = (act_coords*pix_scale + np.array([Y,X])/2).astype(int)
    pix_capsens_coords = (capsens_coords*pix_scale + np.array([Y,X])/2).astype(int)
    
    pix_in = int(r_in*pix_scale)
    pix_out = int(r_out*pix_scale)
    
    CSMAT = np.zeros([len(act_coords),np.sum(1-mask)])
    
    for k, pix_act_coord in enumerate(pix_act_coords):
        act_x, act_y = pix_act_coord
        sens_x, sens_y = pix_capsens_coords[k,:]
        sensor = np.fromfunction(lambda i,j: (d(i-act_y,j-act_x) >= pix_in) * (d(i-sens_y,j-sens_x) < pix_out), [X,Y])
        # img = np.ma.masked_array(sensor,mask), plt.figure(), plt.imshow(img, origin = 'lower')
        
        masked_data = sensor[~mask]
        pix_area = np.sum(masked_data)
        # print(pix_area) # debug
        CSMAT[k,:] = masked_data/pix_area
    
    return CSMAT


def compute_reconstructor(M):
    return np.linalg.pinv(M)


def simulate_influence_functions(act_coords, local_mask, pix_scale, amp = 1.0):
    """ Simulate the influence functions by 
    imposing 'perfect' zonal commands """
    
    n_acts = len(act_coords)
    
    max_x, max_y = np.shape(local_mask)
    
    pix_coords = np.zeros([max_x*max_y,2])
    pix_coords[:,0] = np.repeat(np.arange(max_x),max_y)
    pix_coords[:,1] = np.tile(np.arange(max_y),max_x)
    
    act_pix_coords = np.zeros([n_acts,2])
    act_pix_coords[:,0] = (act_coords[:,1] * pix_scale + max_x/2).astype(int)
    act_pix_coords[:,1] = (act_coords[:,0] * pix_scale + max_y/2).astype(int)
    
    img_cube = np.zeros([max_x,max_y,n_acts])

    for k in range(n_acts):
        act_data = np.zeros(n_acts)
        act_data[k] = amp
        tps = ThinPlateSpline(alpha=0.0)
        tps.fit(act_pix_coords, act_data)
        flat_img = tps.transform(pix_coords)
        img_cube[:,:,k] = np.reshape(flat_img, [max_x,max_y])

    # Masked array
    cube_mask = np.tile(local_mask,n_acts)
    cube_mask = np.reshape(cube_mask, np.shape(img_cube), order = 'F')
    cube = np.ma.masked_array(img_cube,cube_mask)
    
    return cube


def calculate_influence_functions(act_coords, local_mask, normalized_act_radius):
    """ Project the actuator influence functions 
    on the mask via grid interpolation """
    
    pps = 9
    mesh = _define_mesh(act_coords, pps, normalized_act_radius)
    iff_data = _compute_iff_data(mesh)
    
    n_acts = len(act_coords)
    
    # Grid interpolation
    x, y = act_coords[:, 0], act_coords[:, 1]
    npix_x, npix_y = np.shape(local_mask)  
    new_x = np.linspace(min(x), max(x), npix_y) # x coordinate is the img column!
    new_y = np.linspace(min(y), max(y), npix_x) # y coordinate is the img row!
    gx, gy = np.meshgrid(new_x, new_y)
    img_cube = griddata((x, y), iff_data, (gx, gy), fill_value = 0., method = 'linear')
    
    # Masked array
    cube_mask = np.tile(local_mask,n_acts)
    cube_mask = np.reshape(cube_mask, np.shape(img_cube), order = 'F')
    cube = np.ma.masked_array(img_cube,cube_mask)
    
    return cube


def cube2mat(cube):
    """ Get influence functions matrix 
    from the image cube """
    
    n_acts = np.shape(cube)[2]
    valid_len = np.sum(1-cube[:,:,0].mask)
    
    # Save valid data to IFF (full) matrix
    flat_cube = cube.data[~cube.mask]
    local_IFF = np.reshape(flat_cube, [valid_len, n_acts])
    
    IFF = np.array(local_IFF)
    
    return IFF
    


def compute_zernike_matrix(mask, n_modes):
    """ Computes the zernike matrix: 
        [n_pixels,n_modes] """
    
    noll_ids = np.arange(n_modes) + 1
    mat = assemble_zern_mat(noll_ids, mask)
    
    return mat   
    

def _define_mesh(act_coords, points_per_side, normalized_act_radius):
    """ Defines the local mesh point coordinates
    on the segment starting from a semi-structured
    point cloud and adding the act locations to it"""
    
    point_cloud = semi_structured_point_cloud(points_per_side)
    
    # Add points on the actuator locations
    n_acts = len(act_coords)
    flat_act_coords = np.tile(act_coords,5).flatten()
    up_down_left_right = np.array([[0,0],[0,1],[0,-1],[-1,0],
                                   [1,0]]).flatten()
    UDLR = np.tile(up_down_left_right,n_acts)
    act_points = flat_act_coords + UDLR*normalized_act_radius
    act_points = np.reshape(act_points,[n_acts*5,2])
    
    mesh_points = np.concatenate((point_cloud, act_points))
    
    local_mesh_coords = mesh_points
    
    return local_mesh_coords


def _compute_iff_data(mesh):
    raise NotImplementedError("Function not yet implemented! Use the spline interpolation for IFF")


    

    