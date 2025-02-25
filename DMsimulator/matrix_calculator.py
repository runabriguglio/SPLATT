import numpy as np
from tps import ThinPlateSpline # for the simulated IFF
# from scipy.interpolate import griddata

from zernike_polynomials import generate_zernike_matrix as assemble_zern_mat


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

def compute_reconstructor(M, thr:float = 0.):
    """
    Computes the reconstructor (pseudo-inverse) 
    for the interaction matrix M.

    Parameters
    ----------
    M : ndarray(float) [Npix,N]
        Interaction matrix to be inverted.
        
    thr : float, optional
        Threshold for the inverse eigenvalues. 
        The eigenvalues v  s. t. 1/v < thr 
        are discarded when computing the inverse.
        By default, no eigenvalues are discarded.

    Returns
    -------
    Rec : ndarray(float) [N,Npix]
        The pseudo-inverse of M.

    """
    
    U,S,Vh = np.linalg.svd(M, full_matrices=False)
    Sinv = 1/S
    Sinv[Sinv < thr] = 0
    Rec = (Vh.T * Sinv) @ U.T
    return Rec


def simulate_influence_functions(act_coords, local_mask, pix_scale, amps = 1.0):
    """ Simulate the influence functions by 
    imposing 'perfect' zonal commands """
    
    n_acts = len(act_coords)
    
    H,W = np.shape(local_mask)
    
    # pix_coords = np.zeros([max_x*max_y,2])
    # pix_coords[:,0] = np.repeat(np.arange(max_x),max_y)
    # pix_coords[:,1] = np.tile(np.arange(max_y),max_x)
    
    # act_pix_coords = np.zeros([n_acts,2])
    # act_pix_coords[:,0] = (act_coords[:,1] * pix_scale + max_x/2).astype(int)
    # act_pix_coords[:,1] = (act_coords[:,0] * pix_scale + max_y/2).astype(int)
    
    pix_coords = get_pixel_coords(local_mask, pix_scale)
    act_pix_coords = get_pixel_coords(local_mask, pix_scale, act_coords)
    
    img_cube = np.zeros([H,W,n_acts])
    
    if isinstance(amps,float):
        amps *= np.ones(n_acts)

    for k in range(n_acts):
        act_data = np.zeros(n_acts)
        act_data[k] = amps[k]
        tps = ThinPlateSpline(alpha=0.0)
        tps.fit(act_pix_coords, act_data)
        flat_img = tps.transform(pix_coords)
        uint8_img = scale2uint8(flat_img) # scale to uint8
        img_cube[:,:,k] = np.reshape(uint8_img, [H,W])

    # Masked array
    cube_mask = np.tile(local_mask,n_acts)
    cube_mask = np.reshape(cube_mask, np.shape(img_cube), order = 'F')
    cube = np.ma.masked_array(img_cube, cube_mask, dtype=np.uint8)
    
    return cube


def calculate_influence_functions(act_coords, local_mask, normalized_act_radius):
    """ Project the actuator influence functions 
    on the mask via grid interpolation """
    
    raise NotImplementedError()


def cube2mat(cube):
    """ Get influence functions matrix 
    from the image cube """
    
    n_acts = np.shape(cube)[2]
    valid_len = np.sum(1-cube[:,:,0].mask)
    
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


def scale2uint8(img):
    """ Scales an ndarray to uint8 """
    
    img -= min(img)
    img = img/max(img) * (2**8-1)
    img = np.round(img)
    
    return (img).astype(np.uint8)


def get_pixel_coords(mask, pix_scale:float, coords = None):
    """ Convert x,y coordinates in coords to pixel coordinates
    or get the pixel coordinates of a mask

    Parameters
    ----------
    mask : ndarray(bool)
        The image mask where the pixels are.
        
    pix_scale : float
        The number of pixels per meter.
        
    coords : ndarray(float) [N,2], optional
        The N coordinates to convert in pixel coordinates.
        Defaults to all pixels on the mask.

    Returns
    -------
    pix_coords : ndarray(int)
        The obtained pixel coordinates.

    """
    
    H,W = np.shape(mask)
    
    if coords is not None:
        pix_coords = np.zeros([len(coords),2],dtype=int)
        pix_coords[:,0] = np.round(coords[:,1]*pix_scale + H/2)
        pix_coords[:,1] = np.round(coords[:,0]*pix_scale + W/2)
    else:
        pix_coords = np.zeros([H*W,2])
        pix_coords[:,0] = np.repeat(np.arange(H),W)
        pix_coords[:,1] = np.tile(np.arange(W),H)
    
    return pix_coords


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

    pix_act_coords = get_pixel_coords(mask, pix_scale, act_coords)
    pix_capsens_coords = get_pixel_coords(mask, pix_scale, capsens_coords)
    
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
        CSMAT[k,:] = masked_data/pix_area
    
    return CSMAT



# def rescale_img(img, n_pix: int):
    
#     len_x,len_y = np.shape(img)
#     rows = np.repeat(np.arange(len_x),len_y) + max(0,(len_y-len_x)/2)
#     cols = np.tile(np.arange(len_y),len_x) + max(0,(len_x-len_y)/2)
#     xy = np.linspace(0, max(len_x,len_y), n_pix)
#     gx,gy = np.meshgrid(xy,xy)
    
#     rescaled_img = griddata((cols,rows),img.flatten(),(gx,gy))
#     rescaled_img = np.reshape(rescaled_img,[n_pix,n_pix])
    
#     if hasattr(img, 'mask'):
#         rescaled_mask = griddata((cols,rows),img.mask.flatten(),(gx,gy))
#         rescaled_mask = np.reshape(rescaled_mask,[n_pix,n_pix])
#         rescaled_img = np.ma.masked_array(rescaled_img, rescaled_mask)
    
#     return rescaled_img
    




    