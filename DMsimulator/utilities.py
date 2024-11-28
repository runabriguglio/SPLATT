import numpy as np
import matplotlib.pyplot as plt

# from scipy.sparse import csr_matrix
from segmented_deformable_mirror import SegmentedMirror
from hexagonal_geometry import HexGeometry
from matrix_calculator import matmul
import read_and_write_fits as myfits


def dm_system_setup(TN:str, n_global_zern:int = 11, n_local_zern:int = 11):
    """
    Performs the basic operations to define a segmented
    mirror object from the data in the configuration file

    Parameters
    ----------
    TN : string
        Configuration file tracking number.
    n_global_zern : int, optional
        Number of Zernike modes for the IM on the global mask. The default is 7.
    n_local_zern : int, optional
        Number of Zernike modes for the IM on the single segment. The default is 15.

    Returns
    -------
    sdm : segmented deformable mirror class
        Segmented deformable mirror object.

    """
    print('Defining mirror geometry ...')
    geom = HexGeometry(TN)
    
    print('Initializing segmented mirror class ...')
    sdm = SegmentedMirror(geom)
    
    # Global Zernike matrix
    print('Computing ' + str(n_global_zern) + ' modes global Zernike interaction matrix ...')
    sdm.compute_global_zern_matrix(n_global_zern)
    tiptilt = np.zeros(n_global_zern)
    tiptilt[1] = 1
    tiptilt[2] = 1
    wf = matmul(sdm.glob_ZM,tiptilt)
    sdm.surface(wf, 'Global Tip/Tilt')
    
    # Local Zernike matrix
    print('Computing ' + str(n_local_zern) + ' modes local Zernike interaction matrix ...')
    sdm.compute_local_zern_matrix(n_local_zern)
    INTMAT = sdm.ZM
    n_hex = np.shape(INTMAT)[0]
    cmd_ids = np.arange(n_local_zern-1)+1
    cmd_ids = np.tile(cmd_ids,int(np.ceil(n_hex/(n_local_zern-1))))
    cmd_ids = cmd_ids[0:n_hex]
    cmd_ids = cmd_ids + n_local_zern*np.arange(n_hex)
    modal_cmd = np.zeros(n_hex*n_local_zern)
    modal_cmd[cmd_ids] = 1
    flat_shape = matmul(INTMAT,modal_cmd)
    sdm.surface(flat_shape, 'Zernike modes')
    
    # Global influence functions and global reconstructor
    print('Initializing influence functions and reconstructor matrices ...')
    sdm.initialize_IFF_and_R_matrices(simulate = True)
    IFF = sdm.IFF
    R = sdm.R
    cmd_for_zern = matmul(R,flat_shape)
    flat_img = matmul(IFF,cmd_for_zern)
    sdm.surface(flat_img, 'Reconstructed Zernike modes')
    
    n_acts = np.shape(IFF)[2]
    cmd = np.zeros(n_hex*n_acts)
    cmd_ids = np.arange(n_acts)
    cmd_ids = np.tile(cmd_ids,int(np.ceil(n_hex/(n_acts))))
    cmd_ids = cmd_ids[0:n_hex]
    cmd_ids = cmd_ids + n_acts*np.arange(n_hex)
    cmd[cmd_ids] = 1 #np.ones(n_hex)
    flat_img = matmul(IFF,cmd)
    sdm.surface(flat_img, 'Actuators influence functions')
    
    return sdm


def fitting_error_plots(mask, IM, IFF, R, pix_ids = None):
    """
    Computes the fitting error for the reconstructor on a mask,
    plotting the N = np.shape(IM)[1] residual images

    Parameters
    ----------
    mask : bool ndarray [Npix,]
        Boolean mask of the flattened image.
    IM : float ndarray [Npix,Nmodes]
        Interaction matrix of the images to fit.
    IFF : float ndarray [Npix,Nacts]
        Influence functions matrix from actuators data.
    R : float ndarray [Nacts,Npix]
        Reconstructor matrix (IFF pseudo-inverse).
    pix_ids : ndarray(int), optional
        The indices of the valid pixel indices on the mask. The default is None.

    Returns
    -------
    RMS_vec : float ndarray [Nmodes,]
        Vector of ratios of the STDs of the reconstructed
        image and the desired Zernike mode from the IM.

    """
    
    N = np.shape(IM)[-1]
    RMS_vec = np.zeros(N)
    
    flat_img = np.zeros(np.size(mask))
    flat_mask = mask.flatten()
    # cube = np.zeros([N_modes,np.shape(mask)[0],np.shape(mask)[1]])
    
    for k in range(N):
        des_shape = IM[:,k]
        act_cmd = matmul(R,des_shape)
        act_shape = matmul(IFF,act_cmd)
        shape_err = des_shape-act_shape
        
        # if isinstance(shape_err, csr_matrix):
        #     shape_err = (shape_err).toarray()
        #     shape_err = shape_err[:,0]
            
        shape_RMS = np.std(des_shape)
        if shape_RMS < 1e-15: # avoid division by zero
            shape_RMS = 1
        RMS_vec[k] = np.std(shape_err)/shape_RMS
        
        if pix_ids is None:
            pix_ids = ~flat_mask
            
        flat_img[pix_ids] = shape_err
        img = np.reshape(flat_img, np.shape(mask))
        img = np.ma.masked_array(img, mask)
        # cube[k] = img
        plt.figure()
        plt.imshow(img, origin = 'lower', cmap = 'inferno')
        plt.colorbar()
        plt.title('Mode ' + str(k+1) + ' shape error\n RMS: ' + str(RMS_vec[k]) )
        
    plt.figure()
    plt.plot(RMS_vec,'o')
    plt.xlabel('Mode index')
    plt.ylabel('Shape RMS')
    plt.title('Fitting error')
    plt.grid('on')
    
    return RMS_vec


def fitting_error(IM, IFF, R):
    """
    Computes the fitting error for the reconstructor

    Parameters
    ----------
    IM : float ndarray [Npix,Nmodes]
        Interaction matrix of the images to fit.
    IFF : float ndarray [Npix,Nacts]
        Influence functions matrix from actuators data.
    R : float ndarray [Nacts,Npix]
        Reconstructor matrix (IFF pseudo-inverse).

    Returns
    -------
    RMS_vec : float ndarray [Nmodes,]
        Vector of ratios of the STDs of the reconstructed
        image and the desired image from the IM.

    """
    
    act_cmds = matmul(R,IM)
    act_shapes = matmul(IFF,act_cmds)
    res_shapes = IM - act_shapes
    
    img_rms = np.std(IM,axis = 0)
    img_rms[img_rms < 1e-15] = np.ones(len(img_rms[img_rms < 1e-15])) # avoid division by 0
    RMS_vec = np.std(res_shapes,axis = 0)/img_rms
    
    return RMS_vec


def segment_scramble(sdm, mode_amp = 10e-6, apply_shape:bool = False):
    """
    Applies a random shape to all segments using
    a random linear combination of Zernike modes,
    scaled by the inverse of the Noll number

    Parameters
    ----------
    sdm : 
        egmented deformable mirror class.
    mode_amp : float, optional
        Amplitude of the segment scramble. The default is 10e-6.
    apply_shape : bool, optional
        Wheter to apply the shape to the sdm or not. The default is False.

    Returns
    -------
    None.

    """
    file_name = sdm.geom.savepath + 'initial_segment_scramble.fits'
    
    try:
        flat_img = myfits.read_fits(file_name)
    except FileNotFoundError:
        Nsegments = sdm.geom.n_hex
        ZMat = sdm.ZM
        
        # Retrieve number of modes from the interaction matrix
        n_modes = np.shape(ZMat)[-1]
        
        # Generate random mode coefficients
        mode_vec = np.random.randn(Nsegments*n_modes)
        
        # Probability inversely proportional to spatial frequency
        m = int(np.ceil((np.sqrt(8*n_modes)-1.)/2.))
        freq_vec = np.repeat(np.arange(m)+1,np.arange(m)+1)
        prob_vec = 1./freq_vec[:n_modes]
        prob_vec_rep = np.tile(prob_vec,Nsegments)
        
        # Modulate on the probability
        mode_vec = mode_vec * prob_vec_rep
        
        # Amplitude
        mode_vec *= mode_amp
        
        # Matrix product
        flat_img = matmul(ZMat,mode_vec)
        
        # # Global modes
        # n_glob_modes = np.shape(sdm.glob_int_mat)[1]
        # glob_mode_vec = np.random.randn(n_glob_modes)
        # flat_img += matmul(sdm.glob_ZM,glob_mode_vec)
    
    sdm.surface(flat_img,'Segment scramble')
    
    if apply_shape:
        sdm.shape += flat_img
        myfits.write_to_fits(flat_img, file_name)
    
    



