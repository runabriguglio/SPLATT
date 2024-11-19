import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from segmented_mirror import SegmentedMirror
# from read_and_write_fits import write_to_fits


def dm_system_setup(TN):
    """
    Performs the basic operations to define a segmented
    mirror object from the data in the configuration file

    Parameters
    ----------
    TN : string
        Configuration file tracking number.

    Returns
    -------
    sdm : segmented deformable mirror class
        Segmented deformable mirror object.

    """
    
    # Build segmented deformable mirror
    sdm = SegmentedMirror(TN)
    
    # Global interaction matrix
    N_global_modes = 11
    sdm.compute_global_interaction_matrix(N_global_modes)
    glob_INTMAT = sdm.glob_IM
    tiptilt = np.zeros(N_global_modes)
    tiptilt[1] = 1
    tiptilt[2] = 1
    wf = glob_INTMAT * tiptilt
    sdm.plot_wavefront(wf, 'Global Tip/Tilt')
    
    # Interaction matrix
    N_modes = 11
    sdm.compute_interaction_matrix(N_modes)
    INTMAT = sdm.IM
    n_hex = int(np.shape(INTMAT)[1]/N_modes)
    cmd_ids = np.arange(N_modes-1)+1
    cmd_ids = np.tile(cmd_ids,int(np.ceil(n_hex/(N_modes-1))))
    cmd_ids = cmd_ids[0:n_hex]
    cmd_ids = cmd_ids + N_modes*np.arange(n_hex)
    modal_cmd = np.zeros(n_hex*N_modes)
    modal_cmd[cmd_ids] = 1
    flat_shape = INTMAT * modal_cmd
    sdm.plot_wavefront(flat_shape, 'Zernike modes')
    
    # Global influence functions and global reconstructor
    sdm.assemble_IFF_and_R_matrices(simulated_IFFs = True)
    IFF = sdm.IFF
    R = sdm.R
    cmd_for_zern = R * flat_shape
    flat_img = IFF * cmd_for_zern
    sdm.plot_wavefront(flat_img, 'Reconstructed Zernike modes')
    
    cmd = np.zeros(np.shape(IFF)[1])
    cmd_ids = np.arange(n_hex)*n_hex + np.arange(n_hex)
    cmd[cmd_ids] = np.ones(n_hex)
    flat_img = IFF * cmd
    sdm.plot_wavefront(flat_img, 'Actuators influence functions')

    return sdm


def fitting_error(mask, IM, IFF, R):
    """
    Computes the fitting error for the reconstructor on a mask

    Parameters
    ----------
    mask : bool ndarray [Npix,]
        Boolean mask of the flattened image.
    IM : float ndarray [Npix,Nmodes]
        Interaction matrix from the Zernike modes.
    IFF : float ndarray [Npix,Nacts]
        Influence functions matrix from actuators data.
    R : float ndarray [Nacts,Npix]
        Reconstructor matrix (IFF pseudo-inverse).

    Returns
    -------
    RMS_vec : float ndarray [Nmodes,]
        Vector of ratios of the STDs of the reconstructed
        image and the desired Zernike mode from the IM.

    """
    
    N_modes = np.shape(IM)[1]
    RMS_vec = np.zeros(N_modes)
    
    flat_img = np.zeros(np.size(mask))
    flat_mask = mask.flatten()
    # cube = np.zeros([N_modes,np.shape(mask)[0],np.shape(mask)[1]])
    
    for k in range(N_modes):
        des_shape = IM[:,k]
        act_cmd = matmul(R,des_shape)
        act_shape = matmul(IFF,act_cmd)
        shape_err = des_shape-act_shape
        
        if isinstance(shape_err, csr_matrix):
            shape_err = (shape_err).toarray()
            shape_err = shape_err[:,0]
        # if isinstance(des_shape, csr_matrix):
        #     des_shape = (des_shape).toarray()
        #     des_shape = des_shape[:,0]
        RMS_vec[k] = np.std(shape_err) #/np.std(des_shape)
        
        flat_img[~flat_mask] = shape_err
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


def segment_scramble(sdm, mode_amp = 10e-6):
    """
    Applies a random shape to all segments using
    a random linear combination of Zernike modes,
    scaled by the inverse of the Noll number

    Parameters
    ----------
    sdm : segmented deformable mirror class object
        The segmented mirror to apply the sgment scramble to.
    mode_amp : float, optional
        The modal scale amplitude. Default is 10e-6.

    Returns
    -------
    None.

    """
    
    n_hex = len(sdm.hex_centers)
    
    # Retrieve number of modes from the interaction matrix
    n_modes = int(np.shape(sdm.IM)[1]/n_hex)
    
    # Generate random mode coefficients
    mode_vec = np.random.randn(n_hex*n_modes)
    
    # Probability inversely proportional to spatial frequency
    m = int(np.ceil((np.sqrt(8*n_modes)-1.)/2.))
    freq_vec = np.repeat(np.arange(m)+1,np.arange(m)+1)
    prob_vec = 1./freq_vec[:n_modes]
    prob_vec_rep = np.tile(prob_vec,n_hex)
    
    # Modulate on the probability
    mode_vec = mode_vec * prob_vec_rep
    
    # Amplitude
    mode_vec *= mode_amp
    
    # Matrix product
    flat_img = sdm.IM*mode_vec
    
    # # Global modes
    # n_glob_modes = np.shape(sdm.glob_int_mat)[1]
    # glob_mode_vec = np.random.randn(n_glob_modes)
    # flat_img += sdm.glob_int_mat*glob_mode_vec
    
    # Save values in segments' shape
    for k in range(n_hex):
        row_ids = sdm.valid_ids[k]
        sdm.segment[k].shape = flat_img[row_ids]
    
    # Savr values in globla shape
    sdm.shape = flat_img
    
    # Plot the result
    sdm.plot_wavefront(flat_img, 'Segment scramble')
    
    


def matmul(matrix, vector):
    """
    Simple function to perform matrix multiplication
    for both sparse and regular matrices

    Parameters
    ----------
    matrix : TYPE
        DESCRIPTION.
    vector : TYPE
        DESCRIPTION.

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    """
    
    if isinstance(matrix, csr_matrix):
       res = matrix * vector
    else:
        res = matrix @ vector
        
    return res
