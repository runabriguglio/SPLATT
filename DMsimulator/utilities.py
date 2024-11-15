import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from segmented_mirror import SegmentedMirror


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



def dm_system_setup(TN):
    
    # Build segmented deformable mirror
    dm = SegmentedMirror(TN)
    
    # Plot global mask
    plt.figure()
    plt.imshow(dm.global_mask, origin = 'lower', cmap='gray')
    plt.title('Global Mask')
    
    # Global interaction matrix
    N_global_modes = 11
    dm.compute_global_interaction_matrix(N_global_modes)
    glob_INTMAT = dm.glob_IM
    tiptilt = np.zeros(N_global_modes)
    tiptilt[1] = 1
    tiptilt[2] = 1
    wf = glob_INTMAT * tiptilt
    dm.plot_wavefront(wf, 'Global Tip/Tilt')
    
    # Interaction matrix
    N_modes = 11
    dm.compute_interaction_matrix(N_modes)
    INTMAT = dm.IM
    n_hex = int(np.shape(INTMAT)[1]/N_modes)
    cmd_ids = np.arange(N_modes-1)+1
    cmd_ids = np.tile(cmd_ids,int(np.ceil(n_hex/(N_modes-1))))
    cmd_ids = cmd_ids[0:n_hex]
    cmd_ids = cmd_ids + N_modes*np.arange(n_hex)
    modal_cmd = np.zeros(n_hex*N_modes)
    modal_cmd[cmd_ids] = 1
    flat_img = INTMAT * modal_cmd
    dm.plot_wavefront(flat_img, 'Zernike modes')
    
    # Initial segment scramble
    dm.segment_scramble()
    w_scramble = dm.shape.copy()
    dm.plot_wavefront(w_scramble, 'Segment scramble')
    
    # Actuator coordinates
    plt.figure()
    plt.scatter(dm.act_coords[:,0],dm.act_coords[:,1],s=1)
    plt.axis('equal')
    plt.title('Global actuator locations')
    
    # Global influence functions and global reconstructor
    dm.assemble_IFF_and_R_matrices()
    
    return dm


def fitting_error(mask, IM, IFF, R):
    
    N_modes = np.shape(IM)[1]
    RMS_vec = np.zeros(N_modes)
    
    for k in range(N_modes):
        des_shape = IM[:,k]
        act_cmd = matmul(R,des_shape)
        act_shape = matmul(IFF,act_cmd)
        shape_err = des_shape-act_shape
        if isinstance(shape_err, csr_matrix):
            shape_err = (shape_err).toarray()
            shape_err = shape_err[:,0]
        RMS_vec[k] = np.std(shape_err)
        
        img = np.zeros(np.size(mask))
        flat_mask = mask.flatten()
        img[~flat_mask] = shape_err
        img = np.reshape(img, np.shape(mask))
        img = np.ma.masked_array(img, mask)
        plt.figure()
        plt.imshow(img, origin = 'lower', cmap = 'inferno')
        plt.colorbar()
        plt.title('Mode ' + str(k) + ' shape error\n RMS: ' + str(RMS_vec[k]) )
        
    plt.figure()
    plt.plot(RMS_vec,'o')
    plt.xlabel('Mode index')
    plt.ylabel('Shape RMS')
    plt.title('Fitting error')
    plt.grid('on')
    
    return RMS_vec



# # Plot IFF data on segments
# n_hex = len(dm.hex_valid_ids)
# glob_img = np.zeros(np.size(dm.global_mask))
# point_glob_img = np.zeros(np.size(dm.global_mask))
# for kk in range(n_hex):
#     glob_img[dm.hex_valid_ids[kk]] = hexA.sim_IFF[:,kk]#loc_IFF[:,kk]
#     # point_glob_img[dm.hex_valid_ids[kk]] = hexA.IFF[:,kk]
    
# # dm.plot_wavefront(glob_img, 'Actuator Influence Functions')
# full_img = np.reshape(glob_img, np.shape(dm.global_mask))
# img = np.ma.masked_array(full_img, dm.global_mask)
# plt.figure()
# plt.imshow(img, origin = 'lower', cmap='inferno')
# # plt.axis([1215,1620,1220,1600])
# plt.colorbar()


