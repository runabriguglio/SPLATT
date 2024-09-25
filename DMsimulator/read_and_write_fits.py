from astropy.io import fits
import numpy as np
from scipy.sparse import csr_matrix


def write_to_fits(data, file_path):
    
    if isinstance(data,list): # list (e.g. sparse matrix)
        fits.writeto(file_path, data[0], overwrite=True)
        for vec in data[1:]:
            fits.append(file_path, vec)  
            
    else:
        # Print to a .fits file
        fits.writeto(file_path, data, overwrite=True)
        
        if hasattr(data,'mask'): # masked array
            fits.append(file_path, (data.mask).astype(np.uint8))
            
            
def read_fits_file(file_path, is_bool = False, is_ma = False,
                   list_len = 1, sparse_shape = None):
    
    with fits.open(file_path) as hdu:
        
        if is_bool: # Boolean array
            data_out = np.array(hdu[0].data).astype(bool)
            
        elif list_len > 1: # List of arrays
            data_out = []
            for i in range(list_len):
                data_out.append(np.array(hdu[0].data))
            
        elif sparse_shape is not None: # sparse matrix
            mat_data = hdu[0].data
            indices = hdu[1].data
            indptr = hdu[2].data
            data_out = csr_matrix((mat_data,indices,indptr), sparse_shape)
     
        elif is_ma: # masked array
            img = np.array(hdu[0].data)
            img_mask = np.array(hdu[1].data).astype(bool)
            data_out = np.ma.masked_array(img,mask=img_mask)
                
        else: # default
            data_out = np.array(hdu[0].data)
                
        
    return data_out