from astropy.io import fits
import numpy as np
from scipy.sparse import csr_matrix


# Writing files
def write_to_fits(data, file_path):
    
    if isinstance(data,list): # list (e.g. sparse matrix)
    
        fits.writeto(file_path, data[0], overwrite=True)
        
        if hasattr(data[0],'mask'): # masked array
            fits.append(file_path, (data[0].mask).astype(np.uint8))
        
        for vec in data[1:]:

            if hasattr(vec,'mask'): # masked array
                fits.append(file_path, vec.data) 
                fits.append(file_path, (vec.mask).astype(np.uint8))
                
            else:
                fits.append(file_path, vec) 
            
    else:
        # Print to a .fits file
        fits.writeto(file_path, data, overwrite=True)
        
        if hasattr(data,'mask'): # masked array
            fits.append(file_path, (data.mask).astype(np.uint8))
            
            
            
# Reading files           
def read_fits_file(file_path, idx, is_bool, is_ma, sparse_shape):
    
    with fits.open(file_path) as hdu:
        
        if is_bool: # Boolean array
            data_out = np.array(hdu[idx].data).astype(bool)
            
        elif sparse_shape is not None: # sparse matrix
            mat_data = hdu[idx*3].data
            indices = hdu[idx*3+1].data
            indptr = hdu[idx*3+2].data
            data_out = csr_matrix((mat_data,indices,indptr), sparse_shape)
     
        elif is_ma: # masked array
            img = np.array(hdu[idx*2].data)
            img_mask = np.array(hdu[idx*2+1].data).astype(bool)
            # img_mask = np.array(hdu[idx].mask).astype(bool)
            data_out = np.ma.masked_array(img,mask=img_mask)
                
        else: # default
            data_out = np.array(hdu[idx].data)
            
    return data_out
                    
            
    
def read_fits(file_path, list_len = 1, is_bool = False,
                   is_ma = False, sparse_shape = None):
    
    data_out = []
    
    for k in range(list_len):
        out = read_fits_file(file_path, k, is_bool, is_ma, sparse_shape)
        data_out.append(out)
        
    if list_len == 1:
        data_out = data_out[0]
                
        
    return data_out
            
            
# def read_fits_file(file_path, is_bool = False, is_ma = False,
#                    list_len = 1, sparse_shape = None):
    
#     with fits.open(file_path) as hdu:
        
#         if is_bool: # Boolean array
#             data_out = np.array(hdu[0].data).astype(bool)
            
#         elif list_len > 1: # List of arrays
#             data_out = []
#             for i in range(list_len):
#                 data_out.append(np.array(hdu[0].data))
            
#         elif sparse_shape is not None: # sparse matrix
#             mat_data = hdu[0].data
#             indices = hdu[1].data
#             indptr = hdu[2].data
#             data_out = csr_matrix((mat_data,indices,indptr), sparse_shape)
     
#         elif is_ma: # masked array
#             img = np.array(hdu[0].data)
#             img_mask = np.array(hdu[1].data).astype(bool)
#             data_out = np.ma.masked_array(img,mask=img_mask)
                
#         else: # default
#             data_out = np.array(hdu[0].data)
                
        
#     return data_out