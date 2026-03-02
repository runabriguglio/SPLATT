import numpy as np
import time
from datetime import datetime

from devices.webDAQ import WebDAQ as wbdq
import acceleration_analysis as sp
from devices.wavegenerators import WaveGenerator

import opticalib
import matplotlib.pyplot as plt
import os

from opticalib.ground import osutils as osu


# Connect to WebDAQ
webdaq = wbdq()
webdaq.connect()

# Connect to wavegenerator
wg = WaveGenerator()

def perform_sweep(fmin,fmax,ch=1,ampCh=1,duration=10):

    wg.set_wave(ch,ampl=ampCh,offs=0,freq=fmin,wave_form='SIN')
    wg.set_sweep(ch,fmin,fmax,duration,amp=ampCh,return_time=0)

    # Start acquisition and sweep
    webdaq.start_schedule()
    time.sleep(2)
    wg.trigger_sweep(ch)

    webdaq.stop_schedule_when_job_ends(job_id = 0)

    tnnow = datetime.now()
    tn = tnnow.strftime('%Y%m%d_%H%M%S')
    sp.wdsync()
    wdf = sp.last_wdfile(ext='PI')
    sp.savefile(wdf, tn)
    return tn

#tn = perform_sweep(ch=1,fmin=1500,fmax=8000)

################ CAMERA ########################
path = '/mnt/libero/S331sl/'

def show_frame(frame,title:str=''):
    plt.figure()
    plt.imshow(frame,origin='lower',cmap='RdBu')
    plt.colorbar()
    plt.title(title)
    plt.show()

def save_data(frames):
    tn = osu.newtn()
    osu.save_fits(os.path.join(path,tn)+'.fits',frames,overwrite=True)
    print(tn)
    return tn
	
def load_data(tn,show:bool=True):
    data = osu.load_fits(os.path.join(path,tn)+'.fits')
    if show:
        show_frame(data)
    return data

def stability_measurement(cam,N:int=20,Nframes:int=1,pause:float=50):
    tn_list=[]
    for j in range(N):
        fr=cam.acquire_frames(Nframes)
        tn=save_data(fr)
        tn_list.append(tn)
        time.sleep(pause)
    return tn_list

def connect_to_camera():
    cam = opticalib.devices.cameras.AVTCamera(name='GT2460')
    return cam


##################### IM ########################
def best_circle_fit(image, thr:float=0.5, show:bool=False):
    arr = np.asarray(image, dtype=np.float64)
    maxv = arr.max()
    thresh = thr * maxv
    ys, xs = np.nonzero(arr >= thresh)
    x = xs.astype(np.float64)
    y = ys.astype(np.float64)
    A = np.column_stack([2 * x, 2 * y, np.ones_like(x)])
    b = x * x + y * y
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy, c0 = c
    radius = np.sqrt(cx * cx + cy * cy + c0)
    if show:
        plt.figure()
        plt.imshow(arr, origin='lower', cmap='gray')
        theta = np.linspace(0, 2*np.pi, 200)
        circle_x = cx + radius * np.cos(theta)
        circle_y = cy + radius * np.sin(theta)
        plt.plot(circle_x, circle_y, '--', c='C0')
        plt.plot(cx, cy, '+', c='C0')
        plt.xlim([cx-2.4*radius,cx+2.4*radius])
        plt.ylim([cy-2.4*radius,cy+2.4*radius])
    return cx, cy, radius


def normalize_frame(frame):
    arr = np.asarray(frame, dtype=np.float64)
    mean_val = np.mean(arr)
    return arr / mean_val


def circularize_ellipse(cam, freq:float, Nits:int=1, dV:float=0.2, dphi:float=2.0, 
                        gain:float=1.0, Nframes:int=2, debug:bool=False,
                        amp1=1.0, amp2=1.0, phi2=90):

    # Set input waves to correct frequency and align phase
    wg.set_wave(ch=1,ampl=amp1,offs=0,freq=freq,wave_form='SIN')
    time.sleep(0.4)
    wg.set_wave(ch=2,ampl=amp2,offs=0,freq=freq,wave_form='SIN')    
    time.sleep(0.4)
    wg.set_phase(ch=2,phase=phi2)
    time.sleep(0.4)
    wg.phase_align()

    # Acquire reference frame
    print('Acquiring reference frame ...')
    start_frame = cam.acquire_frames(Nframes)
    cx,cy,radius = best_circle_fit(start_frame)

    ymin = max(int(np.floor(cy - 1.4*radius)), 0)
    ymax = min(int(np.ceil(cy + 1.4*radius)), start_frame.shape[0])
    xmin = max(int(np.floor(cx - 1.4*radius)), 0)
    xmax = min(int(np.ceil(cx + 1.4*radius)), start_frame.shape[1])
    arr_crop = normalize_frame(start_frame[ymin:ymax, xmin:xmax])
    if debug:
        show_frame(arr_crop,title='Cropped frame')
    
    # Reference frame
    cx_local = cx - xmin
    cy_local = cy - ymin
    yy, xx = np.indices(arr_crop.shape)
    rr = np.sqrt((xx - cx_local) ** 2 + (yy - cy_local) ** 2)
    r_idx = np.rint(rr).astype(np.int64)
    radial_sum = np.bincount(r_idx.ravel(), weights=arr_crop.ravel())
    radial_count = np.bincount(r_idx.ravel())
    radial_profile = np.divide(
        radial_sum,
        radial_count,
        out=np.zeros_like(radial_sum),
        where=radial_count > 0,
    )
    ref_frame = radial_profile[r_idx]
    if debug:
        show_frame(ref_frame, title='Circular reference frame')

    IM = np.zeros([arr_crop.size,3])

    # Change wave 1 amplitude
    print('Computing interation matrix: ch1 amplitude ...')
    wg.set_wave(ch=1,ampl=amp1+dV,offs=0,freq=freq,wave_form='SIN')
    push_frame = normalize_frame(cam.acquire_frames(Nframes)[ymin:ymax, xmin:xmax])
    wg.set_wave(ch=1,ampl=amp1-dV,offs=0,freq=freq,wave_form='SIN')
    pull_frame = normalize_frame(cam.acquire_frames(Nframes)[ymin:ymax, xmin:xmax])
    wg.set_wave(ch=1,ampl=amp1,offs=0,freq=freq,wave_form='SIN')
    dframe = push_frame - pull_frame
    dframe -= np.mean(dframe)
    IM[:,0] = dframe.flatten()/(2*dV)
    if debug:
        show_frame(dframe/(2*dV), title='Wave 1 Amplitude Change')

    # Change wave 2 amplitude    
    print('Computing interation matrix: ch2 amplitude ...')
    wg.set_wave(ch=2,ampl=amp2+dV,offs=0,freq=freq,wave_form='SIN')
    push_frame = normalize_frame(cam.acquire_frames(Nframes)[ymin:ymax, xmin:xmax])
    wg.set_wave(ch=2,ampl=amp2-dV,offs=0,freq=freq,wave_form='SIN')
    pull_frame = normalize_frame(cam.acquire_frames(Nframes)[ymin:ymax, xmin:xmax])
    wg.set_wave(ch=2,ampl=amp2,offs=0,freq=freq,wave_form='SIN')
    dframe = push_frame - pull_frame
    dframe -= np.mean(dframe)
    IM[:,1] = dframe.flatten()/(2*dV)
    if debug:
        show_frame(dframe/(2*dV), title='Wave 2 Amplitude Change')

    # Change wave 2 phase
    print('Computing interation matrix: ch2 phase ...')
    wg.set_phase(ch=2,phase=phi2+dphi)
    wg.phase_align()
    push_frame = normalize_frame(cam.acquire_frames(Nframes)[ymin:ymax, xmin:xmax])
    wg.set_phase(ch=2,phase=phi2-dphi)
    wg.phase_align()
    pull_frame = normalize_frame(cam.acquire_frames(Nframes)[ymin:ymax, xmin:xmax])
    wg.set_phase(ch=2,phase=phi2)
    wg.phase_align()
    dframe = push_frame - pull_frame
    dframe -= np.mean(dframe)
    IM[:,2] = dframe.flatten()/(2*dphi)
    if debug:
        show_frame(dframe/(2*dphi), title='Wave 2 Phase Change')

    Rec = np.linalg.pinv(IM)
    end_crop = arr_crop.copy()
    start_err_rms = np.sqrt(np.mean((arr_crop.flatten() - ref_frame.flatten())**2))
    print(f'Reconstructing (start error RMS: {start_err_rms:.4g})...')
    for i in range(Nits):
        arr_crop = end_crop.copy()
        error = ref_frame.flatten() - arr_crop.flatten()
        # error -= np.mean(error)
        if debug:
            show_frame(error.reshape(arr_crop.shape), title='Detected difference from circle')
        dcmd = Rec @ error
        dcmd[:2] = np.clip(dcmd[:2], -dV, dV)
        dcmd[2] = np.clip(dcmd[2], -dphi, dphi)
        dcmd[:2] /= 2 # split correction between two channels
        dcmd *= gain
        new_circle = arr_crop.flatten() + IM @ dcmd
        if debug:
            show_frame(new_circle.reshape(arr_crop.shape), title='Predicted circle after correction')

        amp1 += dcmd[0]
        amp2 += dcmd[1]
        phi2 += dcmd[2]
        wg.set_wave(ch=1,ampl=amp1,offs=0,freq=freq,wave_form='SIN')
        wg.set_wave(ch=2,ampl=amp2,offs=0,freq=freq,wave_form='SIN')
        wg.set_phase(ch=2,phase=phi2)
        wg.phase_align()

        end_frame = cam.acquire_frames(Nframes)
        end_crop = normalize_frame(end_frame[ymin:ymax, xmin:xmax])

        if debug:
            show_frame(new_circle.reshape(arr_crop.shape)-end_crop, title='Predicted - obtained circle after correction')
        end_err_rms = np.sqrt(np.mean((end_crop.flatten() - ref_frame.flatten())**2))
        print(f'\tIteration {i+1:1.0f}: error RMS = {end_err_rms:.4g}, ch1+={dcmd[0]*1e+3:1.0f}mV, ch2+={dcmd[1]*1e+3:1.0f}mV, {dcmd[2]:1.3f}°')
        print(f'\t\tParameters: ch1={amp1*1e+3:1.0f}mV, ch2={amp2*1e+3:1.0f}mV, phi2={phi2:1.3f}°\n')

    if debug:
        show_frame(end_crop, title='Final frame after correction')
        show_frame(end_crop - arr_crop, title='Difference between final and reference frame')
    

    return start_frame, end_frame, amp1, amp2, phi2