
from splattsw.devices.deformable_mirror import SPLATTEngine
from splattsw.devices.wavegenerators import WaveGenerator as wg

from threading import Thread
import time
import numpy as np
import queue


if __name__ == "__main__":

    dm = SPLATTEngine()
    wavegen = wg()

    q_pos = queue.Queue()
    q_cur = queue.Queue()
    q_tn = queue.Queue()

    # Prepare sweep parameters
    wavegen.set_wave(ch=1,ampl=2,offs=0,freq=30,wave_form='SIN')
    time.sleep(1) # wait for steady state
    wavegen.set_sweep(1,30,100,10,amp=2)

    def trigger():
        print('Now starting trigger ...')
        wavegen.trigger_sweep(1)
        print('Trigger sent!')
        time.sleep(2)

    def buffer(Nsamples, dec):
        print('Sending buffer ...')
        mean_pos,mean_cur,tn = dm.read_buffers(external=True, n_samples=Nsamples, decimation=dec)
        print('Buffer completed!')
        q_pos.put(mean_pos)
        q_cur.put(mean_cur)
        q_tn.put(q_tn)
        
    
    wg_thread = Thread(target=trigger)
    buf_thread = Thread(target=buffer, args=(300,0))

    wg_thread.start()
    buf_thread.start()

    wg_thread.join()
    buf_thread.join()

    tn = q_tn.get()
    print(tn)

    pos = q_pos.get()
    print(pos)

    cur = q_cur.get()
    print(cur)

