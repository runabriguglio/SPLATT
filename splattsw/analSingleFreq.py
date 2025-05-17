import pandas as pd
import seaborn as sns
import numpy as np

fvec = np.array([3,5,7,10,15,20,25,30])
pvec = np.array([550,550,550,550,400,400,250,250,102,102,-2,-2,-2,250,550])*-1e-3
Kpvec = np.array([500,400,300,200,500,200,500,200,500,200,500,200,0,0,0])

L = len(fvec)
N = len(Kpvec)

fvec = np.tile(fvec,N)
pvec = np.repeat(pvec,L)
Kpvec = np.repeat(Kpvec,L)

tn_list = [
['20250514_111143', '20250514_111216', '20250514_111249', '20250514_111322',
'20250514_111355', '20250514_111428', '20250514_111501', '20250514_111532'],
['20250514_111908', '20250514_111940', '20250514_112013', '20250514_112046',
'20250514_112119', '20250514_112152', '20250514_112225', '20250514_112256'],
['20250514_112528', '20250514_112601', '20250514_112634', '20250514_112706',
'20250514_112739', '20250514_112812', '20250514_112844', '20250514_112916'],
['20250514_113329', '20250514_113401', '20250514_113434', '20250514_113507',
'20250514_113540', '20250514_113612', '20250514_113645', '20250514_113716'],

['20250514_114307', '20250514_114339', '20250514_114412', '20250514_114444',
'20250514_114517', '20250514_114549', '20250514_114622', '20250514_114653'],
['20250514_115029', '20250514_115101', '20250514_115134', '20250514_115207',
'20250514_115239', '20250514_115312', '20250514_115344', '20250514_115415'],

['20250514_115901', '20250514_115934', '20250514_120006', '20250514_120039',
'20250514_120111', '20250514_120143', '20250514_120216', '20250514_120247'],
['20250514_120832', '20250514_120904', '20250514_120937', '20250514_121009',
'20250514_121041', '20250514_121113', '20250514_121145', '20250514_121216'],

['20250514_121725', '20250514_121757', '20250514_121829', '20250514_121902',
'20250514_121934', '20250514_122006', '20250514_122038', '20250514_122109'],
['20250514_122441', '20250514_122513', '20250514_122546', '20250514_122618',
'20250514_122650', '20250514_122722', '20250514_122753', '20250514_122824'],

['20250514_123340', '20250514_123412', '20250514_123444', '20250514_123517',
'20250514_123549', '20250514_123621', '20250514_123653', '20250514_123724'],
['20250514_124112', '20250514_124145', '20250514_124217', '20250514_124249',
'20250514_124321', '20250514_124353', '20250514_124425', '20250514_124456'],

['20250514_130017', '20250514_130052', '20250514_130126', '20250514_130200',
'20250514_130234', '20250514_130307', '20250514_130341', '20250514_130414'],
['20250514_131205', '20250514_131240', '20250514_131314', '20250514_131348',
'20250514_131422', '20250514_131456', '20250514_131530', '20250514_131602'],
['20250514_135150', '20250514_135225', '20250514_135259', '20250514_135333',
'20250514_135407', '20250514_135440', '20250514_135515', '20250514_135547']]


import acceleration_analysis as sp
import splatt_utilities as buf

max_displ = np.zeros([N*L,4])
max_spe = np.zeros([N*L,4])
phi_max_spe = np.zeros([N*L,4])

freq_WebDAQ =  1651.6129 #Hz; minimum sampling frequency
g0 = 9.807

def analyse_tn(tn, ex_freq, bins: int = 2):
    max_spe = np.zeros(4)
    phi_max_spe = np.zeros(4)
    max_displ = np.zeros(4)

    data = sp.openfile(tn)
    spe, fvec, phis = sp.get_spectrum(data, 1/freq_WebDAQ, phase = True)

    om = 2*np.pi*ex_freq

    f_id = np.argmin(np.abs(fvec-ex_freq))
    print(spe[:,f_id])
    min_id = int( np.max((0,f_id-bins)))
    max_id = int (np.min((f_id+bins+1,len(fvec))))
    max_spe = np.sum(spe[:,min_id:max_id],axis=1)*g0/om**2
    ph = phis[:,f_id]
    ph += np.pi*(ph<0)
    phi_max_spe = ph*180/np.pi

    bufdata, _, fvec = buf.analyse_buffer(tn)

    zspe = bufdata['pos_zernike_spectrum']
    f_id = np.argmin(np.abs(fvec-ex_freq))
    max_displ[:3] = zspe[f_id,:3]

    max_id = np.argmax(zspe[f_id,3:])
    max_displ[-1] = zspe[f_id,max_id]

    return max_spe, phi_max_spe, max_displ

if __name__ == "__main__":
    # Cycle on all tns
    for ii in range(N*L):
        spe, phis, displ = analyse_tn(tn_list[int((ii-(ii%L))/L)][ii%L], fvec[ii])
        max_spe[ii] = spe
        phi_max_spe[ii] = phis
        max_displ[ii] = displ

    # Fill dataframe
    ph3 = phi_max_spe[:,3]
    dictframe = {'Pressure [bar]': pvec, 'Kp gain': Kpvec, 'Frequency [Hz]': fvec,
                'OBB front displacement [m]': max_spe[:,0], 'OBB rear displacement [m]': max_spe[:,1], 'Stand displacement [m]': max_spe[:,2], 'Stand arm displacement [m]': max_spe[:,3],
                'OBB front phase wrt stand arm [deg]': phi_max_spe[:,0]-ph3, 'OBB rear phase wrt stand arm [deg]': phi_max_spe[:,1]-ph3, 'Stand phase wrt stand arm [deg]': phi_max_spe[:,2]-ph3, 'Stand arm phase [deg]': ph3,
                'Piston amplitude [m]': max_displ[:,0], 'Tip amplitude [m]': max_displ[:,1], 'Tilt amplitude [m]': max_displ[:,2], 'Other modes amplitude [m]': max_displ[:,3]}

    dframe = pd.DataFrame(dictframe)

    # Save dataframe
    filename = 'SingleFreqAccDataFrame.csv'
    dframe.to_csv(filename,index=False)
