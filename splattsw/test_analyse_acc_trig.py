import numpy as np
import matplotlib.pyplot as plt
import splattsw.acceleration_analysis as anal
import pandas as pd


def analyse_dataframe(dataframe):

    tn_list = dataframe['wdf_tns']

    peak = []
    peakf = []
    trig = []
    trigf = []

    for tn in tn_list:
        spe_sig1, spe_sig2, f_sig, spe_trig1, spe_trig2, f_trig = get_acc_and_trig_spectra(tn)

        spe_sig = np.sqrt(spe_sig1**2+spe_sig2**2)/2
        spe_trig = np.sqrt(spe_trig1**2+spe_trig2**2)/2

        peak.append(np.max(spe_trig1))
        peakf.append(f_sig[np.argmax(spe_sig)])
        trig.append(np.max(spe_trig))
        trigf.append(f_trig[np.argmax(spe_trig)])

    dataframe['max_sig_acc'] = np.array(peak)
    dataframe['max_sig_freq'] = np.array(peakf)
    dataframe['max_trig_acc'] = np.array(trig)
    dataframe['max_trig_freq'] = np.array(trigf)

    return dataframe




def get_acc_and_trig_spectra(wdfile, f_thr = 200):

    N = 6
    max6ids = np.zeros(N)
    min6ids = np.zeros(N)

    data = anal.openfile(wdfile)

    acc = data[1]
    ids = np.argsort(acc)

    # Find the N largest and smallest values
    max_ctr = -1
    min_ctr = 0

    max6ids[0] = ids[max_ctr]
    min6ids[0] = ids[min_ctr]

    thr = 500
    for j in np.arange(1,N):

        max_ctr -= 1
        min_ctr += 1

        while np.min(np.abs(max6ids[:j]-ids[max_ctr])) <= thr:
            max_ctr -= 1

        while np.min(np.abs(min6ids[:j]-ids[min_ctr])) <= thr:
            min_ctr += 1

        max6ids[j] = ids[max_ctr]
        min6ids[j] = ids[min_ctr]

    # Find signal start
    ids = (np.minimum(np.sort(min6ids),np.sort(max6ids))).astype(int)

    sig_len = np.arange(0,1200,dtype=int)
    sig1 = acc[ids[2]+sig_len]
    sig2 = acc[ids[3]+sig_len]

    spe1,f_sig = anal.acc_spectrum(sig1)
    spe2,f_sig = anal.acc_spectrum(sig2)

    # Find accelerometer signal
    acc_len = np.arange(0,500,dtype=int)
    acc11 = acc[ids[0]+acc_len]
    acc12 = acc[ids[1]+acc_len]
    acc21 = acc[ids[4]+acc_len]
    acc22 = acc[ids[5]+acc_len]

    spacc11,f_acc = anal.acc_spectrum(acc11)
    spacc12,f_acc = anal.acc_spectrum(acc12)
    spacc21,f_acc = anal.acc_spectrum(acc21)
    spacc22,f_acc = anal.acc_spectrum(acc22)

    # Perform mean in energy
    spacc1 = np.sqrt(spacc11**2 + spacc21**2)/2
    spacc2 = np.sqrt(spacc12**2 + spacc22**2)/2

    # Remove high frequencies
    spe_sig1 = spe1[f_sig<f_thr]
    spe_sig2 = spe2[f_sig<f_thr]
    spe_trig1 = spacc1[f_acc<f_thr]
    spe_trig2 = spacc2[f_acc<f_thr]
    f_sig = f_sig[f_sig<f_thr]
    f_trig = f_acc[f_acc<f_thr]

    # ref_spe,f_ref = anal.acc_spectrum(acc[0:1000])
    # ref_spe = ref_spe[f_ref<f_thr]
    # f_ref = f_ref[f_ref < f_thr]
    #
    # plt.figure()
    # plt.scatter(f_ref,ref_spe)
    # plt.scatter(f_sig,spe_sig1)
    # plt.grid('on')
    # plt.show()

    return spe_sig1, spe_sig2, f_sig, spe_trig1, spe_trig2, f_trig

# Gap1, Vac2
wdflist = ['OBB-Vibration_2025-02-18T17-03-39-102.wdd',
 'OBB-Vibration_2025-02-18T17-04-35-988.wdd',
 'OBB-Vibration_2025-02-18T17-05-32-785.wdd',
 'OBB-Vibration_2025-02-18T17-06-29-635.wdd',
 'OBB-Vibration_2025-02-18T17-07-26-502.wdd',
 'OBB-Vibration_2025-02-18T17-08-23-331.wdd',

# Gap2, Vac2
'OBB-Vibration_2025-02-18T17-09-26-405.wdd',
'OBB-Vibration_2025-02-18T17-10-23-179.wdd',
'OBB-Vibration_2025-02-18T17-11-19-966.wdd',
'OBB-Vibration_2025-02-18T17-12-16-822.wdd',
'OBB-Vibration_2025-02-18T17-13-13-576.wdd',
'OBB-Vibration_2025-02-18T17-14-10-387.wdd',

# Gap3, Vac2
'OBB-Vibration_2025-02-18T17-15-13-463.wdd',
'OBB-Vibration_2025-02-18T17-16-10-278.wdd',
'OBB-Vibration_2025-02-18T17-17-07-141.wdd',
'OBB-Vibration_2025-02-18T17-18-03-958.wdd',
'OBB-Vibration_2025-02-18T17-19-00-840.wdd',
'OBB-Vibration_2025-02-18T17-19-57-762.wdd',

# Gap1, Vac1
'OBB-Vibration_2025-02-19T16-44-25-595.wdd',
'OBB-Vibration_2025-02-19T16-45-22-171.wdd',
'OBB-Vibration_2025-02-19T16-46-17-195.wdd',
'OBB-Vibration_2025-02-19T16-47-12-137.wdd',
'OBB-Vibration_2025-02-19T16-48-07-156.wdd',
'OBB-Vibration_2025-02-19T16-49-02-115.wdd',

# Gap2, Vac1
'OBB-Vibration_2025-02-19T16-50-03-197.wdd',
'OBB-Vibration_2025-02-19T16-50-58-170.wdd',
'OBB-Vibration_2025-02-19T16-51-53-062.wdd',
'OBB-Vibration_2025-02-19T16-52-47-997.wdd',
'OBB-Vibration_2025-02-19T16-53-43-069.wdd',
'OBB-Vibration_2025-02-19T16-54-38-121.wdd',

# Gap3, Vac1
'OBB-Vibration_2025-02-19T16-55-39-260.wdd',
'OBB-Vibration_2025-02-19T16-56-34-318.wdd',
'OBB-Vibration_2025-02-19T16-57-29-240.wdd',
'OBB-Vibration_2025-02-19T16-58-24-146.wdd',
'OBB-Vibration_2025-02-19T16-59-19-186.wdd',
'OBB-Vibration_2025-02-19T17-00-14-279.wdd',

# Gap1, Vac0
'OBB-Vibration_2025-02-19T17-12-30-930.wdd',
'OBB-Vibration_2025-02-19T17-13-28-856.wdd',
'OBB-Vibration_2025-02-19T17-14-23-680.wdd',
'OBB-Vibration_2025-02-19T17-15-18-422.wdd',
'OBB-Vibration_2025-02-19T17-16-13-135.wdd',
'OBB-Vibration_2025-02-19T17-17-08-012.wdd',

# Gap2, Vac0
'OBB-Vibration_2025-02-19T17-18-09-177.wdd',
'OBB-Vibration_2025-02-19T17-19-04-012.wdd',
'OBB-Vibration_2025-02-19T17-19-58-780.wdd',
'OBB-Vibration_2025-02-19T17-20-53-573.wdd',
'OBB-Vibration_2025-02-19T17-21-48-374.wdd',
'OBB-Vibration_2025-02-19T17-22-43-204.wdd',

# Gap3, Vac0
'OBB-Vibration_2025-02-19T17-23-44-112.wdd',
'OBB-Vibration_2025-02-19T17-24-38-761.wdd',
'OBB-Vibration_2025-02-19T17-25-33-692.wdd',
'OBB-Vibration_2025-02-19T17-26-28-370.wdd',
'OBB-Vibration_2025-02-19T17-27-23-218.wdd',
'OBB-Vibration_2025-02-19T17-28-17-934.wdd']

gap_vec = np.tile(np.repeat([1,2,3],6),3)
vac_vec = np.repeat([2,1,0],18)
mode_vec = np.tile([1,2,3],18)

gap_list = []
gap_list += ("Gap%s" %g for g in gap_vec)

vac_list = []
vac_list += ("Vac%s" %v for v in vac_vec)

mode_list = []
mode_list += ("Mode%s" %m for m in mode_vec)


dataset = {
    'gaps' : gap_list,
    'vacuum' : vac_list,
    'modes' : mode_list,
    'wdf_tns' : wdflist
}

pdset = pd.DataFrame(dataset)