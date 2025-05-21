import opticalib
pyconf = '/mnt/libero/SPLATTData/Data/SysConfigurations/configuration.yaml'
opticalib.load_configuration_file(pyconf)

from splattsw import splatt_analysis as sp


tnrip  = '20250424_155319'
tnlist= ['20250424_162052','20250424_163213','20250424_163920','20250424_164635']
ntn = len(tnlist)
accrip = 'OBB-Vibration_2025-04-24T15-52-38-996.wdd'
acclist = ['OBB-Vibration_2025-04-24T16-20-12-361.wdd','OBB-Vibration_2025-04-24T16-31-32-937.wdd','OBB-Vibration_2025-04-24T16-38-40-115.wdd','OBB-Vibration_2025-04-24T16-45-54-526.wdd']

lab = ['RIP','P=1b','P=0.8b','P=0.7b','P=0.55b']
ttr = sp.tiltvec(tnrip)
sper, f = sp.tt_spectrum(tnrip, tt=ttr)

tt = []
spec = []
for i in tnlist:
    tv = sp.tiltvec(i)
    tt.append(tv)
    spe, f = sp.tt_spectrum(i, tt=tv)
    spec.append(spe)
    sp.plot_ttspectra(spe,f)

figure()
plot(f, sper[1,:])
for i in range(ntn):
    spv = spec[i]
    plot(f, spv[1,:])
    
xlim(30,120); yscale('log')
ylabel('Amp Spectrum [nm]')
xlabel('Freq [Hz]')
legend(lab)

sp.plot_multispec(spec,f, tnlist)
plot(f, sper[1,:],'-.')
xlim(40,110)


## accelerometers
wrip,fa0 = sp.mech_spectrum(accrip)
plot(fa0, wrip[0])
xlim(40,110)

wset = []
wall = []
ff = []
for i in acclist:
    ww, fa = sp.mech_spectrum(i)
    ff.append(fa)
    wset.append(ww[0])
    wall.append(ww)

figure()
plot(fa0, wrip[0])
for i in range(ntn):
    plot(ff[i], wset[i])

yscale('log')
legend(lab)
xlim(20,120)
xlabel('Freq [Hz]')

figure()
plot(fa0, wrip[0])
for i in range(2):
    plot(ff[i], wset[i])

yscale('log')
legend(lab)
xlim(20,120)
xlabel('Freq [Hz]')


#plot di una sola serie per verifica
figure()
wpl=wall[0]
for i in range(4):
    plot(ff[0], wpl[i,:])
yscale('log')
legend(['RBf','RBb','Stand','Piezo'])
xlim(20,120)
xlabel('Freq [Hz]')
ylabel('Oscill amp, by Acc integration [m]')

plot(fa0, wrip[0])
xlim(40,110)


gg = sp.openwdfile(accrip)




#### debug sweep analysis

tnconf= '20250428_100000.meas'

tninfo,measinfo = sp.read_analysisconf(tnconf)


#re-saving acc files
import acceleration_analysis as ac

tnl  = ['20250424_155319','20250424_162052','20250424_163213','20250424_163920','20250424_164635']
ntn = len(tnl)
accl = ['OBB-Vibration_2025-04-24T15-52-38-996.wdd','OBB-Vibration_2025-04-24T16-20-12-361.wdd','OBB-Vibration_2025-04-24T16-31-32-937.wdd','OBB-Vibration_2025-04-24T16-38-40-115.wdd','OBB-Vibration_2025-04-24T16-45-54-526.wdd']
for i in range(ntn):
    print(accl[i]+'  --> '+tnl[i])
    ac.savefile(accl[i],tnl[i])


freq = [3,5,7,10,15,20,25,30]
tnlist = ['20250514_111143', '20250514_111216', '20250514_111249', '20250514_111322','20250514_111355', '20250514_111428', '20250514_111501', '20250514_111532']
ntn = len(tnlist)
modesid = [1,2,4,5];nmodes = len(modesid)
pvec = np.zeros([ntn, nmodes])

for i in range(ntn):
    pp = singlefreq(tnlist[i], freq[i], peakw=1)
    pvec[i,:]= pp


    zv, sv = sp.zvec(tnlist[i])
    freq4d = sp.read_4dfreq(tnlist[i])
    spe, f = th.spectrum(zv,dt=1/freq4d)
    for m in range(nmodes):
        peakt, pf = sp.find_peak(spe[modesid[m],:],f,bound = [freq[i]-1,freq[i]+1])
        pvec[i,m]= peakt
        print(pf)


def dataprocess(tn, meastype = 'sweep', freq=None, nbins=1):
    if meastype == 'single' and freq == None:
        print('Single freq measurement, frequency is requested')
        return 
    zv, sv = sp.zvec(tn)
    freq4d = sp.read_4dfreq(tn)#implementare la lettura del file da TestConfig
    spe, f = th.spectrum(zv,dt=1/freq4d)
    spes,f = th.spectrum(sv, dt=1/freq4d)
    if meastype == 'sweep':
        frunn    = runningMean(f, nbins, dim=0)
        specz    = runningMean(spe, nbins)
        specs    = runningMean(spes, nbins,dim=0)
        return frunn, specz, specs

    if meastype == 'single':
        peaks, pf = sp.find_peak(spes,f,bound = [freq-nbins,freq+nbins])
        pvec = []
        for m in range(len(zv)):
            peakz, pf = sp.find_peak(spe[m,:],f,bound = [freq-nbins,freq+nbins])
            pvec.append(peakz)
            print(pf)
        return freq, pvec, peaks


def singlefreq(tn, freq, peakw=1):
    zv, sv = sp.zvec(tn)
    freq4d = sp.read_4dfreq(tn)
    spe, f = th.spectrum(zv,dt=1/freq4d)
    pvec = []
    for m in range(nmodes):
        peakt, pf = sp.find_peak(spe[modesid[m],:],f,bound = [freq-peakw,freq+peakw])
        pvec.append(peakt)
        print(pf)
    return pvec

