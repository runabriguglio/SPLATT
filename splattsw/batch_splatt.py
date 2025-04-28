import aoptics
pyconf = '/mnt/libero/SPLATTData/Data/SysConfigurations/configuration.yaml'
aoptics.load_configuration_file(pyconf)

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
