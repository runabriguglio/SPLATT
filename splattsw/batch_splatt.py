from SPLATT.splattsw import splatt_analysis as sp
from astropy.io import fits as pyfits
#from M4.m4.mOTT_analysis import timehistory as th
from M4.m4.mini_OTT import timehistory as th

from importlib  import reload #for reload
wdlist0 = ['SPLATT_Test_2021-12-23T12-34-16-379.wdd','SPLATT_Test_2021-12-23T12-38-01-620.wdd','SPLATT_Test_2021-12-23T12-40-30-542.wdd','SPLATT_Test_2021-12-23T12-43-56-167.wdd','SPLATT_Test_2021-12-23T12-45-46-445.wdd','SPLATT_Test_2021-12-23T12-47-48-308.wdd','SPLATT_Test_2021-12-23T12-49-47-559.wdd']  #redpitaya
wdlist1 = ['SPLATT_Test_2021-12-23T13-01-45-656.wdd','SPLATT_Test_2021-12-23T13-03-43-910.wdd','SPLATT_Test_2021-12-23T13-05-53-760.wdd','SPLATT_Test_2021-12-23T13-07-59-066.wdd','SPLATT_Test_2021-12-23T13-09-49-420.wdd','SPLATT_Test_2021-12-23T13-11-57-179.wdd','SPLATT_Test_2021-12-23T13-13-56-266.wdd'] #redpitaya
wdlist2 = ['SPLATT_AIV_2021-12-22T09-48-29-847.wdd','SPLATT_AIV_2021-12-22T09-49-29-631.wdd','SPLATT_AIV_2021-12-22T09-49-49-699.wdd','SPLATT_AIV_2021-12-22T09-50-29-202.wdd','SPLATT_AIV_2021-12-22T09-51-10-432.wdd','SPLATT_AIV_2021-12-22T09-51-30-096.wdd','SPLATT_AIV_2021-12-22T09-51-48-677.wdd'] #tektronix
a0 = []
b0 = []
a1 = []
b1 = []
a2 = []
b2 = []
p0 = '-x'
p1 = '-y'
p2 = '-o'
ct1500 = []
ct1501 = []
f150 = 454
nn = size(wdlist0)
for i in range(nn):
    wd = sp.openfile(wdlist0[i])
    w0 = wd[:,0:4999]
    w1 = wd[:,5000:]
    spe,f=sp.acc_spectrum(w0)
    spe1,f1=sp.acc_spectrum(w1)
    a0.append(np.mean(spe[0,0:f150+10]))
    b0.append(np.mean(spe1[0,0:f150+10]))
    
    wd = sp.openfile(wdlist1[i])
    w0 = wd[:,0:4999]
    w1 = wd[:,5000:]
    spe,f=sp.acc_spectrum(w0)
    spe1,f1=sp.acc_spectrum(w1)
    a1.append(np.mean(spe[0,0:f150+10]))
    b1.append(np.mean(spe1[0,0:f150+10]))

    wd = sp.openfile(wdlist2[i])
    w0 = wd[:,0:4999]
    w1 = wd[:,5000:]
    spe,f=sp.acc_spectrum(w0)
    spe1,f1=sp.acc_spectrum(w1)
    a2.append(np.mean(spe[0,0:f150+10]))
    b2.append(np.mean(spe1[0,0:f150+10]))

plot(a0,p0)
plot(a1, p1)
plot(a2, p2)
plot(b0,p0)
plot(b1,p1)
plot(b2,p2)

wdlist = wdlist0
for i in wdlist:
    wd = sp.openfile(wdlist[i])
    plot(wd[0,:])
    xlim(0,250)
    i = i+1
   
figure(5)
wdlist = wdlist0
i=0
wd = sp.openfile(wdlist[i])
spe,f=sp.acc_spectrum(wd)
plot(f, spe[0,:])
xlim(0,200)
yscale('log')



tn0='20211231_125510'
tn1='20211229_153550'

tt0 = sp.tiltvec(tn0)
tt1 = sp.tiltvec(tn1)
spe, f=sp.tt_spectrum(tt, tn)



rp.pulse_train(1,0.1,0,50,0.3,12)
tn0='20211231_142416'#set
tn1='20211231_142859' #rip
tn11 = '20220103_154858' #RIP, repetition
tn111= '20220103_160017'#RIP, repetition

tt0 = sp.tiltvec(tn0)
tt1 = sp.tiltvec(tn1)
tt11 = sp.tiltvec(tn11)
tt111 = sp.tiltvec(tn111)

spe0, f0=sp.tt_spectrum(tt0, tn0)
spe1, f1=sp.tt_spectrum(tt1, tn1)
spe11, f11=sp.tt_spectrum(tt11, tn11)
spe111, f111=sp.tt_spectrum(tt111, tn111)

clf()
plot(f0, spe0[1,:])
plot(f1, spe1[1,:])
plot(f11, spe11[1,:])
plot(f111, spe111[1,:])


yscale('log')


from M4.m4.mOTT_analysis import timehistory as th
npo=10
spe00=th.runningMean(spe0[1,:],npo)
spe11=th.runningMean(spe1[1,:],npo)
f00 = th.runningMean(f0, npo)
f11 = th.runningMean(f1, npo)
clf()
plot(f00, spe00)
plot(f11, spe11)
yscale('log')

clf()
rr= spe00/spe11
#rr=spe0[1,:]/spe1[1,:]
plot(f00, rr)
yscale('log')
title('Ratio: TTset/TTrip')
nn=size(spe00)
plot(f00, np.ones(nn)

#comparison with single freq data
fs = np.array([0.5,1,5,10,15,20,25,30,40,50,60,72,80,88,100,110,120])
g=np.array([1.08226923321263,1.07776156553086,1.12374546318523,1.22240756373459,1.24354131147671,1.12448759738125,0.985361477418602,0.896473933063444,0.600898529159274,0.614510197258442,0.775761696341308,0.341217885375494,0.257842441447835,0.722265201918614,0.933762500506093,0.460261206274201,0.438239082074698])
figure(2)
plot(fs,g,'X')
yscale('log')



def rms(v):
    v1=np.sum(v**2)
    v1 = np.sqrt(v1)
    return v1




rp.pulse_train(1,0.05,0,5,0.1,3)
wd0='SPLATT_Test_2022-01-03T14-18-32-607.wdd'

rp.pulse_train(1,0.05,0,5,0.1,10)
wdfile0='SPLATT_Test_2022-01-03T14-26-28-062.wdd'

rp.pulse_train(1,0.05,0,20,0.1,10)
wdfile1='SPLATT_Test_2022-01-03T14-27-25-914.wdd'

rp.pulse_train(1,0.05,0,1,0.1,10)
wdfile2='SPLATT_Test_2022-01-03T14-33-57-477.wdd'

rp.pulse_train(1,0.05,0,0.5,0.05,10)
wdfile3='SPLATT_Test_2022-01-03T14-36-42-438.wdd'

rp.pulse_train(1,0.05,0,1,0.5,10)
wdfile4='SPLATT_Test_2022-01-03T14-56-15-938.wdd'

sp0,f0=sp.acc_spectrum(wdfile0)
sp1,f1=sp.acc_spectrum(wdfile1)
sp2,f2=sp.acc_spectrum(wdfile2)
sp3,f3=sp.acc_spectrum(wdfile3)
sp4,f4=sp.acc_spectrum(wdfile4)


#repeatability of pulses with tektronix
w0='SPLATT_Test_2022-01-03T16-36-02-096.wdd'
w1='SPLATT_Test_2022-01-03T16-36-10-961.wdd'
w2='SPLATT_Test_2022-01-03T16-36-18-620.wdd'
s0,f0=sp.acc_spectrum(w0)
s1,f1=sp.acc_spectrum(w1)
s2,f2=sp.acc_spectrum(w2)
clf()
plot(f0, s0[0,:],'.')
xlim(0,100)
yscale('log')
plot(f1,s1[0,:],'.')
plot(f2,s2[0,:],'.')



#optical meas pulses with tektronix
#RIP
r0='20220103_164252'
r1='20220103_164457'
t0='20220103_170408'
t1='20220103_170611'

ttr0 = sp.tiltvec(r0)
ttr1 = sp.tiltvec(r1)
ttt0 = sp.tiltvec(t0)
ttt1 = sp.tiltvec(t1)

spr0, f0=sp.tt_spectrum(ttr0, r0)
spr1, f0=sp.tt_spectrum(ttr1, r1)
spt0, f0=sp.tt_spectrum(ttt0, t0)
spt1, f0=sp.tt_spectrum(ttt1, t1)

clf()
plot(f0, spr0[1,:])
plot(f0, spr1[1,:])
plot(f0, spt0[1,:])
plot(f0, spt1[1,:])
yscale('log')
r=spr0[1,:]+spr1[1,:]
t=spt0[1,:]+spt1[1,:]
clf()
npp=10
tr=th.runningMean(t/r, npp)
plot(th.runningMean(f0,npp), tr,'.')
yscale('log')
plot(f0, ones(size(f0)))

#10000frames
t0='20220103_173059' #set
r0='20220103_173903' # rip

ttr0 = sp.tiltvec(r0)
ttt0 = sp.tiltvec(t0)

spr0, f0=sp.tt_spectrum(ttr0, r0)
spt0, f0=sp.tt_spectrum(ttt0, t0)

clf()
v=spt0[1,:]/spr0[1,:]
npp=10
ff=th.runningMean(f0[1:],npp)
vv=th.runningMean(v[1:],npp)
plot(ff, vv, '.')
plot(ff, ones(size(ff)))

yscale('log')
clf()


#comb test, analysis
#shell rip
tnr0='20220104_173528'
tnr1='20220104_174056'
tnr2='20220104_174357'
tnr3='20220104_174509'
tns0='20220104_175105'
tns1='20220104_175503'
tns2='20220104_175622'
tns3='20220104_175736'

ttr0 = sp.tiltvec(tnr0)
spr0, f=sp.tt_spectrum(ttr0, tnr0)
fc, vr0 = sp.comb_analysis(f,spr0[1,:],10, freqbin=4)

ttr1 = sp.tiltvec(tnr1)
spr1, f=sp.tt_spectrum(ttr1, tnr1)
fc, vr1 = sp.comb_analysis(f,spr1[1,:],10, freqbin=4)

ttr2 = sp.tiltvec(tnr2)
spr2, f=sp.tt_spectrum(ttr2, tnr2)
fc, vr2 = sp.comb_analysis(f,spr2[1,:],10, freqbin=4)

ttr3 = sp.tiltvec(tnr3)
spr3, f=sp.tt_spectrum(ttr3, tnr3)
fc, vr3 = sp.comb_analysis(f,spr3[1,:],10, freqbin=4)


tts0 = sp.tiltvec(tns0)
sps0, f=sp.tt_spectrum(tts0, tns0)
fc, vs0 = sp.comb_analysis(f,sps0[1,:],10, freqbin=4)

tts1 = sp.tiltvec(tns1)
sps1, f=sp.tt_spectrum(tts1, tns1)
fc, vs1 = sp.comb_analysis(f,sps1[1,:],10, freqbin=4)

tts2 = sp.tiltvec(tns2)
sps2, f=sp.tt_spectrum(tts0, tns2)
fc, vs2 = sp.comb_analysis(f,sps2[1,:],10, freqbin=4)

tts3 = sp.tiltvec(tns3)
sps3, f=sp.tt_spectrum(tts3, tns3)
fc, vs3 = sp.comb_analysis(f,sps3[1,:],10, freqbin=4)

vr = (vr0+vr1+vr2+vr3)/4
vs = (vs0+vs1+vs2+vs3)/4
clf()
yscale('log')
plot(fc, vr)
plot(fc, vs,'x')

clf()
plot(fc, vs/vr)
yscale('log')
plot(fc, ones(np.size(fc))

fs = np.array([0.5,1,5,10,15,20,25,30,40,50,60,72,80,88,100,110,120])
g=np.array([1.08226923321263,1.07776156553086,1.12374546318523,1.22240756373459,1.24354131147671,1.12448759738125,0.985361477418602,0.896473933063444,0.600898529159274,0.614510197258442,0.775761696341308,0.341217885375494,0.257842441447835,0.722265201918614,0.933762500506093,0.460261206274201,0.438239082074698])
figure(2)
plot(fs,g,'-X')
yscale('log')

#sweep analysis
tn0=['20220112_174748','20220112_175137','20220112_175525']
tn1 = ['20220112_180426','20220112_180813','20220112_181157']
nn=2501
nreb = 20

f, v00 = sp.sweep_analysis(tn0[0])
f, v01 =sp.sweep_analysis(tn0[1])
f, v02 = sp.sweep_analysis(tn0[2])

f, v10 = sp.sweep_analysis(tn1[0])
f, v11 =sp.sweep_analysis(tn1[1])
f, v12 = sp.sweep_analysis(tn1[2])

v0 = (v00+v01+v02)/3
ff = sp.runningMean(f, nreb)
v0 = sp.runningMean(v0, nreb)

v1 = (v10+v11+v12)/3
ff = sp.runningMean(f, nreb)
v1 = sp.runningMean(v1, nreb)

rr = v1/v0
nf = size(ff)
ref = ones(nf)
plot(ff, ref)
yscale('log')
plot(ff, rr)
title('SPLATT attenuation')



'''v=np.full([nn, 3],0)
for i in np.arange(0,3):
    f,vi = sp.sweep_analysis(tn0[i])
    v[:,i]=vi

v0=np.mean(v,axis=1)
ff = sp.runningMean(f, nreb)
v00 = sp.runningMean(v0, nreb)
v=np.full([nn, 3],0)
for i in np.arange(0,3):
    f,vi = sp.sweep_analysis(tn1[i])
    v[:,i]=vi
'''
v1=np.mean(v,axis=1)
v11 = sp.runningMean(v1, nreb)
rr = v11/v00
nf = size(ff)
ref = ones(nf)
plot(ff, ref)
yscale('log')
plot(ff, rr)
title('SPLATT attenuation')



#sweep at different loop
tn0='20220117_104821'
tn=['20220117_114233','20220117_110423','20220117_112421']
kp=[250,500,1250]
f,v0=sp.sweep_analysis(tn0, tiltsum=1)
plot(f,v0,'k'); yscale('log');xlim(65,115); title('OBB Response vs Kp')
nn = size(tn)
nv=size(v0)
lb = []
figure(3)

vk = np.zeros([nn, nv])
for i in np.arange(nn):
    f,v=sp.sweep_analysis(tn[i], tiltsum=1)
    vk[i,:] = v
    lb.append('Kp: %s' %kp[i])
#    plot(f,v, cols[i], label=lb[i])
for i in np.arange(nn):
    plot(f,vk[i,:],cols[i], label=lb[i])
legend()




#sweep at different gap
gap = np.array([250,200,150,100,75,40,20])
tn0='20220117_104821'
tn=['20220117_140144','20220117_112421','20220117_115657','20220117_120640','20220117_121257','20220117_121931','20220117_122856']
wdlist = ['SPLATT_Test_2022-01-17T14-01-44-319.wdd','SPLATT_Test_2022-01-17T11-24-21-547.wdd','SPLATT_Test_2022-01-17T11-56-57-406.wdd','SPLATT_Test_2022-01-17T12-06-40-732.wdd','SPLATT_Test_2022-01-17T12-12-57-272.wdd','SPLATT_Test_2022-01-17T12-19-31-500.wdd','SPLATT_Test_2022-01-17T12-28-56-065.wdd']
norm4d = 1e-9
normacc =  0.002
cols = ['y','b','c','g','r','m','y']
nn = size(tn)
figure(3)

spea, fa=sp.acc_spectrum(wdlist[0])
nacc=size(fa)
accrb = np.zeros([nn, nacc])
accstand = np.zeros([nn,nacc])
for i in np.arange(nn):
    spea, fa=sp.acc_spectrum(wdlist[i])
    accrb[i,:]=spea[0,:]
    accstand[i,:]=spea[1,:]


#plot
f,v0=sp.sweep_analysis(tn0, tiltsum=1)
plot(f,v0,'k'); yscale('log');xlim(65,115); title('OBB Response vs gap')
nn = size(tn)
nv=size(v0)
lb = []
vgap = np.zeros([nn, nv])
for i in np.arange(nn):
    f,v=sp.sweep_analysis(tn[i], tiltsum=1)
    vgap[i,:] = v
    lb.append('Gap: %sum' %gap[i])
#    plot(f,v, cols[i], label=lb[i])
vgapd = np.zeros([nn, nv])

for i in np.arange(nn):
    vgapd[i,:] = vgap[i,:]/v0

for i in np.arange(nn):
    plot(f,vgapd[i,:],cols[i], label=lb[i])

    legend()


fname = '/mnt/jumbo/SPLATT/20220117_gap_TT.fits'
pyfits.writeto(fname,vgap)
fname = '/mnt/jumbo/SPLATT/20220117_gap_acc.fits'
pyfits.writeto(fname,accrb)


for i in np.arange(nn):
    plot(f,vgap[i,:]/norm4d,cols[i], label=lb[i])
legend()

for i in np.arange(nn):
    plot(fa,accrb[i,:]/normacc,cols[i], label=lb[i])
legend()



# reading of dpos and current
p = [20,45,90,125,172,218]
c = [1402,1408,1433,1450,1477,1505]
figure(4)
plot(p,c); xlabel('Delta Position [um]']; ylabel('Tot current [mA]')

#plot curr vs pos
p = [20,45,90,125,172,218]
     ...: c = [1402,1408,1433,1450,1477,1505]
     ...: figure(4)
     ...: plot(p,c,'-x')


#mean shape at resonance
#delay test
delay = np.arange(11)
tn = ['20220118_151622','20220118_151645',	 '20220118_151713',	 '20220118_151751',	 '20220118_151822',	 '20220118_151854',	 '20220118_151924',	 '20220118_151951',	 '20220118_152033',	 '20220118_152056',	 '20220118_152129']

ntn = size(tn)
fl0=sp.fileList(tn[0])
nn = size(fl)
qmm=[]
for i in np.arange(ntn):
    qm0=[]
    fli = sp.fileList(tn[i])
    for j in np.arange(nn):
        qm0.append(th.frame(j, fli))
    m0=np.ma.mean(qm0, axis=0)
    qmm.append(m0)

v = []
v1= []
zz0 = np.array([1,2,3])
zz = np.array([1,2,3,4,5,6])

for i in np.arange(ntn):
    q1=qmm[i]
    q1 = th.removeZernike(q1,zz0)
    v1.append(np.ma.std(q1))
    v.append(np.ma.std(qmm[i]))


#pixel wide FFT and time delay
tn='20220118_151622'
fl = sp.fileList(tn)
nn = size(fl)
mm=[]
for i in np.arange(nn):
    v= th.frame(i,fl)
    mm.append(v-np.mean(v))

mma = np.asarray(mm)
 
spe  = np.fft.rfft(mm, axis=0, norm='ortho')
npo   = np.sqrt(spe.shape[0])   #modRB 
spe  = (np.abs(spe)) / npo
spe = spe**2
freq = np.fft.rfftfreq((shape(mm))[0], d=dt)




# comparison of mean shapes at diff freq. 200 um, kp=1250
tn0='20211229_175305' #72 Hz
tn1='20211229_180014' #100 Hz
tn='20211229_154247' #shell rip
tn=tn0
fl=sp.fileList(tn)
mm = []
nn = size(fl)
for i in np.arange(nn):
    mm.append(th.frame(i,fl))
    
mm = np.ma.mean(mm, axis=0)
mm = th.removeZernike(mm, np.array([1,2,3]))

# pixel-wide spectrum
tn='20220117_112421' #shell set 200 um kp=1250
tn0='20220117_104821' #shell RIP
freq = 250
dt = 1/freq
fl=sp.fileList(tn)
mm = []
nn = size(fl)
for i in np.arange(nn):
    v= th.frame(i,fl)
    mm.append(v-np.mean(v))

mm = np.asarray(mm)

spe  = np.fft.rfft(mm, axis=0, norm='ortho')
npo   = np.sqrt(spe.shape[0])   #modRB 
spe  = (np.abs(spe)) / npo
spe = spe**2
freq = np.fft.rfftfreq((shape(mm))[0], d=dt)


f2check = [600, 950]
imshow(spe[f2check[0],:,:]); title('Freq: %sHz' %freq[f2check[0]]);


spe  = np.fft.rfft(mm, axis=0, norm='ortho')

freq = np.fft.rfftfreq((shape(mm))[0], d=dt)

freqmask = np.zeros(shape(spe))
freqmask[984,:,:]=1
spe_filtered = spe*freqmask
mminv  = np.real(np.fft.irfft(spe_filtered, axis=0, norm='ortho'))


tn='20220202_154506'
v,f= wf_spectrum(tn, freq=tn:

#20220202
#analysis sweep at different gap
gap = [0,30,70,105,140,189,225,260]
tnlist=['20220202_142211','20220202_154506','20220202_164103','20220202_164515','20220202_164929','20220202_165336','20220202_165906','20220202_170335']
wdlist = ['SPLATT_Test_2022-02-02T14-22-11-621.wdd','SPLATT_Test_2022-02-02T15-45-06-886.wdd','SPLATT_Test_2022-02-02T16-41-03-073.wdd','SPLATT_Test_2022-02-02T16-45-15-294.wdd','SPLATT_Test_2022-02-02T16-49-29-567.wdd','SPLATT_Test_2022-02-02T16-53-36-217.wdd','SPLATT_Test_2022-02-02T16-59-06-452.wdd','SPLATT_Test_2022-02-02T17-03-35-908.wdd']

lab = []
ntn = size(tnlist)
cols = ['k','y','b','c','g','r','m','y']
fl=th.fileList(tn, fold=sp.basepath4d)
npo=int(size(fl)/2+1)
ff = np.zeros(npo)
vv = np.zeros([npo,ntn])

for i in range(ntn):
    f,v = sp.sweep_analysis(tnlist[i],tiltsum=1)
    lab_i = ("Gap: %s" %gap[i])
    lab.append(lab_i)
    ff=f
    vv[:,i] = v
fname = '/mnt/jumbo/SPLATT/OPTData/20220202_sweep_freq.fits'
pyfits.writeto(fname,ff)
fname = '/mnt/jumbo/SPLATT/OPTData/20220202_sweep_spec_tilt.fits'
pyfits.writeto(fname,vv)


vast = np.zeros([npo,ntn])

sp.z2fit=np.array([1,2,3,4,5,6])
for i in range(ntn):
    f,v = sp.sweep_analysis(tnlist[i],tiltsum=np.array([5,6]))
    lab_i = ("Gap: %s" %gap[i])
    lab.append(lab_i)
    ff=f
    vast[:,i] = v
fname = '/mnt/jumbo/SPLATT/OPTData/20220202_sweep_astigm.fits'
pyfits.writeto(fname,vast)


for i in range(ntn):
    plot(ff,vv[:,i],cols[i], label=lab[i])

legend()

wd = sp.openfile(wdlist[0])
npacc = np.shape(wd)[1]
nsacc = int(npacc/2+1)
fa = np.zeros(nsacc)
va = np.zeros([4,nsacc, ntn])
for i in range(ntn):
    spa, facc = sp.acc_spectrum(wdlist[i]) 
    fa=facc
    va[:,:,i]=spa


for i in range(ntn):
    plot(ff,vast[:,i],cols[i], label=lab[i])

legend()


theacc = 0
for i in range(ntn):
    plot(fa,va[theacc,:,i],cols[i], label=lab[i])

legend()
   


g = vv*0
for i in range(ntn):
    g[:,i] = vv[:,i]/vv[:,0]

for i in range(ntn):
    plot(sp.runningMean(ff,10),sp.runningMean(g[:,i],10),cols[i], label=lab[i])

legend()



#helium test
#air
tnlist0 = ['20220211_153030','20220211_153407','20220211_153635','20220211_153949','20220211_154116','20220211_154244']
gaplist0= np.array([0,38,62,409,145,181])
f0,v0=sp.sweep_analysis_sequence(tnlist0)
ntn0=size(tnlist0)
lab0 = []
for i in gaplist0:
    lab0.append("Gap: %s" %i)
plot(f0,v0[:,0], 'k', label=lab0[0])

for i in range(ntn0-1):
    plot(f0,v0[:,i+1], label=lab0[i+1])
legend()

tnlist1=['20220211_153030','20220211_163059','20220211_163355','20220211_163602','20220211_163826','20220211_164048']
gaplist1=np.array([0,37,60,105,141,186])
f1,v1=sp.sweep_analysis_sequence(tnlist1)
ntn1=size(tnlist1)
lab1 = []
for i in gaplist1:
    lab1.append("Gap: %s" %i)
plot(f1,v1[:,0], 'k', label=lab1[0])

for i in range(ntn1-1):
    plot(f1,v1[:,i+1], label=lab1[i+1])
legend()

a0=v0*0
a1=v1*0
f = f0[1:]
for i in range(ntn0):
    a0[1:,i]=v0[1:,i]/v0[1:,0]
for i in range(ntn1):
    a1[1:,i]=v1[1:,i]/v1[1:,0]

nr = 10
figure(1)
clf()
plot(sp.runningMean(f,nr),sp.runningMean(a0[1:,0],nr), 'k', label=lab0[0])
for i in range(ntn0-1):
    plot(sp.runningMean(f,nr),sp.runningMean(a0[1:,i+1],nr), label=lab0[i+1])
    
legend() 
yscale('log')
xlim(50,110)
title('Test in Air')
figure(2)
clf()
plot(sp.runningMean(f,nr),sp.runningMean(a1[1:,0],nr), 'k', label=lab1[0])
for i in range(ntn1-1):
    plot(sp.runningMean(f,nr),sp.runningMean(a1[1:,i+1],nr), label=lab1[i+1])

legend()
yscale('log')
xlim(50,110)
title('Test in Helium')

#comparing air and helium result at small gap regime
gapid =5
clf()
yscale('log')
xlim(40,120)

plot(sp.runningMean(f,nr),sp.runningMean(a0[1:,0],nr), 'k', label=lab0[0])
plot(sp.runningMean(f,nr),sp.runningMean(a0[1:,gapid],nr), label='Air:'+lab0[gapid])
plot(sp.runningMean(f,nr),sp.runningMean(a1[1:,gapid],nr), label='He:'+lab1[gapid])
legend()
yscale('log')
xlim(40,120)



#analysis of single frequency data
flist = np.array([5,10,15,20,25,30,40,50,60,10,80,90,100,110])
tnlist0 = ['20220211_155826','20220211_155757','20220211_155736','20220211_155717','20220211_155656','20220211_155602','20220211_155430','20220211_155411','20220211_155350','20220211_155332','20220211_155312','20220211_155250','20220211_155137','20220211_155106']

tnlist1 =['20220211_164525','20220211_164541','20220211_164558','20220211_164615','20220211_164631','20220211_164648','20220211_164254','20220211_164308','20220211_164323','20220211_164337','20220211_164353','20220211_164407','20220211_164422','20220211_164436']

nf = size(flist)
s0 = []
s1 = []
for i in range(nf):
    bd = np.array([flist[i]-1,flist[i]+1])
    sp.peak_analysis(tnlist0[i], bound=bd, sumtilt=1)
    sp.peak_analysis(tnlist1[i], bound=bd, sumtilt=1)


tn=['20220211_154140','20220211_163804']
gapl=['','']
ff0,rr0=sp.sweep_analysis_sequence(tn, gapl)
ff0=ff0[1:]
rr1=rr0[1:,0]/v1[1:,0]
rr2=rr0[1:,1]/v1[1:,0]


plot(sp.runningMean(ff0,nr),sp.runningMean(rr1,nr))
plot(sp.runningMean(ff0,nr),sp.runningMean(rr2,nr))




tnlist=['20220211_164048','20220211_164110','20220211_164919','20220211_165013','20220211_165036','20220211_165122','20220211_165147','20220211_165316','20220211_165419','20220211_165539','20220211_171734']
jj,kk=sp.sweep_analysis_sequence(tnlist)
ntn = size(tnlist)

jj1=jj[1:]
kk1=kk*0
for i in range(ntn):
    kk1[1:,i]=kk[1:,i]/v0[1:,0]

for i in range(ntn):
    rr1=rr0[1:,0]/v1[1:,0]
rr2=rr0[1:,1]/v1[1:,0]

lab = ''
plot(ones(150), col=('k'))
for i in range(ntn):
    plot(sp.runningMean(jj1,nr),sp.runningMean(kk1[1:,i],nr),label=tnlist[i] )


legend()
plot(ones(150), color=('k'))
yscale('log');xlim(55,110)
i=9
plot(sp.runningMean(jj1,nr),sp.runningMean(kk1[1:,i],nr),label =tnlist[i] );i=i+1
plot(sp.runningMean(jj1,nr),sp.runningMean(kk1[1:,i],nr),label =tnlist[i] );i=i+1


i=4

plot(sp.runningMean(f,nr),sp.runningMean(a1[1:,i+1],nr), label=lab1[i+1]+' He')
plot(sp.runningMean(f,nr),sp.runningMean(a0[1:,i+1],nr),'-x',label=lab0[i+1]+' Air')



#Helium test, inspection of Acc data
wd0 = 'SPLATT_Test_2022-02-11T15-42-44-986.wdd' #air, 180 um gap
wd1 = 'SPLATT_Test_2022-02-11T16-40-48-194.wdd' #helium 
wd11= 'SPLATT_Test_2022-02-11T16-55-39-600.wdd'
wd12 = 'SPLATT_Test_2022-02-11T17-17-34-044.wdd'
spa0, facc = sp.acc_spectrum(wd0)
spa1, facc = sp.acc_spectrum(wd1)
spa11, facc = sp.acc_spectrum(wd11)
spa12, facc = sp.acc_spectrum(wd12)

yscale('log');xlim(55,110)

plot(facc, spa1[0,:])
plot(facc, spa11[0,:])
plot(facc, spa12[0,:])
plot(facc, spa1[1,:])


wdlist0 = ['SPLATT_Test_2022-02-11T16-40-48-194.wdd','SPLATT_Test_2022-02-11T16-41-10-686.wdd','SPLATT_Test_2022-02-11T16-49-19-145.wdd','SPLATT_Test_2022-02-11T16-50-37-010.wdd','SPLATT_Test_2022-02-11T16-51-22-774.wdd','SPLATT_Test_2022-02-11T16-53-16-451.wdd','SPLATT_Test_2022-02-11T16-54-19-426.wdd','SPLATT_Test_2022-02-11T16-55-39-600.wdd','SPLATT_Test_2022-02-11T17-17-34-044.wdd']


ntn = size(wdlist0)
spa0, facc = sp.acc_spectrum(wdlist0[0])
nnp = shape(spa0)
aac = np.zeros([ntn,nnp[0], nnp[1]])
for i in range(ntn):
    spa, facc = sp.acc_spectrum(wdlist0[i])
    fa=facc
    aac[i,:,:]=spa





#autoinduction test
tnlist = ['SPLATT_CoilTest_2022-03-18T17-54-53-708.wdd','SPLATT_CoilTest_2022-03-18T17-55-15-708.wdd' ,'SPLATT_CoilTest_2022-03-18T17-55-39-716.wdd' ,'SPLATT_CoilTest_2022-03-18T17-56-03-075.wdd' ,'SPLATT_CoilTest_2022-03-18T17-56-31-278.wdd' ,'SPLATT_CoilTest_2022-03-18T17-56-49-625.wdd' ,'SPLATT_CoilTest_2022-03-18T17-57-10-124.wdd' ,'SPLATT_CoilTest_2022-03-18T17-57-39-906.wdd' ,'SPLATT_CoilTest_2022-03-18T17-58-02-835.wdd' ,'SPLATT_CoilTest_2022-03-18T17-58-21-483.wdd' ,'SPLATT_CoilTest_2022-03-18T17-58-42-380.wdd']
ntn=size(tnlist)
freq=[10,20,30,40,50,60,70,80,90,100,110]
a = []
b = []
for i in range(ntn):
    spa, f=sp.acc_spectrum(tnlist[i])
    p0,p1 = sp.find_peak(spa[2,:], freq=f, bound=[freq[i]-1,freq[i]+1])
    a.append(sp.acc_integrate(p0,p1))
    b0,b1 = sp.find_peak(spa[3,:], freq=f, bound=[freq[i]-1,freq[i]+1])
    b.append(b0)

a = np.array(a)
b = np.array(b)

plot(freq, b,'-x')
xtitle('Freq [Hz]')
ytitle('Coil voltage [V]')
displ=mean(a)
obbdisp = 100e-9
rr=obbdisp/displ
coilres = 10.
curr = b[9]/coilres*rr

plot(freq, curr,'-x')
xtitle('Freq [Hz]')
ytitle('Current [Amp]')
title('Auto-induced current on OBB')


#test with different amplitudes
tnlist=['SPLATT_CoilTest_2022-03-21T10-37-48-568.wdd',
'SPLATT_CoilTest_2022-03-21T10-38-09-972.wdd','SPLATT_CoilTest_2022-03-21T10-38-31-076.wdd', 'SPLATT_CoilTest_2022-03-21T10-39-12-055.wdd','SPLATT_CoilTest_2022-03-21T10-39-29-947.wdd','SPLATT_CoilTest_2022-03-21T10-39-48-439.wdd']
ntn=size(tnlist)
amp=[0.1,0.2,0.3,0.5,0.8,1]
a = []
b = []
freq = 90
for i in range(ntn):
    spa, f=sp.acc_spectrum(tnlist[i])
    p0,p1 = sp.find_peak(spa[2,:], freq=f, bound=[88,92])
    a.append(sp.acc_integrate(p0,p1))
    b0,b1 = sp.find_peak(spa[3,:], freq=f, bound=[88,92])
    b.append(b0)

a = np.array(a)
b = np.array(b)
plot(amp, b,'-x')
xlabel('Magnet disp. command [V]')
ylabel('Coil voltage [V]')
displ=mean(a)
obbdisp = 100e-9
coilres = 10.

curr = b/coilres
plot(a, curr,'-x')
xlabel('Magnet oscill. amplitude [um]')
ylabel('Coil current [A]')

