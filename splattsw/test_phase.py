#this script is intended to check if there is a dependancy of identified spectrum peak vs signal phase and poor sampling (5, 6 samples per cycle)
#check also the HW verification with the RedPitaya in test_peak_phaseRP.py
nn=10000
t=np.arange(nn)/nn
f=100
ph=90
ph=ph/180*np.pi
v1=sin(2*np.pi*f*t-ph)
v0 = sin(2*np.pi*f*t)
clf()
plot(v0)
plot(v1)
xlim(0,1000)
dt=1/nn
spe1  = np.fft.rfft(v1,  norm='ortho')
freq = np.fft.rfftfreq(v1.shape[0], d=dt)
spe0  = np.fft.rfft(v0,  norm='ortho')
freq = np.fft.rfftfreq(v0.shape[0], d=dt)
normv   = np.sqrt(spe1.shape[0])
spe1  = (np.abs(spe1)) / normv
spe0  = (np.abs(spe0)) / normv


print(max(spe0))
print(max(spe1))
nx   = np.sqrt(spe.shape[0])
spe  = (np.abs(spe)) / nx
#the spectrum is correct
sp = np.fft.fft(v)
npi = int(np.size(sp)/2)
sp = sp[0:npi+1]
ph = np.angle(sp)*180/np.pi

idp = np.argmax(spe)#20
phi=ph[idp]
print(phi)

#decimated sampling
dd=120
st=0.4
st=int(st*500)
v11= v1[st::dd]
v00= v0[0::dd]
dt1=dt*dd
spe11  = np.fft.rfft(v11,  norm='ortho')
freq = np.fft.rfftfreq(v11.shape[0], d=dt1)
spe00  = np.fft.rfft(v00,  norm='ortho')
freq = np.fft.rfftfreq(v00.shape[0], d=dt1)
normv   = np.sqrt(spe11.shape[0])
spe11  = (np.abs(spe11)) / normv
spe00  = (np.abs(spe00)) / normv

print(max(spe00))
print(max(spe11))


