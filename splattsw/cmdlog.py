

freqM=np.array([0.5,25])
ampPI = 2
nfr = 2000
rp.splatt_trigger(freqM[0], freqM[1], ampPI)
webdaq.start_schedule()
tn=comm4d.capture(nfr);print(tn); comm4d.produce(tn);wdfile=sp.lastwdfile(); print(wdfile); sp.freq4d=freqM[1]:q

from SPLATT.splattsw import splatt_acq as acq

acq.acq(0.5,50, bound=np.array([0.2,0.8]))
acq.acq(1,100, bound=np.array([0.8,2]))
acq.acq(5,100)
acq.acq(10,100)
acq.acq(15,150)
acq.acq(20,100)
acq.acq(25,250)
acq.acq(30,240)
acq.acq(40,240)
acq.acq(50,250)
acq.acq(60,240)
acq.acq(72,252,  ampPI=3)
acq.acq(80,260,  ampPI=3)
acq.acq(88,264,  ampPI=3)
acq.acq(100,250, ampPI=3)
acq.acq(110,264, ampPI=3)
acq.acq(120,264, ampPI=3)


acq.acq(40,432)
acq.acq(50,432)
acq.acq(60,432)
acq.acq(72,432,  ampPI=3)
acq.acq(80,432,  ampPI=3)
acq.acq(88,432,  ampPI=3)
acq.acq(100,432, ampPI=3)
acq.acq(110,432, ampPI=3)
acq.acq(120,432, ampPI=3)
acq.acq(150,432, ampPI=3)

acq.acq(40,240)
acq.acq(50,250)
acq.acq(60,240)
acq.acq(72,252,  ampPI=3)
acq.acq(80,260,  ampPI=3)
acq.acq(88,264,  ampPI=3)
acq.acq(100,250, ampPI=3)
acq.acq(110,264, ampPI=3)
acq.acq(120,264, ampPI=3)


acq.acq_pulse(0.1, 50)

acq.acq_bench(2,10)
acq.acq_bench(2,15)
acq.acq_bench(2,20)
acq.acq_bench(2,25)
acq.acq_bench(2,30)
acq.acq_bench(2,40)
acq.acq_bench(2,50)
acq.acq_bench(2,60)
acq.acq_bench(3,72)
acq.acq_bench(3,80)
acq.acq_bench(3,88)
acq.acq_bench(3,100)
acq.acq_bench(3,110)
acq.acq_bench(3,120)

acq.acq_sweep(2500, 250)

f0=80
for i in range(80,120,4):
    acq.acq(i,250,ampPI=10)

tn0=['20220114_182245','20220114_182439','20220114_182634','20220114_182829','20220114_183024','20220114_183218','20220114_183414','20220114_183608','20220114_183804','20220114_183959']
tn1=['20220114_175913','20220114_180108','20220114_180303','20220114_180458','20220114_180652','20220114_180850','20220114_181045','20220114_181242','20220114_181437','20220114_181633']
fvec = np.array([80,84, 88, 92, 96, 100, 104, 108, 112,116]
nn = size(fvec)
for i in np.arange(nn):
    sp.peak_analysis(tn0[i],bound=np.array([fvec[i]-2,fvec[i]+2]),sumtilt=1)
    


#measurement 20220202
acq.acq_sweep(4200, 300) #14 s

acq.nframes = 1000
acq.acq(1,100, bound=np.array([0.8,2]))
acq.acq(5,100)
acq.acq(10,100)
acq.acq(15,150)
acq.acq(20,100)
acq.acq(25,250)
acq.acq(30,300)




from SPLATT.splattsw import splatt_acq_he as acq_he
from SPLATT.splattsw.devices import redpitaya as rp

freq4d = 250.03
duration_sweep= 12
duration_4d  =10
nframes = int(duration_4d*freq4d) #2500

#dry run, misura di prova
#impostare sweep
#impostare durata webdaq 10s
#far termalizzare il RefBody
rp.wave_on(2)
#configurare sweep
#verificare che sweep parte
#impostare Freq 4D a 250 Hz
# 2 repetitions for each gap
#Shell RIP
acq_he.acq_sweep(nframes, freq4d)
acq_he.acq_sweep(nframes, freq4d)
#Shell set, gap xx





#measurement in standard air 
#Shell RIP
acq_he.acq_sweep(nframes, freq4d)
acq_he.acq_sweep(nframes, freq4d)
#Shell set, gap xx
acq_he.acq_sweep(nframes, freq4d)
acq_he.acq_sweep(nframes, freq4d)


 
#gap = 220 um
acq_he.acq_sweep(nframes, freq4d)
acq_he.acq_sweep(nframes, freq4d)


#measurement in He
freq4d=150.04
acq_he.acq(5,freq4d)
acq_he.acq(10,freq4d)
acq_he.acq(15,freq4d)
acq_he.acq(20,freq4d)
acq_he.acq(25,freq4d)
acq_he.acq(30,freq4d)
freq4d=250.03
acq_he.acq(40,freq4d)
acq_he.acq(50,freq4d)
acq_he.acq(60,freq4d)
acq_he.acq(70,freq4d,  ampPI=3)
acq_he.acq(80,freq4d,  ampPI=3)
acq_he.acq(90,freq4d,  ampPI=3)
acq_he.acq(100,freq4d, ampPI=3)
acq_he.acq(110,freq4d, ampPI=3)

#processing
comm4d.produce(tn)
comm4d.frames_transfer(tn)


