"""
Convert input muon and track data files to same format as PF group

Tracks: 96 bits

N_BITS_TRACK_PHI   = 19
N_BITS_TRACK_INVPT = 15
N_BITS_TRACK_ETA   = 14
N_BITS_TRACK_Z0    = 11
N_BITS_TRACK_SECTOR = 5

Muons: 32 bits
N_BITS_MUON_PT  = 9;
N_BITS_MUON_ETA = 9;
N_BITS_MUON_PHI = 10;
"""

def writeHeader():
    header = []
    lines = '='*120+"\n"
    links = 'Input '
    for i in range(67):
        links += 'LINK_{0:02d} '.format(i)
    links += "\n"

    header.append( lines )
    header.append( links )
    header.append( lines )

    return header

def file2list(filename,comment="#"):
    """Load text file and dump contents into a list"""
    listOfFiles = open( filename,'r').readlines()
    listOfFiles = [i.rstrip('\n') for i in listOfFiles if not i.startswith(comment)]
    return listOfFiles


class TkMu(object):
    def __init__(self):
        self.pt      = -1
        self.rinv    = -1
        self.sinhEta = -1
        self.eta     = -1
        self.phi     = -1
        self.charge  = -1
        self.z0      = -1
        self.sector  = -1

class Event(object):
    def __init__(self):
        self.tracks = []
        self.muons  = []


tfile = file2list("../hls-modules/data/CleanTracksAll_binary.dat")
mfile = file2list("../hls-modules/data/hw_muon_data.dat")
ofile = open("outfile.txt","w")


N_MAX    = 38
N_TRACKS = 11  # 33 LINKS
N_MUONS  = 5   #  5 LINKS
events   = [Event() for _ in xrange(10000)]  # list of Event() objects


## read tracks
for line in tfile:
    # track metadata
    values = line.split(' ')

    if line.startswith("BX"):
        eventNumber = int(values[3])-1
        phisector   = "{0:b}".format(int(values[7]))
    else:
        # track properties
        values2 = line.split(' ')
        rinv    = values[1]
        phi     = values[2]
        sinhEta = values[3]
        z0      = values[4]
    
        event = events[eventNumber]
    
        tk = TkMu()
        tk.rinv    = rinv
        tk.sinhEta = sinhEta
        tk.phi     = phi
        tk.z0      = z0
        tk.sector  = phisector
        event.tracks.append(tk)
    
        events[eventNumber] = event



## read muons
eventNumber = -1
for line in mfile:
    if line.startswith("\n"): continue
 
    values = line.split(' ')
    if line.startswith("Begin"):
        eventNumber = int(values[2])-1
    elif any( [line.startswith(i) for i in ["MuonB","MuonO","MuonE","Empty"]] ):
        frame3 = values[1][:32]
        frame2 = values[1][32:64]
        frame1 = values[1][64:]

        pt_str      = frame1[23:32]
        quality_str = frame3[19:23]
        eta_str     = frame1[10:19]
        phi_str     = frame3[1:11]
        q_str       = frame2[31]

        if eventNumber>=0:
            event = events[eventNumber]
            
            mu = TkMu()
            mu.pt      = pt_str
            mu.eta     = eta_str
            mu.phi     = phi_str
            mu.charge  = q_str
            event.muons.append(mu)
            
            events[eventNumber] = event
    elif line.startswith("End"):
        eventNumber = -1

# write to new data file
empty  = "0x00000000"
HEADER = writeHeader()

for h in HEADER:
    ofile.write(h)

for e,event in enumerate(events):
    data = "{0:#0{1}x}".format(e,7)

    n_tk = len(event.tracks)
    tk_fillEmpty = N_TRACKS-n_tk

    for nt,track in enumerate(event.tracks):
        t1 = "{0:#0{1}x}".format(int(track.phi,base=2),10)
        t2 = "{0:#0{1}x}".format(int(track.rinv+track.sinhEta,base=2),10)
        t3 = "{0:#0{1}x}".format(int(track.z0+track.sector,base=2),10)
        tdata = t1+" "+t2+" "+t3
        data += " "+tdata
        if nt==N_TRACKS-1: break
    for _ in range(tk_fillEmpty):
        data += " "+empty

    n_mu = len(event.muons)
    mu_fillEmpty = N_MUONS-n_mu

    for muon in event.muons:
        m1 = "{0:#0{1}x}".format(int(muon.pt+muon.eta+muon.phi,base=2),10)
        data += " "+m1

    for _ in range(mu_fillEmpty):
        data += " "+empty

    ofile.write(data+"\n")




## THE END
