import sys,pandas as pd, matplotlib , matplotlib.pyplot as plt, matplotlib.lines , numpy as np, math, pylab,time, multiprocessing

import ROOT
#import uproot
import root_pandas
masses = {211: 0.13957018, -211: 0.13957018, 321: 0.493677, -321: 0.493677, 2212: 0.93827208816, -2212: 0.93827208816}
def mixed_quantities(E, e_px, e_py, e_pz, h_px, h_py, h_pz,h_pid, h2_px, h2_py, h2_pz,h2_pid):
    beam = ROOT.TLorentzVector()
    beam.SetPxPyPzE(0,0,E,E)
    electron = ROOT.TLorentzVector()
    electron.SetPxPyPzE(e_px,e_py,e_pz,np.sqrt(0.000511**2+e_px**2+e_py**2+e_pz**2))
    target = ROOT.TLorentzVector()
    target.SetPxPyPzE(0,0,0,masses[2212])
    virtual_photon = beam-electron
    cm = virtual_photon+target
    hadron = ROOT.TLorentzVector(h_px,h_py,h_pz, np.sqrt(masses[h_pid]**2+h_px**2+h_py**2+h_pz**2))
    hadron2 = ROOT.TLorentzVector(h2_px,h2_py,h2_pz, np.sqrt(masses[h2_pid]**2+h2_px**2+h2_py**2+h2_pz**2))                                                                                                  
    
    
    
    di = {"mx_eh1x":(beam+target-electron-hadron).M(), 
          "mx_eh2x":(beam+target-electron-hadron2).M(),
          "mx_eh1h2x":(beam+target-electron-hadron-hadron2).M(),
          "pair_mass":(hadron+hadron2).M()}
    
    virtual_photon = beam-electron
    cm = virtual_photon+target
    
    hadron2.RotateZ(-cm.Phi());
    hadron2.RotateY(-cm.Theta());
    hadron2.Boost(0,0,-cm.Beta());
    
    di['h2_cm_pt'] = hadron2.Pt()
    di['h2_cm_rap'] = hadron2.Rapidity()
    di['h2_cm_ph'] = hadron2.Phi()
    
    
    return di
    

#new version
def mix_from_singles(df, binvars=''.split(), nbins=1,maxEvents=None, nAssocPerTrigger=1,j0=0):
    start = time.perf_counter()
    df['h_z'] = df.z
    
    Q2 = df['Q2']
    h_p = df['h_p']
    #try:
    N = len(df)
    if maxEvents != None:
        N = min(N, maxEvents)
    h_pid = df['h_pid']
    partitions = {}
    for var in binvars:
        partitions[var] = [df[var].quantile(i/3) for i in range(nbins+1)]
    df_out = pd.DataFrame()
    
    fields = []
    
    d = {}
    for name in df.columns:
        fields.append(name)
        if name[:2] == 'h_':
            d["h1_" + name[2:]] = []
            d["h2_" + name[2:]] = []
        else :
            d[name] = []
        
    
    d['diff_phi_cm'] = []
    d['diff_phi_lab'] = []
    d['diff_rap_cm'] = []
    d['diff_eta_cm'] = []
    
    def same_bin(i,j):
        #print("checking ", i, j)
        if Q2[i] == Q2[j] or Q2[i] == 0 or Q2[j] == 0 or h_pid[i] == 0 or h_pid[j] == 0: #don't mix with the same event
            return False
        if h_p[i] < h_p[j]:
            return False
        for var in binvars:
            good = False
            xi = df[var][i]
            xj = df[var][j]
            pvar = partitions[var]
            for b in range(nbins):
                if xi >= pvar[b] and xi < pvar[b+1] and\
                        xj >= pvar[b] and xj < pvar[b+1]:
                    good=True
                    break
            if not good:
                return False
        return True
            
    
    #i = trigger
    #j = associated
    j = j0
    for i in range(N):
        
        if i % 10000 == 0:
            duration = time.perf_counter()-start
            print("%.1f"%(i/N*100),"% complete, time so far: ",duration//3600,"hours", 
                  (duration//60)%60, "minutes", int(duration % 60), "seconds")
        
        found_mix= 0
        if df['h_z'][i] < 0.4 or abs(df['h_pid'][i])!= 211:
            continue
        #print('check1')
        for ii in range(nAssocPerTrigger):
            for k in range(100):
                #print("check2")
                j += 1
                if j >= N:
                    j = 0
                #if df['h_z'][j] > df['h_z'][i]:
                #    continue;
                if same_bin(i,j):
                    found_mix = k+1
                    break
            
            #print("check3")
            if found_mix == 0:
                continue
            #print("check4")
            for name in fields:
                if name[:2] == 'h_':
                    d["h1_" + name[2:]].append(df[name][i])
                    if not 'cm' in name and not '_z' in name:
                        d["h2_" + name[2:]].append(df[name][j])
                else :
                    d[name].append(df[name][i])

            d['diff_phi_lab'].append(df.h_ph[i]-df.h_ph[j] 
                                    + 2*np.pi*(df.h_ph[i]-df.h_ph[j]<-np.pi)
                                    -2*np.pi*(df.h_ph[i]-df.h_ph[j]>np.pi))
            #for the cm variables, use the trigger particle's electron info to get the cm frame
            di = mixed_quantities(df.E[i], df.e_px[i], df.e_py[i], df.e_pz[i],df.h_px[i], df.h_py[i], df.h_pz[i],df.h_pid[i],
                                                df.h_px[j], df.h_py[j], df.h_pz[j],df.h_pid[j])
            for name in di.keys():
                if not name in d.keys():
                    d[name] = []
                d[name].append(di[name])
            
            #d["h2_cm_pt"].append(pt)
            #d["h2_cm_rap"].append(rap)
            #d["h2_cm_ph"].append(phi)
            d["h2_z"].append(df.h_z[j]*df.nu[j]/df.nu[i])
            
            d['diff_phi_cm'].append(df.h_cm_ph[i]-di['h2_cm_ph'] 
                                    + 2*np.pi*(df.h_cm_ph[i]-di['h2_cm_ph']<-np.pi)
                                    -2*np.pi*(df.h_cm_ph[i]-di['h2_cm_ph']>np.pi))

            d['diff_rap_cm'].append(df.h_cm_rap[i]-di['h2_cm_rap'])
            #d['diff_eta_cm'].append(df.h_cm_eta[i]-df.h_cm_eta[j])

        

    
    
    duration = time.perf_counter()-start
    print("total time: ",duration//3600,"hours", (duration//60)%60, "minutes", int(duration % 60), "seconds")
    #os.system('say "your program has finished"')
    #print("%.5f"%(100*np.mean(df['mix_found']!=0)),"% of the events have been mixed")
    for key in list(d.keys()):
        print(key, len(d[key]))
        if len(d[key])==0:
            del d[key]
    ret = pd.DataFrame.from_dict(d)
    ret['diff_phi_lab'] = ret.h1_ph - ret.h2_ph + 2*np.pi*(ret.h1_ph-ret.h2_ph<-np.pi)-2*np.pi*(ret.h1_ph-ret.h2_ph>np.pi)
    return ret

if __name__ == '__main__':
    infile = sys.argv[1]
    outfile=sys.argv[2]
    print("infile is", infile, ";  outfile is", outfile)
    if ".root" in infile:
        df = root_pandas.read_root(infile,'hadrons')
        print("opened file ", infile, ";  type(df) is ", type(df))
        #df = root_pandas.read_root(infile,'hadrons')
    if ".pkl" in infile:
        df = pd.read_pickle(infile)
        
    maxEvents = None
    nAssocPerTrigger=1
    j0=0
    singleThread=False
    for arg in sys.argv[2:]:
        if '-N=' in arg:
            maxEvents=int(arg[3:])
        elif '-n=' in arg:
            nAssocPerTrigger=int(arg[3:])
        elif arg == '-s':
            singleThread=True
    def process(j0,i):
        df_mixed = mix_from_singles(df, binvars=''.split(), nbins=1,maxEvents=maxEvents, nAssocPerTrigger=1,j0=j0)
        filei=outfile +("%s.pkl"%i)                                                   
        pd.to_pickle(df_mixed,filei)
        print("wrote to file "+filei)
    
    processes = []
    if not singleThread:
        for i in range(0,10):
            j0 = i*500
            if not singleThread:
                p = multiprocessing.Process(target=process, args=(j0,i))
                processes.append(p)
                p.start()
        for process in processes:
            process.join()

    else :
        df_mixed = mix_from_singles(df, binvars=''.split(), nbins=1,maxEvents=maxEvents, nAssocPerTrigger=nAssocPerTrigger)
        df_mixed.to_root(outfile,key='dihadrons')
        
        
    
'''    import threading
    class myThread (threading.Thread):
        def __init__(self, j0):
            threading.Thread.__init__(self)
            self.j0 = j0
            
        def run(self):
            df_mixed = mix_from_singles(df, binvars=''.split(), nbins=1,maxEvents=maxEvents, nAssocPerTrigger=1,j0=self.j0)
            filei=outfile +("%s.pkl"%self.j0)
            pd.to_pickle(df_mixed,filei)
            print("wrote to file "+filei)
    for j0 in range(nAssocPerTrigger):
        a = myThread(j0)
        a.start()
'''
