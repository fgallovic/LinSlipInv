from numpy import log10,sin,rad2deg,deg2rad,arctan2,sqrt,loadtxt
import pyproj

NAME="2009 L'Aquila"
DATE='04/06/2009'
REF='Gallovic et al. (2015)'
EventTAG='s2009LAquilaGALLOVIC'
origlat=42.339
origlon=13.381

out=open('srcmod.fsp', 'w')
out.write('% ----------------------------------  FINITE-SOURCE RUPTURE MODEL  --------------------------------\n')
out.write('% \n')
out.write('% Event: '+NAME+'    '+DATE+'    ['+REF+']\n')
out.write('% EventTAG: '+EventTAG+'\n')
out.write('% \n')

# Read input.dat
inputdat=open('input.dat','r')
dum=inputdat.readline()
dum=inputdat.readline()
dum=inputdat.readline()
T, TS=inputdat.readline().split()[0:2]
T,TS=float(T),float(TS)
dum=inputdat.readline()
dum,NSeg=inputdat.readline().split()
NSeg=int(NSeg)
dum=inputdat.readline()
NRseis,NRgps=inputdat.readline().split()
dum=inputdat.readline()
NLW=inputdat.readline().split()
NL=[int(x) for x in NLW[0::2]]
NW=[int(x) for x in NLW[1::2]]
dum=inputdat.readline()
M0=inputdat.readline().strip()
Mw=log10(float(M0))/1.5-6.07
dum=inputdat.readline()
SDR=inputdat.readline().split()
strike=[float(x) for x in SDR[0::3]]
dip=[float(x) for x in SDR[1::3]]
rake=[float(x) for x in SDR[2::3]]
dum=inputdat.readline()
hypodepth=[float(x)/1.e3 for x in inputdat.readline().split()]
dum=inputdat.readline()
lengwidt=inputdat.readline().split()
leng=[float(x)/1.e3 for x in lengwidt[0::2]]
widt=[float(x)/1.e3 for x in lengwidt[1::2]]
Dx=[leng[i]/float(NL[i]) for i in range(NSeg)]
Dz=[widt[i]/float(NW[i]) for i in range(NSeg)]
dum=inputdat.readline()
epicLW=inputdat.readline().split()
epicL=[float(x)/1.e3 for x in epicLW[0::2]]
epicW=[float(x)/1.e3 for x in epicLW[1::2]]
Htop=[hypodepth[i]-sin(deg2rad(dip[i]))*(widt[i]-epicW[i]) for i in range(NSeg)]
dum=inputdat.readline()
np=int(inputdat.readline())
dt=T/float(np)
Ssvd=int(TS/dt+1.)
dum=inputdat.readline()
dum=inputdat.readline()
dum=inputdat.readline()
dum=inputdat.readline()
f1,f2=inputdat.readline().split()

seg=0   #so far saving only the first segment!!!!
out.write('% Loc  : LAT  = '+str(origlat)+'\t LON = '+str(origlon)+'\t DEP = '+str(hypodepth[seg])+' km\n')
out.write('% Size : LEN  = '+str(leng[seg])+' km\t WID = '+str(widt[seg])+' km\t Mw = '+str(Mw)+'\t Mo = '+M0+' Nm\n')
out.write('% Mech : STRK = '+str(strike[seg])+' \t DIP = '+str(dip[seg])+' \t RAKE = '+str(rake[seg])+' \t Htop = '+str(Htop[seg])+' km\n')
out.write('% Rupt : HypX = '+str(epicL[seg])+' km\t HypZ = '+str(epicW[seg])+' km\t avTr = -99.0 s\t avVr = -99.0 km/s')

out.writelines("""
% 
% ----------------------------------  inversion-related parameters  --------------------------------
% 
""")

out.write('% Invs :  Nx  =  '+str(NL[seg])+' \t Nz = '+str(NW[seg])+'\t Fmin = '+f1+' Hz\t Fmax = '+f2+' Hz\n')
out.write('% Invs :  Dx  =  '+str(Dx[seg])+' km\t Dz  = '+str(Dz[seg])+' km\n')
out.write('% Invs :  Ntw =  '+str(Ssvd)+' \t Nsg =  '+str(NSeg)+'\t \t \t (# of time-windows,# of fault segments)\n')
out.write('% Invs :  LEN =  '+str(dt)+' s \t SHF =  '+str(dt)+' s\t \t \t (time-window length and time-shift)\n')
out.write('% SVF  :  delta  \t \t \t \t \t (type of slip-velocity function used)\n')
out.write('% \n')
out.write('% Data :  \tSGM \tTELE \tTRIL \tLEVEL \tGPS \tINSAR \tSURF \tOTHER \tOther\n')
out.write('% Data :  \t'+NRseis+' \t0 \t0 \t0 \t'+NRgps+' \t0 \t0 \t0 \t0\n')
out.write('% PHImx:  	999 	0 	0 	0 	0 	0 	0 	0 	0\n')
out.write('% Rmin : 	999.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0\n')

out.writelines("""% 
% --------------------------------------------------------------------------------------------------
% 
% VELOCITY-DENSITY STRUCTURE 
""")
crustal=open('crustal.dat','r')
dum=crustal.readline()
dum=crustal.readline()
Nlay=int(crustal.readline())
out.write('% No.  of layers =  '+str(Nlay)+' \n')
out.write('% \n')
out.write('%\tDEPTH\tP-VEL\tS-VEL\tDENS\tQp\tQs\n')
out.write('%	[km]	[km/s]	[km/s]	[g/cm^3]\n')
dum=crustal.readline()
dum=crustal.readline()
for i in range(Nlay):
  out.write('% '+crustal.readline())

out.writelines("""% 
% --------------------------------------------------------------------------------------------------
% % 11-Mar-2016 Frantisek Gallovic (gallovic@karel.troja.mff.cuni.cz)
% --------------------------------------------------------------------------------------------------
% 
% SOURCE MODEL PARAMETERS
""")
out.write('% \tNsbfs = '+str(NL[seg]*NW[seg])+' subfaults\n')
out.writelines("""% 
% 		X,Y,Z coordinates in km; SLIP in m 
% 		if applicable: RAKE in deg, RISE in s, TRUP in s, slip in each TW in m 
% 
%   Coordinates are given for top-center of each subfault or segment: |'| 
%   Origin of local coordinate system at epicenter: X (EW) = 0,  Y (NS) = 0  
%    LAT       LON      X==EW    Y==NS     Z       SLIP     """)
for i in range(Ssvd):
  out.write('  TW'+str(i+1)+'  rakeTW'+str(i+1))
out.write('\n% -------------------------------------------------------------------------------------------------- \n')
sources=loadtxt('sources.dat')
mtilde=loadtxt('mtilde.dat')
g=pyproj.Geod(ellps='WGS84')
#staN=12.25
#staE=-0.5
dist=sqrt(sources[:,1]**2+sources[:,2]**2)*1.e3
az=rad2deg(arctan2(sources[:,2],sources[:,1]))
for i in range(NL[seg]*NW[seg]):
  (sublon, sublat, backaz) = g.fwd(origlon,origlat,az[i],dist[i])
  sr=mtilde[i*Ssvd:(i+1)*Ssvd]*dt
  dum=[sublat,sublon,sources[i,2],sources[i,1],sources[i,3],sr.sum()]
  for j in range(Ssvd): dum.append(sr[j]);dum.append(rake[seg])
  out.write(' '.join(format(x, "8.4f") for x in dum)+'\n')
