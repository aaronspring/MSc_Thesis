""" 
   Usage: python dic_alk_plot.py 

   Depends on: numpy, netCDF4, socket, matplotlib, warnings

   module load python
"""  

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
#########################################################################
def chemcon(ptiestu,psao,ptho):
  """
  Set chemical constants 
 
  Takes:
  ptiestu water depth of level: ke+1
  psao salt : (ie,je,ke)
  ptho temp : (ie,je,ke)

  Returns:
  pakw3  KW
  pakb3  KB
  pak13  K1
  pak23  K2
  prrrcl 
  pak0   K0

  """
  #     -----------------------------------------------------------------
  #         1. SET HALF PRECISION CONSTANTS
  #           --- ---- --------- ---------

  perc = 1.e-2
  smicr = 1.e-6
  #  3. SET CONVERSION FACTOR SALINITY -> CHLORINITY
  #               ------ ------- -- ---- ------- -- ----------
  #               (AFTER WOOSTER ET AL., 1969)
  salchl = 1./1.80655

  #     -----------------------------------------------------------------
  #*         8. SET COEFFICIENTS FOR SEAWATER PRESSURE CORRECTION OF
  #             (AFTER ZEEBE and WOLF-GLADROW (2001)
  #             ------- -------- --- ---------- ----- --- -------- -
  #             in general the pressure dependency for dissoziation contants is given by
  #             ln (K_p / K_0) = - delta V/RT *P + 0.5 delta kappa/RT *P**2
  #             with delta V corresponds to molal volume  := delta V = dev + devt* TC + devt2*TC**2
  #              and delta kappa to compressibility; here neglected
  #             Thus  K_p = K_0 * exp (-(dekv + devkt*TC + devkt2*TC**2)*CP)
  #             with CP = P/RT
  #             K_p = K_0 * exp( -devk - devkt*TC - devkt2*TC**2)*CP
  #             Note that in table A.11.1 all orginal devk (1. column) are negative; herefore, the sign
  #             already is included in the definition of the constant
  #             devt2 for carbon, calcite and aragonit is equal 0.0
  #
  devk1  = 25.50
  devk1t = 0.1271
  devk2  = 15.82
  devk2t = -0.0219
  devkb  = 29.48
  devkbt = 0.1622
  devkbt2= 2.6080*1.E-3 
  devkw  = 25.60
  devkwt = 0.2324       
  devkwt2= -3.6246*1.E-3

  #     -----------------------------------------------------------------
  #*        11. SET GAS CONSTANT  for pressure correction after Millero(1995)
  #             --- --------- --- --------
  #             in cm**3 bar / mol K (from Zeebe and Wolf-Gladrow)
  #
  rgas = 83.131

  #     -----------------------------------------------------------------
  #*        12. SET BORON CONCENTRATION IN SEA WATER
  #*            IN G/KG PER O/OO CL ACCORDING
  #             TO RILEY AND SKIRROW, 1965 (P.250)
  #             -- ----- --- -------- ---- ------- -
  #
  bor1 = 0.00023

  #    -----------------------------------------------------------------
  #*        13. SET INVERSE OF ATOMIC WEIGHT OF BORON [G**-1]
  #             (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
  #             ----- -- ------- -------- ----- ----- ---- ---------------
  #
  bor2 = 1./10.82



  prrrcl = salchl * 1.025 * bor1 * bor2

  #     -----------------------------------------------------------------
  #*        15. SET COEFF. FOR 1. DISSOC. OF CARBONIC ACID
  #             DOE (1994)
  #             ------------------------------------------
  #
  c1 = np.array([2.83655, -2307.1266,-1.5529413, 0.207608410,4.0484, \
                 0.0846834,-0.00654208,0.001005])


  #     -----------------------------------------------------------------
  #*        16. SET COEFF. FOR 2. DISSOC. OF CARBONIC ACID
  #             DOE (1994)
  #             ------------------------------------------
  #
  c2 =np.array([ -9.226508,-3351.6106,-0.2005743,0.106901773,23.9722, \
                  0.1130822,-0.00846934,0.001005])

  #
  #     -----------------------------------------------------------------
  #*        17. SET COEFF. FOR 1. DISSOC. OF BORIC ACID
  #             DOE (1994) (based on by Dickson 1990)
  #             ---------------------------------------
  #
  cb = np.array([-8966.90,-2890.53,-77.942, 1.728,-0.0996,148.0248,\
                 137.1942,1.62142,24.4344,25.085,0.2474,0.053105])
  #
  #     -----------------------------------------------------------------
  #*        18. SET COEFF. FOR ION PRODUCT OF WATER
  #             DOE (1994)
  #             ---------------------------------------
  #
  cw =np.array([ 148.96502,-13847.26,-23.6521,118.67,-5.977,1.0495,-0.01615])

 #     -----------------------------------------------------------------
  #*        19. SET VOLUMETRIC SOLUBILITY CONSTANTS FOR CO2 IN ML/L
  #             WEISS, R. F. (1974)
  #             CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A
  #             NON IDEAL GAS. MARINE CHEMISTRY, VOL. 2, 203-215.
  #     -----------------------------------------------------------------
  c0=np.array([9345.17,-60.2409, 23.3585, 0.023517, -0.00023656, 0.0047036])
  # 


  #*     22.1 APPROX. SEAWATER PRESSURE AT U-POINT DEPTH (BAR)
  #  ----------------------------------------------------------------

  p = np.multiply(1.025e-1 , ptiestu)


  #
  #*    22.1 SET ABSOLUTE TEMPERATURE
  # ----------------------------------------------------------------

  t = np.add(ptho, 273.15)
  ti = np.divide(1.0, t)
  log_t = np.log(t)
  q = np.multiply(t,perc)
  q2 =np.square(q)
  log_q = np.log(q)
  qi = np.divide(1.0,  q)
  
  s = np.maximum(25. , psao)
  sqrt_s = np.sqrt(s)

  #
  #        21.3 SOLUBILITY OF GASES
  #       -------------------------------------------------
  #
  #      CO2 CARBON DIOXIDE, SEA SURFACE
  #      --------------------------------------
  #         LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
  #
  cek0 = np.multiply(c0[0],ti) + c0[1] + np.multiply(c0[2] , log_q) + \
         np.multiply(s , (c0[3] + np.multiply(c0[4],t) + np.multiply(c0[5],q2)))
  ak0 = np.multiply(np.exp(cek0),smicr)

  #*    22.5 PK1, PK2 OF CARBONIC ACID, PKB OF BORIC ACID, PKW OF WATER
  # ----------------------------------------------------------------
  ck1 = np.add(c1[0],np.multiply(c1[1],ti))  + np.multiply(c1[2],log_t)  \
                - np.multiply(np.add(c1[3],np.multiply(c1[4],ti)),sqrt_s) \
                + np.multiply(c1[5],s) + np.multiply(c1[6],np.multiply(s,sqrt_s))\
                + np.log(np.add(1.,np.multiply(-c1[7],s)))


  ck2 = np.add(c2[0]  , np.multiply(c2[1],ti))  + np.multiply(c2[2],log_t)\
               -np.multiply(np.add(c2[3],np.multiply(c2[4],ti)),sqrt_s) \
                + np.multiply(c2[5],s) + np.multiply(c2[6],np.multiply(s,sqrt_s)) \
                + np.log(np.add(1.,np.multiply(-c2[7],s)))
 
  ckb= np.multiply(np.add(cb[0],np.multiply(cb[1],sqrt_s)) + np.multiply(cb[2],s) \
                  + np.multiply(cb[3],np.multiply(s,sqrt_s)) + np.multiply(cb[4],np.square(s)),ti) \
                  + np.add(cb[5],np.multiply(cb[6],sqrt_s)) + np.multiply(cb[7],s)      \
                  -np.multiply(np.add(cb[8],np.multiply(cb[9],sqrt_s)) + np.multiply(cb[10],s),log_t )  \
                  + np.multiply(cb[11],np.multiply(sqrt_s,t))

  ckw = np.add(cw[0], np.multiply(cw[1],ti)) + np.multiply(cw[2],log_t)                        \
                    +np.multiply(sqrt_s,np.add(np.multiply(cw[3],ti),cw[4]) +np.multiply(cw[5],log_t)) \
                    + np.multiply(cw[6],s)


  #
  #*    22.6 K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O)
  # ----------------------------------------------------------------
  ak1 = np.exp(ck1)
  ak2 = np.exp(ck2)
  akb = np.exp(ckb)
  akw = np.exp(ckw)



  #       FORMULA FOR CP
  #       ----------------------------------------------

  cp = np.divide(p ,np.multiply(rgas,t))

  
  tc = ptho
  tc2 = np.square(tc)



  #
  #       PRESSURE DEPENDENT KB, K1,K2 AND KW
  #       ----------------------------------------------

  pakb3 = np.multiply(akb,np.exp(np.multiply(cp,(np.add(devkb,(-np.multiply(devkbt,tc)-np.multiply(devkbt2,tc2)))))))
  pak13 = np.multiply(ak1,np.exp(np.multiply(cp,np.add(devk1,-np.multiply(devk1t,tc)))))
  pak23 = np.multiply(ak2,np.exp(np.multiply(cp,np.add(devk2,-np.multiply(devk2t,tc)))))
  pakw3 = np.multiply(akw,np.exp(np.multiply(cp,np.add(devkw,np.multiply(-devkwt,tc))- np.multiply(devkwt2,tc2))))  
   
  return pakw3,pakb3,pak13,pak23,prrrcl,ak0


#########################################################################


ens='178'
years='1985_1992'
#ens='143'
#years='1995_2002'

area='global_-60--50'
#area='global_-50--40'

######CONFIGURATION
# salinity
sao=34.04
tho=277.81-273.15#-0.25

# thermal and salt effect
thermal=-8.4
salt=0.7

filename='best.'+ens+'_'+years+"."+area+"."

print(filename)

f1=Dataset('/work/mh0727/m300524/post/hamocc/trends_ymjm/ens'+ens+'_'+years+'_positive_trend8/co2_separation_Lovenduski/lkm0'+ens+'_mm_'+years+'_co2flux_separation_lovenduski.nc.fldmean_'+area+'.nc.ymjm.nc')
actalk=f1.variables['surf_alkali'][0,0,0,0]*1.e3
actdic=f1.variables['surf_sco212'][0,0,0,0]*1.e3 # concentration units: m mol C m-3 
print("actalk=",actalk,"; actdic=",actdic)
mld=f1.variables['zmld'][0,0,0,0]
### BIO
# coex90
#actcoex90=f1.variables['intpp'][0,0,0,0]*3600.*24.*30./mld*1000.*12.
actcoex90=f1.variables['coex90'][0,0,0,0]*3600.*24.*30./mld*1000.*12.
#actcaex90=f1.variables['intdelcar'][0,0,0,0]*3600.*24.*30./mld*1000.*12.
actcaex90=f1.variables['caex90'][0,0,0,0]*3600.*24.*30./mld*1000.*12.
diccoex90=actdic-actcoex90
alkcoex90=actalk+16./122.*actcoex90
diccaex90=actdic-actcaex90
alkcaex90=actalk-2.*actcaex90
dicbio=actdic-actcoex90-actcaex90
alkbio=actalk+16./122.*actcoex90-2.*actcaex90
print("BIO: alkcoex90=",alkcoex90," diccoex90=",diccoex90,";  alkcaex90=",alkcaex90," diccaex90=",diccaex90,";  alkbio=",alkbio," dicbio=",dicbio)

### exchange
actco2flux=f1.variables['co2flux'][0,0,0,0]/12.*1000.*1000.*3600.*24.*30./mld*12.
dicco2flux=actdic-actco2flux
alkco2flux=actalk

### freshwater
actpem=f1.variables['pem'][0,0,0,0]*3600.*24.*30.*12. # org units m/s to m/yr
print("actpem=",actpem,"; meters/month water flux into ocean")
dicpem=actdic*(mld+actpem)/mld
alkpem=actalk*(mld+actpem)/mld
print("FW: dicpem=",dicpem,"; alkpem=",alkpem)


# rsidual
resdic=actdic-(dicco2flux-actdic)-(dicbio-actdic)-(dicpem-actdic)
resalk=actalk-(alkco2flux-actalk)-(alkbio-actalk)-(alkpem-actalk)
print("resdic = ",resdic, " resalk = ",resalk)


### END OF PERIOD

actalke=f1.variables['surf_alkali'][7,0,0,0]*1.e3
actdice=f1.variables['surf_sco212'][7,0,0,0]*1.e3
mlde=f1.variables['zmld'][7,0,0,0]
### BIO
# coex90
#actcoex90e=f1.variables['intpp'][7,0,0,0]*3600.*24.*30./mlde*1000.*12.
actcoex90e=f1.variables['coex90'][7,0,0,0]*3600.*24.*30./mlde*1000.*12.
actcaex90e=f1.variables['caex90'][7,0,0,0]*3600.*24.*30./mlde*1000.*12.
#actcaex90e=f1.variables['intdelcar'][7,0,0,0]*3600.*24.*30./mlde*1000.*12.
diccoex90e=actdice-actcoex90e
alkcoex90e=actalke+16./122.*actcoex90e
diccaex90e=actdice-actcaex90e
alkcaex90e=actalke-2.*actcaex90e
dicbioe=actdice-actcaex90e-actcoex90e
alkbioe=actalke+16./122.*actcoex90e-2.*actcaex90e
print("BIO: alkcoex90e=",alkcoex90e," diccoex90e=",diccoex90e,";  alkcaex90e=",alkcaex90e," diccaex90e=",diccaex90e,";  alkbioe=",alkbioe," dicbioe=",dicbioe)

### exchange
actco2fluxe=f1.variables['co2flux'][7,0,0,0]/12.*1000.*1000.*3600.*24.*30./mlde*12.
dicco2fluxe=actdice-actco2fluxe
alkco2fluxe=actalke

### freshwater
actpeme=f1.variables['pem'][7,0,0,0]*3600.*24.*30.*12. # org units m/s
print("actpeme=",actpeme,"; meters/month water flux into ocean")
dicpeme=actdice*(mlde+actpeme)/mlde
alkpeme=actalke*(mlde+actpeme)/mlde
print("FW: dicpeme=",dicpeme,"; alkpeme=",alkpeme)

#res																																			#thermal #salt
resdice=actdice-(dicco2fluxe-actdice)-(dicbioe-actdice)-(dicpeme-actdice)
resalke=actalke-(alkco2fluxe-actalke)-(alkbioe-actalke)-(alkpeme-actalke)
print("resdice = ",resdice, " resalke = ",resalke)


f1.close()








# DIC, & alkalinity ranges
step=0.000001
start=0.002240
ende=0.002280
alks= np.arange(start,ende+step,step)

 #step=0.000005
start=0.00202
ende=0.00212
dic= np.arange(start,ende+step,step)



ph=np.zeros((np.shape(alks)[0],np.shape(dic)[0]))
logco3=np.zeros((np.shape(alks)[0],np.shape(dic)[0]))
pco2=np.zeros((np.shape(alks)[0],np.shape(dic)[0]))

tiestu=1.
#Chemical constants:
akw3,akb3,ak13,ak23,rrrcl,ak0=chemcon(tiestu,sao,tho)


alki=0
for alkali in alks:
 dici=0
 for sco212 in dic: 

  hi=6.5e-9
    
  # Get hi, co3 corresponding to dic, alk
  for iter in range(0,9): # usually converges already after 2-3 iterations 

     h = hi
     c = sco212
     t1  = np.divide(h,ak13)
     t2  = np.divide(h,ak23)
     ak2 = ak23
     akw = akw3
     bt  = np.multiply(sao,rrrcl)
     akb = akb3
     alk = alkali
     # a = c*(2+t2)/(1+t2+t1*t2)+akw/h-h+bt/(1+h/akb)-alk
     a   = np.multiply(c,np.divide(np.add(2.,t2),np.add(np.add(1.,t2),np.multiply(t2,t1))))+np.divide(akw,h)-h+np.divide(bt,np.add(1.,np.divide(h,akb)))-alk
     # dadh = c*(1/(ak2*(1+t2+t1*t2)) - c*(2+t2)(1/ak2+2*t1/ak2)/((1+t2+t1*t2)**2) - akw/h**2 - 1 - (bt/akb)/((1+h/akb)**2)
     dadh=np.multiply(c,np.divide(1.,np.multiply(ak2,1.+np.add(t2,np.multiply(t2,t1)))))  - np.multiply(c,np.multiply(np.add(2,t2),np.divide(np.divide(1.,ak2)+ np.multiply(2.,np.divide(t1,ak2)),np.square(np.add(1.,t2+np.multiply(t2,t1)))))) -np.add(np.divide(akw,np.square(h)),1.) -np.divide(np.divide(bt,akb),np.square(np.add(1.,np.divide(h,akb))))

     # dddhhh=a/dadh
     dddhhh=np.divide(a,dadh)
     h = h-dddhhh
     hi = np.maximum(hi-dddhhh,1.e-10)
     #co3 = c/(1+hi*(1+hi/ak13)/ak23) 
     co3 =   np.divide(c,np.add(1., np.multiply(hi , np.divide(np.add(1.,np.divide(hi,ak13)) , ak23))))
		 #print("iteration = ",iter)

  #pH=-log10(hi)  
  ph[alki,dici]=-np.log10(hi)
   #logco3=-log10(co3)
  logco3[alki,dici]=-np.log10(co3)
 # co2 = DIC/(1.+K1/h1+K1*K2/h1**2) 
  co2=np.divide(sco212,1. + np.divide(ak13,hi) + np.divide(np.multiply(ak13,ak23),np.square(hi)))
 # pco2 = CO_2/k0 
  pco2[alki,dici]=np.divide(co2,ak0)
  dici+=1
 alki+=1

#print("alks",alks*1.e6)
#print("actalk",actalk)
#print("dic",dic*1.e6)
#print("actdic",actdic)
ialk=np.argmin(np.abs(np.subtract(alks*1.e6,actalk)))
idic=np.argmin(np.abs(np.subtract(dic*1.e6,actdic)))
#print(pco2) 
print("RESULTS Start of Period")
pco2_actual=pco2[ialk,idic]
pco2_bio_caex90=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkcaex90))),np.argmin(np.abs(np.subtract(dic*1.e6,diccaex90)))]
pco2_bio_coex90=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkcoex90))),np.argmin(np.abs(np.subtract(dic*1.e6,diccoex90)))]
pco2_bio=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkbio))),np.argmin(np.abs(np.subtract(dic*1.e6,dicbio)))]
pco2_ex=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkco2flux))),np.argmin(np.abs(np.subtract(dic*1.e6,dicco2flux)))]
pco2_fw=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkpem))),np.argmin(np.abs(np.subtract(dic*1.e6,dicpem)))]
pco2_res=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,resalk))),np.argmin(np.abs(np.subtract(dic*1.e6,resdic)))]#-thermal-salt


print("RESULTS End of Period")
pco2_actuale=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,actalke))),np.argmin(np.abs(np.subtract(dic*1.e6,actdice)))]
pco2_bio_caex90e=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkcaex90e))),np.argmin(np.abs(np.subtract(dic*1.e6,diccaex90e)))]
pco2_bio_coex90e=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkcoex90e))),np.argmin(np.abs(np.subtract(dic*1.e6,diccoex90e)))]
pco2_bioe=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkbioe))),np.argmin(np.abs(np.subtract(dic*1.e6,dicbioe)))]
pco2_exe=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkco2fluxe))),np.argmin(np.abs(np.subtract(dic*1.e6,dicco2fluxe)))]
pco2_fwe=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,alkpeme))),np.argmin(np.abs(np.subtract(dic*1.e6,dicpeme)))]
pco2_rese=pco2[np.argmin(np.abs(np.subtract(alks*1.e6,resalke))),np.argmin(np.abs(np.subtract(dic*1.e6,resdice)))]-thermal-salt


diff=pco2_actuale-pco2_actual
print("RESULTS")
print("Start End Difference in [ppm]")
print("pco2: ",pco2_actual," ",pco2_actuale," ",pco2_actuale-pco2_actual)
print("pco2_bio_caex90: ",pco2_bio_caex90," ",pco2_bio_caex90e," ",pco2_bio_caex90e-pco2_bio_caex90-diff)
print("pco2_bio_coex90: ",pco2_bio_coex90," ",pco2_bio_coex90e," ",pco2_bio_coex90e-pco2_bio_coex90-diff)
print("pco2_bio: ",pco2_bio," ",pco2_bioe," ",pco2_bioe-pco2_bio-diff)
print("pco2_ex: ",pco2_ex," ",pco2_exe," ",pco2_exe-pco2_ex-diff)
print("pco2_fw: ",pco2_fw," ",pco2_fwe," ",pco2_fwe-pco2_fw-diff)
print("pco2_res: ",pco2_res," ",pco2_rese," ",pco2_rese-pco2_res)
print("pco2 thermal trend:				",thermal)
print("pco2 salt trend: 				",salt)




# PLOT
fac=1.e6
#fig1=plt.figure(1)
# contour line spacing for LOG(co3) 
#levs=np.arange(3.1,6.8,0.1)
#cs=plt.contour(dic*fac,alks*fac,logco3,colors='k',levels=levs)
#plt.clabel(cs,inline=1,fontsize=10,fmt="%.1f")
#plt.xlabel(r'DIC ($\mu$mol L$^{-1}$)')
#plt.ylabel(r'Alkalinity ($\mu$mol L$^{-1}$)')
#plt.title('-log10(co3) S='+str(sao)+' T='+str(tho)+'degC')
#plt.savefig('dic_alk_logco3.png',dpi=fig1.dpi)

#fig2=plt.figure(2)
# contour line spacing for pH 
#levs=np.arange(6.5,9.5,0.1)
#cs=plt.contour(dic*fac,alks*fac,ph,colors='k',levels=levs)
#plt.clabel(cs,inline=1,fontsize=10,fmt="%.1f")
#plt.xlabel(r'DIC ($\mu$mol L$^{-1}$)')
#plt.ylabel(r'Alkalinity ($\mu$mol L$^{-1}$)')
#plt.title('pH S='+str(sao)+' T='+str(tho)+'degC')
#plt.savefig('dic_alk_pH.png',dpi=fig2.dpi)

fig3=plt.figure(3)
# contour line spacing for pCO2
#l1=np.arange(0.,1000.,5.)
levs=np.arange(0.,1000.,5.)
#l2=np.arange(1000.,10000.,200.)
#levs=np.concatenate((l1,l2))
cs=plt.contour(dic*fac,alks*fac,pco2,colors='k',levels=levs)
#plt.plot(actdicy,actalky,'ro')
#plt.plot(actdice,actalke,'bx')
#plt.plot(actdicc5,actalkc5,'go')
#plt.plot(actdicc5e,actalkc5e,'gx')
#plt.plot(actdicu,actalku,'mo')

#plt start of period
plt.plot(diccoex90,alkcoex90,'co',label="coex")
plt.plot(diccaex90,alkcaex90,'mo',label="caex")
plt.plot(dicbio,alkbio,'go',label="bio")
plt.plot(dicco2flux,alkco2flux,'ro',label="co2flux")
plt.plot(dicpem,alkpem,'yo',label="fw")
plt.plot(resdic,resalk,'bo',label="res")

plt.plot(actdic,actalk,'ko',label="total")


#plt end of period
plt.plot(diccoex90e,alkcoex90e,'cs')
plt.plot(diccaex90e,alkcaex90e,'ms')
plt.plot(dicbioe,alkbioe,'gs')
plt.plot(dicco2fluxe,alkco2fluxe,'rs')
plt.plot(dicpeme,alkpeme,'ys')

plt.plot(actdice,actalke,'ks')
plt.plot(resdice,resalke,'bs')

plt.clabel(cs,inline=1,fontsize=10,fmt="%.2f")
plt.xlabel(r'DIC ($\mu$mol L$^{-1}$)')
plt.ylabel(r'Alkalinity ($\mu$mol L$^{-1}$)')

plt.legend(bbox_to_anchor=(0., 0.5), loc=2, borderaxespad=0.)

plt.title('pco2 (ppm) S='+str(sao)+' T='+str(tho)+'degC')
plt.savefig(filename+'dic_alk_pco2.png',dpi=fig3.dpi)

plt.show()

