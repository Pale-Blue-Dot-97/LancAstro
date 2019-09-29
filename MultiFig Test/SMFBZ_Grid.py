import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('text', usetex = True)
import matplotlib.cm as cm
import math
import pandas as pd
import numpy as np
from astropy.io import fits
from scipy import optimize
from labellines import labelLine, labelLines


# FIGURE AND AXIS SIZES ======================================================
mpl.rcParams['xtick.labelsize'] = 14        # Size of x-tick labels
mpl.rcParams['ytick.labelsize'] = 14        # Size of y-tick labels
mpl.rcParams['axes.labelsize'] = 14         # Size of axes labels
mpl.rcParams['image.origin'] = 'lower'      # 

mpl.rcParams['axes.linewidth'] = 2.5        # Size of axes line width
mpl.rcParams['xtick.major.size'] = 14       # Size of major x-ticks 
mpl.rcParams['xtick.minor.size'] = 6        # Size of minor x-ticks
mpl.rcParams['xtick.major.width'] = 2.5     # Width of major x-ticks  
mpl.rcParams['xtick.minor.width'] = 1.5     # Width of minor x-ticks
mpl.rcParams['xtick.direction'] = 'in'    # Sets ticks to be inside
mpl.rcParams['ytick.major.size'] = 14       # Size of major y-ticks
mpl.rcParams['ytick.minor.size'] = 6        # Size of minor y-ticks
mpl.rcParams['ytick.major.width'] = 2.5     # Width of major y-ticks
mpl.rcParams['ytick.minor.width'] = 1.5     # Width of minor y-ticks
mpl.rcParams['ytick.direction'] = 'in'    # Sets ticks to be inside

# Sets font style and size
mpl.rcParams.update({'font.size': 12, 'font.weight': 'bold'}) 
############################################################################






















# Create arrays of Band values, and redshift values for use 
NBs=np.array([392, 501, 711, 816]) # NB filter band numbers
IAs=np.array([427, 464, 484, 505, 527, 574, 624, 679, 709, 711, 738, 767, 827]) # IA filter band numbers
Zs=[2.22,2.5,3.1, 3.9, 4.7, 5.4] # All Redshift values
zErrors=[0.02,0.1,0.4,0.3,0.2,0.5] # Redshift error
newZs=[2.5,3.1, 3.9, 4.7, 5.4] # reduced Redshift values, neglecting 2.22 as it causes an error when the NB filters are not used
newzErrors=[0.1,0.4,0.3,0.2,0.5] # Redshift range
uniID=np.concatenate((NBs,IAs),axis=None) #Create an array of both NB and IA band values
IA_VOLUMES=np.array([4.0e6,4.2e6,4.3e6,4.3e6,4.5e6,4.9e6,5.2e6,5.5e6,5.1e6,1.2e6,5.1e6,5.5e6,4.9e6]) # Surveyed volume in IA filters
NB_VOLUMES=np.array([0.6e6,0.9e6,1.2e6,1.8e6]) # Surveyed volume in NB volume
uniVol=np.concatenate((NB_VOLUMES,IA_VOLUMES),axis=None) # #Create an array of both NB and IA surveyed volume values
Mmins_orig=[9.08,9.08,9.08,9.41,9.41,9.41,9.42,9.42,9.42,9.57,9.57,9.58,9.92]
Mmins=[9.1,9.1,9.1,9.4,9.4,9.4,9.4,9.4,9.4,9.6,9.6,9.6,9.9]
Mminsbins=[9.1,9.4,9.4,9.6,9.9]
Mminsall=Mmins+Mminsbins
All_Zerrs=[]







# Open Relevant fits file, and assign column arrays to variables.
hdul =fits.open('JOSH4_SC4K_LAEs_catalogue_Reduced_sorted.fits')
data = hdul[1].data
ID = np.array(data['ID_SC4K_1'])
BAND_ARRAY=[int(i[7:10]) for i in ID] #ID=SC4K_IA427_1234 becomes BAND=427
BAND_ARRAY=np.array(BAND_ARRAY)
uniqueID=np.unique(BAND_ARRAY)
Comments = np.array(data['Comments'])
Res = np.array(data['Response'])
AllMass = np.array(data['Mstar'])
Redshift = np.array(data['LAE_Redshift'])
Xray= np.array(data['XRAY_DET'])
Radio=np.array(data['RADIO_DETECTED'])
Flux=np.array(data['Flux'])

Flux=[math.log(x,10) for x in Flux]
Flux=np.array(Flux)
uniqueZ=np.unique(Redshift)
print(uniqueZ)
# Initialise arrays to hold arrays of values obtained in initial for loop 
zMass=[] 
Phis=[]
Phinocos=[]
PhiErrors=[]
CoPhiErrors=[]
Masxs=[]
Masx0s=[]
Masx0=[]
Masx=[]
zID=[]
zF=[]
zRed=[]
zlen=[]
edg=[]
Phi0s=[]
fit_logPs=[]
Param1s=[]
Param2s=[]
Param1errs=[]
Param2errs=[]


# Inititalise arrays to hold column values
Colcountnocomps = []
Colcountcomps=[]
Colphinocomps = []
Colphis=[]
Colphierrnocomps =[]
ColCofas = []
Colphierrs=[]
ColMasxs=[]
ColVols = []
ColWids = []
#print(uniID)
#print(uniqueZ)


# Define Schechter function parameters
ALPHA=-1.4
logM_star=10.6
######################################################################################################
# Define Schechter function

def schfunct(x,a,b):
    #b=logM_star
    x=10.**x
    a=10.**a
    b=10.**b
    
    return np.log10(((a/b)*(x/b)**(ALPHA)*np.exp(-x/b)*x*np.log(10.)))

# Define Schechter Phi only function

def schfunct1(x,a):
    b=logM_star
    x=10.**x
    a=10.**a
    b=10.**b
    
    return np.log10(((a/b)*(x/b)**(ALPHA)*np.exp(-x/b)*x*np.log(10.)))
######################################################################################################
def make_bins(a,N,ID,F,Mmin):
    
    #Set max and in values of the array
    MAX=11.
    MIN=9.
    # Calculates the width of the bins
    WID=(MAX-MIN)/(N)
    
    # Find the maximum and minimum Mass axis points
    
    CENT_MIN=MIN+WID/2.
    CENT_MAX=MAX-WID/2.
    #Print statement for each bin completion
    print("Constructing bins")
    # Constructs the centre of the bins from the defined mid, max and width
    bins = np.arange(CENT_MIN,CENT_MAX+WID,WID)
    
    #Create MASK to eliminate values outside of defined max and min values
    MASK=(bins<=MAX)*(bins>=MIN)
    
    bins=bins[MASK]
    
    # Obtain array of bin edges for output
    edg = np.arange(MIN,MAX+WID,WID)
    #MASK2=(len(edg)<=len(bins)+1)
    #edg=edg[MASK2]
    
    if len(edg)>N+1:
       edg=edg[:-1]
       
    #obtain array of unique band arrays
    uID=np.unique(ID)
    
    #Initialise empty arrays for use by counts and individual filter volumes
    countsFilt=[]
    FiltsVol=[]
    
    #run a loop for each individual ID, to obtain an array of all filter volumes
    for u in uID:
        #print(u)
        MASK3=(ID==u)  #Define MASK to obtain indexes for respective ID in total ID array in this make bin iteration  
        MASK4= (IAs==u) #Define MASK to obtain indexes for respective ID in all the measured IDs, defined at the start of the script
        FiltID=ID[MASK3] #Filtered array of IDs in this iteration for this specific ID
        FiltVol=IA_VOLUMES[MASK4] #Obtain respective volume
        countFilt=len(FiltID) #obtain total count of this array (NOT USED, but might be useful)
        countsFilt.append(countFilt) #Form Array of total count of this array (NOT USED, but might be useful)
        FiltsVol.append(FiltVol[0]) # Form Array of filter volumes      
    
    #Vol=np.sum([x*y for x,y in zip(countsFilt,FiltsVol)])/sum(countsFilt)
    
    #obtain total volume surveyed in this make bin iteration, this is an output variable
    Vol=np.sum(np.unique(FiltsVol))
    print(Vol)
    
    
    # Creates arrays of the same length as bins to hold the counts per bin
    counts = []  
    counts = bins - bins
    countsComp = []
    countsComp = bins - bins
    
    
    # Loops through all data points, for each one cycling through the bins to
    # see if it resides within it. If so, it is added to counts for that bin
    for j in range(len(a)):
        for k in range(len(bins)):
                if a[j] <= bins[k] + WID/2.0:
                    if a[j] >= bins[k] - WID/2.0:
                        iD=ID[j] #Obtain individual iD value in order to determine which completness function to open
                        flu=F[j] #Obtain individual flux value in order to determine where to interpolate completeness function to, to determine completness
                        hdul2 =fits.open('COMPLETENESS_functions/COMPLETENESS_IA%d.fits'%iD) #Open relevant completeness function
                        data2 = hdul2[1].data
                        flux = np.array(data2['lineflux']) #Open flux column of completness function 
                        comp = np.array(data2['completeness']) #Open completness column of completness function 
                        comp= np.interp(flu, flux, comp) # interpolate completeness function to relevant flux value
                        comp=comp/100 # convert completeness from percentage
                        countsComp[k] = countsComp[k] + 1/comp # add individual source to count and apply completeness effect
                        counts[k] = counts[k] + 1 # add individual source to count with no completeness 
    # Normalises the data
    NORM = sum(counts) # (NOT USED) 

    
    co = [] # Array of y-axis plot points (count values for each bin) without completeness
    cocomp = [] # Array of y-axis plot points (count values for each bin) with completeness
    for j in range(len(bins)):
        co.append(counts[j])
        cocomp.append(countsComp[j])
    return bins, co, cocomp, edg, WID, NORM, MIN, MAX, Vol
#########################################################################################


# Run a for loop to obtain mass, flux values for all sources in catalogue with desired redshift/band, for each individual redshift
for j in uniqueZ:
    print(j)
    k=0.01  # sets an 'error' for the masking redshift as the convertion to float can add an error which means the mask does not work
    #Create Mask to remove all non desired values of mass/ID/flux.
    MASK =(Radio==False)*(Xray==False)*(Redshift<=j+k)*(Redshift>=j-k)*(BAND_ARRAY!=392)*(BAND_ARRAY!=501)*(BAND_ARRAY!=816)
    
    New_mass=AllMass[MASK] # Masked mass array
    New_ID=BAND_ARRAY[MASK] # Masked ID array
    New_Red=Redshift[MASK] # Masked redshift
    uni_New_Redshift=np.unique(Redshift[MASK]) # Unique redshifts masked redshift array
    New_flux=Flux[MASK] # Masked redshift
    
    if len(New_ID)>0: # Condition to avoid the addition of empty arrays to tensors
        zRed.append(uni_New_Redshift[0])# Masked redshift tensor
        zMass.append(New_mass) # Masked mass tensor
        zID.append(New_ID) # Masked ID tensor
        zF.append(New_flux) # Masked flux tensor
        #zlen.append(len(New_mass))
        
        
# Run a for loop to obtain mass, flux values for all sources in catalogue with desired redshift/band, for each redshift bin
for j,k in zip(newZs,newzErrors):
    print(j)
    #Create Mask to remove all non desired values of mass, ID, and flux.
    MASK =(Radio==False)*(Xray==False)*(Redshift<=j+k)*(Redshift>=j-k)*(BAND_ARRAY!=392)*(BAND_ARRAY!=501)*(BAND_ARRAY!=816)
    
    New_mass=AllMass[MASK] # Masked mass array
    New_ID=BAND_ARRAY[MASK] # Masked ID array
    New_Red=Redshift[MASK] # Masked redshift
    New_flux=Flux[MASK] # Masked redshift
    New_New_mass=[] # Initialise Mass tensor
    New_New_ID=[] #Initialise ID tensor
    New_New_flux=[] #Initialise flux tensor
    
    # for loop runs through each individual redshift in each bin and then appends their property values to overaqll tensor, Not actually needed anymore, but may be useful thus left in
    for m in zRed:
        #
        l=0.01
        #MASK7=(New_Red<=m+l)*(New_Red>=m-l)*(New_mass>=Mmins[i])
        MASK7=(New_Red<=m+l)*(New_Red>=m-l) # Add mask for each 
        New_New_mass=np.append(New_New_mass,New_mass[MASK7]) # Create tensor for each ID mass in redshift bin
        New_New_ID=np.append(New_New_ID,New_ID[MASK7]) # Create tensor for each ID  in redshift bin
        New_New_flux=np.append(New_New_flux,New_flux[MASK7]) # Create tensor for each ID flux in redshift bin

    print(New_New_mass)
    if len(New_ID)>0:
        
        zMass.append(New_New_mass)
        zID.append(New_New_ID)
        zF.append(New_New_flux)
        
for M,I,F,i in zip(zMass,zID,zF,np.arange(0,len(zMass),1)):
    
    Mask=(str(M)!='nan')  # to remove all nan values 
    
    #'Clean' respective mass/ID/Flux arrays
    CleanM=M[Mask]
    CleanI=I[Mask]
    CleanF=F[Mask]
    
    #Since mask produces 1d tensor, obtain useful array using [0] index check
    useI=CleanI[0]
    useF=CleanF[0]
    useM=CleanM[0]
    nbins = 15     #Number of bins to divide data into
    
    #Use make_bins method to obtain necessary values
    Masx, Phinoco,Phi,edg, Wid, Norm, Min, Max, Vol=make_bins(useM, nbins,useI,useF,Mminsall[i]) 
    
    zlen.append(sum(Phinoco)) #obtain total number of counts for each redshift/redshift bin
    print(sum(Phinoco))
    
    # Convert to numpy arrays/create arrays of width and volume of the same length as phi values
    Phi=np.array(Phi)
    Volarr = np.array(np.full((1,len(Phi)),Vol))
    Widarr = np.array(np.full((1,len(Phi)),Wid))
    Phinoco=np.array(Phinoco)
    
    # create bin data columns to be sent to a fits file later
    Colcountnocomp = fits.Column(name=('N_sources'),array=Phinoco,format='D')
    Colcountcomp=fits.Column(name=('N_sources_completeness'),array=Phi,format='D')
    ColVol=fits.Column(name=('Volume_Surveyed'),array=Volarr[0],format='D')
    ColWid=fits.Column(name=('Delta_Log_Mass'),array=Widarr[0],format='D')
    
    # obtain Phi error by taking sqrt of non completeness count values, and then applying mass function equation
    PhiError=np.array([math.sqrt(x) for x in Phinoco])
    PhiError= PhiError / Vol / Wid  
    
    
    Phi   = Phi / Vol / Wid   #Apply mass function equation to completeness values
    Phinoco   = Phinoco / Vol / Wid #Apply mass function equation to non-completeness values
    
    #Make more bin data columns to write to fits file
    Colphinocomp = fits.Column(name=('Phi_no_completeness'),array=Phinoco,format='D')
    Colphi=fits.Column(name=('Phi_completeness'),array=Phi,format='D')
    Colphierrnocomp = fits.Column(name=('Phi_Error_no_completeness'),array=PhiError,format='D')
    
    #Obtain Cofactor values (the factor between phi with compleness and without
    Cofa= np.array([x/y for x,y in zip(Phi,Phinoco)])
    #print(Cofa)
    MASK11=(str(Cofa)=='nan')
    Cofa[MASK11]= 0.
    ColCofa = fits.Column(name=('Completeness_Factor'),array=Cofa,format='D') #write to column   
    
    #Create completeness errors by combining cofactor and completenessless error
    CoPhiError= np.array([x*y for x,y in zip(PhiError,Cofa)])
    print(CoPhiError)
    Colphierr=fits.Column(name=('Phi_Error_completeness'),array=CoPhiError,format='D')
    ColMasx=fits.Column(name=('Logmass_center'),array=Masx,format='D')
    
    
    ##########################################
    #Apply MASK2 to output arrays to eliminate 0 values.
    MASK2=(Phi > 0.)
    PhiError = PhiError[MASK2]
    Masxa = Masx[MASK2]
    Phia= Phi[MASK2]
    Phinoco = Phinoco[MASK2]
    CoPhiError = CoPhiError[MASK2]
    #print(PhiError)
    ##########################################
    MASK5=(Phi==0.)
    #print(Phi)
    #print(Masx)
    Masx0=Masx[MASK5]
    #print(Masx0)
    Phi0= 1 / Vol / Wid
    Phi0 = np.array(np.full((1,len(Masx0)),Phi0))
    #print(Phi0)
    #################################
    #Create ArrayArray for all bin data columns
    
    ColWids.append(ColWid)
    ColVols.append(ColVol)
    Colcountnocomps.append(Colcountnocomp)
    Colcountcomps.append(Colcountcomp)
    Colphinocomps.append(Colphinocomp)
    Colphis.append(Colphi)
    Colphierrnocomps.append(Colphierrnocomp)
    ColCofas.append(ColCofa)
    Colphierrs.append(Colphierr)
    ColMasxs.append(ColMasx)
    
    
    
    
    
    
    # Create tensors for plotting data to be used in later for loop
    Phis.append(Phia)
    Phi0s.append(Phi0[0])
    PhiErrors.append(PhiError)
    CoPhiErrors.append(CoPhiError)
    Phinocos.append(Phinoco)
    Masxs.append(Masxa)
    Masx0s.append(Masx0)
###################################################

 #Create the Respective labels arrays
labels=[]
labels2=[]
New_unique_Zs=[]

for d in np.unique(zRed):
    
    New_unique_Zs.append('{:.2f}'.format(d))
    All_Zerrs.append(0.)
#print(New_unique_Zs)
for y in New_unique_Zs:
    z='z=' + y 
    a='z=' + y + ' (Without Completeness)'
    labels.append(z)
    labels2.append(a)
    
for y,x in zip(newZs,newzErrors):
    z='z=' + str(y) + '$\pm$' + str(x)
    a='z=' + str(y) + '$\pm$' + str(x) + ' (Without Completeness)'
    labels.append(z)
    labels2.append(a)    
########################################
#create colour arrays for use in plots
colors= cm.rainbow(np.linspace(0.5,1, len(Phis)))
colors2= cm.rainbow(np.linspace(0,0.5, len(Phis)))


#####################################################
All_Zs=np.append(New_unique_Zs,newZs)
All_Zerrs=All_Zerrs+newzErrors
LabMin=[labels,Mminsall]
pd.DataFrame(labels).to_csv("file.csv")
pd.DataFrame(LabMin).to_csv("file.csv")
plt.figure(figsize=(20,20))
for i,c,c2,B,No in zip(np.arange(0, len(labels),1),colors,colors2,All_Zs,np.arange(0, len(labels),1)):
    print(i)
    a = plt.subplot(6, 3,i+1)
    c='#FF0000'
    c2='#FFFF00'
    #call in bin data for respective band/filter
    ColVol=ColVols[i]
    ColWid=ColWids[i]
    Colcountnocomp = Colcountnocomps[i]
    Colcountcomp = Colcountcomps[i]
    Colphinocomp = Colphinocomps[i]
    Colphi = Colphis[i]
    Colphierrnocomp = Colphierrnocomps[i]
    ColCofa = ColCofas[i]
    Colphierr = Colphierrs[i]
    ColMasx = ColMasxs[i]
    #print(ColVol)
    t = fits.BinTableHDU.from_columns([ColVol,ColWid,Colcountnocomp,Colcountcomp,Colphinocomp,Colphi,Colphierrnocomp,ColCofa,Colphierr,ColMasx])
    t.writeto('Real_Bin/Bin_data_%s.fits'%B,overwrite=True)
    
    Ms=Masxs[i]
    #print(Ms)
    
    #print(Masx)
    M0s=Masx0s[i]
    #print(M0s)
    Ps=Phis[i]
    P0s=Phi0s[i]
    PErrs=PhiErrors[i] 
    coPErrs=CoPhiErrors[i]            
    Ls=labels[i]
    Ls2=labels2[i]
    Phinc=Phinocos[i]
    Mmin=Mminsall[i]
    
    check=np.array([x-y for x,y in zip(Ps,coPErrs)])
    #MASK4=(check==0.)
    #print(check)
    #print(check)
    #print(Ps)
    #print(coPErrs)
    #print(Ps-coPErrs)
    logP=np.array([math.log(x,10) for x in Ps])
    for i,j in zip(check,np.arange(0,len(check),1)):
        if i <= 1e-15:
            check[j]=10**(logP[j]-1)
            
    #print(check)
    logPhinc=[math.log(x,10) for x in Phinc]
    cologPp=[math.log(x+y,10) for x,y in zip(Ps,coPErrs)]
    cologPm=[math.log(x,10) for x in check]
    minco=np.array([x-y for x,y in zip(logP,cologPm)])
    #print(check)
    Plusco=np.array([x-y for x,y in zip(cologPp,logP)])
    coasy_err=[minco,Plusco]
    #print(coasy_err)
    check=np.array([x-y for x,y in zip(Phinc,PErrs)])
    MASK4=(check==0.)
    
    logPp=[math.log(x+y,10) for x,y in zip(Phinc,PErrs)]
    for i,j in zip(check,np.arange(0,len(check),1)):
        if i <= 1e-15:
            check[j]=10**(logPhinc[j]-1)
    logPm=[math.log(x,10) for x in check]
    minnoco=[x-y for x,y in zip(logPhinc,logPm)]
    Plusnoco=[x-y for x,y in zip(logPp,logPhinc)]
    asy_err=[minnoco,Plusnoco]
    
    
    
    logPhi0=np.array([math.log(x,10) for x in P0s])
    #print(logPhi0)
    
    
    #np.concatenate(np.array(logP),np.array(logPhi0))
    #np.append(Ms,M0s)
    #print(Ms)
    #Output plot data for each band/filter
    ColLogphinocomp = fits.Column(name=('Log_Phi_no_completeness'),array=logPhinc,format='D')
    ColLogphi=fits.Column(name=('Log_Phi_completeness'),array=logP,format='D')
    ColLogphierrnocompup = fits.Column(name=('Log_Phi_Error_no_completeness_plus'),array=Plusnoco,format='D')
    ColLogphierrnocompdown = fits.Column(name=('Log_Phi_Error_no_completeness_minus'),array=minnoco,format='D')
    ColLogphierrdown = fits.Column(name=('Log_Phi_Error_completeness_minus'),array=minco,format='D')
    ColLogphierrup = fits.Column(name=('Log_Phi_Error_completeness_plus'),array=Plusco,format='D')
    ColMass = fits.Column(name=('Log_Mass_Center'),array=Ms,format='D')
    t = fits.BinTableHDU.from_columns([ColLogphinocomp,ColLogphierrnocompup,ColLogphierrnocompdown,ColLogphi,ColLogphierrup,ColLogphierrdown,ColMass])
    t.writeto('Real_Plot/Plot_data%s.fits'%B,overwrite=True)
    
    
    # Create 'invisible' points at low values, with error bars that extend to the upper limit so the schechter fit will account for the upper limit when it fits
    fit_logPhi0=logPhi0-4
    fit_logPhi0_dow_err=np.array([0]*len(fit_logPhi0))
    fit_logPhi0_up_err=np.array([4]*len(fit_logPhi0))
    fit_logPhi0_errs=[fit_logPhi0_dow_err,fit_logPhi0_up_err]
    coerr=np.array((Plusco+minco)/2)
    All_errs=np.append(coerr,fit_logPhi0_up_err)
    #print(np.append(logP,fit_logPhi0))
    All_Ms=np.append(Ms,M0s)
    #print(All_Ms)
    fit_Ms=np.unique(All_Ms)
    MASK10=(fit_Ms>=Mmin)
    fit_Ms=fit_Ms[MASK10]
    #print(fit_Ms)
    All_logPs=np.append(logP,fit_logPhi0)
    fit_logPs=[]
    fit_errs=[]
    for m in fit_Ms:
        MASK6=(All_Ms==m)
        add_logP=All_logPs[MASK6]
        add_errs=All_errs[MASK6]
        fit_logPs.append(add_logP[0])
        fit_errs.append(add_errs[0])
    #print(All_logPs)
    #print(fit_logPs)
    #print(fit_errs)
    #print(fit_Ms)
    #print(fit_logPs)
    #print(fit_errs)
    #Bound=np.array([[-np.inf, -np.inf, -1.5],[np.inf,np.inf,-1.3]])
    
        
    popt,pcov = optimize.curve_fit(schfunct,fit_Ms,fit_logPs,sigma=fit_errs,p0=[-4.0,10.5])
    #popt,pcov = optimize.curve_fit(schfunct,Ms,logP,sigma=coerr,p0=[-4.0,10.5])
    #print(pcov)
    
    perr= np.sqrt(np.diag(pcov))
    if perr[1] < 1.:
        Param1s.append(popt[0])
        Param2s.append(popt[1])
        Param1errs.append(perr[0])
        Param2errs.append(perr[1])
        
        #print fit parameters and 1-sigma estimates
        print('fit parameter 1-sigma error')
        print('———————————–')
        for i in range(len(popt)):
            print(str(popt[i])+' +- '+str(perr[i]))
    
        # prepare confidence level curves
        nstd = 1. # to draw 1-sigma intervals
        nstds = [1] # to draw 3-sigma intervals
        
        Colours=['#247AB6','#B62724','#39B624']
        
        txt2='%s'%Ls
        txt3=r'$\phi_0 = $' + '{:.2f}'.format(popt[0]) + r'$\pm$' + '{:.2f}'.format(perr[0])
        #fig, a = plt.subplots(1)
        a.set_ylim([-7,-3])
        for n,c3 in zip(nstds,Colours):
            
            popt_up = popt + n * perr
            popt_dw = popt - n * perr
            
        
            #####################################
            # create fit y values using original masx values, in order to reach mass values that were previously only covered by limits, which are not currently account for
            fit = schfunct(fit_Ms, *popt)
            fit_up = schfunct(fit_Ms, *popt_up)
            fit_dw = schfunct(fit_Ms, *popt_dw)
            #fig.text(.5, .0, txt, ha='center')
            #fig.text(0.18,0.4,s=txt2, horizontalalignment='right',verticalalignment='top', transform=ax.transAxes, fontsize='xx-large')
            #fig.text(0.38,0.32,s=txt3, horizontalalignment='right',verticalalignment='top', transform=ax.transAxes, fontsize='xx-large')
            plt.plot(fit_Ms, fit_up, label=r'%s $\sigma$' %n,color=c3, linewidth=0.7,zorder=-2)
            plt.plot(fit_Ms, fit_dw, label='_nolegend_',color=c3, linewidth=0.7,zorder=-2)
            #labelLines(plt.gca().get_lines(),zorder=2.5)
        a.fill_between(fit_Ms, fit_up, fit_dw, alpha=.35, label='_nolegend_' ,color='#72B8D5',zorder=-2)
        #plt.show()    
        plt.plot(fit_Ms, fit, label='Schechter Fit',zorder=-2)    
        plt.errorbar(Ms,logP,yerr=coasy_err,fmt='o',markeredgecolor='#010101',color=c,label=Ls,alpha=1,zorder=-1)
        plt.errorbar(Ms,logPhinc,yerr=asy_err,fmt='o',color=c2,label=Ls2,alpha=1,zorder=-1,markeredgecolor='#010101')
        
        plt.plot([Mmin,Mmin+0.000001],[-10,10],'k--')
        plt.errorbar(M0s,logPhi0,fmt='v',color=c,label='_nolegend_',alpha=1,markersize=10,zorder=-1,markeredgecolor='#010101')
        
        # Plot the 'invisible' points to check they are as desired
        #plt.errorbar(M0s,fit_logPhi0,fmt='*',yerr=fit_logPhi0_errs,color=c,label='_nolegend_',alpha=1,markersize=10,zorder=-1,markeredgecolor='#010101')
        
        
        
        #plt.xlabel(r'$\log_{10}(M_\star\,/\,M_\odot)$')
        #plt.ylabel(r'$\log_{10}(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
    
        #plt.xticks(np.arange(9,11.5,0.5),fontsize=14)
        plt.minorticks_on()
        #plt.yticks(fontsize=14)
        #plt.ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
        #plt.ylim(1.e-10,1.e-5)
        #plt.autoscale(enable=True,axis='y',tight=True)
        plt.legend()
        plt.tight_layout()
        #a.tick_params(which='both',direction='in',top=True,right=True,labelbottom=False,labelleft=False)
        a.tick_params(which='both',direction='in',top=True,right=True)
        a.set_yticks(np.arange(-6.5,-3.0,0.5))
        if (No == 0 or No == 3 or No == 6 or No == 9 or No ==12):
            #a.set_yticklabels([])
            #a.set_xticklabels([])
            #plt.xlabel(r'$\log_{10}(M_\star\,/\,M_\odot)$')
            plt.ylabel(r'$\log_{10}(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')     
        if (No == 15):
            plt.xlabel(r'$\log_{10}(M_\star\,/\,M_\odot)$')
            plt.ylabel(r'$\log_{10}(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
        if (No==16 or No==17):
            plt.xlabel(r'$\log_{10}(M_\star\,/\,M_\odot)$')
            plt.ylabel(r'$\log_{10}(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
            a.set_yticklabels([])
        if (No==1 or No==2 or No==4 or No==5 or No==7 or No==8 or No==10 or No==11 or No==13 or No==14):
            a.set_yticklabels([])
            a.set_xticklabels([])
            #plt.xlabel(r'$\log_{10}(M_\star\,/\,M_\odot)$') 
             #          {'color' : 'k','fontsize' : 30 })
            
            #plt.ylabel(r'log$_{10}$\,($\Phi$/Mpc$^{-3}$)\,dlogL$^{-1}$', {'color' : 'k','fontsize' : 26 })
                       
        #elif(i == 10 or i == 11):
            #a.set_yticks([])
            
           # plt.xlabel(r'log$_{10}$\,(L$_{\rm Ly\alpha}$/erg\,s$^{-1}$)', 
                    #   {'color' : 'k','fontsize' : 30 })
        #else:
           # a.set_xticklabels([])
           # a.set_yticklabels([])
        #a.tick_params(which='both',direction='in',top=True,right=True)
        #plt.show()
        #plt.savefig('RealSMFs/ZBSMF%s_1sig_General.png'%B)
        #plt.clf()
    else:
        popt,pcov = optimize.curve_fit(schfunct1,fit_Ms,fit_logPs,sigma=fit_errs,p0=[-4.0])
    #popt,pcov = optimize.curve_fit(schfunct,Ms,logP,sigma=coerr,p0=[-4.0,10.5])
    #print(pcov)
    
        perr= np.sqrt(np.diag(pcov))
        Param1s.append(popt[0])
        Param2s.append(logM_star)
        Param1errs.append(perr[0])
        Param2errs.append(0)
        
        #print fit parameters and 1-sigma estimates
        print('fit parameter 1-sigma error')
        print('———————————–')
        for i in range(len(popt)):
            print(str(popt[i])+' +- '+str(perr[i]))
    
        # prepare confidence level curves
        nstd = 1. # to draw 1-sigma intervals
        nstds = [1] # to draw 3-sigma intervals
        
        Colours=['#247AB6','#B62724','#39B624']
        
        txt2='%s'%Ls
        txt3=r'$\phi_0 = $' + '{:.2f}'.format(popt[0]) + r'$\pm$' + '{:.2f}'.format(perr[0])
        #fig, ax = plt.subplots(1)
        a.set_ylim([-7,-3])
        for n,c3 in zip(nstds,Colours):
            
            popt_up = popt + n * perr
            popt_dw = popt - n * perr
            
        
            #####################################
            # create fit y values using original masx values, in order to reach mass values that were previously only covered by limits, which are not currently account for
            fit = schfunct1(fit_Ms, *popt)
            fit_up = schfunct1(fit_Ms, *popt_up)
            fit_dw = schfunct1(fit_Ms, *popt_dw)
            #fig.text(.5, .0, txt, ha='center')
            #fig.text(0.18,0.4,s=txt2, horizontalalignment='right',verticalalignment='top', transform=ax.transAxes, fontsize='xx-large')
            #fig.text(0.38,0.32,s=txt3, horizontalalignment='right',verticalalignment='top', transform=ax.transAxes, fontsize='xx-large')
            plt.plot(fit_Ms, fit_up, label=r'%s $\sigma$' %n,color=c3, linewidth=0.7,zorder=-2)
            plt.plot(fit_Ms, fit_dw, label='_nolegend_',color=c3, linewidth=0.7,zorder=-2)
            #labelLines(plt.gca().get_lines(),zorder=2.5)
        a.fill_between(fit_Ms, fit_up, fit_dw, alpha=.35, label='_nolegend_' ,color='#72B8D5',zorder=-2)
        #plt.show()    
        plt.plot(fit_Ms, fit, label='Schechter Fit',zorder=-2)    
        plt.errorbar(Ms,logP,yerr=coasy_err,fmt='o',markeredgecolor='#010101',color=c,label=Ls,alpha=1,zorder=-1)
        plt.errorbar(Ms,logPhinc,yerr=asy_err,fmt='o',color=c2,label=Ls2,alpha=1,zorder=-1,markeredgecolor='#010101')
        
        plt.plot([Mmin,Mmin+0.000001],[-10,10],'k--')
        plt.errorbar(M0s,logPhi0,fmt='v',color=c,label='_nolegend_',alpha=1,markersize=10,zorder=-1,markeredgecolor='#010101')
        
        # Plot the 'invisible' points to check they are as desired
        #plt.errorbar(M0s,fit_logPhi0,fmt='*',yerr=fit_logPhi0_errs,color=c,label='_nolegend_',alpha=1,markersize=10,zorder=-1,markeredgecolor='#010101')
        
        
        
       # plt.xlabel(r'$\log_{10}(M_\star\,/\,M_\odot)$')
       # plt.ylabel(r'$\log_{10}(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
    
        #plt.xticks(np.arange(9,11.5,0.5),fontsize=14)
        plt.minorticks_on()
        #plt.yticks(fontsize=14)
        #plt.ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
        #plt.ylim(1.e-10,1.e-5)
        #plt.autoscale(enable=True,axis='y',tight=True)
        plt.legend()
        plt.tight_layout()
        #plt.show()
        #plt.savefig('RealSMFs/ZBSMF%s_1sig_General.png'%B)
        #plt.clf()
        a.tick_params(which='both',direction='in',top=True,right=True)
        a.set_yticks(np.arange(-6.5,-3.0,0.5))
        if (No == 0 or No == 3 or No == 6 or No == 9 or No ==12):
            #a.set_yticklabels([])
            #a.set_xticklabels([])
            #plt.xlabel(r'$\log_{10}(M_\star\,/\,M_\odot)$')
            plt.ylabel(r'$\log_{10}(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
            #plt.yticks(-6.5,-3.5,0.5)     
        if (No == 15):
            plt.xlabel(r'$\log_{10}(M_\star\,/\,M_\odot)$')
            plt.ylabel(r'$\log_{10}(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
        if (No==16 or No==17):
            plt.xlabel(r'$\log_{10}(M_\star\,/\,M_\odot)$')
            #plt.ylabel(r'$\log_{10}(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
            a.set_yticklabels([])
        if (No==1 or No==2 or No==4 or No==5 or No==7 or No==8 or No==10 or No==11 or No==13 or No==14):
            a.set_yticklabels([])
            a.set_xticklabels([])
        #if i == 9:
            #plt.xlabel(r'log$_{10}$\,(L$_{\rm Ly\alpha}$/erg\,s$^{-1}$)', 
             #          {'color' : 'k','fontsize' : 30 })
            
            #plt.ylabel(r'log$_{10}$\,($\Phi$/Mpc$^{-3}$)\,dlogL$^{-1}$', {'color' : 'k','fontsize' : 26 })
                       
        #elif(i == 10 or i == 11):
          #  a.set_yticks([])
            
           # plt.xlabel(r'log$_{10}$\,(L$_{\rm Ly\alpha}$/erg\,s$^{-1}$)', 
                    #   {'color' : 'k','fontsize' : 30 })
        #else:
           # a.set_xticklabels([])
           # a.set_yticklabels([])
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('Gridtest.png')    
##########################################################################
# designate and output all the relevant schechter fit data for use in the stellar mass density later.
ColRed = fits.Column(name=('Redshift'),array=All_Zs,format='D')
ColRederrs = fits.Column(name=('Redshift_Range'),array=All_Zerrs,format='D')
ColLab = fits.Column(name=('Labels'),array=labels,format='20A')
Colphipara = fits.Column(name=('LogPhi'),array=Param1s,format='D')
CollogMpara = fits.Column(name=('LogM_star'),array=Param2s,format='D')
Colphiparaerr = fits.Column(name=('logPhierr'),array=Param1errs,format='D')
CollogMparaerr = fits.Column(name=('logM_starerr'),array=Param2errs,format='D')
ColZlen=fits.Column(name=('No.of Sources'),array=zlen,format='D')
t = fits.BinTableHDU.from_columns([ColRed,ColRederrs,ColLab,Colphipara,Colphiparaerr,CollogMpara,CollogMparaerr,ColZlen])
t.writeto('Real_Para/Para_data_General.fits',overwrite=True)
