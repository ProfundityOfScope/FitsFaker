import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.time import Time
import pandas as pd
import os

def gaussian( x, y, xm, ym, amp, std):
    return amp*np.exp(-(x-xm)**2/(2*std**2)) * np.exp(-(y-ym)**2/(2*std**2))

colnames = ['mjd','blank1','blank2','ra','dec','apmag','sbrt','delta',
            'deltadot','SOT','symbol','STO']
system = {}
for file in os.listdir(os.getcwd()):
    if '.txt' in file:
        body = file.replace('.txt','')
        data = pd.read_csv(file, names=colnames, header=None, index_col=False, 
                           usecols=[0,3,4,5,6])
        data.columns = ['mjd'] + [ body+'_'+x for x in data.columns[1:]]
        system[body] = data
    

system_df = pd.DataFrame({'mjd':system['jupiter']['mjd']} )
for body in system.keys():
    system_df = pd.merge(system_df, system[body], on='mjd')

##################################

rinds = np.random.randint(0,len(system_df),40)

for rind in rinds:
    
    date_mjd = system_df['mjd'][rind]
    tod = date_mjd%1
    if 0.25<tod<0.75:
        continue
        
    
    xx = np.linspace(-0.2,0.2,1000)
    XX,YY = np.meshgrid(xx,xx)
    
    gauss = lambda amp, x, y : amp*np.exp(-x**2/10**2)*np.exp(-y**2/10**2)
    
    ZZ = np.random.uniform(1,0.1,XX.shape)
    
    fj = 1e4
    mj = np.mean(system_df['jupiter_apmag'])
    
    plt.figure(figsize=(6,6))
    plt.scatter(0,0,label='jupiter')
    for moon in ['io', 'europa', 'ganymede', 'callisto']:
        xm = system_df[moon+'_ra'][rind] - system_df['jupiter_ra'][rind]
        ym = system_df[moon+'_dec'][rind] - system_df['jupiter_dec'][rind]
        
        plt.scatter(xm,ym,label=moon)
        
        mm = system_df[moon+'_apmag'][rind]
        fm = fj*10**(-0.4*(mm-mj))
        ZZ += gaussian(XX,YY,xm,ym,fm,0.001)
    plt.legend()
    plt.title(system_df['mjd'][rind])    
    plt.xlim(-0.2,0.2)
    plt.ylim(-0.2,0.2)
    plt.savefig(str(rind)+'.pdf')
    
    for i in range(len(xx)):
        for j in range(len(xx)):
            r = np.sqrt(xx[i]**2+xx[j]**2)
            if r < 0.015:
                ZZ[i,j] = 10
            else:
                ZZ[i,j] += gaussian(xx[i],xx[j],0,0,20,0.0127)
    
    
    plt.figure(figsize=(6,6))
    plt.imshow(ZZ, cmap='Greys')
    plt.show()
#    
#    hdr = fits.Header()
#    hdr['OBSERVER'] = 'Seth Bruzewski'
#    hdr['DATE'] = round(Time.now().mjd,3)
#    hdr['DATE-OBS'] = date_mjd
#    hdr['EQUINOX'] = 2000.00
#    hdr['TELESCOP'] = 'Greenwich Totally Real Telescope (GTRT)'
#    hdr['FOC-LEN'] = '5000 mm'
#    hdr['PIX-SIZE'] = '35 um'
#    hdr['INSTRUME'] = 'Greenwich Totally Real Observatory (GTRO)'
#    hdr['OBJECT'] = 'Jupiter'
#    hdu = fits.PrimaryHDU(ZZ, header=hdr)
#    name = 'SB18_{:03d}.fits'.format(rind)
#    hdu.writeto(name, overwrite=True)
