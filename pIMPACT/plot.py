# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:48:25 2016

@author: kilean
"""

import numpy as np
import matplotlib.pyplot as plt
import os

#%%############################################################################
##            internel class and function to layout lattice plot             ##
###############################################################################
class _element:
    def __init__(self, lattice, elem_id):
        self.number=0
        self.elem_id=elem_id
        
        for i in range(len(lattice)):
            if lattice[i,1]==self.elem_id:
                self.number=self.number+1;
        self.positions=np.zeros([4*self.number+2,2])
        count=0
        z=0.0
        for i in range(0,len(lattice)):
            if lattice[i,1]==self.elem_id:
                self.positions[4*count+1,0]=z
                self.positions[4*count+2,0]=z+1E-15;
                self.positions[4*count+3,0]=z+lattice[i,0]-1E-15;
                self.positions[4*count+4,0]=z+lattice[i,0];
                if lattice[i,2]>0:
                    self.positions[4*count+2,1]=1;
                    self.positions[4*count+3,1]=1;
                else :
                    self.positions[4*count+2,1]=-1;
                    self.positions[4*count+3,1]=-1;                    
                count=count+1
            z=z+lattice[i,0]
        self.positions[0,0]=0.0
        self.positions[-1,0]=z               

             
def _lattice(h, lattice_info, fileDir='',):
    
    lattice_offset = lattice_info['offset']
    lattice_scale = lattice_info['scale']
    if 'file' in lattice_info:
        lattice_file = lattice_info['file']
    else:
        lattice_file = 'test.in'
    
    lattice=np.loadtxt(fileDir+lattice_file,
                       comments='!',
                       skiprows=12,
                       usecols=(0,3,4))
                                                  
    RF=_element(lattice,104)
    Bend=_element(lattice,4)
    Quad=_element(lattice,1)
    
    h.fill_between(RF.positions[:,0], lattice_offset, 
                   lattice_offset+lattice_scale*RF.positions[:,1], alpha=0.2, 
                   facecolor='black', color="none", hatch="X")
    h.fill_between(RF.positions[:,0], lattice_offset, 
                   lattice_offset-lattice_scale*RF.positions[:,1], alpha=0.2, 
                   facecolor='black', color="none", hatch="X")
    h.fill_between(Quad.positions[:,0], lattice_offset, 
                   lattice_offset+lattice_scale*Quad.positions[:,1], alpha=0.6, 
                   facecolor='black', linewidth=0.2)
    h.fill_between(Bend.positions[:,0], lattice_offset, 
                   lattice_offset+lattice_scale*Bend.positions[:,1], alpha=0.1, 
                   facecolor='black', color="none", hatch="o")
                   
def _getLatticeInfo(lattice_info=None, minX=0.0, maxX=1.0):
    if lattice_info==None:
        return {'file':'test.in','offset':None, 'scale':None}
  
    lattice_info['offset'] = 0.5*(minX+maxX)
    lattice_info['scale'] = 0.2*( 0.5*(maxX-minX) )
    
    
#%%############################################################################
###############################################################################
###                                  RMS plots                              ###
###############################################################################
############################################################################### 

#%%#####  internel functions  ######
# xy envelope
def _x_envelope(h, fileDir='',flag_maxRadius=False,
                plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    X=np.loadtxt(fileDir+'fort.24')
    X_centroid=X[:,1]*1E3
    X_rms=X[:,2]*1E3
    
    Y=np.loadtxt(fileDir+'fort.25')
    Y_centroid=Y[:,1]*1E3
    Y_rms=Y[:,2]*1E3
    
    S=X[:,0]
    s_max=max(S)
    s_min=min(S) 
    
    h.plot(S,X_centroid+X_rms,'b',S,X_centroid-X_rms,'b')
    h.plot(S,Y_centroid+Y_rms,'r',S,Y_centroid-Y_rms,'r')
    h.set_ylabel(r'$\mathsf{\sigma_{x,y} \, (mm)}$', fontsize=12)    
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min( min(X_centroid -X_rms),  min(Y_centroid -Y_rms) )
            maxX = max( max(X_centroid +X_rms),  max(Y_centroid +Y_rms) )
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 
    
    if flag_maxRadius:
        Rmax=1E3*np.loadtxt(fileDir+'fort.18')[:,-1]
        h.plot(S,-Rmax,'g',S,Rmax,'g')
        
    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
        if 'x' in plotRange:        
            h.set_ylim(plotRange['x'])

# p_xy envelope
def _px_envelope(h, fileDir='', plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    X=np.loadtxt(fileDir+'fort.24')
    Y=np.loadtxt(fileDir+'fort.25')
    
    Px_centroid=X[:,3]*1E3
    Px_rms=X[:,4]*1E3
    Py_centroid=Y[:,3]*1E3
    Py_rms=Y[:,4]*1E3
    
    S=X[:,0]
    s_max=max(S)
    s_min=min(S)   
    
    h.plot(S,Px_centroid+Px_rms,'b',S,Px_centroid-Px_rms,'b')
    h.plot(S,Py_centroid+Py_rms,'r',S,Py_centroid-Py_rms,'r')
    h.set_ylabel(r'$\mathsf{\sigma_{px,py}\, (mrad)}$', fontsize=12)    
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min( min(Px_centroid -Px_rms),  min(Py_centroid -Py_rms) )
            maxX = max( max(Px_centroid +Px_rms),  max(Py_centroid +Py_rms) )
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 
    
    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
        if 'px' in plotRange:        
            h.set_ylim(plotRange['px'])

def _x_emittance(h, halo=None, fileDir='', plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    X=np.loadtxt(fileDir+'fort.24')
    X_emittance=X[:,6]*1E6
    
    Y=np.loadtxt(fileDir+'fort.25')
    Y_emittance=Y[:,6]*1E6
    
    S=X[:,0]
    s_max=max(S)
    s_min=min(S)   
    
    h.plot(S,X_emittance,'b')
    h.plot(S,Y_emittance,'r')
    h.set_ylabel(r'$\mathsf{norm.\, \epsilon \, (mm-mrad)}$', fontsize=11)   
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min( min(X_emittance),  min(Y_emittance) )
            maxX = max( max(X_emittance),  max(Y_emittance) )
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir)    
    
    if isinstance(halo,int) and X.shape[1] == 10:
        h2 = h.twinx()
        if halo==99:
            X_halo=X[:,9]*1E6
            Y_halo=Y[:,9]*1E6
        elif halo==90:
            X_halo=X[:,7]*1E6
            Y_halo=Y[:,7]*1E6
        else :
            X_halo=X[:,8]*1E6
            Y_halo=Y[:,8]*1E6
        
        h2.plot(S,1-X_halo/X_emittance,'b--',alpha=0.5)
        h2.plot(S,1-Y_halo/Y_emittance,'r--',alpha=0.5)
        h2.set_ylabel(r'$\mathsf{ \Delta \epsilon / \epsilon}$', fontsize=12)           
    elif isinstance(halo,int):
        print '90,95,99% rms emittance is not calculate. IMPACT standard output flagg must be set to 2'

    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
            h2.set_xlim(plotRange['s'])  
        if 'ex' in plotRange:        
            h.set_ylim(plotRange['ex'])  
    

#%% Plot longitudinal RMS statistics
def _z_envelope(h, fileDir='', plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    Z=np.loadtxt(fileDir+'fort.26')
    Z_centroid=Z[:,1]
    Z_rms=Z[:,2]
    
    S=Z[:,0]
    s_max=max(S)
    s_min=min(S) 
    
    h.plot(S,Z_centroid+Z_rms,'b',S,Z_centroid-Z_rms,'b')
    h.set_ylabel(r'$\mathsf{\sigma_z \, (deg)}$', fontsize=12)
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min(Z_centroid -Z_rms)
            maxX = max(Z_centroid +Z_rms)
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 

    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
        if 'z' in plotRange:        
            h.set_ylim(plotRange['z'])
            
    
def _pz_envelope(h, fileDir='', plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    Z=np.loadtxt(fileDir+'fort.26')
    Pz_centroid=Z[:,3]
    Pz_rms=Z[:,4]
    S=Z[:,0]
    s_max=max(S)
    s_min=min(S)   
    
    h.plot(S,Pz_centroid+Pz_rms,'b',S,Pz_centroid-Pz_rms,'b')
    h.set_ylabel(r'$\mathsf{\sigma_E \, (MeV)}$', fontsize=12)
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min(Pz_centroid -Pz_rms)
            maxX = max(Pz_centroid +Pz_rms)
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 
    
    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
        if 'px' in plotRange:        
            h.set_ylim(plotRange['px'])         

     
def _z_emittance(h, halo=None, fileDir='', plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    Z=np.loadtxt(fileDir+'fort.26')
    Z_emittance=Z[:,6]   
    S=Z[:,0]
    s_max=max(S)
    s_min=min(S)   
    
    h.plot(S,Z_emittance,'b')
    h.set_ylabel(r'$\mathsf{\epsilon \, (deg-MeV)}$', fontsize=11)
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min(Z_emittance)
            maxX = max(Z_emittance)
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 
     
    
    if isinstance(halo,int) and Z.shape[1] == 10:
        h2 = h.twinx()
        if halo==99:
            Z_halo=Z[:,9]
        elif halo==90:
            Z_halo=Z[:,7]
        else :
            Z_halo=Z[:,8]           
        h2.plot(S,1-Z_halo/Z_emittance,'b--',alpha=0.5)
        h2.set_ylabel(r'$\mathsf{ \Delta \epsilon / \epsilon}$', fontsize=12)  
    elif isinstance(halo,int):
        print '90,95,99% rms emittance is not calculate. IMPACT standard output flagg must be set to 2' 
        
    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
            h2.set_xlim(plotRange['s'])  
        if 'ez' in plotRange:        
            h.set_ylim(plotRange['ez'])  

    
    
#%% Plot RMS statistics
def _getLatticeInfo_4rmsPlot(x=True, px=True, ex=False, z=True, pz=True, ez=False):
    lattice_info = {}
    if x:
        lattice_info['x'] = {'file':'test.in','offset':None, 'scale':None}
    if px:
        lattice_info['px'] = {'file':'test.in','offset':None, 'scale':None}
    if ex:
        lattice_info['ex'] = {'file':'test.in','offset':None, 'scale':None}
    if z:
        lattice_info['z'] = {'file':'test.in','offset':None, 'scale':None}
    if pz:
        lattice_info['pz'] = {'file':'test.in','offset':None, 'scale':None}
    if ez:
        lattice_info['ez'] = {'file':'test.in','offset':None, 'scale':None}
    return lattice_info    
    
def rms(savefileID=0, fileDir='', lattice_info=_getLatticeInfo_4rmsPlot(), flag_maxRadius=False,
        plot_momentum=True, halo=None, plotRange=None):
            
    plt.rcParams.update({'font.size': 9})  
    # transverse
    fig = plt.figure() 
    
    if plot_momentum:
        plotid = 311
    else:
        plotid = 211
    sub1=fig.add_subplot(plotid)
    if 'x' in lattice_info:
        info = lattice_info['x']
    else :
        info = None
    _x_envelope(sub1,fileDir=fileDir, flag_maxRadius=flag_maxRadius,
                 plotRange=plotRange, lattice_info=info)


    if plot_momentum:
        plotid=plotid+1
        if 'px' in lattice_info:
            info = lattice_info['px']
        else :
            info = None
        sub2=fig.add_subplot(plotid)
        _px_envelope(sub2,fileDir=fileDir, plotRange=plotRange, lattice_info=info)


    plotid=plotid+1
    sub3=fig.add_subplot(plotid)
    if 'ex' in lattice_info:
        info = lattice_info['ex']
    else :
        info = None
    _x_emittance(sub3,halo=halo,fileDir=fileDir, plotRange=plotRange, lattice_info=info)
    sub3.set_xlabel('z (m)')
    
    
    plt.tight_layout() 
    plt.savefig(fileDir+'x'+str(savefileID)+'.png', dpi=240)

    
    # longitudinal
    fig = plt.figure() 
    if plot_momentum:
        plotid = 311
    else:
        plotid = 211
    sub1=fig.add_subplot(plotid)
    if 'z' in lattice_info:
        info = lattice_info['z']
    else :
        info = None
    _z_envelope(sub1,fileDir=fileDir, 
                 plotRange=plotRange, lattice_info=info)


    if plot_momentum:
        plotid=plotid+1
        sub2=fig.add_subplot(plotid)
        if 'pz' in lattice_info:
            info = lattice_info['pz']
        else :
            info = None
        _pz_envelope(sub2,fileDir=fileDir, plotRange=plotRange, lattice_info=info)


    plotid=plotid+1
    sub3=fig.add_subplot(plotid)
    if 'ez' in lattice_info:
        info = lattice_info['ez']
    else :
        info = None    
    _z_emittance(sub3,halo=halo,fileDir=fileDir, plotRange=plotRange, lattice_info=info)
    sub3.set_xlabel('z (m)')


    plt.tight_layout() 
    plt.savefig(fileDir+'z'+str(savefileID)+'.png', dpi=240)


#%% Plot max Amplitude
def _getLatticeInfo_4maxPlot(x=True, px=True, z=True, pz=True):
    lattice_info = {}
    if x:
        lattice_info['x'] = {'file':'test.in','offset':None, 'scale':None}
    if px:
        lattice_info['px'] = {'file':'test.in','offset':None, 'scale':None}
    if z:
        lattice_info['z'] = {'file':'test.in','offset':None, 'scale':None}
    if pz:
        lattice_info['pz'] = {'file':'test.in','offset':None, 'scale':None}
    return lattice_info    

def maxAmplitude(savefileID=0, fileDir='', 
        plot_momentum=True, flag_radius=True,
        plotRange=None, lattice_info=_getLatticeInfo_4maxPlot() ):
    
    plt.rcParams.update({'font.size': 9})            

    Data=np.loadtxt(fileDir+'fort.27')*1E3
    S  = Data[:,0]*1E-3
    X  = Data[:,1]
    Px = Data[:,2]
    Y  = Data[:,3]
    Py = Data[:,4]
    R  = np.loadtxt(fileDir+'fort.18')[:,-1]*1E3
    Z  = Data[:,5]*1E-3
    Pz = Data[:,6] 
    s_max=max(S)
    s_min=min(S)   
    
    plotid=111
    fig = plt.figure() 
    sub1=fig.add_subplot(plotid)
    if flag_radius:
        sub1.plot(S,R,'b')
        sub1.set_ylabel(r'$\mathsf{ R_{max} \, (mm)}$', fontsize=12)
    else:
        sub1.plot(S,X,'b',S,Y,'b--',alpha=0.7)
        sub1.set_ylabel(r'$\mathsf{ x,y_{max} \, (mm)}$', fontsize=12)
    sub1.set_xlabel('z (m)')
        
    if lattice_info['x'] != None:
        if lattice_info['x']['scale']==None:
            minX = min(min(X),min(Y)) 
            maxX = max(max(X),max(Y)) 
            _getLatticeInfo(lattice_info['x'], minX, maxX)
        _lattice(sub1, lattice_info['x'], fileDir=fileDir) 
        
    if plotRange==None:
        sub1.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            sub1.set_xlim(plotRange['s'])        
        if 'x' in plotRange:        
            sub1.set_ylim(plotRange['x'])        
    
    sub2 = sub1.twinx()  
    sub2.plot(S,Z,'r')
    sub2.set_ylabel(r'$\mathsf{ z_{max} \, (deg)}$', fontsize=12)       
    if plotRange !=None and 'z' in plotRange:       
        sub2.set_ylim(plotRange['z'])    
        
    plt.tight_layout() 
    plt.savefig(fileDir+'max_radius'+str(savefileID)+'.png', dpi=240)



    fig = plt.figure() 
    sub1=fig.add_subplot(plotid)
    sub1.plot(S,Pz,'b')
    sub1.set_ylabel(r'$\mathsf{ \Delta E_{max} \, (MeV)}$', fontsize=12)
    sub1.set_xlabel('z (m)')
    
    if lattice_info['pz'] != None:
        if lattice_info['pz']['scale']==None:
            minX = min(Pz)
            maxX = max(Pz)  
            _getLatticeInfo(lattice_info['pz'], minX, maxX)
            _lattice(sub1, lattice_info['pz'], fileDir=fileDir) 
            
    if plotRange==None:
        sub1.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            sub2.set_xlim(plotRange['s'])        
        if 'pz' in plotRange:        
            sub2.set_ylim(plotRange['pz'])
 
    if plot_momentum:
        sub2 = sub1.twinx()  
        sub2.plot(S,Px,'r',S,Py,'r--',alpha=0.7)
        sub2.set_ylabel(r'$\mathsf{ p_{x,y,max} \, (mrad)}$', fontsize=12)
        if plotRange !=None and 'px' in plotRange:
            sub2.set_ylim(plotRange['px'])        

    plt.tight_layout() 
    plt.savefig(fileDir+'max_energy'+str(savefileID)+'.png', dpi=240)
        

        



    
    
#%% Plot particle loss
from impact import readIMPACT
def powerLoss(savefileID=0, fileDir='', lattice_info=None, plotRange=None):
            
    
    dummy = np.loadtxt(fileDir+'fort.18')
    S = dummy [:,0]
    Energy = dummy[:,3]    
    Particles = np.loadtxt(fileDir+'fort.32')[:,1]
    n_total = Particles[0]
     
    beam,lattice=readIMPACT('test.in')
    PowerLoss = (1.0-Particles/n_total)*Energy*beam['current'] #assume all buckets are filled (Watt)
    
    dP = np.gradient(PowerLoss)
    dS = np.gradient(S)
    PowerLossPerMeter = dP/dS
    
    plt.rcParams.update({'font.size': 9})  
    fig = plt.figure() 
    sub1=fig.add_subplot(111)
    sub1.plot(S,PowerLossPerMeter,'r')
    ###
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min(PowerLossPerMeter)
            maxX = max(PowerLossPerMeter) 
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(sub1, lattice_info, fileDir=fileDir) 

    sub1.set_xlabel('z (m)')
    sub1.set_ylabel('power loss (W/m)')

    plt.tight_layout() 
    plt.savefig(fileDir+'loss'+str(savefileID)+'.png', dpi=240)

#%%############################################################################
###############################################################################
###                              phase-space plots                          ###
###############################################################################
############################################################################### 
#def phase_space(n,freq,ke,mass,saveDir='',z=None,plotRange=None):
#    
#    data=np.loadtxt('fort.'+str(n))
#    sigmax=np.std(data[:,0])
#    sigmapx=np.std(data[:,1])
#    mux=-np.mean(data[:,0]*data[:,1])/sigmax/sigmapx;
#    
#    sigmay=np.std(data[:,2])
#    sigmapy=np.std(data[:,3])
#    muy=-np.mean(data[:,2]*data[:,3])/sigmay/sigmapy;
#    
#    sigmaz=np.std(data[:,4])
#    sigmapz=np.std(data[:,5])
#    muz=-np.mean(data[:,4]*data[:,5])/sigmaz/sigmapz;
#    
#    data[:,6]=data[:,0]*data[:,0]/(sigmax*sigmax) \
#              +data[:,1]*data[:,1]/(sigmapx*sigmapx) \
#              +2*mux*data[:,0]*data[:,1]/(sigmax*sigmapx)
#    data[:,7]=data[:,2]*data[:,2]/(sigmay*sigmay) \
#              +data[:,3]*data[:,3]/(sigmapy*sigmapy) \
#              +2*muy*data[:,2]*data[:,3]/(sigmay*sigmapy)
#    data[:,8]=data[:,4]*data[:,4]/(sigmaz*sigmaz) \
#              +data[:,5]*data[:,5]/(sigmapz*sigmapz) \
#              +2*muz*data[:,4]*data[:,5]/(sigmaz*sigmapz)
#    
#
#    gamma = ke/mass+1.0
#    beta = np.sqrt(1.0-1.0/(gamma*gamma))
#    
#    dimX=2*np.pi*freq/299792458
#    dimPx=gamma*beta
#    #dimPx=1;
#    
#    fig = plt.figure()
#    plt.rcParams.update({'font.size': 9})
#    ax=fig.add_subplot(131)
#    ax.scatter(1E3*data[:,0]/dimX, 1E3*data[:,1]/dimPx,c=-np.sqrt(data[:,6]),marker='.',edgecolors='none' )
#    ax.grid()          
#    ax.xlabel("x (m)")
#    ax.ylabel("p_x (mrad) ")
#    #ax.set_xlim([-6,6])
#    #ax.set_ylim([-2.0,2.0])
#    #print([1E3*sigmax/dimX,1E3*sigmapx/dimPx])    
#    
#    
#    ay=fig.add_subplot(132)
#    ay.scatter(1E3*data[:,2]/dimX, 1E3*data[:,3]/dimPx,c=-np.sqrt(data[:,7]),marker='.',edgecolors='none' )
#    ay.grid()
#    ay.set_title("y[mm]-y'[mrad]")
#    #ay.set_xlim([-6,6])
#    #ay.set_ylim([-2.0,2.0])
#    #print([1E3*sigmay/dimX,1E3*sigmapy/dimPx]) 
#    
#    az=fig.add_subplot(133)
#    az.scatter(180/np.pi*data[:,4], 1E-6*mass*data[:,5],
#               c=-np.sqrt(data[:,8]), marker='.', edgecolors='none' )
#    az.grid()
#    az.set_title("ph[deg@"+str(int(freq*1.0E-6))+"]-W[MeV]")
#    #print([180/np.pi*sigmaz, 1E-6*mass*sigmapz]) 
#    #az.set_xlim([-13,13])
#    #az.set_ylim([-1.0,1.0])
#    if plotRange != None:
#        ax.set_xlim(plotRange[0])
#        ax.set_ylim(plotRange[1])
#        ay.set_xlim(plotRange[2])
#        ay.set_ylim(plotRange[3])
#        az.set_xlim(plotRange[4])
#        az.set_ylim(plotRange[5])
#    if z!=None:
#        plt.title()
#    
#    plt.tight_layout()
#    
#    if saveDir == '':
#        plt.savefig('poincare'+str(n)+'.png', dpi=1200)
#    else:        
#        if not os.path.exists(saveDir):
#            os.makedirs(saveDir)
#        plt.savefig('./'+saveDir+'/poincare'+str(n)+'.png', dpi=240)

#%% plot poincare of sliced bunch
from impact import readParticleDataSliced
def phase_space(fileID, ke, mass, freq, zSliced=True, nSlice=1,
                plotRange = None, fileDir='', saveDir='', showAll=False):
    pData = readParticleDataSliced(nSlice, fileID, ke, mass, freq, zSliced)
    for i in range(nSlice):
        data = pData[i]
        n_particles = len(data)
        if n_particles < 10:
            continue
        
        sigmax=np.std(data[:,0])
        sigmapx=np.std(data[:,1])
        dx=data[:,0]-np.mean(data[:,0])
        dpx=data[:,1]-np.mean(data[:,1])
        mux=-np.mean(dx*dpx)
        data[:,6]= np.sqrt( dx*dx*(sigmapx*sigmapx) 
                            +2.0*mux*dx*dpx+
                            +dpx*dpx*(sigmax*sigmax) )
                
        dy=data[:,2]-np.mean(data[:,2])
        dpy=data[:,3]-np.mean(data[:,3])
        sigmay=np.std(data[:,2])
        sigmapy=np.std(data[:,3])
        muy=-np.mean(dy*dpy)                 
        data[:,7]= np.sqrt( dy*dy*(sigmapy*sigmapy) 
                            +2.0*muy*dy*dpy+
                            +dpy*dpy*(sigmay*sigmay) )
                          
        muxy=-np.mean(dx*dy);         
        data[:,8]= np.sqrt( dx*dx*(sigmay*sigmay) 
                            +2.0*muxy*dx*dy+
                            +dy*dy*(sigmax*sigmax) )
        
        if showAll :
            row_col_i = 141
            fig = plt.figure(figsize=(9.6, 2.4))
        else:
            row_col_i = 131   
            fig = plt.figure(figsize=(7.2, 2.4))
        
        plt.rcParams.update({'font.size': 9})
        ax=fig.add_subplot(row_col_i)
        ax.scatter(1E3*data[:,0], 1E3*data[:,1],
                   c=-data[:,6],marker='.',edgecolors='none' )
        ax.grid()          
        ax.set_xlabel(r"$\mathsf{x}$"+' (mm)',fontsize=10) 
        ax.set_ylabel(r"$\mathsf{p_x}$"+' (mrad)',fontsize=10) 
        
        
        ay=fig.add_subplot(row_col_i+1)
        ay.scatter(1E3*data[:,2], 1E3*data[:,3],
                   c=-data[:,7],marker='.',edgecolors='none' )
        ay.grid() 
        ay.set_xlabel(r"$\mathsf{y}$"+' (mm)',fontsize=10) 
        ay.set_ylabel(r"$\mathsf{p_y}$"+' (mrad)',fontsize=10) 

        
        axy=fig.add_subplot(row_col_i+2)
        axy.scatter(1E3*data[:,0], 1E3*data[:,2],
                   c=-data[:,8], marker='.', edgecolors='none' )
        axy.grid()
        axy.set_xlabel(r"$\mathsf{x}$"+' (mm)',fontsize=10) 
        axy.set_ylabel(r"$\mathsf{y}$"+' (mm)',fontsize=10) 

        if showAll :
            sigmaz=np.std(data[:,4])
            sigmapz=np.std(data[:,5])
            dz=data[:,4]-np.mean(data[:,4])
            dpz=data[:,5]-np.mean(data[:,5])
            muz=-np.mean(dz*dpz);        

            data[:,8]= np.sqrt( dz*dz*(sigmapz*sigmapz) 
                                +2.0*muz*dz*dpz+
                                +dpz*dpz*(sigmaz*sigmaz) )                  
                      
            az=fig.add_subplot(row_col_i+3)
            az.scatter(data[:,4], data[:,5],
                       c=-data[:,8], marker='.', edgecolors='none' )
            az.grid()
            az.set_xlabel(r"$\mathsf{\phi}$"+' (degree)',fontsize=10) 
            az.set_ylabel(r"$\mathsf{\Delta E}$"+' (MeV)',fontsize=10) 

        if plotRange != None:
            ax.set_xlim(plotRange[0])
            ax.set_ylim(plotRange[1])
            ay.set_xlim(plotRange[2])
            ay.set_ylim(plotRange[3])
            axy.set_xlim(plotRange[4])
            axy.set_ylim(plotRange[5])
            
        plt.tight_layout()
        
        if zSliced==True:
            if nSlice !=1:
                fig.suptitle( str(n_particles)+" particles of bunch slice at "+
                              "{0:.3f}".format(np.mean(data[:,4]))+" degree", 
                              y=1.05, fontsize=13)
            if saveDir == '':
                plt.savefig('poincare'+str(fileID)+'_z'+str(i)+'.png', dpi=240)
            else:        
                if not os.path.exists(saveDir):
                    os.makedirs(saveDir)
                plt.savefig('./'+saveDir+'/poincare'+str(fileID)+'_z'+str(i)+'.png', dpi=240)                          
        else:
            if nSlice !=1:
                fig.suptitle( str(n_particles)+" particles of bunch slice at "+
                              "{0:.3f}".format(np.mean(data[:,5]))+"MeV", 
                               y=1.05, fontsize=13)
            if saveDir == '':
                plt.savefig('poincare'+str(fileID)+'_E'+str(i)+'.png', dpi=240)
            else:        
                if not os.path.exists(saveDir):
                    os.makedirs(saveDir)
                plt.savefig('./'+saveDir+'/poincare'+str(fileID)+'_E'+str(i)+'.png', dpi=240)                                                                 
                         
        
        

        
    
        
        
