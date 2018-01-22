# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 21:23:03 2016

@author: kilean
"""
import numpy as np
import os
# IMPACT input generator, output reader and plot 
# key data structure are 'beam' and 'element' dictionary


def run(nCore=None):
    impact_path = os.path.abspath(os.path.dirname(__file__))
    if nCore == None:
        return os.system(impact_path+'/ImpactZ > log')
    else:
        return os.system('mpirun -n '+str(nCore) + ' '+impact_path+'/ImpactZexe > log')
###############################################################################
###############################################################################
###                      IMPACT INPUT GENERATOR                             ###
###############################################################################
###############################################################################
#%%============================================================================
#                                    beam                               
#==============================================================================
def getBeam() :
    """
    Beam = getBeam()
    get a template of a beam dictionary. 
    units : 
        mass : eV
        energy = kinetic energy : eV
        current : Ampere
        x00, x11, x01 : IMPACT internel dimension
        frequency : Hz
        phase : radian
        charge per mass = charge number / mass : 1/eV
    """
    beam = {'nCore_x':1,'nCore_y':1,
            'mass': 938.27231e6, # eV/c^2
            'energy': 150.0e6, # eV
            'n_particles': 1,
            'distribution id':3, # water bag distribution
            'error study flag':0,
            'restart flag':0,
            'standard output':1,
            'current' : 0.0, # ampere
            'x00': 0.0, 'x11': 0.0, 'x01': 0.0, # IMPACT internel dimension
            'y00': 0.0, 'y11': 0.0, 'y01': 0.0,
            'z00': 0.0, 'z11': 0.0, 'z01': 0.0,
            'frequency': 650.0e6, # Hz
            'phase': 0.0,   #radian
            'charge per mass' : 0,
            'mesh_x' : 16, 'mesh_y' : 16, 'mesh_z' : 16}
      
    beam['charge per mass']=1.0/beam['mass']
    return beam

def twiss2beam(beam,betx=0.0,alfx=0.0,norm_ex=0.0,
                    bety=0.0,alfy=0.0,norm_ey=0.0,
                    betz=0.0,alfz=0.0,norm_ez=0.0):
    """
    twiss2beam(beam,betx=0.0,alfx=0.0,norm_ex=0.0,
                    bety=0.0,alfy=0.0,norm_ey=0.0,
                    betz=0.0,alfz=0.0,norm_ez=0.0)
    update the beam distribution using twiss parameters 
    !! IMPORTANT : energy and frequency must be updated before hand    
    input 
        beam = (dict) beam dictionary to be updated
        betx = beta-function in x-direction
        alfx = alpha-function in x-direction
        norm_ex = normalized emittance in x-direction    
        betz = beta-function in z-direction
        alfz = alpha-function in z-direction
        norm_ez = normalized emittance in z-direction      
    """
    clight = 299792458  # m/s
    rel_gamma = 1.0 + beam['energy']/beam['mass']
    rel_beta = (1.0 - 1.0/rel_gamma/rel_gamma)**0.5
    dummy = 2*3.14159265359*beam['frequency']/clight
    
    if betx==0:
        beam['x00']=0.0
        beam['x11']=0.0
        beam['x01']=0.0        
    else:  
        beam['x00']=dummy*( betx* norm_ex/rel_gamma/rel_beta /(1.0+alfx*alfx) )**0.5
        beam['x11']=rel_gamma*rel_beta*( norm_ex/rel_gamma/rel_beta /betx )**0.5
        beam['x01']=alfx / (1.0+alfx*alfx)**0.5
        
    if bety==0:
        beam['y00']=0.0
        beam['y11']=0.0
        beam['y01']=0.0        
    else:  
        beam['y00']=dummy*( bety* norm_ey/rel_gamma/rel_beta /(1.0+alfy*alfy) )**0.5
        beam['y11']=rel_gamma*rel_beta*( norm_ey/rel_gamma/rel_beta /bety )**0.5
        beam['y01']=alfy / (1.0+alfy*alfy)**0.5          

    if betz==0:
        beam['z00']=0.0
        beam['z11']=0.0
        beam['z01']=0.0 
    else:  
        beam['z00']=3.14159265359/180*( betz*norm_ez/(1.0+alfz*alfz) )**0.5
        beam['z11']=( norm_ez/betz )**0.5 /beam['mass']*1.0E6
        beam['z01']=alfz / (1.0+alfz*alfz)**0.5    
        
        
def beam2str(beam):
    """
    beamStrList = beam2str(beam)
    from beam to string list of IMPACT format
    input 
       x = (dict) beam 
    output 
        beamStrList = (list) list of string of IMPACT format of beam
    """
    beamStrList=[str(beam['nCore_x'])+' '+str(beam['nCore_y'])+' \n',
                 '6 '+str(beam['n_particles'])+' 2 '+\
                 str(beam['error study flag'])+' '+\
                 str(beam['standard output'])+' \n',
                 str(beam['mesh_x'])+' '+str(beam['mesh_y'])+' '+\
                 str(beam['mesh_z'])+' 1 0.1 0.1 0.1 \n',
                 str(beam['distribution id'])+' '+\
                 str(beam['restart flag'])+' 0 1 \n',
                 str(beam['n_particles'])+' \n',
                 str(beam['current'])+' \n',
                 str(beam['charge per mass'])+' \n',
                 str(beam['x00'])+' '+str(beam['x11'])+' '+str(beam['x01'])+' '+'1 1 0 0 \n',
                 str(beam['y00'])+' '+str(beam['y11'])+' '+str(beam['y01'])+' '+'1 1 0 0 \n',
                 str(beam['z00'])+' '+str(beam['z11'])+' '+str(beam['z01'])+' '+'1 1 0 0 \n',
                 str(beam['current'])+' '+str(beam['energy'])+' '+str(beam['mass'])+' 1.0 '+\
                 str(beam['frequency'])+' '+str(beam['phase'])+' \n',
                 '!==lattice=======================================\n' 
                 ]
    return beamStrList

     
def str2beam(beamStr):
    """
    beam = str2beam(beamStr)
    from string list of IMPACT format to beam
    input 
        beamStr = (list) list of string of IMPACT format
    output 
        beam = (dict) beam dictionary
    """
    beam=getBeam()
    beam['nCore_x'],beam['nCore_y']=[int(beamStr[0].split()[0]) for i in range(2)]
    beam['n_particles']=int(beamStr[1].split()[1])
    beam['mesh_x'],beam['mesh_y'],beam['mesh_z']=[int(beamStr[2].split()[0]) for i in range(3)]
    #beam['distribution id']=int(beamStr[3].split()[0])
    #beam['restart flag']=int(beamStr[3].split()[1])
    beam['current']=float(beamStr[5].split()[0])
    beam['charge per mass']=float(beamStr[6].split()[0])
    beam['x00'],beam['x11'],beam['x01']=[float(beamStr[7].split()[i]) for i in range(3)]
    beam['y00'],beam['y11'],beam['y01']=[float(beamStr[8].split()[i]) for i in range(3)]
    beam['z00'],beam['z11'],beam['z01']=[float(beamStr[9].split()[i]) for i in range(3)]
    beam['current'],beam['energy'],beam['mass']=[float(beamStr[10].split()[i]) for i in range(3)]
    beam['frequency'],beam['phase']=[float(beamStr[10].split()[i]) for i in range(4,6)]
    return beam         

#%%============================================================================
#                                   lattice                                  
#==============================================================================

#%%#================================element====================================
def getElem(elemType) : 
    """
    f = getElem(elemType)
    get a template of an element dictionary.  
    input 
        elemType = (str) element type. one of the following
                   'drift', 'quad', 'bend', 'scrf', 
                   'kick', 'write full', 'restart', 'halt'
    output 
        f = (dict) element dictionary
    """
    if elemType=='drift' :
        return {'type':'drift', 'length': 0.1, 'n_sckick': 1, 
                'n_map': 1, 'radius': 1.0}   
    elif elemType=='quad' :
        return {'type':'quad', 'length': 0.1, 'n_sckick': 2, 'n_map': 1, 
                'B1': 17.0, 'input file id' : 0, 'radius': 1.0,    # B1 [T/m]
                'dx': 0.0, 'dy': 0.0 }
    elif elemType=='bend' :
        return {'type':'bend', 'length': 1.0, 'n_sckick': 25, 'n_map': 1, 
                'angle': 0.9, 'k1': 0.0, 'input switch' : 150,     #angle [rad]
                'radius': 1.0, 'entrance edge':0.0, 'exit edge':0.0, #edge[rad]
                'entrance curvature':0.0, 'exit curvature':0.0, 'FINT':0.0}   
    elif elemType=='scrf' :
        return {'type':'scrf', 'length': 0.948049, 'n_sckick': 100, 'n_map': 1,
                'scale': 34.0e6, 'frequency': 650e6 , 'phase': 0.0, 
                'input file id' : 1, 'radius': 1.0}   # phase [degree]
    elif elemType=='kick' :
        return {'type':'kick', 'length': 0.0,'dx': 0.0 , 'dpx': 0.0, 'dy' : 0.0, 'dpy': 0.0, 
                'dz' : 0.0, 'dpz': 0.0}   
                # dx,dy in meter, dpx,dpy in radian, dz in degree, dpz in MeV   
    elif elemType=='write full' :
        return {'type':'write full', 'length': 0.0, 'file id': 1000}
    elif elemType=='restart' :
        return {'type':'restart', 'length': 0.0}      
    elif elemType=='halt' :
        return {'type':'halt', 'length': 0.0}              
        
def elem2str(elemDict): 
    """
    f = elem2str(elemDict)
    from element to (IMPACT format) string
    input 
        x = (dict) element dictionary
    output 
        f = (str) element string in IMPACT format
    """
    if elemDict['type']=='drift':
        return str(elemDict['length'])+' '+str(elemDict['n_sckick'])+' '+\
               str(elemDict['n_map'])+' 0 '+str(elemDict['radius'])+' /\n'
    elif elemDict['type']=='quad' :
        f = str(elemDict['length'])+' '+str(elemDict['n_sckick'])+' '+\
               str(elemDict['n_map'])+' 1 '+str(elemDict['B1'])+' '+\
               str(elemDict['input file id'])+' '+str(elemDict['radius'])
        if 'dx' in elemDict :
            f = f+ ' '+str(elemDict['dx'])
            if 'dy' in elemDict :
                f = f+ ' '+str(elemDict['dy'])
        return f+' /\n' 
    elif elemDict['type']=='bend' :
        return str(elemDict['length'])+' '+str(elemDict['n_sckick'])+' '+\
               str(elemDict['n_map'])+' 4 '+str(elemDict['angle'])+' '+\
               str(elemDict['k1'])+' '+str(elemDict['input switch'])+' '+\
               str(elemDict['radius'])+' '+str(elemDict['entrance edge'])+' '+\
               str(elemDict['exit edge'])+' '+\
               str(elemDict['entrance curvature'])+' '+\
               str(elemDict['exit curvature'])+' '+str(elemDict['FINT'])+' /\n' 
    elif elemDict['type']=='scrf' :
        return str(elemDict['length'])+' '+str(elemDict['n_sckick'])+' '+\
               str(elemDict['n_map'])+' 104 '+str(elemDict['scale'])+' '+\
               str(elemDict['frequency'])+' '+str(elemDict['phase']) +' '+\
              str(elemDict['input file id'])+' '+str(elemDict['radius'])+' /\n' 
    elif elemDict['type']=='kick' : 
        return '0.0 0 0 -21 1.0 '+str(elemDict['dx'])+' '+\
               str(elemDict['dpx'])+' '+str(elemDict['dy'])+' '+\
               str(elemDict['dpy'])+' '+str(elemDict['dz'])+' '+\
               str(elemDict['dpz'])+' /\n' 
    elif elemDict['type']=='write full' : 
        return '0.0 0 '+str(elemDict['file id'])+' -2 0.0 1 /\n'
    elif elemDict['type']=='restart' : 
        return '0.0 0 0 -7 0 /\n'
    elif elemDict['type']=='halt' : 
        return '0.0 1 1 -99 0 /\n'        

    
def str2elem(elemStr): 
    """
    elemtDict = str2elem(elemStr)
    from (IMPACT format) string to element  
    input 
        elemStr = (str) string of a IMPACT lattice line
    output 
        elemtDict = (dict) element dictionary 
    """
    elemStr = elemStr.split()
    elemID=int(float(elemStr[3]))
    if elemID == 0:
        elemtDict = {'type':'drift',
                     'length': float(elemStr[0]),
                     'n_sckick': int(elemStr[1]), 
                     'n_map': int(elemStr[2]), 
                     'radius': float(elemStr[4])
                    }   
    elif elemID == 1:
        elemtDict = {'type':'quad',
                     'length': float(elemStr[0]),
                     'n_sckick': int(elemStr[1]), 
                     'n_map': int(elemStr[2]), 
                     'B1': float(elemStr[4]), 
                     'input file id': int(elemStr[5]), 
                     'radius': float(elemStr[6]) }
        if len(elemStr)>8:
                     elemtDict['dx']=float(elemStr[7])
        if len(elemStr)>9:
                     elemtDict['dy']=float(elemStr[8]) 
    elif elemID == 4:
        elemtDict = {'type':'bend',
                     'length': float(elemStr[0]),
                     'n_sckick': int(elemStr[1]), 
                     'n_map': int(elemStr[2]), 
                     'angle': float(elemStr[4]), 
                     'k1': float(elemStr[5]), 
                     'input switch': int(float(elemStr[6])),
                     'radius': float(elemStr[7]),
                     'entrance edge': float(elemStr[8]),
                     'exit edge': float(elemStr[9]),
                     'entrance curvature': float(elemStr[10]),
                     'exit curvature': float(elemStr[11]),
                     'FINT': float(elemStr[12]),
                    }                     
    elif elemID == 104:
        elemtDict= {'type':'scrf',
                    'length': float(elemStr[0]),
                    'n_sckick': int(elemStr[1]), 
                    'n_map': int(elemStr[2]), 
                    'scale': float(elemStr[4]), 
                    'frequency': float(elemStr[5]), 
                    'phase': float(elemStr[6]), 
                    'input file id': int(elemStr[7]), 
                    'radius': float(elemStr[8])
                   }
    elif elemID == -2:
        elemtDict= {'type':'write full',
                    'file id': int(elemStr[2])}                   
    elif elemID == -7:
        elemtDict= {'type':'restart'}
    elif elemID == -99:
        elemtDict= {'type':'halt'} 
    else :
        elemtDict= {}            
    return elemtDict
  
#%%=================================lattice====================================
#lattice is a list of element dictionaries
def lattice2str(lattice):
    """
    f = lattice2str(lattice)
    from lattice to string list of IMPACT format
    input 
        x = (list) lattice 
    output 
        f = (list) list of string of of IMPACT format
    """
    latticeStr = []
    for i in range(len(lattice)):
        latticeStr.append(elem2str(lattice[i]))
    return latticeStr
    
def str2lattice(latticeStr):
    """
    lattice = str2lattice(latticeStr)
    from string list of IMPACT format to lattice
    input 
        latticeStr = (list) list of string of of IMPACT format
    output 
        lattice = (list) list of element dictionaries
    """
    lattice = []
    for i in range(len(latticeStr)):
        if latticeStr[i] in ['\n', '\r\n'] or latticeStr[i][0]=='!':
            continue
        elem = str2elem(latticeStr[i])
        if elem : #check if elem is not empty
            lattice.append(elem)
    return lattice    


def getElemIndex(lattice,typename):
    """
    oneElemLattice = getElemIndex(lattice,typename)
    from lattice list of dictionary to single element lattice list of dictionary
    input 
        lattice = (list) list of elements dictionaries
    output 
        oneElemLattice = (list) list of one type(=typename) element dictionaries
    """
    f=[]
    for i in range(len(lattice)):
        if lattice[i]['type']==typename:
            f.append(i)
    return f
#%%============================================================================
#                               IMPACT test.in I/O                                
#==============================================================================  
def writeIMPACT(filename,beam,lattice=[]):
    """
    write a IMPACT input file
    input 
        beam = (dict) beam dictionary
        lattce = (list) list of element dictionaries
    """
    beamStrList=beam2str(beam)                 
    latticeStrList=lattice2str(lattice)
    
        
    f=open(filename,'w') 
    f.writelines(beamStrList)
    f.writelines(latticeStrList)
    f.close()


def updateLattice(lattice):
    """
    update lattice such that each begining location of element is saved
    input 
        lattice = (list) list of elemenet dictionaries
    """
    z=0.
    for i in range(len(lattice)):
        lattice[i]['z']=z
        z=z+lattice[i]['length']
        
        
def readIMPACT(filename='test.in'):
    """
    beam, lattice = readIMPACT(filename='test.in')
    read a IMPACT input file 
    output : (list) element dictionaries
    """
    file = open(filename,'r')
    lines = file.readlines()
    file.close()
    row_end = len(lines)
    row_start = 11
    for i in range(0,row_end):
        if lines[i][:12] == '!==lattice==':
            row_start=i
            break
    lattice=str2lattice(lines[row_start+1:])
    updateLattice(lattice)
    beam=str2beam(lines[0:row_start])
    return beam, lattice
        
 
#%%############################################################################
###############################################################################
###                          IMPACT OUTPUT Reader                           ###
###############################################################################
###############################################################################   
def readReferenceOrbit(fileloc=''):
    file = open(fileloc+'fort.18','r')
    lines = file.readlines()
    file.close()
    f=[]
    for j in range(len(lines)) :
        f.append( [ float(lines[j].split()[i]) for i in range(5) ] )
    return f

def readReferenceOrbitAt(zIndex,fileloc=''):
    file = open(fileloc+'fort.18','r')
    lines = file.readlines()
    file.close()
    return [float(lines[zIndex].split()[i]) for i in range(5)]

def readReferenceOrbitAtEnd(fileloc=''):
    return readReferenceOrbit(fileloc)[-1]

def readZIndex(z,fileloc='',flagRef=False):
    rf = readReferenceOrbit(fileloc)
    for i in range(1,len(rf)) :
        z1 =  rf[i][0]
        if z-z1 < 0 :
            break
    if z1-z < z-rf[i-1][0]:
        if flagRef:
            return i,rf[i]
        else:
            return i,rf[i][0]
    else:
        if flagRef:
            return i-1,rf[i-1]
        else:
            return i-1,rf[i-1][0]
    
    
def readBeamSize(direction,nSkip=1,fileLoc=''):
    """
    f = readBeamSize(direction,nSkip=1,fileLoc='')
    Read RMS beam size
    input 
        direction = (char) 'x', 'y' or 'z'
        nSkip = (int>0) number of lines to skip when reading output    
    output 
        f = (list) each element of list is a vector of 
                   rms_x,px,y,py,z,E (meter, rad, deg, MeV)
    """
    if direction == 'x':
        file = open(fileLoc+'fort.24','r')
    elif direction == 'z':
        file = open(fileLoc+'fort.26','r')
    else  :
        file = open(fileLoc+'fort.25','r')
    lines = file.readlines()
    file.close()
    f=[]     
    for i in range(0,len(lines),nSkip) :
        f.append(float(lines[i].split()[2]))
    return f


def readBeamSizeAt(zIndex,direction,nSkip=1,fileLoc=''):
    """
    f = readBeamSizeAt(zIndex,direction,nSkip=1,fileLoc=''):
    Read RMS beam size at location corresponds to zIndex
    input 
        direction = (char) 'x', 'y' or 'z'
        nSkip = (int>0) number of lines to skip when reading output    
    output 
        f = (list) each element of list is a vector of 
                   rms_x,px,y,py,z,E (meter, rad, deg, MeV)
    """
    if direction == 'x':
        file = open(fileLoc+'fort.24','r')
    elif direction == 'z':
        file = open(fileLoc+'fort.26','r')
    else  :
        file = open(fileLoc+'fort.25','r')
    lines = file.readlines()
    file.close()
    
    return float(lines[zIndex].split()[2])    
    
def readBeamSizeAtEnd(fileLoc=''):
    return [readBeamSizeAt(-1,'x',fileLoc),
            readBeamSizeAt(-1,'y',fileLoc), readBeamSizeAt(-1,'z',fileLoc) ]
           
          
    
def readOptics(direction,nSkip=1,fileLoc=''):
    """
    f = readOptics(direction,nSkip=1,fileLoc='')
    Read Optics functions ( Optics ftn is calcualted using beam porfile )
    input 
        zIndex = (int) 
        nSkip  = sampling rate from beam distribution output file of IMPACTz
        fileLoc = (string) path
    output 
        f = optics parameters at every (sampling rate of nSkip) 
    """
    if direction == 'x':
        file = open(fileLoc+'fort.24','r')
    elif direction == 'z':
        file = open(fileLoc+'fort.26','r')
    else  :
        file = open(fileLoc+'fort.25','r')
    lines = file.readlines()
    file.close()
    f=[]
    for j in range(0,len(lines),nSkip) :
        sigmax, sigmap, alpha, emittance_norm = [ float(lines[j].split()[i]) for i in [2,4,5,6] ]
        beta=(1+alpha*alpha)**0.5 *sigmax/sigmap
        if direction == 'z':
            f.append( [beta, alpha, emittance_norm] )
        else:
            f.append( [beta, alpha, emittance_norm] )
    return f
    
def readOpticsAt(zIndex, direction, fileLoc=''):
    """
    beta, alpha, emittance_norm = readOpticsAt(zIndex, direction, fileLoc=''):
    Read Optics function at location corresponds to zIndex
    input 
        zIndex = (int) 
        fileLoc = (string) path
    output 
        optics parameters at the location specified by zIndex
    """
    if direction == 'x':
        file = open(fileLoc+'fort.24','r')
    elif direction == 'z':
        file = open(fileLoc+'fort.26','r')
    else  :
        file = open(fileLoc+'fort.25','r')
    lines = file.readlines()
    file.close()
    sigmax, sigmap, alpha, emittance_norm = [ float(lines[zIndex].split()[i]) for i in [2,4,5,6] ]
    beta=(1+alpha*alpha)**0.5 *sigmax/sigmap
    
    if direction == 'z':
        return beta, alpha, emittance_norm
    else  :
        return beta, alpha, emittance_norm
   
def readOpticsAtEnd(fileLoc=''):
    return readOpticsAt(-1,'x',fileLoc) +\
           readOpticsAt(-1,'y',fileLoc) + readOpticsAt(-1,'z',fileLoc)
   
def readCentroid(direction, nSkip=1, fileLoc=''):
    """
    f = readCentroid(direction, nSkip=1, fileLoc='')
    Read RMS beam size
    input 
        direction = (char) 'x', 'y' or 'z'
        nSkip = (int>0) number of lines to skip when reading output    
    output 
        f = (list) each element of list is a vector of 
                   rms_x,px,y,py,z,E (meter, rad, deg, MeV)
    """
    if direction == 'x':
        file = open(fileLoc+'fort.24','r')
    elif direction == 'z':
        file = open(fileLoc+'fort.26','r')
    else  :
        file = open(fileLoc+'fort.25','r')
    lines = file.readlines()
    file.close()
    f=[]     
    for i in range(0,len(lines),nSkip) :
        f.append( [float(lines[i].split()[1]), float(lines[i].split()[3])] )
    return f    

def readCentroidAt(zIndex, direction, fileLoc=''):
    """
    f = readCentroidAt(zIndex, direction, fileLoc='')
    Read RMS beam size
    input 
        direction = (char) 'x', 'y' or 'z'
        nSkip = (int>0) number of lines to skip when reading output    
    output 
        f = (list) each element of list is a vector of 
                   rms_x,px,y,py,z,E (meter, rad, deg, MeV)
    """
    if direction == 'x':
        file = open(fileLoc+'fort.24','r')
    elif direction == 'z':
        file = open(fileLoc+'fort.26','r')
    else  :
        file = open(fileLoc+'fort.25','r')
    lines = file.readlines()
    file.close()
    return float(lines[zIndex].split()[1]), float(lines[zIndex].split()[3])

def readCentroidAtEnd(fileLoc=''):
    return readCentroidAt(-1,'x',fileLoc) +\
           readCentroidAt(-1,'y',fileLoc) + readCentroidAt(-1,'z',fileLoc)    
           
def readLoss(nSkip=1,fileLoc=''):
    file = open(fileLoc+'fort.32','r')
    lines = file.readlines()
    file.close()
    f=[]     
    for i in range(0,len(lines),nSkip) :
        f.append( int(lines[i].split()[1]) )
    return f    
    
def readLossAt(zIndex, fileLoc=''):
    file = open(fileLoc+'fort.32','r')
    lines = file.readlines()
    file.close()
    return int(lines[zIndex].split()[1])
    
def readLossAtEnd(fileLoc=''):
    file = open(fileLoc+'fort.32','r')
    lines = file.readlines()
    file.close()
    return int(lines[-1].split()[1])
    
        
#%%############################################################################
###############################################################################
###                      Particle data Manipulator                          ###
###############################################################################
############################################################################### 

def normalizeParticleData(data, ke, mass, freq):
    gamma = ke/mass+1.0
    beta = np.sqrt(1.0-1.0/(gamma*gamma))
    x_norm = 2*freq*3.141592653589793/299792458
    px_norm = gamma*beta
    data[:,0] = data[:,0]*x_norm
    data[:,1] = data[:,1]*px_norm
    data[:,2] = data[:,2]*x_norm
    data[:,3] = data[:,3]*px_norm
    data[:,4] = np.pi/180*data[:,4]
    data[:,5] = data[:,5]/mass
    return data
    
def unNormalizeParticleData(data, ke, mass, freq):
    gamma = ke/mass+1.0
    beta = np.sqrt(1.0-1.0/(gamma*gamma))
    x_norm = 2*freq*3.141592653589793/299792458
    px_norm = gamma*beta
    data[:,0] = data[:,0]/x_norm
    data[:,1] = data[:,1]/px_norm
    data[:,2] = data[:,2]/x_norm
    data[:,3] = data[:,3]/px_norm
    data[:,4] = 180/np.pi*data[:,4]
    data[:,5] = mass*data[:,5]
    return data

def readParticleData(fileID, ke, mass, freq, fileLoc=''):
    data=np.loadtxt(fileLoc+'fort.'+str(fileID))
    return unNormalizeParticleData(data, ke, mass, freq)
    
def readParticleDataSliced(nSlice, fileID, ke, mass, freq, zSliced=True, fileLoc=''):
    data=np.loadtxt(fileLoc+'fort.'+str(fileID))
    data=unNormalizeParticleData(data, ke, mass, freq)
    
    f=[]    
    if zSliced:
        z_min = min(data[:,4])
        z_max = max(data[:,4])
        dz = (z_max-z_min)/float(nSlice)
        for i in range(nSlice):
            temp = []
            for j in range(len(data[:,4])):
                if z_min + i*dz < data[j,4] < z_min + (i+1)*dz :
                    temp.append(data[j,:])
            f.append(np.array(temp))
        return f
        
    else:
        ke_min = min(data[:,5])
        ke_max = max(data[:,5])
        dke = (ke_max-ke_min)/float(nSlice)
        for i in range(nSlice):
            temp = []
            for j in range(len(data[:,5])):
                if ke_min + i*dke < data[j,5] < ke_min + (i+1)*dke :
                    temp.append(data[j,:])
            f.append(np.array(temp))  
        return f
        
def writeParticleData(data, ke, mass, freq, fileLoc='',filename='partcl.data'):
    data=normalizeParticleData(data, ke, mass, freq)
    np.savetxt(filename,data,header=str(len(data))+' 0. 0.',comments='')
#%%############################################################################
###############################################################################
###                           Lattice Manipulator                           ###
###############################################################################
############################################################################### 

def getZIndex(lattice,z):
    updateLattice(lattice)
    
    for i in range(len(lattice)):
        z1 =  lattice[i]['z']+lattice[i]['length']
        if z-z1 < 0 :
            break
    if z1-z < z-lattice[i]['z'] :
        return i+1,z1
    else:
        return i,lattice[i]['z']
        
    
#def divideLattice(Z,lattice):
#    updateLattice(lattice)
#    z_index=-1
#    while Z(z_index) > lattice[-1]['z']+lattice[-1]['length']:
#        z_index = z_index -1
#        
#    for i in range(-len(lattice), 0):
#        if lattice[i]['z'] < Z[z_index]:
#            if lattice[i]['type'] == 'scrf':
                             
            
#%%############################################################################
###############################################################################
###                   Physical analysis with Impact                         ###
###############################################################################
###############################################################################  

def getTransferMap(lattice,q,mass,ke,freq,
                   epsilon=[1e-8,1e-6,1e-8,1e-6,1e-7,1e-8],
                   #epsilon=[1e-11,1e-9,1e-11,1e-9,1e-11,1e-11],
                   fname='test.in' ):
    """
    M = getTransferMap(lattice,q,mass,ke,freq,
                       epsilon=[1e-9,1e-6,1e-9,1e-6,1e-3,1e-9],
                       fname='test.in' )
    get linear transfer map (without space-charge)  by tracking 6 particles
    whose initial phase-space perturbation given by epsilon
    input                
        lattice = (dict) lattice dictionary whose transvermap to be determined
        q = charge (-1.0 for electron)
        ke = reference particle energy at lattice entrance [MeV]
        mass = particle mass [MeV]
        freq = reference RF frequency [Hz]
        epsilon = 6 dimension array of perturbation for 
                  x,px,y,py, z*360/v/freq, E  in unit of 
                  [m],[rad],[m],[rad],[deg],[MeV]
                  default : epsilon = [1e-e-8,1e-6,1e-8,1e-6,1e-7,1e-8]
        fname = input file name for IMPACTz 
                depending on version of IMPACTz the file name may change
                default is 'test.in'
    """
    beam = getBeam()
    beam['mass'] = mass*1e6
    beam['charge per mass'] = q/beam['mass']
    beam['energy'] = ke*1e6
    beam['n_particles'] = 6
    beam['frequency'] = freq
    beam['distribution id'] = 19

    fPool = []
    for i in range(len(lattice)):
        if 'file id' in lattice[i]:
            fPool.append(lattice[i]['file id'])
    fileID = 5555
    while True:
        if fileID in fPool:
            fileID = np.random.randint(low=1000,high=9999)
        else:
            break
    lattice.append(getElem('write full'))
    lattice[-1]['file id'] = fileID
    writeIMPACT(fname,beam,lattice)
    lattice.pop()
    
    data = np.zeros([6,9])
    for i in range(6):
        data[i,i] = epsilon[i]
        data[i,8] = i+1
    data[:,6] = beam['charge per mass']
    writeParticleData(data, ke, mass, freq)
    
    run()
    
    dataOut = readParticleData(fileID, ke, mass, freq)[:,:6]
    #os.system('rm fort.'+str(fileID))
    m,n = dataOut.shape
    M = np.zeros([6,6])
    if m<6:
        print 'particle lost observed. too large inital perturbation'
    else:
        for i in range(6):
            M[:,i] = dataOut[i,:]/epsilon[i]
    return M
    

    


#%%############################################################################
###############################################################################
###                           OUTPUT Manipulator                            ###
###############################################################################
###############################################################################           
