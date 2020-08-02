import numpy as np
import IPython
import ipywidgets as widgets
import colorsys
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit,QuantumRegister,ClassicalRegister
from qiskit import execute, Aer, BasicAer
from qiskit.visualization import plot_bloch_multivector
from qiskit.tools.jupyter import *
from qiskit.visualization import *
import os
import glob
import moviepy.editor as mpy
import seaborn as sns
sns.set()


'''========State Vector======='''

def getStateVector(qc):
    '''get state vector in row matrix form'''
    backend = BasicAer.get_backend('statevector_simulator')
    job = execute(qc,backend).result()
    vec = job.get_statevector(qc)
    return vec


def vec_in_braket(vec: np.ndarray) -> str:
    '''get bra-ket notation of vector'''
    nqubits = int(np.log2(len(vec)))
    state = ''
    for i in range(len(vec)):
        rounded = round(vec[i], 3)
        if rounded != 0:
            basis = format(i, 'b').zfill(nqubits)
            state += np.str(rounded).replace('-0j', '+0j')
            state += '|' + basis + '\\rangle + '
    state = state.replace("j", "i")
    return state[0:-2].strip()

def vec_in_text_braket(vec):
    return '$$\\text{{State:\n $|\\Psi\\rangle = $}}{}$$'\
                              .format(vec_in_braket(vec))

def writeStateVector(vec):
    return widgets.HTMLMath(vec_in_text_braket(vec))



'''==========Bloch Sphere ========='''

def getBlochSphere(qc):
    '''plot multi qubit bloch sphere'''
    vec = getStateVector(qc)
    return plot_bloch_multivector(vec)


def getBlochSequence(path,figs):
    '''plot block sphere sequence and save it
       to a folder for gif movie creation'''
    try:
        os.mkdir(path)
    except:
        print('Directory already exist')
    for i,fig in enumerate(figs):
        fig.savefig(path+"/rot_"+str(i)+".png")
    return

        
def getBlochGif(figs,path,fname,fps,remove = True):
    '''create gif movie from provided images'''
    file_list = glob.glob(path + "/*.png") 
    list.sort(file_list, key=lambda x: int(x.split('_')[1].split('.png')[0]))
    clip = mpy.ImageSequenceClip(file_list, fps=fps)
    clip.write_gif('{}.gif'.format(fname), fps=fps)
    '''remove all image files after gif creation'''
    if remove:
        for file in file_list:
            os.remove(file)
    return

'''=========Matrix================='''

def getMatrix(qc):
    '''get numpy matrix representing a circuit'''
    backend = BasicAer.get_backend('unitary_simulator')
    job = execute(qc, backend)
    ndArray = job.result().get_unitary(qc, decimals=3)
    Matrix = np.matrix(ndArray)
    return Matrix
    
    
def plotMatrix(M):
    '''visualize a matrix using seaborn heatmap'''
    MD = [["0" for i in range(M.shape[0])] for j in range(M.shape[1])]
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            r = M[i,j].real
            im = M[i,j].imag
            MD[i][j] =  str(r)[0:4]+ " , " +str(im)[0:4]
    plt.figure(figsize = [2*M.shape[1],M.shape[0]])
    sns.heatmap(np.abs(M),\
                annot = np.array(MD),\
                fmt = '',linewidths=.5,\
                cmap='Blues')
    return


'''=========Measurement========'''

def getCount(qc):
    backend= Aer.get_backend('qasm_simulator')
    result = execute(qc,backend).result()
    counts = result.get_counts(qc)
    return counts

def plotCount(counts,figsize):
    plot_histogram(counts)


    
'''========Phase============'''


def getPhaseCircle(vec):
    '''get phase color, angle and radious of phase circir'''
    Phase = []
    for i in range(len(vec)):
        angles = (np.angle(vec[i]) + (np.pi * 4)) % (np.pi * 2)
        rgb = colorsys.hls_to_rgb(angles / (np.pi * 2), 0.5, 0.5)
        mag = np.abs(vec[i])
        Phase.append({"rgb":rgb,"mag": mag,"ang":angles})
    return Phase
    
    
def getPhaseDict(QCs): 
    '''get a dictionary of state vector phase circles for 
       each quantum circuit and populate phaseDict list'''
    phaseDict = []
    for qc in QCs:
        vec = getStateVector(qc)
        Phase = getPhaseCircle(vec)
        phaseDict.append(Phase)  
    return phaseDict


def plotiPhaseCircle(phaseDict,depth,path,show=False,save=False):
    '''plot any quantum circuit phase circle diagram
       from provided phase Dictionary'''
    r = 0.30
    dx = 1.0
    nqubit = len(phaseDict[0])
    fig = plt.figure(figsize = [depth,nqubit])
    for i in range(depth):
        x0 = i
        for j in range(nqubit):
            y0 = j+1
            try:
                mag = phaseDict[i][j]['mag']
                ang = phaseDict[i][j]['ang']
                rgb = phaseDict[i][j]['rgb']
                ax=plt.gca()
                circle1= plt.Circle((dx+x0,y0), radius = r, color = 'white')
                ax.add_patch(circle1)
                circle2= plt.Circle((dx+x0,y0), radius= r*mag, color = rgb)
                ax.add_patch(circle2)
                line = plt.plot((dx+x0,dx+x0+(r*mag*np.cos(ang))),\
                            (y0,y0+(r*mag*np.sin(ang))),color = "black")
                
            except:
                ax=plt.gca()
                circle1= plt.Circle((dx+x0,y0), radius = r, color = 'white')
                ax.add_patch(circle1)
    plt.ylim(nqubit+1,0)
    plt.yticks([y+1 for y in range(nqubit)])
    plt.xticks([x for x in range(depth+2)])
    plt.xlabel("Circuit Depth")
    plt.ylabel("Basis States")
    if show:
        plt.show()
        plt.savefig(path+".png")
        plt.close(fig)
    if save:
        plt.savefig(path +".png")
        plt.close(fig)
    return 

def plotiPhaseCircle_rotated(phaseDict,depth,path,show=False,save=False):
    '''plot any quantum circuit phase circle diagram
       from provided phase Dictionary'''
    r = 0.30
    dy = 1.0
    nqubit = len(phaseDict[0])
    fig = plt.figure(figsize = [nqubit,depth])
    for i in range(depth):
        y0 = i
        for j in range(nqubit):
            x0 = j+1
            try:
                mag = phaseDict[i][j]['mag']
                ang = phaseDict[i][j]['ang']
                rgb = phaseDict[i][j]['rgb']
                ax=plt.gca()
                circle1= plt.Circle((x0,dy+y0), radius = r, color = 'white')
                ax.add_patch(circle1)
                circle2= plt.Circle((x0,dy+y0), radius= r*mag, color = rgb)
                ax.add_patch(circle2)
                line = plt.plot((x0,x0+(r*mag*np.cos(ang))),\
                            (dy+y0,dy+y0+(r*mag*np.sin(ang))),color = "black")
                
            except:
                ax=plt.gca()
                circle1= plt.Circle((x0,dy+y0), radius = r, color = 'white')
                ax.add_patch(circle1)
    plt.ylim(0,depth+1)
    plt.yticks([x+1 for x in range(depth)])
    plt.xticks([y for y in range(nqubit+2)])
    plt.ylabel("Circuit Depth")
    plt.xlabel("Basis States")
    if show:
        plt.show()
        plt.savefig(path+".png")
        plt.close(fig)
    if save:
        plt.savefig(path +".png")
        plt.close(fig)
    return 

def getPhaseSequence(QCs,path,rotated=False):
    '''plot a sequence of phase circle diagram for a given
       sequence of quantum circuits'''
    
    try:
        os.mkdir(path)
    except:
        print("Directory already exist")
    depth = len(QCs)
    phaseDict =[]
    for i,qc in enumerate(QCs):
        vec = getStateVector(qc)
        Phase = getPhaseCircle(vec)
        phaseDict.append(Phase)  
        ipath = path + "phase_" + str(i)
        if rotated:
            plotiPhaseCircle_rotated(phaseDict,depth,ipath,save=True,show=False)
        else:
            plotiPhaseCircle(phaseDict,depth,ipath,save=True,show=False)
            
        
    return
        
    
def getPhaseGif(path,fname,fps,remove = True):
    '''create a gif movie file from phase circle figures'''
    file_list = glob.glob(path+ "/*.png") 
    list.sort(file_list, key=lambda x: int(x.split('_')[1].split('.png')[0]))
    clip = mpy.ImageSequenceClip(file_list, fps=fps)
    clip.write_gif('{}.gif'.format(fname), fps=fps)
    '''remove all image files after gif creation'''
    if remove:
        for file in file_list:
            os.remove(file)
    return  
    
    
    
    
    
    