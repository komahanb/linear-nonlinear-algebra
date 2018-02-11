# Import numpy
import numpy as np

# Import matplotlib
import matplotlib.pylab as plt

# Configure 
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

# Optionally set font to Computer Modern to avoid common missing font errors
params = {
  'axes.labelsize': 24,
  'legend.fontsize': 14,
  'xtick.labelsize': 24,
  'ytick.labelsize': 24,
  'text.usetex': True}
plt.rcParams.update(params)

# Latex math
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}']
plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'courier'
plt.rcParams['font.size'] = 18
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.color'] = 'r'

# Make sure everything is within the frame
plt.rcParams.update({'figure.autolayout': True})

# Set marker size
markerSize = 11.0

# bar chart settings
alpha = 0.8

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

mevery = 90

class PostOpt:
    """
    Class to handle plotting of post optimization data.
    """
    @staticmethod
    def plot_convergence(data1, data2, data3, name):
        plt.figure()
        fig, ax = plt.subplots()

        # remove the upper and right borders
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
       
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.loglog(data1[:,0], data1[:,1], '-', label='jacobi',
                     ms=markerSize, mec='black', color=tableau20[0], 
                     alpha=alpha)
        ax.loglog(data2[:,0], data2[:,1], '-', label='seidel',
                     ms=markerSize, mec='black', color=tableau20[2], 
                     alpha=alpha)
        ax.loglog(data3[:,0], data3[:,1], '-', label='optimal sor',
                     ms=markerSize, mec='black', color=tableau20[4], 
                     alpha=alpha)
        
        ax.set_xlabel('Number of Iterations')
        ax.set_ylabel('Error (Tolerance)')
        #niters_labl = [1, 20, 40, 60, 80, 100]
        ax.set_yticks(np.logspace(-6, 0, 7, endpoint=True))
        #ax.set_xticks(np.linspace(1, 75, 5, endpoint=True, dtype=np.int))
        ax.legend(loc='upper right', framealpha=0.0)
        
        plt.savefig(name, bbox_inches='tight', pad_inches=0.05)
        
        return

    @staticmethod
    def plot_relaxation(data, name):
        plt.figure()
        fig, ax = plt.subplots()

        # remove the upper and right borders
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
       
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.semilogy(data[:,0], data[:,2], '-o', 
                     ms=markerSize, mec='black', color=tableau20[0], markevery=1,
                     alpha=alpha)
        
        ax.set_xlabel('Relaxation Parameter')
        ax.set_ylabel('Number of Iterations')
        #niters_labl = [1, 20, 40, 60, 80, 100]
        ax.set_yticks(np.logspace(3, 5, 3, endpoint=True))
        #ax.set_xticks(np.linspace(1, 75, 5, endpoint=True, dtype=np.int))
        #ax.legend(loc='upper left', framealpha=0.0)
        
        plt.savefig(name, bbox_inches='tight', pad_inches=0.05)
        
        return
    
    @staticmethod
    def plot_solution(data1, data2, data3, name):
        plt.figure()
        fig, ax = plt.subplots()

        # remove the upper and right borders
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
       
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.plot(data1[:,0], data1[:,1], '-D', label='jacobi',
                     ms=markerSize, mec='black', color=tableau20[0], markevery=mevery,
                     alpha=alpha)
        ax.plot(data2[:,0], data2[:,1], '-v', label='seidel',
                     ms=markerSize, mec='black', color=tableau20[2], markevery=mevery,
                     alpha=alpha)
        ax.plot(data3[:,0], data3[:,1], '-*', label='optimal sor',
                     ms=markerSize, mec='black', color=tableau20[4], markevery=mevery,
                     alpha=alpha)
        ax.plot(data3[:,0], data3[:,2], '--', label='exact',
                ms=markerSize, mec='black', color=tableau20[6], markevery=mevery,
                alpha=alpha)
        
        ax.set_ylabel('Solution u(x)')
        ax.set_xlabel('Domain x')
        #niters_labl = [1, 20, 40, 60, 80, 100]
        #ax.set_yticks(np.logspace(-5, 0, 6, endpoint=True))
        #ax.set_xticks(np.linspace(1, 75, 5, endpoint=True, dtype=np.int))
        ax.legend(loc='upper left', framealpha=0.0)
        
        plt.savefig(name, bbox_inches='tight', pad_inches=0.05)
        
        return
 
if __name__ == "__main__":

    inpFile = open("jacobi.dat", "r")
    fjacobi = list(inpFile.readlines())
    inpFile.close()

    inpFile = open("seidel.dat", "r")
    fseidel = list(inpFile.readlines())
    inpFile.close()

    inpFile = open("sor.dat", "r")
    fsor = list(inpFile.readlines())
    inpFile.close()

    soljacobi = []
    solseidel = []
    solsor = []
    
    for line in fjacobi:
        entry = line.split()
        soljacobi.append([float(entry[0]), float(entry[1]), float(entry[2])])
    soljacobi = np.array(soljacobi)

    for line in fseidel:
        entry = line.split()
        solseidel.append([float(entry[0]), float(entry[1]), float(entry[2])])
    solseidel = np.array(solseidel)
    
    for line in fsor:
        entry = line.split()
        solsor.append([float(entry[0]), float(entry[1]), float(entry[2])])
    solsor = np.array(solsor)

    # Plot the solution in the profile
    PostOpt.plot_solution(soljacobi,solseidel, solsor, "profile.pdf")

    # plot omega versus niters
    inpFile = open("omega.dat", "r")
    fomega = list(inpFile.readlines())
    inpFile.close()
    
    solomega = []
    for line in fomega:
        entry = line.split()
        print entry
        solomega.append([float(entry[0]), float(entry[1]), int(entry[2])])
    solomega = np.array(solomega)
   
    PostOpt.plot_relaxation(solomega, "relaxation_study.pdf")

    inpFile = open("jacobi.log", "r")
    fjacobi = list(inpFile.readlines())
    inpFile.close()

    inpFile = open("seidel.log", "r")
    fseidel = list(inpFile.readlines())
    inpFile.close()

    inpFile = open("sor.log", "r")
    fsor = list(inpFile.readlines())
    inpFile.close()

    soljacobi = []
    solseidel = []
    solsor = []
    
    for line in fjacobi:
        entry = line.split()
        soljacobi.append([int(entry[0]), float(entry[1])])
    soljacobi = np.array(soljacobi)

    for line in fseidel:
        entry = line.split()
        solseidel.append([int(entry[0]), float(entry[1])])
    solseidel = np.array(solseidel)
    
    for line in fsor:
        entry = line.split()
        solsor.append([int(entry[0]), float(entry[1])])
    solsor = np.array(solsor)

    # Plot the solution in the profile
    PostOpt.plot_convergence(soljacobi,solseidel, solsor, "convergence.pdf")
