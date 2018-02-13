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
    Class to handle plotting
    """
    @staticmethod
    def plot_solution(data1, data2, name):
        plt.figure()
        fig, ax = plt.subplots()

        # remove the upper and right borders
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
       
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.plot(data1[:,0], data1[:,1], '-D', label='(a) cg',
                     ms=markerSize, mec='black', color=tableau20[0], markevery=mevery,
                     alpha=alpha)
        ax.plot(data1[:,0], data1[:,2], '-v', label='(a) pcg',
                     ms=markerSize, mec='black', color=tableau20[0], markevery=mevery,
                     alpha=alpha)
        ax.plot(data1[:,0], data1[:,3], '-*', label='(a) exact',
                     ms=markerSize, mec='black', color=tableau20[0], markevery=mevery,
                     alpha=alpha)


        ax.plot(data2[:,0], data2[:,1], '-D', label='(b) cg',
                     ms=markerSize, mec='black', color=tableau20[2], markevery=mevery,
                     alpha=alpha)
        ax.plot(data2[:,0], data2[:,2], '-v', label='(b) pcg',
                     ms=markerSize, mec='black', color=tableau20[2], markevery=mevery,
                     alpha=alpha)
        ax.plot(data2[:,0], data2[:,3], '-*', label='(b) exact',
                     ms=markerSize, mec='black', color=tableau20[2], markevery=mevery,
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
    
    inpFile = open("dirichlet.dat", "r")
    fdirichlet = list(inpFile.readlines())
    inpFile.close()

    inpFile = open("mixed.dat", "r")
    fmixed = list(inpFile.readlines())
    inpFile.close()

    soldirichlet = []
    for line in fdirichlet:
        entry = line.split()
        soldirichlet.append([float(entry[0]), float(entry[1]), float(entry[2])])
    soldirichlet = np.array(soldirichlet)

    solmixed = []
    for line in fmixed:
        entry = line.split()
        solmixed.append([float(entry[0]), float(entry[1]), float(entry[2])])
    solmixed = np.array(solmixed)

    # Plot the solution in the profile
    PostOpt.plot_solution(soldirichlet, solmixed, "profile.pdf")
