import SuperQuantModel as SQM
import numpy as np
import matplotlib.pyplot as plt

# module thet contains most common units
import MicrowaveUnits as ut

"""
 _________________________
| DEFINE MODEL PARAMETERS |
|_________________________|
"""
# small parameter. Will be the value of parasitic capacitances
# and will be used to compare collinearity of vectors
SMALL = 1E-4
    
# Code describing the circuit
# pseudographic drawing of the circuit simply helps to recognize the circuit

MyCircuitCode = [\
#                              _____ 
#                                |
#                                | {Terminal 0}
#                       _________|_________
                        [0,1, "C", 10.*ut.fF], # Largee capacitor
#                       ___________________
#                                |
#                                | {Terminal 1}
#                   _____________|___________ 
#                  C                         |
#           Small  C                         |
#        inductor  C                         |
#                  C                         |
             [1,2, "L", 1*ut.pH],  #         |
             [1,2, "C", 10*ut.fF], #         |
#                   |                        |
#                   | {Terminal 2}           |
#                  _|_                      _|_
#                 |\|/|                    |\|/|
#                 |/|\| Jos. junct         |/|\| Jos. junct
             [0,2, "J", 0.3*ut.K],    [0,1, "J", 0.31*ut.K],
             [0,2, "C", 0.5*ut.fF],   [0,1, "C", 0.51*ut.fF]
#                   |                        |
#                   |________________________|
#                             |
#                             | {Terminal 0}
#                           __|__                 
                 ]

# Define mutual inductances (none)
MyMutualInd = []

# Describe number of states for each coordinate
N_of_states = [200,1]
             # 2*200+1 states for cyclic coordinates on the island
             # 1 coordinate in small inductor.
             # coordinate is considered classical (needing only
             # 1 quantum state) if quantum phase fluctuations
             # defined by impedance Z are small.




"""
 ___________________
| START CALCULATION |
|___________________|
| Do not modify this section
| unless you know what you
| are doing.
|___________________ 
"""
EasySQM = SQM.EasyModel(SMALL, MyCircuitCode, MyMutualInd)

EasySQM.DISPLAY = True          # display progress messages on screen
ENFORCE_INT_IN_Z_CYC = True     # Z_cyc matrix will be forced to be integer

# Calculate initial circuit parameters
EasySQM.buildCircuit()
# At this ponint user can also set user-calculated capacitance
# and inductance matrices by directly assigning
# "EasySQM.SQM.CircBuild.Circuit.C" and  "EasySQM.SQM.CircBuild.Circuit.L"

# Perform coordinate transformations
EasySQM.transformCoordinates()

# summarize information about coordinates
EasySQM.summarizeCoord()

# set number of quantum states for each variable
EasySQM.set_N_of_states(N_of_states)
print("    CALCULATING OPERATORS ...")
# get and print total number of quantum states used to model circuit
print("total N_of_states = "+str(EasySQM.SQM.ModelCalc.getTotal_N_of_states()))

# Generate quantum operators
EasySQM.setOperators()


"""
 ____________________________________
|You can freely modify the code below
|____________________________________
"""


"""
 __________________________
| SOLVE EIGENVALUE PROBLEM |
| FOR VARIOUS FLUXES AND   |
| CHARGES                  |
|__________________________|
"""

print("    GENERATING SPECTRA ...")

# define set of flux biases at which to model the circuit
PPhi = np.linspace( 0*ut.Phi_0, 1*ut.Phi_0,  21 )

# Create empty array that will contain eigenvalues
Hamiltonian_size = EasySQM.SQM.ModelCalc.getTotal_N_of_states()
EigVals = np.matrix( np.zeros((Hamiltonian_size,0)) )
# Rows -- eigenenergies. Columns -- iterations of charge bias Q

# Stop displaying progress on screen. Change to True if want to see progress
EasySQM.DISPLAY = False

for Phi in PPhi:                    # scan over flux biases
    
    # charge bias code
    ChargeBiasCode = [ [1, 1.1*ut.e] ]    # put 1.1e charge bias in node 1 
    # flux bias code (vary flux)
    FluxBiasCode = [ [1,2, Phi] ]
        # bias small inductor connecting nodes 1 and 2
    
    # calculate transformed flux and charge bias vectors
    EasySQM.transformBiasVectors(FluxBiasCode,ChargeBiasCode)
    
    # Generate Hamiltonian of the circuit. Save into 
    #                          EasySQM.ModelCalculator.Operators.Hamiltonian
    EasySQM.setHamiltonian()

    # Get eigenvalues (D) and col-eigenvectors (V) of the Hamiltonian
    # use "H" option to diagonalize Hamiltonian as Hermitian       
    D,V = SQM.eigsort(EasySQM.getHamiltonian(), "H")
   
    # Add the eigenenergies to list EigVals
    EigVals = np.hstack((EigVals,D.real))       
    
    # print current step summary
    print( "Phi/Phi_0 = "\
           +str(round(Phi/ut.Phi_0, 2)) \
           +";    E_10/GHz = "\
           +str(round( (D[1,0]-D[0,0]).real/ut.GHz, 3))\
         )

# Resume displaying progress on screen
EasySQM.DISPLAY = True
EasySQM.printMessage("Spectra calculated.")    






"""
 ____________________
| PLOT BANDSTRUCTURE |
|____________________|
"""
# Number of states to plot
N_states_to_plot = 4
# plot bandstructure
plt.plot(PPhi/ut.Phi_0, 
         (EigVals[0:N_states_to_plot,:].T-EigVals[0,0].T)/ut.K, 'o') 
plt.xlabel("$\Phi/\Phi_0$",  size = 15)
plt.ylabel("$E_n/$K", size = 15)
plt.show()
