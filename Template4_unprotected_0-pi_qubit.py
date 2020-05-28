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

# Define circuit parameters
L = 1.*ut.K|ut.E_to_L()
C = 2.*ut.K|ut.E_to_C("2e") # 2K E_C (assuming 2-electron charging energy)
E_J = 2.*ut.K
C_shunt = 70.*ut.fF

MyCircuitCode = [\
#                  _______
#                     |
#                     | {Terminal 0}
#         ____________|____________
#        |                         |
#       C                         C
#       C                         C
#       C                         C
#       C                         C
  [0,2, "L", L],           [0,3, "L", L], 
#        |                        |                
#        | {Terminal 2}           | {Terminal 3}                
#       _|_                      _|_               
#      |\|/|                    |\|/|              
#      |/|\|                    |/|\|                
  [2,1, "J", E_J],         [3,1, "J", E_J], 
  [2,1, "C", C],           [3,1, "C", C], 
#        |                        |                
#        |________________________|
#                    |                
#                    | {Terminal 1}
#                    |                
#         ___________|___________
              [0,1, "C", C_shunt]
#         _______________________
#                    |
#                    | {Terminal 0}
#                 ___|___
                ]

# Define mutual inductances (make small mutual between inductors)
MyMutualInd = [ [0,2,   0,3,   1*ut.pH] ]

# Describe number of states for each coordinate
N_of_states = [10,5,5]
             # 2*10+1=21 states for charge island (Terminal 1)
             # and 5x5 states for oscillator coordinates (Terminals 2,3)




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
PPhi = np.linspace( 0.4*ut.Phi_0, 0.6*ut.Phi_0,  41 )

# Create empty array that will contain eigenvalues
Hamiltonian_size = EasySQM.SQM.ModelCalc.getTotal_N_of_states()
EigVals = np.matrix( np.zeros((Hamiltonian_size,0)) )
# Rows -- eigenenergies. Columns -- iterations of charge bias Q

# Stop displaying progress on screen. Change to True if want to see progress
EasySQM.DISPLAY = False

for Phi in PPhi:                    # scan over flux biases
    
    # charge bias code (none)
    ChargeBiasCode = []    # put 1.1e charge bias in node 1 
    # flux bias code (vary flux -- divide between two inductors)
    FluxBiasCode = [ [0,2, Phi/2], [0,3, -Phi/2] ]
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
         (EigVals[0:N_states_to_plot,:].T-EigVals[0,0].T)/ut.GHz) 
plt.xlabel("$Phi/Phi_0$",  size = 15)
plt.ylabel("$E_n/$GHz", size = 15)
plt.show()
