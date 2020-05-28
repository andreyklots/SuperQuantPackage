"""
  ___________________________________________________________________________
 |                                                                           |
 |   This is a template/tutorial for modeling superconducting circuits.      |
 |                                                                           |
 |   Please modify only the code between markers "START USER-DEFINED MODEL   |
 |   PARAMETERS" and "END USER-DEFINED MODEL PARAMETERS".                    |
 |                                                                           |
 |   First one has to define parameters of a circuit (circuit structure,     |
 |   values of inductors, capacitors, Josephson junctions, mutual            |
 |   inductances, etc.                                                       |
 |   Then one has to define how flux and charge biases in the circuit        |
 |   should vary and what calculated values to display and plot.             |
 |   The code does all the computational code for the user after that.       |
 |                                                                           |
 |   How automatic calculation is performed?                                 |
 |    -- We start with writing the Hamiltonian of a circuit in terms         |
 |       of phases and charges in the nodes of the user-defined circuit      |
 |    -- Then we perform a special linear coordinate transformation that     |
 |       separates cyclic coordinates (ones that obey charge quantization)   |
 |       from oscillatory ones (coordinates that behave like oscillators).   |
 |       This allows to untangle oscillator modes of the circuit and make    |
 |       them independent (at least explicitly) to simplify calculations.    |
 |    -- After that each transformed coordinate is quantized. User has to    |
 |       determine how many quantum states should be taken into account for  |
 |       each transformed coordinate.                                        |
 |    -- Eventually Hamiltonian is written in a matrix form and diagonalized |
 |                                                                           |
 |   Please carefully read the comments to understand how to set up a        |
 |   and calculate its eigenstates                                           |
 |___________________________________________________________________________|
 """

# Import important modules. In fact only modules that are absolutely
# necessary are SuperQuantModel that does all the work and numpy module
# SuperQuantModel only uses numpy and nothing else.
# MicrowaveUnits is just needed for more convenient work with standard units
# commonly used in microwave physics.
# matplotlib is only needed for plotting graphs
import SuperQuantModel as SQM
import numpy as np
import matplotlib.pyplot as plt
# module that contains most common units
import MicrowaveUnits as ut


###############################################################################
###############################################################################
#"""""""""""""""""""""""""""""""""""""""""""""
#"""""""""""""""""""""""""""""""""""""""""""""
#""""                                     """"
#""""             S T A R T               """"
#""""                                     """"
#""""      U S E R - D E F I N E D        """"
#""""                                     """"
#""""   M O D E L   P A R A M E T E R S   """"
#""""                                     """"
#"""""""""""""""""""""""""""""""""""""""""""""
#"""""""""""""""""""""""""""""""""""""""""""""
###############################################################################
###############################################################################


# Define small parameter (SMALL). It serves three functions:
#    (i) this parameter is used to decide which frequency modes
#        are considered effectively zero
#   (ii) this parameter is used to determine collinearity of 
#        vectors when doing linear-algebraic transformations.
#  (iii) when calculating a circuit's capacitance matrix this 
#        parameter is used as value of a smallest capacitance.
#        All nodes of the circuit by default will fave mutual
#        capacitance of the order of this SMALL value.
SMALL = 1E-4


" ========================== "
"|| D E S C R I P T I O N  ||"
"||          O F           ||"
"||     C I R Q U I T      ||"
" ========================== "

# Here starts code describing the circuit
# Please see details of the code syntax in documentation.
# It is advised to use MicrowaveUnits module to describe circuit parameters
# Please also see documentation for MicrowaveUnits module.
# The code can be used without MicrowaveUnits module. But then keep in mind:
#    By default the code works in units where h_bar = 2e = phi_0/2pi = 1.
#    Capacitor charging energy is defined as E_C = (2e)^2/2C = 1/2C
#    (2-electron charging energy). Inductor energy is defined as 
#    E_L = (Phi_0/2pi)^2/L = 1/L.
# Right below is the code describing the circuit.
# Pseudographic drawing of the circuit simply helps to recognize its structure.

# Define parameters of elements in the circuit
L = 1.*ut.K|ut.E_to_L()     # inductance based on inductor energy
C = 2.*ut.K|ut.E_to_C("2e") # capacitance based on 2-electron charging energy
E_J = 2.*ut.K               # Josephson junction energy
C_shunt = 100.*ut.fF        # shunt capacitance


# Define circuit structure
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





# Now we define names of inductors. Inductors are defined by numbers of nodes
# that they connect (written in [...]-brackets.
#                                   Order should e same as in MyCircuitCode).
# Below we assign names (written in quotation marks "...") to those inductors
# that we are going to bias with flux biases
Inductors = { \
               "L1":[0,2],  # inductor connecting nodes 0,2 is now called "L1"
               "L2":[0,3]
            }


# In the same way, assign names to charge islands that we want to bias
# in future. assign name to node numbers that we want to use:
ChargeIslands = { \
                   "I1": 1, # node (island) 1 is now called "I1"
                   "I2": 2
                }

# Define mutual inductances. Syntax is as follows: 
#   Inductors[NAME_1] + Inductors[NAME_2] + [VALUE]. This entry in the
# MyMutualInd list will define mutial inductance VALUE between inductors
# with names NAME_1 and NAME_2.
MyMutualInd = [ \
               Inductors["L1"] + Inductors["L2"] + [1*ut.pH]
               #, add new lines with comma separating code lines
              ]





# Describe number of quantum states that we want to use
# for each coordinate (cyclic or linear) during computation.
# Total number of coordinates (N_tot) is number of nodes in the circuit.
# Number of cyclic coordinates (N_cyc) is a number of isolated charge 
# clusters with quantized charges (i.e. number of segments of the circuit
# that are not connected to ground via any continuous path of inductors).
# Number of linear (oscillatory) coordinates is N_tot - N_cyc.
# The algorithm determines N_cyc and N_osc automatically. Check if these
# numbers make sense upon running the code.
# In the N_of_states array cyclic coordinates go first and oscillatory -- last
# Transformed coordinates with lower frequency modes go first.
# The software displays how transformed coordinates are related to
# original ones in the SUMMARY section once the code is executed.
N_of_states = [\
                 # cyclic coordinates:
                 10,
                 # oscillatory coordinates:
                 5,5
               ]



" ============================= "
"||    S E T T I N G    U P   ||"
"|| F L U X   &   C H A R G E ||"
"||           B I A S         ||"
"||         S W E E P S       ||"
" ============================= "

# Now the system is set up.
# It is time to set up a numerical experiment
# First, define variable t. All flux and charge biases 
# will change as a function of t
t = np.linspace(0,1, 41);


# Define t-dependent flux biases for each inductor as
#     INDUCTOR_NAME: FLUX_BIAS_VALUE as a function of t
# NOTE: FLUX_BIAS_VALUE can also be constant
FluxBiases =   { \
                 "L1":  t*ut.Phi_0/2,
                 "L2": -t*ut.Phi_0/2
               }

# Define t-dependent (or constant) charge biases in the same way
ChargeBiases = { "I1": 0 }



" ===================== "
"||    D E F I N E    ||"
"|| D I S P L A Y E D ||"
"||    V A L U E S    ||"
" ===================== "

# define function that plots the resulting spectra
def plottingFunc(t, EigVals):
    #    In EigVals array (array of energy eigenvalues) 
    # row number corresponds to eigenstate index
    # col number corresponds to each value of t-variable 
    N_bands_to_plot = 5  # number of bands (eigenenstates) to plot
    plt.plot(
             t, EigVals[0:N_bands_to_plot,:].T-np.min(EigVals)
            )   
    plt.xlabel("t-variable (arb.units.)", size = 15)  # set labels
    plt.ylabel("Energy (K)", size = 15)
    plt.show()



# define function for displaying values during computation process
# for each value of t we can print current eigenenergies
def showProgressMsg(t_current, EigVals_at_current_t):
    # t_current -- current element in t during each calculation
    # EigVals_at_current_t -- col-vector of current eigenenergies
    N_bands_to_show = 5
    print( "t = " + str(round(_t,4)))
    print("E_n/K:")
    print( EigVals_at_current_t[0:N_bands_to_show,0]/ut.K )


# define function for displaying chosen values before
# computation of (t-sweep) starts.
def showInitialMsg(EasySQM):
    # EasySQM -- object of class EasyModel containing all the model
    # properties
    print( "_____oOOOo_____" )
    print( "Coeffitients inside cos of JJs (after transform II):" )
    print( np.hstack([EasySQM.SQM.ModelCalc.Matrices.Z_cyc,
                      EasySQM.SQM.ModelCalc.Matrices.M_osc]) )
    print( "_____oOOOo_____" )


###############################################################################
###############################################################################
#"""""""""""""""""""""""""""""""""""""""""""
#"""""""""""""""""""""""""""""""""""""""""""
#"""                                     """
#"""                E N D                """
#"""                                     """
#"""      U S E R - D E F I N E D        """
#"""                                     """
#"""   M O D E L   P A R A M E T E R S   """
#"""                                     """
#"""""""""""""""""""""""""""""""""""""""""""
#"""""""""""""""""""""""""""""""""""""""""""
###############################################################################
###############################################################################






















##############################################################################
##############################################################################
##                                                                          ##
##                                D O   N O T                               ##
##                                                                          ##
##                            M O D I F Y   T H E                           ##
##                                                                          ##
##                             C O D E   B E L O W                          ##
##                                                                          ##
##                     Unless you know what you are doing                   ##
##                                                                          ##
##############################################################################
##############################################################################





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
 __________________________
| SOLVE EIGENVALUE PROBLEM |
| FOR VARIOUS FLUXES AND   |
| CHARGES                  |
|__________________________|
"""

showInitialMsg(EasySQM)

print("    SWEEPING PARAMETER t from 0 to 1 ...")


# Create empty array that will contain eigenvalues
Hamiltonian_size = EasySQM.SQM.ModelCalc.getTotal_N_of_states()
EigVals = np.matrix( np.zeros((Hamiltonian_size,0)) )
# Rows -- eigenenergies. Columns -- iterations of charge bias Q

# Stop displaying progress on screen. Change to True if want to see progress
EasySQM.DISPLAY = False


Step_in_t = 0                   # count elements in t
for _t in t:                    # scan over t
    
    
    # charge bias code
    # start with empty charge bias code
    ChargeBiasCode = []    
    # build Charge bias code from ChargeBiases dictionary for current Step_in_t
    for k in ChargeBiases.keys():
        ChargeBiasCode += [ \
                           [\
                           ChargeIslands[k], (ChargeBiases[k]+0*t)[Step_in_t]
                           ]    # adding 0*t converts constant number to array
                          ]
    FluxBiasCode = []    
    # build Flux bias code from FluxBiases dictionary for current Step_in_t
    for k in FluxBiases.keys():
        FluxBiasCode += [ \
                         [\
                           Inductors[k][0], Inductors[k][1],
                           (FluxBiases[k]+0*t)[Step_in_t]
                         ] 
                        ]


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
    showProgressMsg(t_current = _t, EigVals_at_current_t = D)
    
    Step_in_t += 1              # iterate the step              


# Resume displaying progress on screen
EasySQM.DISPLAY = True
EasySQM.printMessage("Spectra calculated.")    






"""
 ____________________
| PLOT BANDSTRUCTURE |
|____________________|
"""

# plot bandstructure
plottingFunc(t, EigVals)