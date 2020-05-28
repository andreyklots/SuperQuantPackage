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

# Circuit parameters
# transmon:
E_J     = 0.3*ut.K
E_J_coupling = 2*ut.K
C_J     = 0.5*ut.fF
C_shunt = 60*ut.fF
L_wire = 0.01*ut.pH
C_wire = 10*ut.fF
# resonator
C_coupling = 5*ut.fF
f_res = 5*ut.GHz    # desired target resonator frequency
Z     = 50*ut.Ohm   # resonator impedance
L_res = Z/f_res     # calculate resonator inductance
C_res = 1/(Z*f_res) # calculate resonator capacitance

# code for circuit structure
# pseudographic drawing of the circuit simply helps to recognize the circuit
MyCircuitCode = [\
#             _
#            ___  
#           ______
#             |
#             |{Terminal 0: Ground}
#             |
#             |
#  ___________|___________ Large transmon capacitor                 
   [0,2, "C", C_shunt],         [2,4, "J", E_J_coupling],
#  _______________________                 
#             |_____________________\/_____________
#             |  {Terminal 2}       /\            |
#             C                                   |
#    Small    C                                   |  Coupling
# inductor:   C  wire                             |  capacitor
#   Z=1E-3    C                              _____|_____
    [1,2, "L", L_wire],                 [3,4, "C", C_coupling],
    [1,2, "C", C_wire], #                   ___________
#            |                                   | {Terminal 3}
#            | {Terminal 1}    ..................|......................... 
#            |                :           _______|_________               :
#            |                :          C                |  Resonator    :
#            |                :          C                |               :
#           _|_               :          C                |               :
#          |\|/|              :          C                |               :
#          |/|\| Jos. junct   :          C          ______|_____          :
    [0,1, "J", E_J],            [0,3, "L", L_res],     [0,3, "C", C_res],#: 
    [0,1, "C", C_J]      #    :          |          ____________          :
#            |                :          |________________|               :
#            | {Terminal 0}   :                   | {Terminal 0}          :
#            |                :...................|........................
#          __|__                                __|__
                 ]

# Define mutual inductances (none)
MyMutualInd = []  # no mutual inductance



"""
 ___________________
| START CALCULATION |
|___________________|

"""
EasySQM = SQM.EasyModel(SMALL, MyCircuitCode, MyMutualInd)

EasySQM.DISPLAY = True          # display progress messages on screen
ENFORCE_INT_IN_Z_CYC = True     # Z_cyc matrix will be forced to be integer

# Calculate initial circuit parameters based on circuit code above
EasySQM.buildCircuit()
# At this ponint user can also set user-calculated capacitance
# and inductance matrices by directly assigning
# "EasySQM.SQM.CircBuild.Circuit.C" and  "EasySQM.SQM.CircBuild.Circuit.L"

# Perform coordinate transformations
EasySQM.transformCoordinates()

# summarize information about coordinates
EasySQM.summarizeCoord()

"""
 ____________________ 
| EXTRACT PARAMETERS |
|                    |
|      FROM          |
|                    |
|   CALCULATION      |
|____________________|

"""
# Code below displays Hamiltonian written in transformed coordinates

# Choose how tto round numbers when displaying
# Hamiltonian parameters
N_DIGITS_TO_ROUND = 4

# procure all information describing circuit and transformed coordinates
# N of cyclic coordinates
N_cyc = EasySQM.SQM.ModelCalc.Matrices.N_cyc
# N of oscillatory coordinates
N_osc = EasySQM.SQM.ModelCalc.Matrices.N_osc 
# total N of coordinates
N_tot = N_cyc + N_osc                        
# list of Josephson energies
E_jos = EasySQM.SQM.ModelCalc.Circuit.E_jos.round(N_DIGITS_TO_ROUND)  
# N of Josephson junctions
N_jos = len(E_jos)                           
# get parameters of the Hamiltonian
# inv_c_cyc -- inverce capacitance matrix for cyclic coordinates
# sqrt_omega_osc -- square-root oscillator frequencies
# M_III -- coeffitients in front of phases inside Josephson junction cosines
(inv_c_cyc, sqrt_omega_osc, M_III)\
              = EasySQM.SQM.ModelCalc.getHamiltonianParameters3()
inv_c_cyc = inv_c_cyc.round(N_DIGITS_TO_ROUND)
sqrt_omega_osc = sqrt_omega_osc.round(N_DIGITS_TO_ROUND)
M_III = M_III.round(N_DIGITS_TO_ROUND)
# oscillator frequencies
omega_osc = np.square(sqrt_omega_osc).round(N_DIGITS_TO_ROUND)
# coordinate transform matrix
T = EasySQM.SQM.ModelCalc.getTransform3Matrix().round(N_DIGITS_TO_ROUND)



"""
 ____________________ 
| DISPLAY CALCULATED |
|                    |
|       CIRCUIT      |
|                    |
|      PARAMETERS    |
|____________________|

"""


Hamiltonian_expression = """

  H = H    + H    + H
       cyc    osc    jos
"""
Hamiltonian_description = """
Hamiltonian is a sum of
1. cyclic part H   representing variables
                cyc
                
   obeying charge quantization.
2. oscillatory part H    representing 
                     osc
                     
   oscillator-like coordinates associated
   with inductors.
3. Josephson part H    representing
                   jos
                   
   interactions between the modes via
   Josephson junctions.

               tot
  Circuit has N   = """+str(N_tot)+""" terminals.
"""
Hamiltonian_cyc = """
            cyc
           N  
         ====                    
         \     1     cyc  -1
 H    =   >   ---  (c   )     q  q
  cyc    /     2          jk    j  k
         ====
         j,k=1

"""
Hamiltonian_cyc_description = """
  Indexes j,k go over cyclic coordinates.

               cyc
  Circuit has N    = """+str(N_cyc)+""" cyclic
coordinates that follow charge quantization.
  q-variables represent charge operators.
  Inverse capacitance matrix for
cyclic coordinates is

   cyc -1
 (c   )   =
"""+str(inv_c_cyc)+"""
"""

Hamiltonian_osc = """
           tot
          N
         ====                  
         \      1       2     2
 H    =   >    --- w ( q  +  v  )
  osc    /      2   j   j     j
         ====
          tot  cyc
       j=N   -N    
"""
Hamiltonian_osc_description = """
  Index j goes over oscillatory coordinates
of the circuit. The circuit has

 osc
N   = """+str(N_osc)+""" oscillator-like coordinates.
   Here v- and q- variables are phase and charge
operators respectively. Oscillator frequency modes are

w  =
 j
"""+str(omega_osc.T)+"""

"""
Hamiltonian_jos = """

            N            _    tot          _   
             jos        |    N              |              
           ====         |  ====             | 
           \            |  \     (III)      |  
 H    = --  >    E  cos |   >   M       v   |
  jos      /      J     |  /     Jj      j  |
           ====         |  ====             |        
           J=1          |_  j=1            _|
"""
Hamiltonian_jos_description = """
  Index J runs over """+str(N_jos)+""" Josephson
junctions. Josephson junction energies are

E  =
 J
"""+str(E_jos.T)+"""

  Index j runs over all """+str(int(N_tot))+""" coordinates.
Indexes in front of transformed (after transformation III)
phases v are

 (III)
M      =
 Jj
"""+str(M_III)+"""

  Note that first """+str(N_cyc)+""" cols are
indexes in front of cyclic coordinates and are integer.
Remaining """+str(N_osc)+""" indexes correspond to
oscillator-like coordinates and do not obey charge
quantization and hence do not need to be integer.

"""

Coordinate_transform_description1 = """
Transformed coordinates v  are related
                         j
to original coordinates phi  as
                           j
""" 
Coordinate_transform = """
 _    _     _                                _     _    _
|  cyc |   |  cyc    cyc      cyc    cyc      |   |      |
| v    |   | T      T    ... T      T     ... |   |phi   |
|  1   |   |  1,1    1,2      1,j    1,j+1    |   |   1  |
|      |   |                                  |   |      |
|  cyc |   |   cyc   cyc      cyc    cyc      |   |      |
| v    |   |  T     T    ... T      T     ... |   |phi   |
|  2   |   |   2,1   2,2      2,j    2,j+1    |   |   2  |
| .    |   |  .     .    .   .      .         |   | .    |
| .    |   |  .     .     .  .      .         |   | .    |
| .    |   |  .     .      . .      .         |   | .    |
|  osc | = |   osc   osc      osc    osc      | x |      |
| v    |   |  T     T    ... T      T     ... |   |phi   |
|  1   |   |   1,1   1,2      1,j    1,j+1    |   |   j  |
|      |   |                                  |   |      |
|  osc |   |   osc   osc      osc    osc      |   |      |
| v    |   |  T     T    ... T      T     ... |   |phi   |
|  2   |   |   2,1   2,2      2,j    2,j+1    |   |   j+1|
| .    |   |  .     .    .   .      .     .   |   | .    |
| .    |   |  .     .     .  .      .      .  |   | .    |
| .    |   |  .     .      . .      .       . |   | .    |
|_    _|   |_                                _|   |_    _|
"""
Coordinate_transform_description2 = """
Coordinate transformation matrix

 (III)
T      =
 jk
 
"""+str(T)+"""

  First """+str(N_cyc)+""" rows provide extract cyclic
coordinates and remaining """+str(N_osc)+""" rows extract
oscillator-like coordinates from original phase
coordinates \phi .
                j
"""
print("Write down the Hamiltonan:")                 
print(Hamiltonian_expression)
print(Hamiltonian_description)
print("Cyclic part of Hamiltonian:")
print(Hamiltonian_cyc)
print(Hamiltonian_cyc_description)
print("Oscillatory part of Hamiltonian:")
print(Hamiltonian_osc)
print(Hamiltonian_osc_description)
print("Josephson part of the Hamiltonian:")
print(Hamiltonian_jos)
print(Hamiltonian_jos_description)

print("Coordinate transformation:")
print(Coordinate_transform_description1)
print(Coordinate_transform)
print(Coordinate_transform_description2)
