# SuperQuantPackage
Modeling and Analysis of Superconducting Quantum Circuits

This software allows modeling of quantum superconducting circuits.

  Any user of the SuperQuant package must agree with the latest version
of the license agreement published in
github.com/andreyklots/SuperQuantPackage
or in other repository created by author(s) of this code.

  Please look at QuickDemo.pdf to get overall impression.

  Please contact Andrey R Klots (andreyklots@gmail.com)
if you have questions or comments.

 ----------------
|  REQUIREMENTS  |
 ----------------

For its operation the software requires only two modules:
  numpy -- necessary
  matplotlib -- optional: if one wants to plot results of modeling
  
  The software works on python versions 2.7 and above.
  It was not tested on earlier versions of python, but is likely
  to work. It was also tested on mobile platforms (such as
  Python3 on Android OS).

 -------------
|  MAIN IDEA  |
 -------------
Detailed theoretical considerations behind the code are
Presented in

               Theory.pdf

The main idea of the algorithm is:

1.  User determines the structure of superconducting circuit
    consisting of Josephson junctions (JJ), capacitors (C)
    and inductors (L). User defines nodes of the circuit
    and for each element (JJ, L or C) defines which two
    nodes this element connects and what is Josephson
    energy, inductance or capacitance of that element.

2.  Algorithm performs linear coordinate transformations
    that separate cyclic coordinates (corresponding
    to quantized net charge) from oscillatory ones
    (corresponding to oscillator-like degrees of freedom
    that do not obey charge quantization: i.e. in some way are
    connected to ground via inductors). Algorithm also contains
    a method "SuperQuantModel.EasyModel.summarizeCoord()"
    that allows to compare new (i.e. transformed)
    coordinates to original ones (i.e. phase on each node of a 
    circuit).

3.  User defines how many quantum states one needs to model each
    cyclic or oscillatory coordinate.
    If for a particular cyclic coordinate user sets M states,
    then the corresponding charge operator will be computed in
    the basis of |-M>, |-M+1>, ..., |+M> charge states.
    If for a particular oscillatory coordinate user sets N states
    then the corresponding flux operators will be computed in the
    basis of first N harmonic oscillator eigenfunctions. Each oscillatory
    coordinate is characterized by frequency and impedance. These
    values determine the width of a harmonic oscillator function.
    One can often intuitively guess how many states should be needed
    for each transformed coordinate by comparing them to original ones
    using method ".summarizeCoord()" that shows the linear transformation
    between transformed coordinates and phases on circuit nodes.

4.  Algorithm calculates operators in the basis defined by the user.
    Basic operators are:
              a. charge number operators for each cyclic coordinate
              b. photon number operators for each oscillator coordinate
              c. Josephson junction operators that allow interactions
                 between different cyclic and oscillatory coordinates.
              All cyclic and oscillatory coordinates are independent
              and corresponding operators are diagonal.
              Josephson junction operators have off-diagonal elements
              that allow interactions between charges and photons.

5.   User defines charge offsets on nodes and flux offsets in inductors.

6.   Based on previously calculated operators, charge and flux offsets
     the algorithm constructs the Hamiltonian of the circuit.

7.   Algorithm diagonalizes Hamiltonian to find eigenenergies
     and eigenstates of the circuit.



 ---------
|  FILES  |
 ---------

One can find following files in the folder:

0. LICENSE.txt
      License agreement. Any user or developer acknowledges agreement 
      written in LICENSE.txt before working with any related code.

1. SuperQuantModel.py
      This is the main module that describes four classes:
      a. "CircuitBuilder": this class allows to input a user-friendly
         code describing the circuit (see Documentation) and from it
         extract/compute values of Josephson energies, inductances,
         capacitances, parasitic capacitances, mutual inductances, 
         offset flux and charge biases, and structure of a circuit in forms
         of math-friendly vectors and matrices.
     b. "CircuitParameters": this class contains lists/vectors and matrices
        needed mathematically to define a circuit. These parameters can be 
        either set directly or computed using methods of "CircuitBuilder".
     c. "ModelCalculator": this is the main class: it contains methods for:
        (i)  performing coordinate transformations (separating cyclic degrees
             of freedom from oscillatory ones) based on circuit parameters
             and structure (transformations are stored in subclass ".Matrices").
        (ii) computing operators and building Hamiltonian in terms of
             transformed coordinates (based on subclass ".Operators").
     d. "EasyModel": class that contains methods that automatize some
        calculations performed by "CircuitBuilder" and "ModelCalculator"
        classes. Basically User feeds the circuit to the "EasyModel" 
        object and as a return gets calculated coordinate transformations
        and quantum operators. Without this class one needs to sequentially
       call methods from "CircuitBuilder" and
       "ModelCalculator" classes. Class "EasyModel" just simplifies
        the code.

2.  Documentation_SuperQuantModel.txt
       Detailed documentation on SuperQuantModel.py

3.  MicrowaveUnits.py
       Module that contains units commonly used in microwave electronics.
       The code originally works in system of units
                 2e = \h_bar = \Phi_0/2\pi = 1
       leaving only one arbitrary unit: energy. MicrowaveUnits.py sets 
       energy units to Kelvins: ( 1Kelvin = 1 ) and computes all other units
       with respect to this basis. This allows easy conversion of units.

4.   Documentation_MicrowaveUnits.txt
       Documentation for MicrowaveUnits.py

5.   Template_for_SimpleCircuitCalculator.py
       This is a template code that calculates a spectrum for a simple
       superconducting circuit. The template consists of two parts:
       (I)   Part that is meant to be modified by user:
             it has 
             a.  codes describing the circuit,
             b.  how to vary bias charges and fluxes: we introduce a
                 parameter variable t that goes from t=0 to t=1. All
                 charge and flux biases are defined as
                 user-defined functions of t.
             c.  what values to display during calculations and what
                 spectra to plot after calculation is complete.
       (II)  Part that is meant to remain intact:
             it has
             a. code that performs coordinate transformations and computes
                quantum operators. At each step of the calculation
                intermediate results are displayed on the screen. Some of
                these values include frequencies of different circuit modes
                and coordinate transformation matrix.
             b. code that sweeps variable t from 0 to 1 and for each
                corresponding set of t-dependent charge/flux biases
                computes the Hamiltonian, diagonalizes the Hamiltonian and
                records the set of all eigenenergies.
           
6.   Files named "Template1*.py", "Template2*.py", "Template3*.py",...
         Less user-friendly templates that show in greater detail,
         step by step how to use SuperQuantModel module and class
         "EasyModel" in order to solve eigenvalue problems for different
         circuits. These templates allow to better understand the workings
         of Template_for_SimpleCircuitCalculator.py.
       In particular, file named "Template0_Transmon_Cavity_Hamiltonian.py"
         does not perform calculation of eigenstates, unlike other templates.
         Instead it in detail describes a Hamiltonian written in transformed
         system of coordinates. This is a good example to study in the
         beginning.


 ---------------------------------------
| ERRORS, CONTRIBUTIONS AND SUGGESTIONS |
 ---------------------------------------

The current code is a beta-version, so it may contain some errors, bugs
 or discrepancies. Thus, reports of errors or bugs from users are
 greatly appreciated.
