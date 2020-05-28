"""
Copyright (c) 2018 and later, Andrey R. Klots. With supervision
                   of Lev B. Ioffe and with support of
                   United States Army Research Office (ARO).
All rights reserved

This file is a part of "SuperQuant" package.
Terms of use of this code are provided below (until long horizontal line)

Most up-to-date version of this agreement must always be used. Updates
Are available at github.com/andreyklots/SuperQuantPackage.

   Terms used in agreement:
        "modification": for purposes of this agreement, the word "modified"
      means any modification of the current code that includes not only
      changes in original python code, but also translation of
      this code (or parts of it) into another programming or script language
      or compilation of this code and translation into binary form.

   Use of this code:
        If this software (in whole or in part) is used as a part of a
      published research then need to be cited as "SuperQuant package by
      AR Klots" or "github.com/andreyklots/SuperQuantPackage".
        This software is free for use for purely research purposes.
      However, if research using this code or its derivatives is a part of
      development of a commercial or military product then direct permission
      of the designer and owner of this code is required. 
     
   Modification and redistribution:
        Any redistribution of this code, in whole or in part, modified or
      unmodified, must reproduce all the terms of the most up-to-date version
      of the current license agreement.
      Each distribution must explicitly refer to the original
      code.

   Integration and compilation:
        If this code (in whole or in part, modified or unmodified) is
      used as a part of another code or software package, then that code
      or software package must reproduce all the terms of the most
      up-to-date version of the current license agreement
      pertaining to those parts of the new code or
      software that were derived from the current code.
        In order to incorporate the current code (in whole or in part,
      modified or unmodified) as a part of software used for commercial
      or military purposes, then permission of the owner and designer
      of the code is required.
        Name of "SuperQuant" package or names of its contributors
      cannot be used to promote or endorse any other products, even
      if they are derived from this software.
     
   Liability:
        No person or organization that owns this code, participated
      or contributed to its creation, or holds the copyright, shall
      be liable for any damage, loss and all other undesirable and
      unexpected consequences, be they direct or indirect, or caused
      by or related to this software in any other way. Any user of
      the current code or its derivatives acknowledges that this 
      code may contain errors and no party involved in creation of
      this code is not liable for any undesirable consequences.
     
   Exemptions:
        As a supporting organization, ARO is exempt from limitations on
      use, modification, distribution and integration of this code.

--------------------------------------------------------------------------

"""


import numpy as np                    
from numpy import linalg as LA        


LOGO = """ 

    _________________________________________________
  _|                                         |       |
 /                                           |       |
 \_     S U P E R C O N D U C T I N G        X     __|__
 /                                           |     _____
 \_            Q U A N T U M                 X       |
 /                                           |    ___|___
 \_            C I R C U I T                 X   |       |
 /                                           |   X       x
 \_              M O D E L                   X   |_______|
   |     | |                                 |       |
   |_____| |________\/_______________________|_______|
         | |        /\   (c) Andrey Klots
         | |           

""" 
print("""
 ,---CCCCC---X--||------x---------,
 |       SuperQuant Circuit Model |   
_|_              (c) Andrey Klots_|_
""")























# _____________________#    
#! USEFUL FUNCTIONS   !#
#!____________________!#

# same as dir-function, but does not 
#return methods/properties starting with "_"
def dirv(Obj):
    return [ a for a in dir(Obj) if a[0]!='_' ]

# Diagonalize matrix A and sort from smaller to larger eigenvalues
def eigsort(A, METHOD = None):
    if (METHOD==None)or(METHOD=="A")or(METHOD=="a")or(METHOD=="arbitrary"):
        E, u = LA .eig(A)
        #sort eigenvalues by their real part
        idx = np.argsort(E.real)
        E = E[idx]
        u = u[:,idx]                          
        # eigenvectors as matrix type                                              
        return [np.matrix(E).T, np.matrix(u)]
    if (METHOD=="R")or(METHOD=="r")or(METHOD=="real"):
        # extract real-symmetric part of the matrix 
        _A = 0.5*(A+A.H).real
        E, u = LA .eigh(_A)
        return [np.matrix(E).T, np.matrix(u)]
    if (METHOD=="H")or(METHOD=="h")or(METHOD=="hermit")or(METHOD=="Hermit"):
        # extract hermitian part of the matrix 
        _A = 0.5*(A+A.H)
        E, u = LA .eigh(_A)
        return [np.matrix(E).T, np.matrix(u)]
























""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""                                                        """ 
"""                      C L A S S                         """
"""                                                        """ 
"""      T H A T     D E S C R I B E S    B A S I C        """
"""                                                        """
"""        C I R C U I T     P A R A M E T E R S           """
"""                                                        """ 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

class CircuitParameters:
    # show properties of the class
    def __repr__(self):
        return str(dirv(self))
    # square root of inverse Maxwell capacitance matrix $C^{-1/2}$
    Inv_sqrt_C = np.matrix([None])
    # inductor connectivity matrix. Dimensions: (N of inductors) x (N of nodes)
    G = np.matrix([None])    
    # inverse inductance matrix. Dimensions:(N of inductors) x (N of inductors)                       
    Inv_L = np.matrix([None])
    # Josephson connectivity matrix. Same as Inductance 
    # capacitance matrix, but for Josephson junctions
    H = np.matrix([None])               
    # array(vector) of Josephson energies
    E_jos = np.matrix([None])
    # vector of offset fluxes           
    Phi = np.matrix([None])
    # vector of offset charges
    Q = np.matrix([None])




















""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""                                                        """ 
"""                      C L A S S                         """
"""                                                        """ 
"""   T H A T     C A L C U L A T E S    C I R C U I T     """
"""                                                        """
"""   P A R A M E T E R S   A N D   O P E R A T O R S      """
"""                                                        """ 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


class ModelCalculator:
    
    # show properties of the class
    def __repr__(self):
        return str(dirv(self))
    

    SMALL = 1E-4
    
    ##############################################
    #                                            #
    #              C I R C U I T                 #
    #                                            #
    #   Class defining properties of a circuit   #
    #                                            #
    ##############################################
    

    # create a property in SuperCondModel
    Circuit = CircuitParameters()

    ##############################################
    #                                            #
    #             M A T R I C E S                #
    #                                            #
    #   Class and functions defining matrices    #
    #                                            #
    ##############################################

    class __Matrices:
        # show properties of the class
        def __repr__(self):
            return str(dirv(self))
        # unitary matrix for first transformation
        U = np.matrix([None])
        # array of frequency modes
        omega_eigenvalues = np.matrix([None])
        # $N^{cyc}$
        N_cyc = None
        # $N^{osc}$
        N_osc = None
        # matrix defining values inside cosines of Josephson potentials ($M
        # $)
        M = np.matrix([None])
        # cyclic part of M ($\mathbb{M}^{cyc}$)
        M_cyc = np.matrix([None])
        # oscillatory part of M ($\mathbb{M}^{osc}$)
        M_osc = np.matrix([None])
        # largest invertible submatrix of M_cyc ($\mathbb{m}^{cyc}$)
        m_cyc = np.matrix([None])
        # integer form of M_cyc after transformation ($\mathbb{Z}^{cyc}$). 
        # Z_cyc=M_cyc/m_cyc
        Z_cyc = np.matrix([None])
        # transformed flux biases
        F = np.matrix([None])
        # transformed charge biases
        K = np.matrix([None])
    # create a property in SuperCondModel
    Matrices = __Matrices()
    
   # names of functions below are self-evident: calculate different matrices
    def calc_U_omega_eigenvalues(self):
        Mat = self.Circuit.Inv_sqrt_C.H * self.Circuit.G.T \
              * self.Circuit.Inv_L \
              * self.Circuit.G * self.Circuit.Inv_sqrt_C
        [omega_eigenvalues_squared, self.Matrices.U] \
            = eigsort(Mat, METHOD = "real")
        self.Matrices.omega_eigenvalues = \
            np.sqrt(np.abs(omega_eigenvalues_squared))
    # find number of cyclic coordinates and from there get number of oscill
    # atory ones
    def calc_N_cyc_osc_tot(self):
        self.Matrices.N_cyc, self.Matrices.N_osc \
            = self.__calc_N_cyc_osc_tot( 
                    self.Matrices.omega_eigenvalues, self.SMALL )
    def calc_M(self):
        self.Matrices.M \
            = self.Circuit.H * self.Circuit.Inv_sqrt_C * self.Matrices.U
    def calc_M_cyc_osc(self):
        self.Matrices.M_cyc = self.Matrices.M[ :, 0:self.Matrices.N_cyc ] 
        self.Matrices.M_osc \
            = self.Matrices.M[ :, \
                        self.Matrices.N_cyc 
                            : (self.Matrices.N_cyc+self.Matrices.N_osc) ] 
    def calc_m_cyc(self):
        self.Matrices.m_cyc \
            = self.getMaxInvertibleSubmatrix( \
                    self.Matrices.M[:,0:self.Matrices.N_cyc], self.SMALL)
    # calculate Z-matrix. Choose if we want to force round Z-matrix
    def calc_Z_cyc(self, ENFORCE_INT):
        # By default ENFORCE_INT is True
        if ENFORCE_INT==None: ENFORCE_INT = True
        # choose if we want to leave only integer part
        if ENFORCE_INT:
            self.Matrices.Z_cyc \
                = ( self.Matrices.M_cyc*self.Matrices.m_cyc.I  
                   ).round().astype(int) \
                       if np.size(self.Matrices.m_cyc)>0 else np.matrix([[]])
        else:
            self.Matrices.Z_cyc  \
                = ( self.Matrices.M_cyc*self.Matrices.m_cyc.I ) if\
                            np.size(self.Matrices.m_cyc)>0 else np.matrix([[]])
    def calc_Transformed_Flux_Charge_Biases(self):
        self.Matrices.F, self.Matrices.K = self.__calc_Flux_Charge_Biases(
                    self.Circuit.Phi, self.Circuit.Q,\
                    self.Circuit.G, self.Circuit.H,\
                    self.Matrices.m_cyc, self.Matrices.U,\
                    self.Circuit.Inv_sqrt_C)


    # returns transform 1 matrix -- U^T*C^(1/2) -- transforms original coor
    # dinates to new
    def getTransform1Matrix(self):
        return self.Matrices.U.T*self.Circuit.Inv_sqrt_C.I
    # returns transform 2 matrix. Acts on original coordinates.
    def getTransform2Matrix(self):
        m_cyc = self.Matrices.m_cyc
        N_cyc = self.Matrices.N_cyc
        Mat = np.matrix(np.identity(N_cyc+self.Matrices.N_osc))
        Mat[0:N_cyc,0:N_cyc] = m_cyc
        return Mat * self.Matrices.U.T * self.Circuit.Inv_sqrt_C.I
    # returns transform 3 matrix
    def getTransform3Matrix(self):
        m_cyc = self.Matrices.m_cyc
        N_cyc = self.Matrices.N_cyc
        Mat = np.matrix(np.diag(
                    # diagonal matrix of omega eigenvalues
                    np.sqrt(np.array(self.Matrices.omega_eigenvalues.T)[0]) ))
        Mat[0:N_cyc,0:N_cyc] = m_cyc
        return Mat * self.Matrices.U.T * self.Circuit.Inv_sqrt_C.I
    # returns matrices for transformed (3) Hamiltonian: c^-1 and M^III
    def getHamiltonianParameters3(self):
        m_cyc = self.Matrices.m_cyc
        Z_cyc, M_osc = self.Matrices.Z_cyc, self.Matrices.M_osc
        # 1d array of sqrt of omega eigenvalues
        Sqrt_omega = np.array( np.sqrt(self.Matrices.omega_eigenvalues.T) )[0]
        N_cyc, N_osc = self.Matrices.N_cyc, self.Matrices.N_osc
        InvCyclicCapacitanceMatrix = m_cyc*m_cyc.T
        M_osc3 = M_osc*np.matrix(np.diag(  
                        np.reciprocal(Sqrt_omega[N_cyc:N_cyc+N_osc])    ))
        M_3 = np.hstack([Z_cyc, M_osc3])
        return (InvCyclicCapacitanceMatrix,Sqrt_omega[N_cyc:N_cyc+N_osc], M_3)

    ######################################################
    #                                                    #
    #                  O P E R A T O R S                 #
    #                                                    #
    #   Class and functions defining quantum operators   #
    #                                                    #
    ######################################################

    
    class __Operators:
        # show properties of the class
        def __repr__(self):
            return str(dirv(self))
        # list of unity matrices
        I_hat = [None]
        # list of number operators
        N_hat = [None]
        # list of exponential operators
        V_hat = [None]
        # list of numbers of states for each coordinate.
        N_of_states    = [None]
        # note that for cyclic coordinates number N_of_states=M corresponds
        #  to 2M+1 states
        Hamiltonian = np.matrix([None]);
    
    # create a property in SuperCondModel
    Operators = __Operators()
    
    # function that calculates operators. option DISPLAY_PROGRESS determine
    # s if user wants to display current progress of calculating V_hat matr
    # ices
    def obtainOperators(self, DISPLAY_PROGRESS = True):
        self.Operators.I_hat, self.Operators.N_hat, self.Operators.V_hat \
            = self.__obtainOperators(
                    self.Matrices.N_cyc, self.Matrices.N_osc,\
                    self.Operators.N_of_states,\
                    self.Matrices.Z_cyc.round().astype(int),\
                    self.Matrices.M_osc, self.Matrices.omega_eigenvalues,\
                    DISPLAY_PROGRESS)
    def obtainHamiltonian(self, DISPLAY_PROGRESS = True):
        self.Operators.Hamiltonian \
            = self.__obtainHamiltonian(
                    self.Matrices.N_cyc, self.Matrices.N_osc,\
                    self.Matrices.m_cyc, self.Matrices.omega_eigenvalues,\
                    self.Circuit.E_jos, self.Operators.I_hat,\
                    self.Operators.N_hat, self.Operators.V_hat,\
                    self.Matrices.F, self.Matrices.K,   DISPLAY_PROGRESS)
    # calculate total number of quantum states
    def getTotal_N_of_states(self):
        # this will be total number of states 
        Tot_N_states = 1
        # number of current coordinate
        coordNumber = 0
        for N in self.Operators.N_of_states[ \
                                        0:(self.Matrices.N_cyc
                                        +self.Matrices.N_osc) ]:
            Tot_N_states = Tot_N_states * ( (2*N+1) \
                            # multiply all N_of_states (depending if it is 
                            # cyclic or not)
                            if coordNumber<self.Matrices.N_cyc else N )
            coordNumber = coordNumber + 1
        return Tot_N_states

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#

          
    ######################################################
    #                                                    #
    #                  F U N C T I O N S                 #
    #                                                    #
    #         Define more complicated functions          #
    #                                                    #
    ######################################################
    


# ______________________________________________________________#    
#!  FUNCTION TO FIND LARGEST INVERTIBLE SUBMATRIX OF A MATRIX  !#
#!_____________________________________________________________!#
    # function that returns largest invertible submatrix of matrix A with p
    # recision of TOLERANCE
    def getMaxInvertibleSubmatrix(self, A, TOLERANCE=None):
        # note that this is an important function and hence is not hidden
        if TOLERANCE==None: TOLERANCE = 1E-10
        
        if np.size(A)==0:
            return []
        
        if np.max( np.abs( A ) ) <= 2*TOLERANCE:
            # if matrix is too small then return error
            return "ERROR: matrix elements too small"
        
        # get dimensions of A
        N = np.size( A[:,0] )
        M = np.size( A[0,:] )
        # make A1, a copy of A. Will work with A1
        A1 = A
        # check if matrix is horizontal. The algorithm finds linearly indep
        # endent rows. If horizontal, then will need to transpose
        if N < M:
            NeedTranspose = True
            A1 = A.T
            # get dimensions of A transposed
            N = np.size( A1[:,0] )
            M = np.size( A1[0,:] )
        else:
            NeedTranspose = False
        
        # initial max invertible submatrix
        a = np.matrix([])


        # sort rows of A1 from largest absolute value of row-vector to lowest
        # so that norm of invertible submatrix is as large as possible
        # make sure A2 is matrix (and not just 2D array). Will be needed fo
        # r next row
        A2 = np.matrix(A1)
            # make array of absolute values of rows of A
        A1_abs_values = np.array([ (row*row.H)[0,0].real for row in A2 ])
            # sort A1 in the order of decreasing absolute values of row-vectors
        idx = np.argsort(A1_abs_values)
        A1 = A1[idx]
        # now row-vectors are sorted in the order of descending row-vector abs
        
        # add one more non-zero row so that if all rows are zero, then erro
        # r does not occur when n+1'th row is extracted
        A1 = np.vstack( [A1, np.zeros([1,M])] )

        # number of non-zero row
        n_non_zero = 0
            # find first non-zero row
        while A1[n_non_zero,:]*A1[n_non_zero,:].H <= 2*TOLERANCE \
                                                and n_non_zero < N:
            n_non_zero = n_non_zero + 1
        if n_non_zero == N+1:
            # if all rows are zero, return error
            return "ERROR: matrix elements too small"

        else:
            # add first nonzero row to a -- our invertible submatrix
            a = A1[n_non_zero,:]
            # number of rows in array a. Initially 1 
            rows_in_a = 1
            
            # this will be array of indexes later used for inverse sorting 
            # (so that linearly-independent rows are in the same order in a
            #  as in A)
            idx_inverse = [ idx[n_non_zero] ]
            for n in range( n_non_zero+1, N ):
                if np.linalg.matrix_rank(\
                                            np.vstack( [a,A1[n,:]] 
                                        # check if n'th row is linearly ind
                                        # ependent from whatever is already
                                        #  in a
                                        ), TOLERANCE ) > rows_in_a:
                    # add linearly independent row
                    a = np.vstack( [a,A1[n,:]] )
                    # add index of the current row in A
                    idx_inverse = idx_inverse + [ idx[n] ]
                    # increase current rank of a
                    rows_in_a = rows_in_a + 1


        # sort indexes in the order as they appear in original array A
        idx_for_a_inversion = np.argsort(idx_inverse)
        # resort a in such a way that rows in a appear in the same order as
        #  they appear in A
        a = a[idx_for_a_inversion]
        
        
        if NeedTranspose:
            # if needed to transpose original matrix then transpose the res
            # ulting one
            a=a.T
        return a




# ___________________________________#    
#!  CALCULATE N_cyc, N_osc, N_tot   !#
#!__________________________________!#
    def __calc_N_cyc_osc_tot( self, omega_eigenvalues, SMALL ):
        # get total number of generalized coordinates
        N_tot = np.size( omega_eigenvalues )
        # initial number of cyclic coordinates
        N_cyc = 0
        for n in range( 0, N_tot ):
            # if mode frequency eigenvalue is small
            if omega_eigenvalues[n] < np.sqrt(SMALL):
                # then this is a cyclic coordinate and we add it to N_cyc
                N_cyc = N_cyc + 1
        # get N_osc
        N_osc = N_tot - N_cyc
        return [N_cyc, N_osc]
    
    
    

    
# _________________________________________________________#    
#!   FUNCTION TO GET TRANSFORMED FLUX AND CHARGE BIASES   !#
#!________________________________________________________!#
    def __calc_Flux_Charge_Biases(self, Phi, Q,    G, H, m_cyc, U, Inv_sqrt_C):
        # total number of degrees of freedom
        N_tot = np.size(Q)
        # number of cyclic degrees of freedom
        N_cyc = np.size(m_cyc[0,:]) if len(m_cyc)>0 else 0
        # make sure that number of inductors is > 0
        if np.size(G) > 0:
            # calculate equation F = H(G\Phi)
            F = H * ( G.I * Phi )
        else:
            # if no inductors in the circuit, then return empty matrix for 
            # offsets. Transformed offset charge vector is K = diag(m_cyc^T
            # ,1)*U^T*C^(-1/2) Q
            F = np.matrix([])

        # now calculate diag(m_cyc^T,1)
        diag_m_1 = np.matrix(np.identity(N_tot))
        # put m_cyc in the top left corner of the identity matrix
        diag_m_1[0:N_cyc,0:N_cyc] = m_cyc
        # get transformed vector of offset charges     
        K = diag_m_1.I.T * U.T * Inv_sqrt_C * Q

        return [F, K]    


# ________________________________________#    
#!   FUNCTION TO GET QUANTUM OPERATORS   !#
#!_______________________________________!#
    def __obtainOperators( self, N_cyc, N_osc, N_of_states,\
                                      Z_cyc, M_osc,   omega_eigenvalues,\
                                      DISPLAY_PROGRESS ):
        # if DISPLAY_PROGRESS is not set, then it is by default False
        if type(DISPLAY_PROGRESS)==None: DISPLAY_PROGRESS = False

        # get total number of coordinates
        N_tot = N_cyc + N_osc

        # list of identity operators
        I = [0 for _i in range(N_tot)]
        # list of number operators
        N = [0 for _i in range(N_tot)]
        # initialize an empty (N_Jos x N_tot) list for V-operators --  2D l
        # ist of V-matrices $V^{(J,j)}$. J -- first index. j -- second inde
        # x
        V = [[0 for _i in range(N_tot)] for _j in range(np.size(M_osc[:,0]))]

                                                                    
                                                                    
    #  NOTE: here for simplicity we omit "hats" so instead of N_hat,
    #        V_hat we write N, V,...
        
        # Generate identity and number operators
        for n_coord in range(0, N_tot):
            # N of states for n-th coordinate
            _ns = N_of_states[n_coord]
            
            # create identity operators
            if n_coord < N_cyc:
                # create identity matrix for cyclic coordinate
                I[n_coord] = np.matrix(np.identity( 2*_ns+1 ))
                # create charge number operator
                N[n_coord] = np.matrix(np.diag(range(-_ns,_ns+1)))
            else:
                # create identity matrix for oscillatory coordinate
                I[n_coord] = np.matrix(np.identity( _ns ))
                # create oscillator number operator
                N[n_coord] = np.matrix(np.diag(range(0,_ns)))
        
        # this string will contain progress in calculating V_hat- matrices
        _Progress_String = ""
        # get number of Josephson junctions
        N_jos = np.size(M_osc[:,0])
        # scan over Josephson junctions
        for n_jos in range(0,N_jos):
            # scan over coordinates
            for n_coord in range(0, N_tot):
                # N of states for n-th coordinate
                _ns = N_of_states[n_coord]
                # if user chooses to display the progress
                if DISPLAY_PROGRESS:
                   # this string will contain message showing current progress
                    _Progress_String = "...computing V_hat(n_jos = "\
                            +str(n_jos+1)+"/"+str(np.size(M_osc[:,0]))\
                            +", n_coord = "+str(n_coord+1)\
                            +"/"+str(N_tot)+")..."
                    print(_Progress_String)
                # calculate V for cyclic coordinates
                if n_coord < N_cyc:
                    # make empty matrix for V. Elements will be filled late
                    # r
                    currentV = np.matrix(np.zeros([2*_ns+1 ,2*_ns+1]))
                    # fill V- matrix element-by-element
                    for m1 in range(0,2*_ns+1):
                        for m2 in range(0,2*_ns+1):
                            if (m1-m2) == Z_cyc[n_jos,n_coord]:
                                # put 1 into element [m1,m2] if m1-m2=0. Ot
                                # herwise leave 0
                                currentV[m1,m2] = 1
                                
                # calculate V for oscillatory coordinates
                if n_coord >= N_cyc:
                    # make empty matrix for V. Elements will be filled late
                    # r
                    currentV = np.matrix(np.zeros([_ns,_ns]))
                    # fill V- matrix element-by-element
                    for m1 in range(0,_ns):
                        for m2 in range(0,_ns):                    
                            # Now calculate current matrix element
                            Q = M_osc[n_jos, n_coord-N_cyc]/np.sqrt(
                                       omega_eigenvalues[n_coord])
                            Sum = 0
                            for l in range(0,min(m1+1,m2+1)):
                                Sum = Sum + (-1)**l\
                                    * np.sqrt(1.*np.math.factorial(m1) \
                                                     /np.math.factorial(l))\
                                    * np.sqrt(1.*np.math.factorial(m2) \
                                                     /np.math.factorial(l))\
                                    * np.math.factorial(m1-l)**(-1)\
                                    * np.math.factorial(m2-l)**(-1)\
                                    * (Q/np.sqrt(2.)) ** (m1+m2-2*l)
                            currentV[m1,m2] = ((-1)**m2) \
                                            * np.exp(-Q**2 / 4) * Sum    

                # add the calculated V-matrix to the list
                V[n_jos][n_coord] = currentV

        # return obtained matrices
        return [I, N, V]
    
    
# __________________________________#
#!   FUNCTION TO GET HAMILTONIAN   !#                            
#!_________________________________!#
    def __obtainHamiltonian(self, N_cyc, N_osc, m_cyc, \
                                    omega_eigenvalues, E_jos, I, N, V,\
                                    F, K,   DISPLAY_PROGRESS):
        # if DISPLAY_PROGRESS is not set, then it is by default False      
        #                           # First we want to construct identity (
        # I) and number (N) operators that all have greatest dimension (per
        # form Kroneker-multiplication of all identity and number operators
        # )
        if DISPLAY_PROGRESS==None: DISPLAY_PROGRESS = False

        # this will be an array of number operators, but all brought to sam
        # e dimension. While arrays N_hat, I_hat each have dimensions of co
        # rresponding N_of_states. Replace N0,N1,N2,N3,... --> N0 x I1 x I2
        #  x ..., I0 x N1 x I2 x ...,  I0 x I1 x N2 x ..., ....
        _N = [0 for n in range(0,len(N))]
        #  NOTE: here for simplicity we omit "hats" so instead of
        # N_hat, V_hat,.. we write N, V,...

        # now want to get identity operator of total dimension ...
        _I = 1
        # will count identity operators
        current_I_index = -1
        for I_current in I:
            current_I_index=current_I_index+1
            # ... by Kronecker-multiplying all identity operators
            _I = np.kron(_I, I_current)
            # if user chooses to display the progress
            if DISPLAY_PROGRESS:
                # this string will contain message that shows current progress
                # in getting Identity matrix
                _Progress_String = "calc. direct product of identity "\
                                   "matrices: I:=IxI_" \
                                   +str(current_I_index + 1) \
                                   + "/" +str(N_cyc+N_osc) + "."
                print(_Progress_String)
    



        # counts index of the N-operator
        current_N_index = -1
        # scan over N-operators
        for N_current in N:
            current_N_index = current_N_index + 1                     
            # index that will count identity operators
            current_I_index = -1
            # current total operator (i.e. of form I0 x N1 x I2 x ...)
            current_N_tot = 1
            
            # scan over identity operators
            for I_current in I:
                current_I_index = current_I_index + 1

                if current_I_index == current_N_index:
                    # in Kronecker product in space current_N_index put N-o
                    # perator and not I-operator
                    current_N_tot = np.kron(current_N_tot, N_current)
                else:
                    # everywhere else put I-operator of corresponding dimen
                    # sions
                    current_N_tot = np.kron(current_N_tot, I_current)
           
            # add calculated operator to the list _N    
            _N[current_N_index] = current_N_tot

                                # Get cyclic part of the Hamiltonian
        Hamilt_cyc=0*_I
        # effectively (1/2)* Capacitance matrix for cyclic degrees of freed
        # om
        m2 = 0.5*m_cyc*m_cyc.T if len(m_cyc)>0 else 0
        for m in range(0,N_cyc):
            for n in range(0,N_cyc):
                # get cyclic part of the Hamiltonian
                Hamilt_cyc = Hamilt_cyc \
                           + ( _N[m] - np.kron(K[m],_I) ) \
                           * m2[m,n] * ( _N[n] - np.kron(K[n],_I) )
                                # Get oscillatory part of the Hamiltonian
        Hamilt_osc = 0*_I
        for m in range(N_cyc, N_cyc+N_osc):
            # get oscillatory Hamiltonian
            Hamilt_osc = Hamilt_osc \
                       + np.kron( omega_eigenvalues[m], (_N[m] + 0.5*_I) )
                                # Get Josephson Hamiltonian
        Hamilt_Jos_exp=0*_I;
        N_Jos = np.size(E_jos)
        for m in range(0,N_Jos):
            _row_product = 1
            # if user chooses to display the progress
            if DISPLAY_PROGRESS:
                 # this string will contain message that shows current 
                 # progress in calculating Josephson Hamiltonian
                _Progress_String = "calc. Hamilt. of Jos. junction #: "\
                                   +str(m+1)+"/"+str(N_Jos)+"..."
                print(_Progress_String)

            for n in range(0,N_cyc+N_osc):
                _row_product=np.kron(_row_product,V[m][n])
            # get sum of products of V-operators
            Hamilt_Jos_exp = Hamilt_Jos_exp \
                           + np.kron(np.exp(-1j*F[m,0])*E_jos[m],_row_product)

        # get hermitian part of the Hamiltonian
        Hamilt_Jos_exp = ( Hamilt_Jos_exp + Hamilt_Jos_exp.H )/2.
       
        return Hamilt_cyc + Hamilt_osc - Hamilt_Jos_exp























































































""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""                                                        """ 
"""                      C L A S S                         """
"""                                                        """ 
"""  T H A T     D E F I N E S    T H E    C I R C U I T   """
"""                                                        """
"""          A N D    I T S   P A R A M E T E R S          """
"""                                                        """ 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

class CircuitBuilder:
    # show properties of the class
    def __repr__(self):
        return str(dirv(self))
    
    # threshold under which numbers are considered small
    SMALL = 1E-4

    ##################################################
    #                                                #
    #                  C I R C U I T                 #
    #                                                #
    #   Class containing code that defines circuit   #
    #                                                #
    ##################################################

    class __Circuit:
        # show properties of the class
        def __repr__(self):
            return str(dirv(self))

        Code = [None]
        MutualInductanceCode = []
        
        C = np.matrix([None])
        L = np.matrix([None])
        E_jos = np.matrix([None])

        # from CircuitCode extract only lines corresponding to element Elem
        # entType(1-character variable: 'L', 'C' or 'J')
        def extract(self, ElementType):
            return [Line for Line in self.Code if Line[2][0]==ElementType]
        # get number of elements in the circuit (number of lines in Circuit
        # .Code)
        def nElements(self):
            return len(self.Code)
        # number of nodes in the circuit
        def size(self):
            return np.amax( [a[0:2] for a in self.Code] )
        # find number of the element of type "ElementType" connecting nodes
        #  Node1 and Node2
        def findElementNumber(self, Node1,Node2, ElementType):
            # extract subcircuit of inductors 
            SubCircuit = self.extract(ElementType)
            # get number of elements of our type in the circuit
            N_of_Elements_of_Type = len(SubCircuit)
            # default number of the element we are looking for. Initially n
            # ot found returns -1
            ElementNumber = -1
            # scan over all elements of our type
            for n in range(0,N_of_Elements_of_Type):
                # if element with given terminals is found
                if SubCircuit[n][0:2] == [Node1, Node2]:
                    ElementNumber = n
            # then return it
            return ElementNumber
    Circuit = __Circuit()    

    def setCapacitanceMatrix(self):
        # set capacitance matrix in Circuit.C
        self.Circuit.C = self.__obtainCapacitanceMatrix(self.Circuit)
    # set inductance matrix in Circuit.L
    def setInductanceMatrix(self):
        self.Circuit.L = self.__obtainInductanceMatrix(self.Circuit)
        self.Circuit.L = self.__obtainMutualInductances(self.Circuit)

    def setJosephsonEnergies(self):
        self.Circuit.E_jos = np.matrix([Line[3] 
                        # extract Josephson energies from Circuit.Code
                        for Line in self.Circuit.Code if Line[2][0]=='J']).T
    def getInductorConnectivityMatrixG(self):
        return self.obtainConnectivityMatrix(self.Circuit, 'L')
    def getJosephsonConnectivityMatrixH(self):
        return self.obtainConnectivityMatrix(self.Circuit, 'J')

    def getFluxBiasVector(self, FluxBiasCode):
        return self.__getFluxBiasVector(self.Circuit, FluxBiasCode)

    def getChargeBiasVector(self, ChargeBiasCode):
        return self.__getChargeBiasVector(self.Circuit, ChargeBiasCode)


#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#

          
    ######################################################
    #                                                    #
    #                  F U N C T I O N S                 #
    #                                                    #
    #         Define more complicated functions          #
    #                                                    #
    ######################################################
    

# ___________________________________________________________#    
#! get connectivity mattrix for particular type of elements !#
#!__________________________________________________________!#
    # returns connectivity matrix for elements of type "ElementType". Again
    # , this function is important so it is not hidden
    def obtainConnectivityMatrix(self, Circuit, ElementType):
        # extract subcircuit only of elements of type "ElementType"
        SubCircuit = Circuit.extract(ElementType)
        # determine number of nodes in the circuit
        NofNodes = Circuit.size()
        # determine number of elements of type "ElementType"
        NofElements = len(SubCircuit)
        # this will be connectivity matrix
        ConnectivityMatrix = np.matrix(np.zeros([NofElements,NofNodes]))
        
        # this variable counts lines in the connectivity matrix
        CurrentRow = -1
        # by default Factor that fills in connectivity matrix is 1
        Factor = 1
        for Line in SubCircuit:
            # increase the line counter
            CurrentRow = CurrentRow + 1
            
            # if Element is Josephson junction then it can have different o
            # rders (J, J2, J3, ...). Number after 'J' determine order of c
            # osine: E.g.: cos3\phi
            if Line[2][0] == 'J':
                # see if the Josephson command has form "J2", "J3", .... bu
                # t not just "J"(by default equals J1)
                if Line[2]!=Line[2][0]:
                    # extract integer following 'J'. Factor is the integer 
                    # that will be put into the connectivity matrix
                    Factor = int( Line[2][ 1:len(Line[2]) ] )
                else:
                    Factor = 1
            else:
                # for elements that are not "J2,3,4..." integer Factor is 1
                # 
                Factor = 1

            if Line[0] > 0:
                # fill the connectivity matrix elements. Zeroth(ground) ter
                # minal is omitted in connectivity matrix
                ConnectivityMatrix[CurrentRow, Line[0]-1] = Factor
            if Line[1] > 0:
                ConnectivityMatrix[CurrentRow, Line[1]-1] = -Factor

        return ConnectivityMatrix


# __________________________________________________________#    
#! convert between mutual and Maxwell capacitance matrices !#
#!_________________________________________________________!#
    # convert between mutual and Maxwell capacitance matrices
    def convertMutualMaxwellCapMat(self, C):
        # find out number of nodes (size of C)
        N_tot = len(C[:,0])
        
        # extract diagonal part of C
        C_diag = np.matrix(np.diag(np.diag(C)))
        # extract off-diagonal part of C
        C_off_diag = C - C_diag
        
        C_New_diag = np.diag(
                        # get diagonal part of new matrix
                        np.array( np.matrix(np.ones([1,N_tot]))*C )[0]
                            )
        # get transformed capacitance matrix
        C_New = np.matrix(C_New_diag) - C_off_diag
        return C_New




# ___________________________________________________________________________#    
#! get inverse square root of a positively-defined real symmetric matrix    !#
#!__________________________________________________________________________!#
    # calculate an inverse square root of real symmetric (capacitance) matr
    # ix C with tolerance equal to SMALL
    def inv_sqrt(self, C, SMALL):
        # diagonalize C
        D, V = LA.eig(C)
        # check if matrix is symmetric
        if not np.allclose(C.T, C, atol = SMALL):
            return "ERROR: Matrix is not symmetric"
        # check if matrix is real
        if not np.allclose(C.T, C.H, atol = SMALL):
            return "ERROR: Matrix is not symmetric real"
        # check if C is positive
        if min(D)<=0:
            # if not, then try to add small (parasitic) capacitances to sma
            # ll number to fix
            D = D + SMALL
        # if matrix is still not positive after adding parasitic capacitanc
        # es then generate error
        if min(D)<=0:
            return "ERROR: Matrix is non-positive. Couldn't fix"
        else:
            # turn eigenvalues into a diagonal matrix
            D = np.matrix(np.diag(D))
            V = np.matrix(V)
            # take square root of diagonalized matrix
            D_sqrt = np.sqrt(D)
            # inverse it
            D_sqrt_inv = D_sqrt.I
            # return inverse sqrt of original matrix
            return V*D_sqrt_inv*V.T



# ____________________________#    
#! get capacitance matrix    !#
#!___________________________!#
    def __obtainCapacitanceMatrix(self, Circuit):
        # number of degrees of freedom (number of nodes)
        N_tot = Circuit.size()
        # extract capacitor subcircuit
        CapSubCircuit = Circuit.extract('C')
        # this will be initial diagonal part of the capacitance matrix
        Cap_diag = np.matrix(np.zeros([N_tot,N_tot]))
        # this will be initial off-diagonal part of capacitance matrix
        Cap_off_diag = np.matrix(np.zeros([N_tot,N_tot]))
        for Line in CapSubCircuit:
            # check if current node is connected to the ground
            if Line[0]==0 or Line[1]==0:
                # number of a node connected to the ground
                Term = Line[0]+Line[1] - 1
                # add the capacitance to the diagonal of the initial capaci
                # tance matrix. Note that two capacitances in parallel can 
                # be added
                Cap_diag[Term,Term] = Cap_diag[Term,Term]+Line[3]
            else:
                # get terminals
                T1 = Line[0] - 1
                T2 = Line[1] - 1
                # add off-diagonal matrix entries
                Cap_off_diag[T1,T2] = Cap_off_diag[T1,T2] + Line[3]
                # symmetrically
                Cap_off_diag[T2,T1] = Cap_off_diag[T2,T1] + Line[3]
            
        # this will be turned into Maxwell capacitance matrix
        Cap = Cap_diag + Cap_off_diag
        # Now get reduced Maxwell capacitance matrix
        C = self.convertMutualMaxwellCapMat(Cap)
        # make sure C is symmetric
        C = 0.5*(C+C.T)
        # return reduced Maxwell capacitance matrix
        return  C


# ____________________________#    
#! get inductance matrix     !#
#!___________________________!#
    def __obtainInductanceMatrix(self, Circuit):
        Inductances = [Line[3] for Line in Circuit.Code if Line[2]=='L']
        # if no inductances found then should return 0-matrix
        if len(Inductances) == 0:
            Inductances = [0]
        return np.matrix(np.diag(Inductances))




# ______________________________#    
#! get vector of charge biases !#
#!_____________________________!#
    # get vector of charge biases. Input: ChargeArray: contains pairs [Node
    # Nmber, ChargeOffset]
    def __getChargeBiasVector(self, Circuit, ChargeArray):
        # number of nodes in the circuit
        N_Charges = Circuit.size()
        # create empty charge bias vector
        ChargeBiasVector_Q = np.matrix(np.zeros([N_Charges,1]))
        for CurrentBiasCode in ChargeArray:
            # put Charge bias in the row NodeNumber
            ChargeBiasVector_Q[CurrentBiasCode[0]-1,0] = CurrentBiasCode[1]
        return ChargeBiasVector_Q
    


    
# ____________________________#    
#! get vector of flux biases !#
#!___________________________!#
    # get vector of flux biases. Input: FluxArray: contains pairs [Node1, N
    # ode2, FluxOffset]
    def __getFluxBiasVector(self, Circuit, FluxArray):
        # find number of inductors in the circuit
        N_Fluxes = len(self.Circuit.extract('L'))
        # create empty flux bias vector
        FluxBiasVector_Phi = np.matrix(np.zeros([N_Fluxes,1]))

        for CurrentBiasCode in FluxArray:
            InductorNumber = Circuit.findElementNumber( 
                                # find number of current inductor connectin
                                # g Node1 and Node2
                                CurrentBiasCode[0],CurrentBiasCode[1], 'L' )
            if InductorNumber>=0:
                # if inductor found then put flux value in the correspondin
                # g line in FluxBiasVector_Phi
                FluxBiasVector_Phi[InductorNumber,0] = CurrentBiasCode[2]
        return FluxBiasVector_Phi
    
 
# ____________________________#    
#! set mutual inductances    !#
#!___________________________!#
    # add off-diagonal elements to the inductance matrix
    def __obtainMutualInductances(self, Circuit):
        # this will be our inductance matrix. Initially diagonal 
        L = Circuit.L
        for Line in Circuit.MutualInductanceCode:
            # find number of inductor 1
            Inductor1id = Circuit.findElementNumber(Line[0],Line[1], 'L')
            # find number of inductor 2
            Inductor2id = Circuit.findElementNumber(Line[2],Line[3], 'L')
            L[Inductor1id,Inductor2id] = Line[4]                              
            # add mutual inductance to off-diagonal position
            L[Inductor2id,Inductor1id] = Line[4]
        return L
            
       


































































""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""                                                        """ 
"""                      C L A S S                         """
"""                                                        """ 
"""   T H A T     C O N T A I N S   S H O R T - H A N D    """
"""                                                        """
"""    F U N C T I O N S    T H A T   C A L C U L A T E    """
"""                                                        """ 
"""          C I R C U I T    P R O P E R T I E S          """
"""                                                        """ 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""





class EasyModel:
    # show properties of the class
    def __repr__(self):
        return str(dirv(self))


    # this is a small parameter used in the model. Will be the value of par
    # asitic capacitances and will be used to compare collinearity of vecto
    # rs
    SMALL = 1E-4
    # this array will contain messages during calculations
    MESSAGE_ARRAY = []
    # determines if messages will be displayed on the screen
    DISPLAY = True
    # this property determines if we want to force Z_cyc-matrix to be force
    # d to be integer
    ENFORCE_INT_IN_Z_CYC = True
    # this set will contain intermediate results of the calculations
    DEBUG_DATA = {
                      "CircBuild.Circuit.C": "Not obtained",
                      "CircBuild.Circuit.L": "Not obtained",
                      "CircBuild.Circuit.E_jos": "Not obtained",
                      "ModelCalc.Circuit.Inv_sqrt_C": "Not obtained",
                      "ModelCalc.Circuit.G": "Not obtained",
                      "ModelCalc.Circuit.Inv_L": "Not obtained",
                      "ModelCalc.Circuit.H": "Not obtained",
                      "ModelCalc.Circuit.E_jos": "Not obtained",
                      "ModelCalc.Circuit.Phi": "Not obtained",
                      "ModelCalc.Circuit.Q": "Not obtained",
                      "ModelCalc.Matrices.U": "Not obtained",
                      "ModelCalc.Matrices.omega_eigenvalues": "Not obtained",
                      "ModelCalc.Matrices.N_cyc": "Not obtained",
                      "ModelCalc.Matrices.N_osc": "Not obtained",
                      "ModelCalc.Matrices.M": "Not obtained",
                      "ModelCalc.Matrices.M_cyc": "Not obtained",
                      "ModelCalc.Matrices.m_cyc": "Not obtained",
                      "ModelCalc.Matrices.Z_cyc": "Not obtained",
                      "ModelCalc.Matrices.F": "Not obtained",
                      "ModelCalc.Matrices.K": "Not obtained"
                 }
    # Class containing circuit parameters
    class __SQM:
        # show properties of the class
        def __repr__(self):
            return str(dirv(self))
        CircBuild = CircuitBuilder()
        ModelCalc = ModelCalculator()
    SQM = __SQM()


    def __init__(self, SMALL, CircuitCode, MutualInductanceCode):
        self.SMALL = SMALL
        self.SQM = self.__SQM()
        self.Code.CircuitCode = CircuitCode
        self.Code.MutualInductanceCode \
              =MutualInductanceCode if MutualInductanceCode!=None else []


    
    # Class containing code
    class __Code:
        # show properties of the class
        def __repr__(self):
            return str(dirv(self))
        CircuitCode = []
        MutualInductanceCode = []
        FluxBiasCode = []
        ChargeBiasCode = []
    Code = __Code()
    
    # function to display messages during calculations
    def printMessage(self, Message):
        if self.DISPLAY:
            print(Message)
            self.MESSAGE_ARRAY.append(Message)
    # function that builds circuit parameters from code
    def buildCircuit(self):
        self.SQM.CircBuild = self.__buildCircuit(
                                self.Code.CircuitCode, \
                                self.Code.MutualInductanceCode, self.SMALL)
    # function that performes coordinate transformations
    def transformCoordinates(self):
        self.SQM.ModelCalc = self.__transformCoordinates(
                                               self.SQM.CircBuild, self.SMALL)

    def transformBiasVectors(self, FluxBiasCode, ChargeBiasCode):
        self.SQM.ModelCalc.Circuit.Phi \
            = self.SQM.CircBuild.getFluxBiasVector(FluxBiasCode)
        self.Code.FluxBiasCode = FluxBiasCode
        self.SQM.ModelCalc.Circuit.Q \
            = self.SQM.CircBuild.getChargeBiasVector(ChargeBiasCode)
        self.Code.ChargeBiasCode = ChargeBiasCode
        self.SQM.ModelCalc.calc_Transformed_Flux_Charge_Biases()

    # make summary of cyclic and oscillatory coordinates and their correspo
    # nding frequencies
    def summarizeCoord(self):
        self.printMessage( "" )
        self.printMessage( "----SUMMARY----" )
        self.printMessage(  "Frequencies of cyclic coordinates (" \
                            + str(self.SQM.ModelCalc.Matrices.N_cyc) + "): " 
                            + str(self.SQM.ModelCalc.Matrices\
                                  .omega_eigenvalues.T\
                                  [0,0:self.SQM.ModelCalc.Matrices.N_cyc])  )
        self.printMessage("Frequencies of oscillatory coordinates (" 
                        + str(self.SQM.ModelCalc.Matrices.N_osc) + "): " 
                        + str(self.SQM.ModelCalc.Matrices.omega_eigenvalues.T[\
                            0,self.SQM.ModelCalc.Matrices.N_cyc:(
                                    self.SQM.ModelCalc.Matrices.N_cyc
                                   +self.SQM.ModelCalc.Matrices.N_osc)])  )
        self.printMessage("Transformation II matrix Mat: "
                          + "( \\vec{v,z} = Mat\\vec{\phi} ) converts "
                          + "original coordinates to new ones.  Mat:"
                          + "\n"
                          + str( self.SQM.ModelCalc.getTransform2Matrix() )
                         )
        self.printMessage( "---------------" )

    # set n of states from the list
    def set_N_of_states(self, N_of_states):
        self.SQM.ModelCalc.Operators.N_of_states = N_of_states

    # calculate quantum operators
    def setOperators(self):
        self.SQM.ModelCalc.obtainOperators(self.DISPLAY)
    
    def setHamiltonian(self):
        self.SQM.ModelCalc.obtainHamiltonian(self.DISPLAY)
    def getHamiltonian(self):
        return self.SQM.ModelCalc.Operators.Hamiltonian

    
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#

          
    ######################################################
    #                                                    #
    #                  F U N C T I O N S                 #
    #                                                    #
    #         Define more complicated functions          #
    #                                                    #
    ######################################################
    

# _______________________________________________#    
#! call set of functions for building a circuit !#
#!______________________________________________!#
    # build the circuit
    def __buildCircuit(self, CircuitCode, MutualInductanceCode, SMALL):
        self.printMessage("    CALCULATING CIRCUIT PARAMETERS ...")
        # create object that will help us build a circuit
        CircBuild = CircuitBuilder()
        # set small parameter in "CircuitBuilder" that will be used as an i
        # nfinitesmall parasitic capacitance value
        CircBuild.SMALL = SMALL
        # input a code that describes the circuit structure into "CircuitBu
        # ilder.Circuit.Code"
        CircBuild.Circuit.Code = CircuitCode
        # set code to calculate mutual inductances
        CircBuild.Circuit.MutualInductanceCode = MutualInductanceCode
        # from the code describing the circuit, obtain Maxwell capacitance 
        # matrix and put it into "CircuitBuilder.Circuit.C"
        CircBuild.setCapacitanceMatrix()
        # update DEBUG_DATA
        self.DEBUG_DATA["CircBuild.Circuit.C"] = CircBuild.Circuit.C
        self.printMessage("[OK] CircuitBuilder.setCapacitanceMatrix(): "\
                          "Capacitance matrix calculated" if \
                                   np.linalg.det(CircBuild.Circuit.C)>0 else \
                          "[WARNING] CircuitBuilder.setCapacitanceMatrix(): "\
                          "Capacitance matrix not strictly positive. "\
                          "But may be able to fix it. "\
                          # Check if Maxwell capacitance matrix is positive
                          # ly-defined
                          "Still, try adding capacitors to circuit.")
        # Obtain inductance matrix based on "CircuitBuilder.Circuit.Code" a
        # nd "CircuitBuilder.Circuit.MutualInductanceCode". Inductance matr
        # ix is saved in "CircuitBuilder.Circuit.L"
        CircBuild.setInductanceMatrix()
        self.DEBUG_DATA["CircBuild.Circuit.L"] = CircBuild.Circuit.L
        self.printMessage("[OK] CircuitBuilder.setInductanceMatrix(): "\
                          "Inductance matrix set" if \
                                          len(CircBuild.Circuit.L)>0 else \
                                          "[ERROR] CircuitBuilder."\
                                          "setInductanceMatrix(): No "\
                                          # make sure that circuit contains
                                          #  at least one inductor         
                                          #   
                                          "inductors found in circuit")
        # from "CircuitBuilder.Circuit.Code" obtain a vector of Josephson e
        # nergies. The vector is saved in "CircuitBuilder.Circuit.E_jos"
        CircBuild.setJosephsonEnergies()
        self.DEBUG_DATA["CircBuild.Circuit.E_jos"] = CircBuild.Circuit.E_jos
        self.printMessage("[OK] CircuitBuilder.setJosephsonEnergies(): "\
                          "Josephson energies set" if \
                                          len(CircBuild.Circuit.E_jos)>0 else \
                                              "[ERROR] CircuitBuilder."\
                                              "setJosephsonEnergies(): No "\
                                              "Josephson junctions found in "\
                                              # make sure the circuit conta
                                              # ins at least one Josephson 
                                              # junction
                                              "circuit")
        return CircBuild
    

# ______________________________________________________________#    
#! call set of functions to perform cooedinate transformations !#
#!_____________________________________________________________!#
    # perform coordinate transformations
    def __transformCoordinates(self, CircBuild, SMALL):

        # create object that stores circuit parameters most important for o
        # ur calculations
        CircParam = CircuitParameters()
        CircParam.Inv_sqrt_C = CircBuild.inv_sqrt( CircBuild.Circuit.C, 
                                                   # calculate inverse squa
                                                   # re root of the Maxwell
                                                   #  capacitance matrix an
                                                   # d save it in "CircuitP
                                                   # arameters.Inv_sqrt_C"
                                                   CircBuild.SMALL )
        self.DEBUG_DATA["ModelCalc.Circuit.Inv_sqrt_C"]=CircParam.Inv_sqrt_C
        self.printMessage("[OK] CircuitBuilder.inv_sqrt(C): C^(-1/2) "\
                          "calculated" if type(CircParam.Inv_sqrt_C)!=str else\
                          ("[ERROR] CircuitBuilder.inv_sqrt(C): Could not"\
                           "raise capacitance matrix to power -1/2. ErrorMsg"\
                           # Sometimes, Maxwell capacitance matrix can be n
                           # on-positive. The function "CircuitBuilder.inv_
                           # sqrt()" tries to correct that. If correction i
                           # s unsuccessful, return error
                           "<"+CircParam.Inv_sqrt_C)+">")
        # from "CircuitBuilder.Circuit.Code" extract inductance connectivit
        # y matrix (transforms phases on nodes into phase drops across indu
        # ctances). Save it in "CircuitParameters.G"
        CircParam.G = CircBuild.getInductorConnectivityMatrixG()
        self.DEBUG_DATA["ModelCalc.Circuit.G"] = CircParam.G
                           # make sure that inductance matrix is invertible
        self.printMessage( "[ERROR] CircuitBuilder.Circuit.L: Zero "\
                           "inductance(s) in the circuit" if \
                                   np.linalg.det(CircBuild.Circuit.L)==0 else\
                                                                          "" )
        # invert inductance matrix in "CircuitBuilder.Circuit.L" and save i
        # n "CircuitParameters.Inv_L"
        CircParam.Inv_L = CircBuild.Circuit.L.I
        self.DEBUG_DATA["ModelCalc.Circuit.Inv_L"]=CircParam.Inv_L
        # from "CircuitBuilder.Circuit.Code" extract Josephson connectivity
        #  matrix (transforms phases on nodes into phase drops across Josep
        # hson junctions). Save it in "CircuitParameters.H"
        CircParam.H = CircBuild.getJosephsonConnectivityMatrixH()
        self.DEBUG_DATA["ModelCalc.Circuit.H"]=CircParam.H
        # transfer Josephson energies into "CircuitParameters.E_jos"
        CircParam.E_jos = CircBuild.Circuit.E_jos
        self.DEBUG_DATA["ModelCalc.Circuit.E_jos"]=CircParam.E_jos
        # from "CircuitBuilder.Biases.FluxBiasCode" extract vector of flux 
        # biases in inductors. Save into CircuitParameters.CircParam.Phi
        CircParam.Phi = CircBuild.getFluxBiasVector([])
        self.DEBUG_DATA["ModelCalc.Circuit.Phi"]=CircParam.Phi
        # from "CircuitBuilder.Biases.ChargeBiasCode" extract vector of cha
        # rge biases on nodes. Save into CircuitParameters.CircParam.Q
        CircParam.Q = CircBuild.getChargeBiasVector([])
        self.DEBUG_DATA["ModelCalc.Circuit.Q"]=CircParam.Q


        self.printMessage("    COORDINATE TRANSFORM (I) ...")
        # create ModelCalculator -- object that will help us perform proper
        #  coordinate transformations and generate quantum operators
        ModelCalc = ModelCalculator()
        # set small parameter of "ModelCalculator". It will be used to chec
        # k if vectors are collinear
        ModelCalc.SMALL = SMALL
        # import circuit parameters into ModelCalculator
        ModelCalc.Circuit = CircParam
        # perform "coordinate transformation I" that untangles harmonic osc
        # illators inside the circuit
        ModelCalc.calc_U_omega_eigenvalues()
        self.DEBUG_DATA["ModelCalc.Matrices.U"]=ModelCalc.Matrices.U
        self.DEBUG_DATA["ModelCalc.Matrices.omega_eigenvalues"]\
                       = ModelCalc.Matrices.omega_eigenvalues
        self.printMessage("[OK] ModelCalculator.calc_U_omega_eigenvalues(): "\
                         "Harmonic oscillators diagonalized. Frequency modes: "
                         # display frequency modes of harmonic oscillators.
                         #  Save frequency modes in "ModelCalculator.Matric
                         # es.omega_eigenvalues". Save Unitary operator des
                         # cribing the transformation into "ModelCalculator
                         # .Matrices.U"
                         + str(ModelCalc.Matrices.omega_eigenvalues.T)  )
        # determine number of cyclic coordinates (with small frequency mode
        # s below ModelCalculator.SMALL). Then get number of oscillatory co
        # ordinates. These values are stored in "ModelCalculator.Matrices.N
        # _cyc" and "ModelCalculator.Matrices.N_osc" respectively
        ModelCalc.calc_N_cyc_osc_tot()
        self.DEBUG_DATA["ModelCalc.Matrices.N_cyc"]=ModelCalc.Matrices.N_cyc
        self.DEBUG_DATA["ModelCalc.Matrices.N_osc"]=ModelCalc.Matrices.N_osc
        self.printMessage( "[OK] ModelCalculator.calc_N_cyc_osc_tot(): "\
                          "Cyclic and oscillatory coordinates separated. "\
                          "N_cyc: " + str(ModelCalc.Matrices.N_cyc) 
                                    + "; N_osc: "
                                    # display number of cyclic and oscillat
                                    # ory coordinates
                                    + str(ModelCalc.Matrices.N_osc)  )
        # calculate matrix that transforms our new coordinates (after "tran
        # sformation I") into phase drops across Josephson junctions. This 
        # matrix is stored in "ModelCalculator.Matrices.M"
        ModelCalc.calc_M()
        self.DEBUG_DATA["ModelCalc.Matrices.M"]=ModelCalc.Matrices.M
        self.printMessage("    COORDINATE TRANSFORM (II) ...")
        # split M-matrix into cyclic and oscillatory columns. Save it into 
        # "ModelCalculator.Matrices.M_cyc" and "ModelCalculator.Matrices.M_
        # osc"
        ModelCalc.calc_M_cyc_osc()
        self.DEBUG_DATA["ModelCalc.Matrices.M_cyc"]=ModelCalc.Matrices.M_cyc
        self.DEBUG_DATA["ModelCalc.Matrices.M_osc"]=ModelCalc.Matrices.M_osc
        # calculate largest invertible submatrix of M_cyc. This is needed t
        # o restore integer charge quantization for cyclic coordinates. Sav
        # ed in "ModelCalculator.Matrices.m_cyc"
        ModelCalc.calc_m_cyc()
        self.DEBUG_DATA["ModelCalc.Matrices.m_cyc"]=ModelCalc.Matrices.m_cyc
        # if m_cyc is returned as a string, that is an error message
        if type(ModelCalc.Matrices.m_cyc)!=str:
            # if m_cyc is non-empty
            if len(ModelCalc.Matrices.m_cyc)>0:
                # check if matrix m_cyc is square
                if (
                    np.size(ModelCalc.Matrices.m_cyc[0,:]) \
                                   ==np.size(ModelCalc.Matrices.m_cyc[:,0])
                   ) \
                and(
                    # and checking m_cyc has dimensions of N_cyc
                    np.size(ModelCalc.Matrices.m_cyc[0,:]) \
                                    ==ModelCalc.Matrices.N_cyc
                   ):     
                    self.printMessage(  "[OK] ModelCalculator.calc_m_cyc(): "\
                                        "matrix m_cyc obtained: det(m_cyc) = "
                                        + str(np.linalg.det(
                                                    ModelCalc.Matrices.m_cyc))
                                        + " -- larger absolute value of "
                                        # if it is not square, then probabl
                                        # y our circuit contains clusters t
                                        # hat are not connected to ground v
                                        # ia any Josephson junctions
                                        + "determinant is better"  )
                else:
                    self.printMessage("[ERROR] ModelCalculator.calc_m_cyc()"\
                                      ": matrix m_cyc has wrong dimensions. "\
                                      "Circuit "\
                                      "likely contains clusters not connected"\
                                      " to ground via Josephson junctions "\
                                      "(No charge measure is provided to "\
                                      "enforce charge quantization). Try "\
                                      "adding Jos. junctions to circuit."  ) 
        else:
            self.printMessage(  "[ERROR] ModelCalculator.calc_m_cyc(): "\
                                "rows in M_cyc are too small. ErrorMsg <"
                                + ModelCalc.Matrices.m_cyc + ">")
        # calculate integer-valued matrix Z_cyc = M_cyc/m_cyc. Matrix "Mode
        # lCalculator.Matrices.Z_cyc" is forced to be integer-valued
        ModelCalc.calc_Z_cyc(self.ENFORCE_INT_IN_Z_CYC)
        self.DEBUG_DATA["ModelCalc.Matrices.Z_cyc"]=ModelCalc.Matrices.Z_cyc
        # here we want to make sure that M_cyc/m_cyc is indeed close to int
        # eger
        if len(ModelCalc.Matrices.m_cyc)>0:
            Delta = ModelCalc.Matrices.Z_cyc - ModelCalc.Matrices.M_cyc\
                                              *ModelCalc.Matrices.m_cyc.I
            if np.max(np.abs(Delta)) > 2*ModelCalc.SMALL:
                self.printMessage(  "[WARNING] ModelCalculator.calc_Z_cyc():"\
                                    " matrix Z_cyc is far from integer. "\
                                    # if not integer, then make a warning
                                    "Deviations:")
                self.printMessage(  str(Delta)  )
            else:
                self.printMessage(  "[OK] ModelCalculator.calc_Z_cyc(): "\
                                    "matrix Z_cyc calculated:")
        # perform our coordinate transformations on both: charge and flux b
        # ias vectors "ModelCalculator.Circuit.Phi" and "ModelCalculator.Ci
        # rcuit.Q". Transformed charge and flux bias vectors are saved in "
        # ModelCalculator.Matrices.F" and "ModelCalculator.Matrices.K" resp
        # ectively 
        ModelCalc.calc_Transformed_Flux_Charge_Biases()
        self.DEBUG_DATA["ModelCalc.Matrices.F"]=ModelCalc.Matrices.F
        self.DEBUG_DATA["ModelCalc.Matrices.K"]=ModelCalc.Matrices.K

        return ModelCalc

