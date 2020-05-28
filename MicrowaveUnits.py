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


################################################
# This module converts standard units into     #
# dimensionless ones with                      #
#                                              #
#          h_bar = 2e = Kelvin = 1             #
#                                              #
################################################


# pi
pi = 3.14159265358979323846264338327950288419716939937510

# let us measure eveything in units of 1K
K = 1.
mK = 1E-3
uK = 1E-6

# other energy units: frequency units
GHz = 1/20.84
THz = (1E+3)/20.84
MHz = (1E-3)/20.84
kHz = (1E-6)/20.84
# electronVolts
meV = 241.8/20.84
ueV = (1E-3) * 241.8/20.84
eV =  (1E+3) * 241.8/20.84

# time units
ms = 1/kHz
us = 1/MHz
ns = 1/GHz

# capacitance units
fF = 1/(2*3.718)
pF = (1E+3)/(2*3.718)
aF = (1E-3)/(2*3.718)
nF = (1E+6)/(2*3.718)

# inductance units
nH = 1/(7.845)
pH = (1E-3)/(7.845)
uH = (1E+3)/(7.845)

# electron charge (so that 2e = 1)
e = 0.5
Coulomb = 0.5/(1.602176634E-19)
C = Coulomb

# electric current
A = Coulomb*kHz*1E-3  # Coulomb*second
mA = (1E-3)*A
uA = (1E-6)*A
nA = (1E-9)*A

# voltage
V = eV/e
mV = (1E-3)*V
uV = (1E-6)*V
nV = (1E-9)*V

# power
W  = mV*Coulomb*kHz
mW = mV*Coulomb*kHz*1E-3
uW = uV*Coulomb*kHz*1E-3
nW = (1E-3)*uV*Coulomb*kHz*1E-3
pW = (1E-6)*uV*Coulomb*kHz*1E-3

# flux quantum (2\pi)
Phi_0 = 2*pi
mPhi_0 = (1E-3)*2*pi
uPhi_0 = (1E-6)*2*pi

# Planck constant
h_bar = 1
h = 2*pi
   
# Conductance quantum (1/4*pi)
Conductance_quantum = 1/(4*pi)
G_0 = 1/(4*pi)

# Conductance
S = G_0/(7.748091729E-5)
mS = (1E-3)*G_0

# Resistance:
Ohm = 1/S
mOhm = (1E-3)*Ohm
kOhm = (1E+3)*Ohm
MOhm = (1E+6)*Ohm

#----------------------------
# M E T H O D S
#----------------------------


# convert energy to capacitance
def convert_E_to_C(E, CAP_TYPE = None):
    # use one-electron formula
    if (CAP_TYPE==None)or(CAP_TYPE=="e")or(CAP_TYPE=="1e")or(CAP_TYPE==1):
        return 0.25/(2*E)
    # use 2-electron formula
    if (CAP_TYPE=="2e")or(CAP_TYPE==2):
        return 1.0/(2*E)

# way to convert energy to capacitance via command 123*GHz|E_to_C()
#                                               or 123*GHz|E_to_C("2e")
class E_to_C:
    # Capacitance formula type (1e or 2e)
    __CAP_TYPE = None
    def __init__(self, CAP_TYPE = None):
        self.__CAP_TYPE = CAP_TYPE
    def __ror__(self,other):
        return convert_E_to_C(other, self.__CAP_TYPE)
        



# convert energy to inductance    
def convert_E_to_L(E):
    return 1.0/E
  
# way to convert energy to inductance via command 123*GHz|E_to_L
class E_to_L:
    def __ror__(self,other):
        return convert_E_to_L(other)            



# convert capacitance.inductance to energy -- symmetric
convert_C_to_E = convert_E_to_C
convert_L_to_E = convert_E_to_L
C_to_E = E_to_C
L_to_E = E_to_L
