#
#
# Test Analytic Density diagnostic equation(s)
#
# Papers/sources are realllyyyy good at showing how to find upper/lower
# limits to Ne based on collisional and radiative processes in ions,
# but they annoyingly don't explain how to find the in-between densities!!
#
# Here I will try to derive
#
#

import numpy
#import matplotlib.pyplot as plt
import pylab as plt

# Define constants here
C = 8.629E-6			# in collision rate equation
k = 8.61733E-5			# Boltzmann constant
# [OII]: 20 = 3729, 10 = 3726
A20 = 3.6E-5		# spontaneous decay
A10 = 1.8E-4
Om20 = 1.34 		# Collision rate per volume
Om10 = 1.34
g2 = 6.				# statistical weight
g1 = 4.
E20 = 3.324085		# Energy of upper level
E10 = 3.326567


ne = numpy.logspace(0.,9.,num=100)
te = 1.E+4

q20 = (C/te**0.5)*(Om20)#/g2)
q10 = (C/te**0.5)*(Om10)#/g1)



try9 = (Om20/Om10)*(g2/g1)*((1.+(ne*q10/(A10)))/(1.+(ne*q20/(A20))))*numpy.exp(-(E20-E10)/(k*te))*(E10/E20)#*(E20/E10)

fig = plt.figure()
plt.semilogx(ne,try9,label='[OII]')		# YESSSS!!



# [NI]: 20 = 5200, 10 = 5197
A20 = 7.3E-6
A10 = 2.0E-5
Om20 = 0.48
Om10 = 0.48
g2 = 6.
g1 = 4.
E20 = 2.3835297
E10 = 2.384610

q20 = (C/te**0.5)*(Om20)#/g2)
q10 = (C/te**0.5)*(Om10)#/g1)


try9 = (Om20/Om10)*(g2/g1)*((1.+(ne*q10/(A10)))/(1.+(ne*q20/(A20))))*numpy.exp(-(E20-E10)/(k*te))*(E10/E20)	# CORRECT!!!
plt.semilogx(ne,try9,label='[NI]')		# YESSSS!!



# [SII]: 20 = 6716, 10 = 6731
A20 = 2.6E-4
A10 = 8.8E-4
Om20 = 6.98
Om10 = 6.98
g2 = 6.
g1 = 4.
E20 = 1.845471
E10 = 1.841530

q20 = (C/te**0.5)*(Om20)#/g2)
q10 = (C/te**0.5)*(Om10)#/g1)


try9 = (Om20/Om10)*(g2/g1)*((1.+(ne*q10/(A10)))/(1.+(ne*q20/(A20))))*(E10/E20)*numpy.exp(-(E20-E10)/(k*te))		# CORRECT!!!
plt.semilogx(ne,try9,label='[SII] (6700)',linestyle='dashed')		# YESSSS!!



# [SII]: 20 = 4069, 10 = 4076
A20 = 9.1E-2
A10 = 2.2E-1
Om20 = 2.28
Om10 = 2.28
g2 = 4.
g1 = 2.
E20 = 3.046483
E10 = 3.040691

q20 = (C/te**0.5)*(Om20)#/g2)
q10 = (C/te**0.5)*(Om10)#/g1)

# Doesn't quite work.. (not sensitive to density changes until ~1E6)
try9 = (Om20/Om10)*(g2/g1)*((1.+(ne*q10/(A10)))/(1.+(ne*q20/(A20))))*(E10/E20)*numpy.exp(-(E20-E10)/(k*te))		# CORRECT!!!
plt.semilogx(ne,try9,label='[SII] (4000)')		# YESSSS!!



#~ # [MgI]: 20 = 4563, 10 = 4571
#~ A20 = 2.8E-4
#~ A10 = 4.3E+2
#~ Om20 = 2.4
#~ Om10 = 2.4
#~ g2 = 5.
#~ g1 = 3.
#~ E20 = 2.7166398
#~ E10 = 2.7115919

#~ q20 = (C/te**0.5)*(Om20)#/g2)
#~ q10 = (C/te**0.5)*(Om10)#/g1)

#~ # Doesn't quite work.. (not sensitive to density changes until ~1E6)
#~ try9 = (Om20/Om10)*(g2/g1)*((1.+(ne*q10/(A10)))/(1.+(ne*q20/(A20))))*(E10/E20)*numpy.exp(-(E20-E10)/(k*te))		# CORRECT!!!
#~ plt.semilogx(ne,try9,label='[MgI] (4570)')		# YESSSS!!



# [Cl III]: 20 = 5518, 10 = 5538
A20 = 6.2E-4	#8.6E-4
A10 = 4.4E-3	#6.5E-3
Om20 = 3.4	#2.05
Om10 = 3.4	#1.36
g2 = 6.
g1 = 4.
E20 = 2.246399
E10 = 2.238220

q20 = (C/te**0.5)*(Om20)#/g2)
q10 = (C/te**0.5)*(Om10)#/g1)

# Doesn't quite work.. (not sensitive to density changes until ~1E6)
try9 = (Om20/Om10)*(g2/g1)*((1.+(ne*q10/(A10)))/(1.+(ne*q20/(A20))))*(E10/E20)*numpy.exp(-(E20-E10)/(k*te))		# CORRECT!!!
plt.semilogx(ne,try9,label='[Cl III] (5500)')		# YESSSS!!



# [Ar IV]: 20 = 5518, 10 = 5538
A20 = 2.8E-3
A10 = 2.2E-2
Om20 = 3.24
Om10 = 3.24
g2 = 6.
g1 = 4.
E20 = 2.63092
E10 = 2.61490

q20 = (C/te**0.5)*(Om20)#/g2)
q10 = (C/te**0.5)*(Om10)#/g1)

# Doesn't quite work.. (not sensitive to density changes until ~1E6)
try9 = (Om20/Om10)*(g2/g1)*((1.+(ne*q10/(A10)))/(1.+(ne*q20/(A20))))*(E10/E20)*numpy.exp(-(E20-E10)/(k*te))		# CORRECT!!!
plt.semilogx(ne,try9,label='[ArIV] (4700)')		# YESSSS!!


# [Ne III]: 20 = 3967, 10 = 3869
A20 = 5.4E-2
A10 = 1.7E-1
Om20 = 1.65
Om10 = 1.65
g2 = 5.
g1 = 3.
E20 = 3.1241345
E10 = 3.203841

q20 = (C/te**0.5)*(Om20)#/g2)
q10 = (C/te**0.5)*(Om10)#/g1)

# Doesn't quite work.. (not sensitive to density changes until ~1E6)
try9 = (Om20/Om10)*(g2/g1)*((1.+(ne*q10/(A10)))/(1.+(ne*q20/(A20))))*(E10/E20)*numpy.exp(-(E20-E10)/(k*te))		# CORRECT!!!
plt.semilogx(ne,try9,label='[Ne III]')		# YESSSS!!



# [Ne IV]: 20 = 4714, 10 = 4724
A20 = 4.0E-1
A10 = 4.4E-1
Om20 = 0.90
Om10 = 0.51
g2 = 4.
g1 = 4.
E20 = 5.11242
E10 = 5.11793

q20 = (C/te**0.5)*(Om20)#/g2)
q10 = (C/te**0.5)*(Om10)#/g1)

# Doesn't quite work.. (not sensitive to density changes until ~1E6)
try9 = (Om20/Om10)*(g2/g1)*((1.+(ne*q10/(A10)))/(1.+(ne*q20/(A20))))*(E10/E20)*numpy.exp(-(E20-E10)/(k*te))		# CORRECT!!!
plt.semilogx(ne,try9,label='[Ne IV]')		# YESSSS!!





plt.xlim(1,1E6)
plt.ylabel('line ratio')
plt.xlabel(r'N$_{e}$ (for T$_{e}$ = 10$^{4}$) [cm$^{-3}$]')
plt.legend()
plt.show()

