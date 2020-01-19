# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 23:22:36 2020

@author: wesle
"""

import matplotlib.pyplot as plt
AR = [6,12,18]
CLCDS20 = [17.468, 27.583, 35.772] 
CLCDS60 = [17.498, 27.875, 36.636]
CLCDS100 =[17.509, 28.032, 36.91] 

plt.figure()
plt.subplot(121)
plt1 = plt.plot(AR,CLCDS20, label='S = 20 $m^2$', marker='o')
plt2 = plt.plot(AR,CLCDS60, label='S = 60 $m^2$', marker='o')
plt3 = plt.plot(AR,CLCDS100, label='S = 100 $m^2$', marker='o')
plt.legend()
plt.ylabel('$C_L/C_D$')
plt.xlabel('Aspect ratio')
plt.title('Sensitivity $C_L/C_D$ to changes in surface area')

V = [80, 100, 120]
CLCDV = [27.953, 27.875, 27.803]



plt.subplot(122)
plt2 = plt.plot(V,CLCDV, label = 'AR = 12', marker='o')
plt.ylabel('$C_L/C_D$')
plt.xlabel('Apparent windspeed $V_a$ [m/s]')
plt.ylim([27.5,28.5])
plt.title('Sensitivity $C_L/C_D$ to changes in apparent windspeed')
plt.legend()
plt.show()