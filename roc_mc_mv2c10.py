#!/usr/bin/env python
# coding: utf-8
from __future__ import division
from ROOT import TFile
import matplotlib.pyplot as plt
f=TFile('a1_mv2c10_reg.root')
h_0=f.Get('mv2c10_0')
h_2=f.Get('mv2c10_2')
x_vals1 = [h_2.Integral(i,-1)/h_2.Integral(0,-1) for i in range(101)]
y_vals1 = [h_0.Integral(0,-1)/h_0.Integral(i,-1) for i in range(101)]
plt.plot(x_vals,y_vals,"-b", label="tt-bar d0max=5.0mm,z0max=5.0mm,light jet rejection mv2c10")

'''
f1=TFile('mv2c_mod.root')
h_01=f1.Get('w_mv2c10_0')
h_21=f1.Get('w_mv2c10_2')
x_vals1 = [h_21.Integral(i,-1)/h_21.Integral(0,-1) for i in range(101)]
y_vals1 = [h_01.Integral(0,-1)/h_01.Integral(i,-1) for i in range(101)]
plt.plot(x_vals1,y_vals1,"-r", label="tt-bar d0max=7.5mm,z0max=7.5mm,light jet rejection mv2c10")

f3=TFile('mv2c_mod1.root')
h_03=f3.Get('mv2c10_0')
h_31=f3.Get('mv2c10_2')
x_vals3 = [h_31.Integral(i,-1)/h_31.Integral(0,-1) for i in range(101)]
y_vals3 = [h_03.Integral(0,-1)/h_03.Integral(i,-1) for i in range(101)]
plt.plot(x_vals3,y_vals3,"-g", label="d0max=7.5mm,z0max=7.5mm,light jet rejection mv2c10")
'''
plt.xlabel('b-tagging efficiency')
plt.ylabel('light jet rejection')
plt.legend(loc="upper left")
plt.savefig('bte_ljr_reg_mv2c10.png')
plt.show()


plt.figure(1)
h_1=f.Get('mv2c10_1')
x_vals2 = [h_1.Integral(i,-1)/h_1.Integral(0,-1) for i in range(101)]
y_vals2 = [h_0.Integral(0,-1)/h_0.Integral(i,-1) for i in range(101)]
plt.plot(x_vals2,y_vals2,"-b", label="tt-bar d0max=5.0mm,z0max=5.0mm, c jet rejection mv2c10")

'''
h_11=f1.Get('w_mv2c10_1')
x_vals11 = [h_11.Integral(i,-1)/h_11.Integral(0,-1) for i in range(101)]
y_vals11 = [h_01.Integral(0,-1)/h_01.Integral(i,-1) for i in range(101)]
plt.plot(x_vals11,y_vals11,"-r", label="tt-bar d0max=7.5mm,z0max=7.5mm, c jet rejection mv2c10")

h_33=f3.Get('mv2c10_1')
x_vals33 = [h_33.Integral(i,-1)/h_33.Integral(0,-1) for i in range(101)]
y_vals33 = [h_03.Integral(0,-1)/h_03.Integral(i,-1) for i in range(101)]
plt.plot(x_vals33,y_vals33,"-g", label="d0max=7.5mm,z0max=7.5mm, c jet rejection mv2c10")

plt.xlabel('b-tagging efficiency')
plt.ylabel('c jet rejection')
plt.legend(loc="upper left")
plt.savefig('bte_cjr_reg_mv2c10.png')
plt.show()
'''


# In[ ]:




