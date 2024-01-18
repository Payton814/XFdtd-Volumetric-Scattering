#####################################################################################
#
#   The purpose of this code is to read in the S-parameter data
#   Collected in XFdtd for volumetric scattering and calculate
#   The optical depth. This is then compared to the optical
#   depth that is calculated using the scattering cross section
#   from the Ulaby mie scattering matlab code
#
#####################################################################################
#
#   How to use:     You will need to make sure that the df_matlab is reading the correct
#                   CSV file. This may mean you have to change your path to wherever
#                   it was saved.
#
#                   The github as of this writing has data for rocks embedded in regolith
#                   as well as cavities in ice. MAKE SURE you have the file names correct
#                   so that the data being accessed is the correct one.
#                   for instance, for rock/regolith data folder will have: regbk_rock
#
#                   The data and folders are formatted such that you can know some of the
#                   basic XFdtd simulation set up. 
#                   EXAMPLE:    0.5r --> the simulation was populated with scatterers with
#                                        radius = 0.5 m
#
#                               28a14b --> the simulation was placed in a waveguide with 
#                                        width, a = 28 m, and height, b = 14 m
#
#                               8d --> The distance inside the waveguide actually populated
#                                      by scatterers was d = 8 m
#
#                               250f --> The driving frequency for the TE10 mode was f = 250MHz
#
#####################################################################################

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from S_param_df_gen import SP_df_gen




r = 0.5
porosity = 0.005
datalength = 40
a = 22
b = 11
f = 250

df_matlab = pd.read_csv("./VS_regbk_rock_" + str(r) + "r/CS_regbk_rockbk_" + str(f) + "MHz_0.001_step.csv")
iter = 0
while True:
    if df_matlab['Radius'][iter] == r:
        CS_M = df_matlab['Qs'][iter]
        break
    iter += 1
Qs_M = CS_M/(np.pi*r**2)

ddf = []
S2sum = []
ts_err = []
Sf10 = []
numModes = 1
d_xf = np.arange(5, datalength + 5, 1)
d = np.arange(0,datalength + 5,1)
V_s = (4/3)*np.pi*(r)**3


for i in range(datalength):
    print("Getting " + str(5+i))
    Folder = "SP_" + str(porosity) + "p_1m_" + str(a) + "a" + str(b) + "b_" + str(f) + "f/"
    File = str(5 + i) + "d/"
    PATH = "./VS_regbk_rock_" + str(r) + "r/" + Folder + File + "/s_param"
    df, SP = SP_df_gen(PATH, num_modes=1)
    S2sum.append(SP[3][0])
    ddf.append((df, SP))


print(d_xf)
M =[]
S2f = []
optical_depth_XF = []
tau = []
for i in range(datalength):
    S2f.append(ddf[i][1][1][0])

    M = ddf[i][1][1][0]

    ##########################################################
    # Optical depth can be defined as natural log of incident to transmitted power
    # In XF the scattering parameters correspond to the amount of the incident power
    # that scatters into a particular mode inside a waveguide. Therefore any forward
    # scattering will result in transmited power. What we are really interested in
    # is the amount still in the original mode. Then the optical depth is -ln(Sf11^2)

    tau.append(-np.log(M))
print('length of optical:', len(optical_depth_XF))

print("S forward", S2f)
#print(type(optical_depth_XF[:][0]))


ks = 3*porosity*Qs_M/(4*r)
optical_depth_M = ks*d

res = stats.linregress(d_xf, tau)

bf_line = res.intercept + res.slope*d
mu, std = stats.norm.fit(np.divide(tau, d_xf))
print("average value is: ", mu)
print("standard deviation is: ", std)

plt.plot(d, optical_depth_M, label = "MATLAB optical depth, ks = " + str(ks))
plt.plot(d_xf, tau, "D", label = "XFdtd optical depth")

plt.plot(d, res.intercept + res.slope*d, 'r', label='fitted line, slope ks = ' + str(res.slope) )

plt.ylabel("Optical Depth, ts [unitless]")
plt.xlabel("Volume Thickness, d [m]")
plt.title("Optical depth for " + str(r) +"m radius sphere at " + str(f) + "MHz, porosity = " +str(porosity))
plt.legend()

plt.show()

plt.plot(d_xf, np.divide(tau, d_xf), 'o')
plt.axhline(y = mu, xmin = 0, xmax = 1, color = 'r', linewidth = 1.0)
plt.show()

n, bins, patches = plt.hist(x= np.divide(tau, d_xf), bins='auto', color='#0504aa',
                                        alpha=0.7, rwidth=0.85)

xmin, xmax = plt.xlim() 
x = np.linspace(xmin, xmax, 100) 
p = std*stats.norm.pdf(x, mu, std)

#std_g = r/(1.2*a*b*100*porosity)
#p_g = std_g*stats.norm.pdf(x, ks, std_g)
  
plt.plot(x, n.max()*np.sqrt(2*np.pi)*p, 'k', linewidth=2) 
plt.axvline(x = ks, color = 'red', linestyle = '-', label = "real ks = " + str(np.round(ks, 4)))
#plt.plot(x, n.max()*np.sqrt(2*np.pi)*p_g, 'k', linewidth=2)
plt.xlabel("scattering coefficient, ks [1/m]")
plt.title("calculated scattering coeffient, mean = " + str(np.round(mu, 4)) + " std = " + str(np.round(std, 4)))
plt.legend()

plt.show()