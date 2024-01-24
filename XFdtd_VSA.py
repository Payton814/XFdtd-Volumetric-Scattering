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
# TODO:
#  - render plots with LaTeX labels
#  - get rid of unused variables
#  - add plot titles
#  - expand plot range of spread to emphasize the small deviation. current plot doesn't tell us much visually
#  - figure out what ddf[], d_xf[], and d[] are.
# REVISIONS: 
# 22-Jan-2024
#  - deleted declaration of ts_err array and Sf10 array.
#  - deleted calculation of sphere volume
#  - commented out creation of S2sum[]
#  - replaced while loop method of getting CS_M and Qs_M
# 23-Jan-2024
#  - deleted S2sum[] instances
#  - deleted old while loop method to calculate CS_M and Qs_M

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from S_param_df_gen import SP_df_gen

# MODIFIABLE PARAMETERS
r = 0.5             # scatterer radius, meters
porosity = 0.005    # sample porosity (VOID VOLUME / TOTAL VOLUME)
datalength = 33     # TODO: What is this for? Why choose this value?
a = 28              # waveguide width, meters
b = 14              # waveguide height, meters
f = 250             # wave frequency, MHz
step_size = 0.001   # radius increment size in .csv files produced by MATLAB

# read in CSV to np DataFrame df_matlab, calculate CS_M, Qs_M
df_matlab = pd.read_csv("./VS_regbk_rock_" + str(r) + "r/CS_regbk_rockbk_" + str(f) + "MHz_0.001_step.csv")

loc = r / step_size
CS_M = df_matlab['Qs'][loc]
Qs_M = CS_M / (np.pi*r**2)

# TODO: IMPROVE VARIABLE NAMES / DOCUMENTATION. No idea what any of these are currently.
ddf = []
numModes = 1        
d_xf = np.arange(5, datalength + 5, 1)
d = np.arange(0,datalength + 5,1)       

# why do we iterate datalength times?
for i in range(datalength):
    print("Getting " + str(5+i))
    
    # create PATH string to the .csv we want to parse with SP_df_gen(). Path and file name will change dependent on global variables
    Folder = "SP_" + str(porosity) + "p_1m_" + str(a) + "a" + str(b) + "b_" + str(f) + "f/"
    File = str(5 + i) + "d/"
    PATH = "./VS_regbk_rock_" + str(r) + "r/" + Folder + File + "/s_param"
    
    df, SP = SP_df_gen(PATH, num_modes=numModes)
    ddf.append((df, SP))

print(d_xf)

# unsure of what each of these arrays is for / whether they are needed.
M = []
S2f = []
optical_depth_XF = []
tau = []

# I believe we could do this within the first 'for' loop above, as long as we declare the above arrays prior to that loop.
for i in range(datalength):
    S2f.append(ddf[i][1][1][0])

    M = ddf[i][1][1][0]

    ##########################################################
    # Optical depth can be defined as natural log of incident to transmitted power
    # In XF the scattering parameters correspond to the amount of the incident power
    # that scatters into a particular mode inside a waveguide. Therefore any forward
    # scattering will result in transmitted power. What we are really interested in
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

# PLOT RESULTS
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

std_g = r/(1.2*a*b*100*porosity)
p_g = std_g*stats.norm.pdf(x, ks, std_g)
  
#plt.plot(x, n.max()*np.sqrt(2*np.pi)*p, 'k', linewidth=2) 
plt.axvline(x = ks, color = 'red', linestyle = '-', label = "real ks = " + str(np.round(ks, 4)))
#plt.plot(x, n.max()*np.sqrt(2*np.pi)*p_g, 'k', linewidth=2)
plt.xlabel("scattering coefficient, ks [1/m]")
plt.title("calculated scattering coeffient, mean = " + str(np.round(mu, 4)) + " std = " + str(np.round(std, 6)))
plt.legend()

plt.show()


std_l = [0.007405, 0.003457, 0.002695, 0.001981]

plt.plot([14*7, 20*10, 22*11, 28*14], std_l, linestyle = "--", marker = "o", label = "data")
plt.plot([14*7, 20*10, 22*11, 28*14], [r/(14*7*0.5), r/(20*10*0.5), r/(22*11*0.5), r/(28*14*0.5)], label = "r/(ab*porosity)")
plt.xlabel("ab [m^2]")
plt.ylabel("standard deviation")
plt.legend()
plt.show()
