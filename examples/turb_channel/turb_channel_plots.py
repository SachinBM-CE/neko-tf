import numpy as np
import matplotlib.pyplot as plt

# >>>>> INPUTS <<<<<
batch_time = 100
stress_file = 'Re2800_out' # out_DNS' 
csv_file = 'fluid_stats0.csv'

# ***********************************************************************************
# >>>>> Step 1: <<<<<
# Extraction of stress in stream-wise direction
# https://github.com/ExtremeFLOW/neko/blob/develop/examples/turb_channel/README.md
# assuming one saves the log in my_log) and put it into a textfile out
# awk '/forcex/ {print($1,$2,$3)} ' my_log > out
# ***********************************************************************************

# *******************************
# >>>>> Step 2: <<<<<
# Evolution of Re_tau with Time
# *******************************

Re_b = 2800 #Because we are non-dimensional this is ok
dat = np.genfromtxt(stress_file) # Read out into a numpy array

# slicing [start:stop:step]
# Every even-indexed row (starting from index 0) is ::2
# Every odd-indexed row (starting from index 1) is 1::2

# u_tau = sqrt(tau_w/rho) | Re_tau = u_tau*h/nu | Re_b = u_b*h/nu
u_tau_bottom = np.sqrt(dat[::2,2])
u_tau_top = np.sqrt(dat[1::2,2])
Re_tau_bottom = u_tau_bottom * Re_b
Re_tau_top = u_tau_top * Re_b

plt.plot(dat[::2,1],    Re_tau_bottom,    label='Bottom wall')
plt.plot(dat[1::2,1], Re_tau_top, '--', label='Top wall')
plt.ylabel(r'$Re_\tau$')

##plt.plot(dat[::2,1],u_tau_bottom,label='Bottom wall')
##plt.plot(dat[1::2,1],u_tau_top,'--',label='Top wall')
##plt.ylabel(r'$u_\tau$')

plt.xlabel(r'Time, $\delta/U_b$')
plt.legend()
plt.savefig('01_Re2800_Spalding_Re_tau_Time.png', dpi=300, bbox_inches='tight')
plt.show()

# ***********************
# >>>>> Step 3: <<<<<
# Turbulence Statistics
# ***********************

# dat = [output time, coordinate, <p>, <u>, <v>, <w>, <pp>, <uu>, <vv>, <ww>, <uv>, <uw>,...
dat = np.genfromtxt(csv_file, delimiter=',', invalid_raise=False)
dat = dat[~np.isnan(dat).any(axis=1)]  # Remove rows with NaN

#Time of batch of interest (our batchsize is 5 convective time units) and sampling started at 60
#First batch of average between 60 and 65 is therefore written out at T=65 (OBSERVE you need to have run beyond T=65)
# Meaning: in case file under "type": "fluid_stats", "output_control": "simulationtime", "output_value": 5, "start_time": 60
# Lets see how the mean profile looks

#Extract <u>
U_vel = dat[np.abs(dat[:,0]-batch_time)<0.1, 3]
print(f"Length of U_vel     : {len(U_vel)}")
print(f"Length of u_tau_top : {len(u_tau_top)}")

#coordinates
y_coords = dat[np.abs(dat[:,0]-batch_time)<0.1, 1]

#plot the profile
plt.plot(y_coords,U_vel)

#If you run longer you can compute better averages by adding many batches together, for somewhat converged statistics you will need at least 100 time units.
plt.ylabel(r'$u$ ($U_b$)')
plt.xlabel(r'$y$ ($\delta$)')
plt.legend()
plt.savefig('02_Re2800_Spalding_mean_velocity_profile.png', dpi=300, bbox_inches='tight')
plt.show()

# Plotting u+ vs. y+
##nu = 1/Re_b
##y_plus = np.array(y_coords)
