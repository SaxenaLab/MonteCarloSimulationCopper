import numpy as np
import math
import random
from matplotlib import pyplot as plt
from scipy.spatial.transform import Rotation as R

print('What do you want to title this session? We will use the title for the output files. Please don\'t use period (.) in the title.')
title=input()
textfile=open(title+'_inputs.txt','w')

#Set the g (g_perp,g_parallel) and A (A_perp,A_parallel) values of the Cu(II) spin. The A values are in the units of Gauss
g=np.array([2.05704,2.2773])
a=np.array([10.7,161.6])

def B(g_eff):
    h=6.626*10**-34 #m^2*kg/s
    v=34.15*10**9 #1/s
    be=9.274*10**-24 #kg*m^/s^2*T
    return h*v/(g_eff*be)*10000

def ave_g(angle):
    return np.sqrt((g[1]*np.cos(angle))**2+(g[0]*np.sin(angle))**2)
def ave_a(angle):
    return np.sqrt((a[1]*g[1]**2*np.cos(angle))**2+(a[0]*g[1]**2*np.sin(angle))**2)/((g[1]*np.cos(angle))**2+(g[0]*np.sin(angle))**2)

def lorentzian(x, mu, sig):
    return sig**2/((x-mu)**2+sig**2)

def find_index(x_observer,x_points):
    difference=abs(x_observer-x_points)
    return np.argmin(difference)

#Parameters for building the in silico sample. Angle parameters are in units of degrees and distance are in units of nm.

print('Input chi (in degrees). Format: mean,sd. For example: 40,10')
chi_input=input().split(',')
while len(chi_input)<2:
    print('2 inputs are needed')
    print('Input chi (in degrees). Format: mean,sd. For example: 40,10')
    chi_input=input().split(',')
textfile.write('chi_input\t'+chi_input[0]+','+chi_input[1]+'\n')
chi=float(chi_input[0])
chi_sd=float(chi_input[1])
print('Input gamma (in degrees). Format: mean,sd. For example: 40,10')
gamma_input=input().split(',')
while len(gamma_input)<2:
    print('2 inputs are needed')
    print('Input gamma (in degrees). Format: mean,sd. For example: 40,10')
    gamma_input=input().split(',')
textfile.write('gamma_input\t'+gamma_input[0]+','+gamma_input[1]+'\n')
gamma=float(gamma_input[0])
gamma_sd=float(gamma_input[1])
print('Input eta (in degrees). Format: mean,sd. For example: 40,10')
eta_input=input().split(',')
while len(eta_input)<2:
    print('2 inputs are needed')
    print('Input eta (in degrees). Format: mean,sd. For example: 40,10')
    eta_input=input().split(',')
textfile.write('eta_input\t'+eta_input[0]+','+eta_input[1]+'\n')
eta=float(eta_input[0])
eta_sd=float(eta_input[1])
print('Input distance (in nm). Format: mean,sd. For example: 5.3,0.2')
r_input=input().split(',')
while len(r_input)<2:
    print('2 inputs are needed')
    print('Input distance (in nm). Format: mean,sd. For example: 5.3,0.2')
    r_input=input().split(',')
textfile.write('r_input\t'+r_input[0]+','+r_input[1]+'\n')
r_distance=float(r_input[0])
r_sd=float(r_input[1])
print('Input number of molecules you want to simulate. Note, high number of molecules will take require longer simulation time')
n_molecules=int(input())
textfile.write('n_molecules\t'+str(n_molecules)+'\n')

#Parameters for counting spins and simulating field sweep
print('Input broadening parameter (beta)')
beta_input=input()
textfile.write('beta_input\t'+beta_input+'\n')
beta=float(beta_input)



#Build n_molecules number of Spin A oriented to the z-axis
spin_a=np.zeros((n_molecules,3))
spin_a[:,2]=1

#Sample chi, gamma, and eta from gaussian distributions
chi_rot=R.from_euler('y',np.random.normal(chi,chi_sd,size=n_molecules),degrees=True)
gamma_rot=R.from_euler('y',np.random.normal(gamma,gamma_sd,size=n_molecules),degrees=True)
eta_rot=R.from_euler('z',np.random.normal(eta,eta_sd,size=n_molecules),degrees=True)

#Sample the interspin distance from a gaussian distribution
r_mag=np.random.normal(r_distance,r_sd,size=n_molecules)

#Build the interspin vector and Spin B for each Spin A with the sampled chi, gamma, and eta
r=np.zeros((n_molecules,3))
spin_b=np.zeros((n_molecules,3))
for i in range(len(spin_a)):
    r[i]=np.dot(chi_rot[i].as_matrix(),spin_a[i])
    spin_b[i]=np.dot(gamma_rot[i].as_matrix(),spin_a[i])
    spin_b[i]=np.dot(eta_rot[i].as_matrix(),spin_b[i])

### Randomly orienting your spin-pairs
orientations=R.random(n_molecules).as_matrix()
spin_a2=np.zeros((n_molecules,3))
r2=np.zeros((n_molecules,3))
spin_b2=np.zeros((n_molecules,3))
for i in range(len(spin_a)):
    r2[i]=np.dot(orientations[i],r[i])
    spin_a2[i]=np.dot(orientations[i],spin_a[i])
    spin_b2[i]=np.dot(orientations[i],spin_b[i])

### Calculating the effective g, A, and resonant fields
z_angle_a=np.arccos(spin_a2[:,2])
z_angle_b=np.arccos(spin_b2[:,2])
all_g_a=ave_g(z_angle_a)
all_g_b=ave_g(z_angle_b)
all_a_a=ave_a(z_angle_a)
all_a_b=ave_a(z_angle_b)

fields_a=B(all_g_a)
fields_b=B(all_g_b)

x_values=np.linspace(10000,13000,1024)

split_a=[]
for i in range(len(fields_a)):
    peaks=[fields_a[i]-1.5*all_a_a[i],fields_a[i]-0.5*all_a_a[i],fields_a[i]+0.5*all_a_a[i],fields_a[i]+1.5*all_a_a[i]]
    split_a.append(peaks)
split_b=[]
for i in range(len(fields_b)):
    peaks=[fields_b[i]-1.5*all_a_b[i],fields_b[i]-0.5*all_a_b[i],fields_b[i]+0.5*all_a_b[i],fields_b[i]+1.5*all_a_b[i]]
    split_b.append(peaks)
    
components_a=np.zeros((len(split_a),1024))
for i in range(len(split_a)):
    line=np.zeros(1024)
    for k in split_a[i]:
        gauss=lorentzian(x_values,k,beta)
        line+=gauss
    components_a[i]=line
components_b=np.zeros((len(split_b),1024))
for i in range(len(split_b)):
    line=np.zeros(1024)
    for k in split_b[i]:
        gauss=lorentzian(x_values,k,beta)
        line+=gauss
    components_b[i]=line

    
plt.figure(figsize=(12,10))
threshold=np.arange(0,2.5,0.01)
n_counter=np.zeros(len(threshold))
fields=[11000,11200,11400,11600,11800]
for j in fields:
    field_index=find_index(j,x_values)
    for i in range(len(n_counter)):
        check_a=components_a[:,field_index]>threshold[i]
        count_a=len(np.where(check_a==True)[0])
        check_b=components_b[:,field_index]>threshold[i]
        count_b=len(np.where(check_b==True)[0])
        n_counter[i]=count_a+count_b
    plt.plot(threshold,n_counter,linewidth=5,label=str(j)+' G')
plt.ylim((0,20000))
plt.yticks(ticks=[])
plt.xlabel('$\\alpha$',fontsize=30)
plt.xticks(fontsize=25)
plt.ylabel('Number of\nSpins Excited',fontsize=30)
plt.legend(fontsize=25)
plt.savefig(title)