import numpy as np
import random
from matplotlib import pyplot as plt
from scipy.spatial.transform import Rotation as R

print('What do you want to title this session? We will use the title for the output files. Please don\'t use period (.) in the title.')
title=input()
textfile=open(title+'_inputs.txt','w')

#Set the g (g_perp,g_parallel) and A (A_perp,A_parallel) values of the Cu(II) spin. The A values are in the units of Gauss
g=np.array([2.05704,2.2773])
a=np.array([10.7,161.6])
print('Input frequency (in GHz):')
v_input=input()
textfile.write('frequency\t'+v_input+'\n')
v=float(v_input)*10**9 #1/s


def B(g_eff):
    #Calculates the resonant field given a the effective g value
    #Parameter: g_eff=effective g value. See ave_g() for help with obtaining the effective g value of a Cu(II) spin
    #Return: Resonant field in units of Gauss
    h=6.626*10**-34 #m^2*kg/s
    be=9.274*10**-24 #kg*m^/s^2*T
    return h*v/(g_eff*be)*10000

def ave_g(phi):
    #Calculates the effective g value given an angle
    #Parameter: phi=angle between the g_parallel of the spin and the applied magnetic field in units of radians
    #Return: Effective g value
    return np.sqrt((g[1]*np.cos(phi))**2+(g[0]*np.sin(phi))**2)

def ave_a(phi):
    #Calculates the effective A value given an angle
    #Parameter: phi=angle between the g_parallel of the spin and the applied magnetic field in units of radians
    #Return: Effective A value in units of Gauss
    return np.sqrt((a[1]*g[1]**2*np.cos(phi))**2+(a[0]*g[1]**2*np.sin(phi))**2)/((g[1]*np.cos(phi))**2+(g[0]*np.sin(phi))**2)

def lorentzian(x, mu, beta):
    #Obtains the lorentzian lineshape representing the inhomogeneous broadening of the EPR lineshape
    #Parameter:
    #x=an array of size (1,N). x should represent the range of fields of interest
    #mu=an integer or float value. mu should represent the resonant field of the spin
    #beta=an integer or float value. beta represents the broadening parameter of the Lorentzian lineshape
    #Return:
    #an array of size (1,N) that represents the intensity of the Lorentzian signal in arbitrary units
    return beta**2/((x-mu)**2+beta**2)

def find_index(x_observer,x_points):
    #Finds the index for a value x closest to x_observer (integer or float value) within x_points (array of length (1,N))
    difference=abs(x_observer-x_points)
    return np.argmin(difference)

def FS(x,res_fields_a,res_fields_b):
    #Simulates a Field-Swept spectrum from two arrays of resonant fields
    #Parameter:
    #x = an array of (1,N) dimension. x should represent the range of fields of interest with N number of data points.
    #res_fields_a and res_fields_b = arrays of (J,K) dimension where J is the number of Spin A and Spin B and K is 
    #                                the number of resulting resonant fields due to the hyperfine splitting 
    #                                (4 for the case of Cu(II)). Note, both arrays must be equal in dimensions.
    #Return:
    #fs = an array of (1,N) dimension that represents the intensity of the Field-Swept spectrum in arbitrary units
    #counts = an arrays of (1,N) dimension that represents the number of phi excitable at a given field. Note, the number is
    #         dependent on the beta parameter of the lorentzian lineshape and the threshold that distringuishes significant
    #         and non-significant intesity of the lorentzian lineshape
    #components_a = an array of (J,N) dimension that represents the lorentzian lineshape for each of J number of Spin A.
    #components_b = an array of (J,N) dimension that represents the lorentzian lineshape for each of J number of Spin B.
    components_a=np.zeros((len(res_fields_a),len(x)))
    counts=np.zeros(1024)
    fs=np.zeros(1024)
    for i in range(len(res_fields_a)):
        line=np.zeros(1024)
        for k in res_fields_a[i]:
            gauss=lorentzian(x,k,beta)
            line+=gauss
        components_a[i]=line
        counts+=line>alpha
        fs=fs+line
    components_b=np.zeros((len(res_fields_b),len(x)))
    for i in range(len(res_fields_b)):
        line=np.zeros(1024)
        for k in split_b[i]:
            gauss=lorentzian(x_values,k,beta)
            line+=gauss
        components_b[i]=line
        counts+=line>alpha
        fs=fs+line
    return fs,counts,components_a,components_b

def leftover(spin_a,spin_b,r,res_fields_a,res_fields_b,excited_index):
    #Returns a new set of arrays that represents the spin pairs that did not contribute towards the DEER intramolecular
    #signal
    #Parameter:
    #spin_a and spin_b = arrays of (J,3) dimension that represent the three-dimensional orientations of g_parallel
    #                    as vectors for each of J number of Spin A and Spin B. Both arrays must have the same dimensions.
    #r = an array of (J,3) dimension that represents the three-dimensional orientations of the interspin vectors of each
    #    of J number of spin pairs.
    #res_fields_a and res_fields_b = arrays of (J,K) dimension where J is the number of Spin A and Spin B and K is 
    #                                the number of resulting resonant fields due to the hyperfine splitting 
    #                                (4 for the case of Cu(II)). Note, both arrays must be equal in dimensions.
    #excited_index = an array of (1,L) that represents the index on both Spin A and Spin B that contributes to the DEER 
    #                intramolecular signal.
    #Return:
    #spin_a2 and spin_b2 = arrays of (J-L,3) dimension that represent the three-dimensional orientations of g_parallel
    #                      as vectors for each of J-L number of Spin A and Spin B that did not contribute towards the DEER
    #                      signal.
    #r2 = an array of (J-L,3) dimension that represents the three-dimensional orientations of the interspin vectors of each
    #     of J-L number of spin pairs that did not contribute towards the DEER signal.
    #res_fields_a2 and res_fields_b2 = arrays of (J-L,K) dimension where J-L is the number of Spin A and Spin B that did
    #                                  not contribute towards the DEER signal and K is the number of resulting resonant
    #                                  fields due to the hyperfine splitting (4 for the case of Cu(II)).
    spin_a2=np.delete(spin_a,excited_index,0)
    spin_b2=np.delete(spin_b,excited_index,0)
    r2=np.delete(r,excited_index,0)
    res_fields_a2=np.delete(res_fields_a,excited_index,0)
    res_fields_b2=np.delete(res_fields_b,excited_index,0)
    return spin_a2,spin_b2,r2,res_fields_a2,res_fields_b2

def excited_spins(a,b,obs,pump):
    #Identifies the index of Spin A and Spin B that contributes to the DEER intramolecular signal given an observer pulse
    #with a bandwidth of 38 MHz and a pump pulse with a bandwidth of 100 MHz
    #Parameter:
    #a and b = an array of (J,N) dimension that represents the lorentzian lineshape for each of J number of Spin A and 
    #          Spin B. See FS() for more details
    #obs = The center of the observer field in units of Gauss excitable by the observer pulse. 
    #pump = The center of the pump field in units of Gauss excitable by the pump pulse.
    #Return:
    #excited_index = an array of (1,L) that represents the index on both Spin A and Spin B that contributes to 
    #                the DEER intramolecular signal
    relevant_comp_a=set()
    relevant_comp_b=set()
    obs_a=set()
    pump_a=set()
    obs_b=set()
    pump_b=set()
    for i in range(5):
        check=a[:,obs-2+i]>alpha
        index=np.where(check==True)
        set_a=set(index[0])
        relevant_comp_a=set.union(relevant_comp_a,set_a)
        obs_a=set.union(obs_a,set_a)
        check=b[:,obs-2+i]>alpha
        index=np.where(check==True)
        set_b=set(index[0])
        relevant_comp_b=set.union(relevant_comp_b,set_b)
        obs_b=set.union(obs_b,set_b)
    for i in range(13):
        check=b[:,pump-6+i]>alpha
        index=np.where(check==True)
        set_b=set(index[0])
        relevant_comp_b=set.union(relevant_comp_b,set_b)
        pump_b=set.union(pump_b,set_b)
        check=a[:,pump-6+i]>alpha
        index=np.where(check==True)
        set_a=set(index[0])
        relevant_comp_a=set.union(relevant_comp_a,set_a)
        pump_a=set.union(pump_a,set_a)

    intersection1=set(list(pump_a.intersection(obs_b)))
    intersection2=set(list(pump_b.intersection(obs_a)))
    total_intersection=set.union(intersection1,intersection2)
    excited_index=list(total_intersection)
    print('Excited Spin Pairs',len(excited_index))
    return excited_index

def dipolar_freq(r_vec,r_mag):
    #Calculates the dipolar frequency for a given interspin vector
    #Parameter:
    #r_vec = an array of (1,3) dimension that represents the orientation of the interspin vector
    #r_mag = a value that represents the interspin distance
    #Return:
    #v_dipolar = dipolar frequency in units of Hz
    be=9.274*10**-24 #kg*m^/s^2*T
    mu=1.25663706212*10**-6 #kg*m/s^2*A^2
    h=6.626*10**-34 #m^2*kg/s
    return mu/4/np.pi/h*((g[0]*2+g[1])/3)**2*be**2/(r_mag*10**-9)**3*(3*(np.cos(np.arccos(r_vec[2])))**2-1)

def DEER_signal(t,frequency):
    #Obtains the intramolecular DEER signal from the dipolar frequency
    #Parameter:
    #t = an array of (1,M) dimension that represents the dipolar evolution time with M data points
    #frequency = dipolar frequency in units of Hz
    #Return:
    #P(t) = The intramolecular DEER signal
    return np.cos(t*10**-9*frequency*2*np.pi)

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
print('Input broadening parameter (beta) and threshold (alpha). We recommend 40,0.4 if you are unsure. Format: beta,alpha. For example: 40,0.4')
beta_alpha_input=input().split(',')
while len(beta_alpha_input)<2:
    print('2 inputs are needed')
    print('Input broadening parameter (beta) and threshold (alpha). We recommend 40,0.4 if you are unsure. Format: beta,alpha. For example: 40,0.4')
    beta_alpha_input=input().split(',')
textfile.write('beta_alpha_input\t'+beta_alpha_input[0]+','+beta_alpha_input[1]+'\n')
beta=float(beta_alpha_input[0])
alpha=float(beta_alpha_input[1])

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

split_a=np.zeros((n_molecules,4))
for i in range(len(fields_a)):
    peaks=[fields_a[i]-1.5*all_a_a[i],fields_a[i]-0.5*all_a_a[i],fields_a[i]+0.5*all_a_a[i],fields_a[i]+1.5*all_a_a[i]]
    split_a[i]=np.array(peaks)
split_b=np.zeros((n_molecules,4))
for i in range(len(fields_b)):
    peaks=[fields_b[i]-1.5*all_a_b[i],fields_b[i]-0.5*all_a_b[i],fields_b[i]+0.5*all_a_b[i],fields_b[i]+1.5*all_a_b[i]]
    split_b[i]=np.array(peaks)
    
    
### Simulating field sweep
x_values=np.linspace(10000,13000,1024)
fs1,counts1,comp_a,comp_b=FS(x_values,split_a,split_b)


### Setting your observer and pump pulses
print('Input three pump-fields you want to use (in Gauss). For Cu(II)NTA-dHis spin-label, we recommend ~100 G, ~580 G, ~827 G lower than the maximum of your FS-ESE spectrum. Format: field1,field2,field3. For example: 11739,11240,11011')
field_input=input().split(',')
while len(field_input)<3:
    print('3 inputs are needed')
    print('Input three pump-fields you want to use (in Gauss). For Cu(II)NTA-dHis spin-label, we recommend ~100 G, ~580 G, ~827 G lower than the maximum of your FS-ESE spectrum. Format: field1,field2,field3. For example: 11739,11240,11011')
    field_input=input().split(',')
textfile.write('field_input\t'+field_input[0]+','+field_input[1]+','+field_input[2]+'\n')
pump_field=find_index(float(field_input[0]),x_values)
pump_field2=find_index(float(field_input[1]),x_values)
pump_field3=find_index(float(field_input[2]),x_values)

obs_field=find_index(float(field_input[0])-54,x_values)
obs_field2=find_index(float(field_input[1])-54,x_values)
obs_field3=find_index(float(field_input[2])-54,x_values)

### Excite spin-pairs
relevant_index1=excited_spins(comp_a,comp_b,obs_field,pump_field)
relevant_index2=excited_spins(comp_a,comp_b,obs_field2,pump_field2)
relevant_index3=excited_spins(comp_a,comp_b,obs_field3,pump_field3)


### Simulate DEER signal
sample_freq=np.zeros(n_molecules)
print('Input dipolar evolution time of your DEER signal (in ns) and number of points. Format: evolution_time,#_points. For example: 7000,210')
acq_time_input=input().split(',')
while len(acq_time_input)<2:
    print('2 inputs are needed')
    print('Input dipolar evolution time of your DEER signal (in ns) and number of points. Format: evolution_time,#_points. For example: 7000,210')
textfile.write('acq_time\t'+acq_time_input[0]+','+acq_time_input[1]+'\n')
acq_time=np.linspace(0,float(acq_time_input[0]),int(acq_time_input[1]))
sample_signal=np.zeros((n_molecules,len(acq_time)))
for i in range(len(sample_freq)):
    sample_freq[i]=dipolar_freq(r2[i],r_mag[i])
    sample_signal[i]=DEER_signal(acq_time,sample_freq[i])
total_signal=np.sum(sample_signal,0)
total_signal=total_signal/total_signal[0]

excited_freq=np.zeros(len(relevant_index1)+len(relevant_index2)+len(relevant_index3))
excited_signal=np.zeros((len(relevant_index1)+len(relevant_index2)+len(relevant_index3),len(acq_time)))
k=0
for i in range(len(relevant_index1)):
    excited_freq[k]=dipolar_freq(r2[relevant_index1[i]],r_mag[relevant_index1[i]])
    excited_signal[k]=DEER_signal(acq_time,excited_freq[k])
    k+=1
for i in range(len(relevant_index2)):
    excited_freq[k]=dipolar_freq(r2[relevant_index2[i]],r_mag[relevant_index2[i]])
    excited_signal[k]=DEER_signal(acq_time,excited_freq[k])
    k+=1
for i in range(len(relevant_index3)):
    excited_freq[k]=dipolar_freq(r2[relevant_index3[i]],r_mag[relevant_index3[i]])
    excited_signal[k]=DEER_signal(acq_time,excited_freq[k])
    k+=1
    
total_excited_signal=np.sum(excited_signal,0)
total_excited_signal=total_excited_signal/total_excited_signal[0]

### Plot DEER signal
print('Would you like to add noise? (y/n)')
noise_check=input()
textfile.write('noise_check\t'+noise_check+'\n')
if noise_check=='y':
    print('Input level of noise as a value between 0 and 1')
    noise_level=float(input())
    textfile.write('noise_level\t'+str(noise_level)+'\n')
    noise_ideal = np.random.normal(0,noise_level,len(acq_time))
    noise_excited = np.random.normal(0,noise_level,len(acq_time))
    plt.rcParams.update({'font.size': 35})
    plt.figure(figsize=(9,8))
    plt.plot(acq_time,total_signal+noise_ideal,linewidth=5,label='ideal')
    plt.plot(acq_time,total_excited_signal+noise_excited,linewidth=5,label='excited')
    plt.xlabel('time (ns)')
    plt.ylabel('$V/V_0$')
    plt.gca().get_yaxis().set_visible(False)
    plt.legend(frameon=False)
    plt.xticks(fontsize=25)
    plt.savefig(title+'_plot')
    np.savetxt('title'+'_output.txt',np.transpose(np.array([acq_time,(total_excited_signal/total_excited_signal[0])])))
else:    
    plt.rcParams.update({'font.size': 35})
    plt.figure(figsize=(9,8))
    plt.plot(acq_time,total_signal,linewidth=5,label='ideal')
    plt.plot(acq_time,total_excited_signal,linewidth=5,label='excited')
    plt.xlabel('time (ns)')
    plt.ylabel('$V/V_0$')
    plt.gca().get_yaxis().set_visible(False)
    plt.legend(frameon=False)
    plt.xticks(fontsize=25)
    plt.savefig('DEER '+'chi'+str(chi)+'gamma'+str(gamma)+'eta'+str(eta))
    np.savetxt('title'+'_output.txt',np.transpose(np.array([acq_time,(total_excited_signal/total_excited_signal[0])])))

textfile.close()

   

