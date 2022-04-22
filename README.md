# MonteCarloSimulationCopper
This set of scripts is designed to perform Monte-Carlo based simulation of DEER of a Cu(II) system to solve orientational selectivity effects in Q-band DEER. In general, the script builds an in silico sample containing a number of molecules (defined by the user) doubly-labeled by Cu(II) with diverse orientations (defined by user) using Monte-Carlo approach. From this in silico sample, we can obtain important information that helps us obtain a minimal acquisition scheme for DEER.

Run the script by running the following command in a terminal:
$~ python [script.py]
Each script will prompt the user for any required inputs.

Each script provides a different role:
-Monte Carlo Alpha Determination:
  If you are unsure where to start, we recommend start with this script. We want the simulation to replicate the experimental FS-ESE spectrum by inputing
  the frequency of the experiment, the g and A values, and the broadening parameter (beta). However, once beta is determined, a proper alpha must also be determined.
  Alpha is a threshold parameter that dictates whether a spin can be considered to be excitable or not at a given field. This threshold must not be too low, leading
  to overcounting, but also not too high, leading to undercounting. To obtain a proper value of alpha, this script will output a plot of number of excitable spins as
  as function of alpha at different fields across the FS-ESE spectrum. By inspecting the plot, one can determine an optimal alpha that leads to stable counting of spins.
  
-Monte Carlo Phi Curve
  Once alpha and beta have been determined, this script allows you to identify the three fields that can excite the largest number of phi angles. The script will output
  a Phi curve. This curve describes the number of spins that can be excited at a given field. By identifying the maximum of the Phi curve, one can determine a field
  that can excite the largest number of orientations of the spin. For Cu(II), a single field is generally not sufficient to excite all possible angles. Therefore,
  the script will simulate a DEER at the maximum of the Phi curve and return the leftover spins that are not yet excited. From the unexcited spins, another Phi curve
  is plotted to identify another maximum field that can best excite the rest of the spins. This process are done three times in total to provide three fields that
  are most promising for efficiently exciting all orientations of the molecule.
  Note, the fields identified here should be used as the center field of the pump pulse in the DEER experiment. Additionally, the user may one to run this script in
  repeats as the outcome can be variable due to the incomplete of sampling from the Monte Carlo simulation. By rerunning the script, the user should consider fields
  that are consistently identified in multiple simulations.
  
-Monte Carlo DEER time trace simulation
  This script takes three magnetic fields and outputs a DEER time domain signal from performing DEER at the three fields. For comparison, the ouput also contains
  the ideal DEER time domain signal assuming all possible spin-pairs are excited. The ouutput can inform the user if the three identified fields can sufficiently
  minimze the effect of orientational selectivity for Cu(II) Q-band DEER.
  
Details of Input
-g and A values
  The script currently assumes g and A values for dHis-Cu(II)NTA spin-label. If the user wants to change these values to represent a different Cu(II) species, these
  values should be changed in the script itself. The A values specifically are in the units of Gauss.
-frequency
  The frequency of the DEER experiment should be inputted in units of GHz
-number of molecules
  We recommend to start with 10000 molecules.
-chi
  angle in degrees that represents the angle between g_parallel of Spin A and the interspin vector
-gamma
  angle in degrees that represents the angle between g_parallel of Spin A and the g_parallel of Spin B
-eta
  angle in degrees that represents the angle between g_perpendicular of Spin A and g_perpendicular of Spin B
-beta
  broadening parameter in units of Gauss
-alpha
  threshold parameter. If unsure, use the script 'Monte Carlo Alpha Determination.py' to get obtain a correct value

More details can be seen in the script itself. If you have any questions, please refer to xxx. Additionally, if you have any questions, feel free to contact
Zikri Hasanbasri at zih12@pitt.edu
