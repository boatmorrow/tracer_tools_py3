
import numpy as np
import tracer_tools.noble_gas_tools as ng
import pdb
import pandas as pd
import scipy.integrate as intgrt
from scipy.interpolate import interp1d
import importlib.resources as pkg_resources
import os
import scipy.constants as scic
import tracer_tools.He_tools as He
import matplotlib.pyplot as plt
import copy 

#set file path for saving results
parentfolder = 'C:/Users/Fanka Neumann/Documents/Daten/Production_rates'
sample_name = 'Aro_1'
file_path = parentfolder + '/' + sample_name

#----------------------------------------------------

# defining a rock sample

#input_composition = {'Li':3.4e-5,  'C':2.3e-4, 'O':0.494, 'F':1.28e-3, 'Na':0.03, 'Mg':4.72e-3, 'Al':0.0815, \
#                            'Si': 0.318, 'P':  6.11e-4, 'K':0.0370, 'Ca':0.0130, 'Ti': 2.84e-3, 'Fe':0.0185,\
#                            'Sr':4.75e-4, 'Ba': 1.88e-3, 'Th': 2.47e-5, 'U':  1.95e-6 } #Granite from Sramek

input_composition = {'N':8.3e-5,  'C':1.36e-3, 'O':0.48, 'F':5.57e-4, 'Na':0.0243, 'Mg':0.015, 'Al':0.0815, \
                            'Si': 0.311, 'P':  6.55e-4, 'K':0.0232, 'Ca':0.0257, 'Ti': 3.84e-3, 'Fe':0.0385,\
                            'Cl':3.7e-4, 'Mn': 7.74e-4, 'Th': 1.05e-5, 'U':  2.7e-6 } #Average upper crust from Sramek

input_composition['K']=0.0122
input_composition['Ca'] = 0.0082

rock = He.rock_type(composition=input_composition, MC_flag=1)
input_density = 2.6 #g/cm^3
rock.density = input_density
input_avg_depth_cm = 3.5 #cm average sample depth
input_avg_depth = input_density*input_avg_depth_cm # g/cm^2

#give depth range
rock.APR.depth = np.arange(0,input_avg_depth,input_avg_depth/50) #g/cm^2
# enable EXPACS calculation 
rock.APR.excel=1 
#Setting location parameters
input_elev = 2387
rock.APR.elev = input_elev
input_latitude = 46.01824
rock.APR.latitude = input_latitude
input_w_value = 50
rock.APR.w_value=input_w_value
input_soilmoisture = 0.25
rock.APR.soilmoisture=input_soilmoisture
input_shielding = 0.977
rock.APR.top_shielding = input_shielding

# define the function to reset parameters at the end of one MC iteration
def resetting_rock_parameters():
    rock.density = input_density
    rock.APR.depth = np.arange(0,input_avg_depth,input_avg_depth/50) #g/cm^2
    rock.APR.elev = input_elev
    rock.APR.latitude = input_latitude
    rock.APR.w_value=input_w_value  
    rock.APR.soilmoisture=input_soilmoisture
    rock.APR.top_shielding = input_shielding
    rock.composition = input_composition

#---------------------------------------------------------------

# running the calculation for the best estimate
rock.calcAr_prod_rate_whole_rock(solar='avg',calc_flux=1)

# calculating sample average production rate (Ar39 atoms/gram rock / year)
Ar39_avg_prod_rate = intgrt.trapezoid(rock.APR.Ar39_value, rock.APR.depth)/input_avg_depth

#-------------------------------------------------------------

# setting the uncertainties (1 sigma relative or absolute deviation)
avg_depth_cm_error = 0.5 #abs
density_error = 0.3 #abs.
soil_moisture_error = 0.05 #absolute deviation
w_number_error = 40 #abs.
altitude_error = 0.5 #abs
shielding_factor_error = 0.001 #abs
xsec_relative_intensity = 0.28 # basen on Sramek 2017
ar_yield_per_mu_uncy = 0.1 #rel
mu_stop_rate_uncy = 0.3 #rel
n_yield_per_mu_uncy = 0.1 #rel
P_n_UTh_uncy = 0.3 #rel
K_uncy = 0.07 #rel
Ca_uncy = 0.07 #rel
comp_uncy = 0.3 #rel

# flag for varying depth         
depth_flag = 1 
#if flag=0 depth is not varied to calculate the uncertainties for a depth-dependant production rate profile
#if flag not 0 depth is varied to account for the uncertainty of a specific sample geometry. Only sample averaged production rate is saved. 

#----------------------------------

# setting up the file for saving the results

def write_to_file(file_path, line):
    # Open the file in append mode (creates the file if it doesn't exist)
    with open(file_path, 'a') as file:
        # Write the line of text to the file
        file.write(line + '\n')

# Check if the file exists and create a header, if it does not

file_exists = os.path.exists(file_path)

# Open the file in append mode and create the header if necessary
with open(file_path, 'a') as file:
    if not file_exists:
        # If the file doesn't exist, write the header
        #sample name
        file.write("Sample name" + sample_name + "\n" + 'Parameters' + "\n")
        # list of sample parameters
        file.write("density [g/cm3], avg. depth [g/cm2], elev [masl], lat. [degr], w_value [1], soil moist. [m3/m3], shielding factor [1] \n")
        parameter_string = str(rock.density) + ", " + str(input_avg_depth) + ', ' + str(rock.APR.elev) + ', ' + str(rock.APR.latitude)  \
        + ', ' + str(rock.APR.w_value) + ', ' + str(rock.APR.soilmoisture) + ', ' + str(rock.APR.top_shielding) + '\n'
        file.write(parameter_string + 'Composition: \n')
        # Join keys of the composition dictionary into a single string 
        keys_string = ', '.join(input_composition.keys())
        # Join values into a single string separated by commas
        values_string = ', '.join(map(str, input_composition.values()))
        # write down the composition
        file.write(keys_string + '\n' + values_string + '\n')
        # writing down the best estimate of the production rate (not MC iterated)
        file.write('Argon 39 production rate [atoms/g/yr] (avg, depth-dependant) \n')
        depth_str = 'avg, '+ np.array2string(rock.APR.depth, separator=', ', max_line_width=800)+'\n'
        file.write(depth_str)
        prod_rate_str = str(Ar39_avg_prod_rate) + ', ' + np.array2string(rock.APR.Ar39_value, separator=', ', max_line_width=800) + '\n'
        file.write(prod_rate_str)
        # writing down the uncertainties
        file.write('1 sigma - uncertainties, relative value unless specified \n')
        file.write("Depth[cm] (abs), Density (abs), Soil moist. (abs), W-Nr. (abs), Shielding factor (abs)," + \
        " Cross-sec (rel), Ar39 yield/mu, Muon stop. rate, N yield/mu, P_n(U/Th), K conc., Ca conc., Composition \n")
        uncertainty_parameters_str = str(avg_depth_cm_error) + ', ' + str(density_error) + ', ' + str(soil_moisture_error) \
        + ', ' + str(w_number_error) + ', ' + str(altitude_error) + ', ' + str(shielding_factor_error) + ', ' \
        + str(xsec_relative_intensity) + ', ' + str(ar_yield_per_mu_uncy) + ', ' + str(mu_stop_rate_uncy) + ', ' \
        + str(n_yield_per_mu_uncy) + ', ' + str(P_n_UTh_uncy) + ', ' + str(K_uncy) + ', ' + str(Ca_uncy) + ', ' + \
        str(comp_uncy) + '\n'
        file.write(uncertainty_parameters_str)
        file.write('Monte Carlo iterations of the Ar39 production rate [atoms/g/yr] \n')
        if depth_flag == 0:
            file.write(depth_str)
        else:
            file.write ('sample average \n')
#---------------------------------------------

# Monte Carlo Iterations

# read in the the data which stays unchanged:

current_dir = os.path.dirname(os.path.abspath(__file__)) #path to current folder (example)
parent_dir = os.path.dirname(current_dir) #path to parent directory (tracer_tools_py3)
# path to the xsec-files in the parallel folder
# read in K39-cross-section: 
inp_file = os.path.join(parent_dir, 'data', 'K39_XS.csv')
K39xsec = pd.read_csv(inp_file)
print(K39xsec.shape) #should be (660,2)
#read in Ca42 cross-section
inp_file = os.path.join(parent_dir, 'data', 'Ca42_XS.csv')
Ca42xsec = pd.read_csv(inp_file)

#read in spectral data: normalized differential flux from Felsenbergkeller (Grieger et al. 2021)
inp_file = os.path.join(parent_dir, 'data', 'Differential_neutron_flux_normalized.csv')
diff_flux_normal = pd.read_csv(inp_file) #neutrons/cm^2/s/MeV
mu_shape_deviation = abs(diff_flux_normal['norm_phi_n_mu_avg']-diff_flux_normal['norm_phi_n_mu_MK1'])
alpha_shape_deviation = abs(diff_flux_normal['norm_phi_n_alpha_avg']-diff_flux_normal['norm_phi_n_alpha_MK1'])

def gaussian_noise_weights(df, column_name, relative_intensity):
    column_values = df[column_name]
    noise = np.random.normal(1, relative_intensity, len(column_values)) #mean 1, standard deviation, length of array
    noise[noise < 0] = 0 # set negativ values to zero
    return noise

def mc_iteration():
    # vary sample parameters
    rock.density = abs(np.random.normal(rock.density, density_error))
    if depth_flag == 1:
        avg_depth_cm = abs(np.random.normal(input_avg_depth_cm, avg_depth_cm_error))
        avg_depth = avg_depth_cm*rock.density
        rock.APR.depth = np.arange(0,avg_depth,avg_depth/50) #g/cm^2
    rock.APR.soilmoisture = abs(np.random.normal(rock.APR.soilmoisture, soil_moisture_error))
    rock.APR.w_value = abs(np.random.normal(rock.APR.w_value, w_number_error))
    rock.APR.elev = abs(np.random.normal(rock.APR.elev, altitude_error))
    rock.APR.top_shielding = abs(np.random.normal(rock.APR.top_shielding, shielding_factor_error))

    # vary composition
    exclude_entries = ['K', 'Ca', 'O']
    varied_composition = copy.deepcopy(input_composition)
    for key, value in varied_composition.items():
        if key not in exclude_entries:
            value = abs(value*np.random.normal(1,comp_uncy))
    varied_composition['K']=abs(varied_composition['K']*np.random.normal(1, K_uncy))
    varied_composition['Ca'] = abs(varied_composition['Ca']*np.random.normal(1, Ca_uncy))
    varied_composition['O'] = 1 - sum(varied_composition.values()) + varied_composition['O'] #Setting oxygen so total composition is 1
    rock.composition = varied_composition

    # vary general production rate parameters 

    # cross-sections 
    #add variation to cross-section
    rock.APR.K39xsec_noise = gaussian_noise_weights(K39xsec, 'Cross-Section [barns]', xsec_relative_intensity) # uncertainty set above
    rock.APR.Ca42xsec_noise = gaussian_noise_weights(Ca42xsec, 'Cross-Section [barns]', xsec_relative_intensity)

    # muon parameters
    rock.APR.Ar_yield_noise = abs(np.random.normal(1,ar_yield_per_mu_uncy))
    rock.APR.n_yield_noise = abs(np.random.normal(1,n_yield_per_mu_uncy))
    rock.APR.stop_rate_noise = abs(np.random.normal(1,mu_stop_rate_uncy))

    # underground neutron production uncertainty
    rock.APR.P_U_Th_noise = abs(np.random.normal(1,P_n_UTh_uncy))

    # spectral shape of alpha- and muon induced neutrons
    rock.APR.n_mu_shape_noise = np.random.normal(0, mu_shape_deviation)
    rock.APR.n_alpha_shape_noise = np.random.normal(0, alpha_shape_deviation)

    # running the calculation
    rock.calcAr_prod_rate_whole_rock(solar='avg',calc_flux=1)
    # calculating sample average production rate (Ar39 atoms/gram rock / year)
    Ar39_avg_prod_rate = intgrt.trapezoid(rock.APR.Ar39_value, rock.APR.depth)/avg_depth

    if depth_flag == 0:
        write_to_file(file_path, str(Ar39_avg_prod_rate) + ', ' +  np.array2string(rock.APR.Ar39_value, separator=', ', max_line_width=800))
    else:
        write_to_file(file_path, str(Ar39_avg_prod_rate))

    resetting_rock_parameters()

# run the iteration
iteration_number = 500
for i in range(0,iteration_number):
    mc_iteration()
    i += 1    
print('Calculation finished')
rock.close_excel()

""" plt.figure(figsize=(16,10))
plt.plot(K39xsec['Energy [MeV]'], K39xsec['Cross-Section [barns]'], color='black', linewidth=2.5, label = 'K39(n,p)Ar39 x-section')
plt.plot(Ca42xsec['Energy [MeV]'], Ca42xsec['Cross-Section [barns]'], color='red', linewidth=2.5, label='Ca42(n,alpha)Ar39 x-section')
# add the noise 
i = 0
for i in range(5):
    K39xsec_var = gaussian_noise_weights(K39xsec, column_name, relative_intensity)*K39xsec['Cross-Section [barns]']
    plt.plot(K39xsec['Energy [MeV]'], K39xsec_var, color='blue', linestyle='dashed', alpha=0.5)
    Ca42xsec_var = gaussian_noise_weights(Ca42xsec, column_name, relative_intensity)*Ca42xsec['Cross-Section [barns]']
    plt.plot(Ca42xsec['Energy [MeV]'], Ca42xsec_var, color='magenta', linestyle='dashed', alpha=0.5)
    i +=1
plt.title('Reaction cross-sections for Ar39 production, with 5 iterations of added gaussian white noise, relative intensity 0.28')
plt.xlabel('Energy [MeV]')
plt.xscale('log')
plt.ylabel('Cross-section [barns]')
plt.yscale('log')
plt.legend()
plt.show() """

""" variables  
#variablen: (multiplied)

# self.APR.stop_rate_noise  
# self.APR.Ar_yield_noise = 1
# self.APR.n_yield_noise = 1
# self.APR.P_U_Th_noise = 1

# self.composition

# variablen: added
# self.APR.n_mu_shape_noise
# self.APR.n_alpha_shape_noise"""

# reset parameters to original before new setting