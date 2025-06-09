#!/usr/bin/python
#
# H3O_Grotthuss.py
# Written by: Lars Konermann, November 2021
# konerman@uwo.ca
#
# Program to simulate Grotthuss shuttling of H+, i.e. from H3O+ to adjacent H2O
#
############################################################################################################
# When using this file, please cite                                                                        #
# "Grotthuss Molecular Dynamics Simulations for Modeling Proton Hopping in Electrosprayed Water Droplets"  #
# L. Konermann and S. Kim, J. Chem. Theory Comput. 18, 3781–3794 (2022)                                    #
############################################################################################################
#
# April 2022: code is now updated for Python3, multiprocessor operations have been cleaned up (no more transfer files)
#
# calling program from Python directory:  $ python          H3O_Grotthuss.py input_file.gro <n_hop> <hop_cutoff> <out_control> <log file name>
# calling program from another directory: $ python ~/Python/H3O_Grotthuss.py input_file.gro <n_hop> <hop_cutoff> <out_control> <log file name>
# 
# Command line arguments:
# - input file is a gro file; output files will be named input_file_grott<i>.gro
# - n_hop is the number of proton hopping events allowed for each H+
# - hop_cutoff is the maximum distance (nm) between H3O-HWx and H2O_OW for which H+ transfer is allowed
# - out_control is either 'out-all' or 'out-last'; for writing all Grotthuss steps to .gro files or just the last
# - log file name: Keeps track of what happens over multiple iterations. A new file is created, or data are appended to an existing file.

import sys
import numpy as np
import os.path
from os import path
import multiprocessing as mp
import time
import shutil

#------------------------------------------------------------------------------------------
# This function writes operational statistics for each round of H+ hopping to a log file.

#global log_file

def write_log_file():

  global log_file

  # Log file is opened after the first round of H+ hopping
  if (i_hop == 0):
    path.exists(log_file_name)
    # Log file already exists; the current information is appended to that file. 
    if (os.path.isfile(log_file_name) == True):
      log_file  = open(log_file_name,'a+')

    # Log file does not yet exist; a new file is generated and the current information is written to that file 
    if (os.path.isfile(log_file_name) == False):
      log_file  = open(log_file_name,'w+')
      output_line = '  H3O_Grotthuss.py      n_hop = ' + format(n_hop, '3') + '   hop_cutoff = ' + format(hop_cutoff, '5.3f') + ' nm'  + '\n'
      log_file.write(output_line)
  
  # Converting list to string and cleaning it up
  string_H3O_hop_success = str(H3O_hop_success).replace("'", "").replace("[", "").replace("]", "").replace(",", "")
  output_line = gro_name_ini + '  '  + string_H3O_hop_success + '\n'
  log_file.write(output_line)

  # Log file is closed after the last round of H+ hopping  
  if ( i_hop == (n_hop-1) ):
    output_line = '\n'
    log_file.write(output_line)
    log_file.close()
    print("  H+ hopping summary written to: ",log_file_name)
    print("  ")

# End of function write_log_file():
#------------------------------------------------------------------------------------------
# This function assigns charges to all atoms; it provides the groundwork for function electrostatic_H2O_selection().
def  assign_charges():

  global atom_x, atom_y,atom_z,atom_q

  charge_tot = 0.

  atom_x      = []   # Numbering of these four float arrays is such that it coincides with atom numbering in gro file
  atom_y      = []   # (i_atom = 1 corresponds to atom 1)
  atom_z      = []
  atom_q      = []

  i_atom      = 1
 
  atom_x.append(0.)  # Filling element no. 0 with dummy information that is never used,
  atom_y.append(0.)  # such that the first actual atom has index variable i_atom = 1
  atom_z.append(0.)
  atom_q.append(0.)

  for i_line in range(2,n_atom+2):
    help_line = gro_curr_content[i_line]
    atom_x.append(float(help_line[20:28]))
    atom_y.append(float(help_line[28:36]))
    atom_z.append(float(help_line[36:44]))

    # Assigning charges for each atom (will have to do this via input file starting at some future program version)
    recognition_tag =     help_line[5:15]
    molecule_number = int(help_line[0:5])
    if (recognition_tag == 'SOL     OW'):    # Assigning charges for TIP4P/2005 water
      charge = 0.
    if (recognition_tag == 'SOL    HW1'):
      charge = 0.5564
    if (recognition_tag == 'SOL    HW2'):
      charge = 0.5564
    if (recognition_tag == 'SOL     MW'):
      charge = -1.1128

    if (recognition_tag == 'H3O     OW'):   # Assigning charges for TIP4P H3O+
      charge = 0.                           
    if (recognition_tag == 'H3O    HW1'):
      charge = 0.416
    if (recognition_tag == 'H3O    HW2'):
      charge = 0.416
    if (recognition_tag == 'H3O    HW3'):
      charge = 0.416
    if (recognition_tag == 'H3O     MW'):
      charge = -0.248   

    if (recognition_tag == 'NA      NA'):   # Assigning charges for other types of atoms
      charge = 1.
    if (recognition_tag == 'CL      CL'):
      charge = -1.

    if (recognition_tag == 'HCL    CL1'):   # Assigning charges for other types of atoms
      charge = -0.182
    if (recognition_tag == 'HCL   HCL1'):
      charge = 0.182

    atom_q.append(charge)

    i_atom = i_atom + 1


    if (i_atom == 99999):
      sys.exit('   ERROR: i_atom cannot exceed 99999') 

    if (
        (recognition_tag != 'SOL     OW') and    #checking for unassignable atom types
        (recognition_tag != 'SOL    HW1') and
        (recognition_tag != 'SOL    HW2') and
        (recognition_tag != 'SOL     MW') and
        (recognition_tag != 'H3O     OW') and
        (recognition_tag != 'H3O    HW1') and
        (recognition_tag != 'H3O    HW2') and
        (recognition_tag != 'H3O    HW3') and
        (recognition_tag != 'H3O     MW') and
        (recognition_tag != 'NA      NA') and
        (recognition_tag != 'CL      CL') and
        (recognition_tag != 'HCL    CL1') and
        (recognition_tag != 'HCL   HCL1')
        ):
      print("  "),recognition_tag
      sys.exit('   ERROR:this  atom type is not defined')

    charge_tot = charge_tot + charge

  if (i_hop == 0):
    print("  ")
    print("  total system charge = ","{:.2f}".format(charge_tot))
    print("  ")


 
# End of function assign_charges():
#------------------------------------------------------------------------------------------
# This function tests the current H3O with its three closests neighbors (H-bonded to the three H3O hydrogens).
# to determine the most likely H+ acceptor.
# This is done by temporarily removing the H3O and the three H2O from the droplet, i.e. all their atomic charges are set to zero.
# We then place a dummy charge of +1e into the former MW positions of the three H2O, to determine which of the three charge 
# positions favorable (or least unfavorable); it is the one with the smallest electrostatic energy.
# We also calculate the electrostatic energy of the system when the dummy charge is located in the original H3O_MW position.
# The corresponding energy difference is delta_energy = energy_i - energy_orig.
# The chosen water molecule has the lowest delta_energy (least positive or most negative delta_energy).
# 
# Another way of describing what this function does is calculate the electrostatic potential at the four sites:
# (this avoids the "test charge" concept.
# Note that      electrostatic potential = SUM(qi/rij) vs.   electrostatic energy = SUM(q_test * qi/rij)
# In MD units, el. potential and el. energy are numerically the same, because q_test = 1.


def  electrostatic_H2O_selection(i_H3O,n_atom,atom_x,atom_y,atom_z,atom_q,gro_curr_content,
                              water_number_1,water_number_2,water_number_3,H3O,H3O_hop_success_prev):

  # FIRST STEP: using charges from assign_charges(), but:
  # For calculations, we use this list of temporary charge values, because some charges have to be temporarily set to zero
  atom_q_temp = []
  atom_q_temp.append(0.)
  for i in range(n_atom):
    i_atom = i + 1
    atom_q_temp.append(atom_q[i_atom])

  i_atom      = 1
  for i_line in range(2,n_atom+2):
    help_line = gro_curr_content[i_line]
    recognition_tag =     help_line[5:15]
    molecule_number = int(help_line[0:5])
    if (recognition_tag == 'SOL     MW'):
      # Store MW coordinates for the three possible proton acceptor H2O molecules
      if (molecule_number == water_number_1[i_H3O]):
        H2O_1_charge_center = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])]) 
      if (molecule_number == water_number_2[i_H3O]):
        H2O_2_charge_center = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])
      if (molecule_number == water_number_3[i_H3O]):
        H2O_3_charge_center = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])]) 

    if (recognition_tag == 'H3O     MW'):
      if (molecule_number == H3O[i_H3O]):
      # Store MW coordinates for the proton donor H3O
        H3O_charge_center = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])
       
    if ((molecule_number == water_number_1[i_H3O]) or
        (molecule_number == water_number_2[i_H3O]) or
        (molecule_number == water_number_3[i_H3O]) or 
        (molecule_number == H3O[i_H3O]) ):           # We temporarily assign a zero charge to all atoms of the current H3O
      atom_q_temp[i_atom] = 0.                       # and all atoms of the three possible H2O proton acceptor molecules.
    
    i_atom = i_atom + 1
  

  # SECOND STEP: electrostatic calculations

  energy_1    = 0. #electrostatic energy of the system when H+ is moved to H2O_1 position
  energy_2    = 0. #electrostatic energy of the system when H+ is moved to H2O_2 position
  energy_3    = 0. #electrostatic energy of the system when H+ is moved to H2O_3 position
  energy_orig = 0. #electrostatic energy of the system when H+ is left in its original position

#  if(H3O_hop_success_prev[i_H3O] == 'y'):  #this is now tested in the main program

  for i_atom in range(1,n_atom+1):
    distance_1 = np.sqrt( (H2O_1_charge_center[0]-atom_x[i_atom])**2 +
                          (H2O_1_charge_center[1]-atom_y[i_atom])**2 +
                          (H2O_1_charge_center[2]-atom_z[i_atom])**2 )
    if (distance_1 != 0):
      energy_1 = energy_1 + atom_q_temp[i_atom] / distance_1  # energy = SUM(q_test*qi / distance), with q_test = +1e

    distance_2 = np.sqrt( (H2O_2_charge_center[0]-atom_x[i_atom])**2 +
                          (H2O_2_charge_center[1]-atom_y[i_atom])**2 +
                          (H2O_2_charge_center[2]-atom_z[i_atom])**2 )
    if (distance_2 != 0):
      energy_2 = energy_2 + atom_q_temp[i_atom] / distance_2  # energy = SUM(q_test*qi / distance), with q_test = +1e

    distance_3 = np.sqrt( (H2O_3_charge_center[0]-atom_x[i_atom])**2 +
                          (H2O_3_charge_center[1]-atom_y[i_atom])**2 +
                          (H2O_3_charge_center[2]-atom_z[i_atom])**2 )
    if (distance_3 != 0):
      energy_3 = energy_3 + atom_q_temp[i_atom] / distance_3  # energy = SUM(q_test*qi / distance), with q_test = +1e

##   To speed up calculations, deactivate distance_orig calculations. This will not affect the electrostatic selection of H2O acceptor.
############ START DEACTIVATE (OPTIONAL)
#    distance_orig = np.sqrt( (H3O_charge_center[0]-atom_x[i_atom])**2 +
#                             (H3O_charge_center[1]-atom_y[i_atom])**2 +
#                             (H3O_charge_center[2]-atom_z[i_atom])**2 )
#    if (distance_orig != 0):
#      energy_orig = energy_orig + atom_q_temp[i_atom] / distance_orig  # energy = SUM(q_test*qi / distance), with q_test = +1e
############ END DEACTIVATE (OPTIONAL)

    #--- end of 'for i_atom in range(1,n_atom+1):'

  energy_1    = energy_1    * 138.935         # The factor 138.935 converts the result to kJ/mol
  energy_2    = energy_2    * 138.935
  energy_3    = energy_3    * 138.935
  energy_orig = energy_orig * 138.935


  # THIRD STEP: assigning acceptor H2O that has the most favorable energy upon protonation

  if(energy_1 <= energy_2) and (energy_1 <= energy_3):
    H2O_acceptor_number = water_number_1[i_H3O]
    delta_energy = energy_1 - energy_orig

  if(energy_2 <= energy_1) and (energy_2 <= energy_3):
    H2O_acceptor_number = water_number_2[i_H3O]
    delta_energy = energy_2 - energy_orig

  if(energy_3 <= energy_1) and (energy_3 <= energy_2):
    H2O_acceptor_number = water_number_3[i_H3O]   
    delta_energy = energy_3 - energy_orig


  print("  H3O+",i_H3O,H3O[i_H3O]," H2Os","{:5}".format(water_number_1[i_H3O]),"{:5}".format(water_number_2[i_H3O]),"{:5}".format(water_number_3[i_H3O]),
        "  V:","{:6.1f}".format(energy_orig),"{:6.1f}".format(energy_1),"{:6.1f}".format(energy_2),"{:6.1f}".format(energy_3),"kJ/mol - selecting H2O",
           "{:5}".format(H2O_acceptor_number),"delta_V","{:6.1f}".format(delta_energy),"kJ/mol") 

  return H2O_acceptor_number

# End of function electrostatic_H2O_selection():
#------------------------------------------------------------------------------------------
# This function writes data to output gro file after proton hopping

def  write_output():

  for i_H3O in range(n_H3O):
    for i_line in range(2,n_atom+2):

      # First write H2Os that are former H3Os to gro_out_content ...
      # If hopping has taken place, all H2O atoms are given the same velocity as the former H3O OW atom

      recognition_line = gro_curr_content[i_line]
      molecule_number = int(recognition_line[0:5])
      recognition_tag = recognition_line[5:15]
      if ((molecule_number == H2O[i_H3O]) and (recognition_tag == 'SOL     OW')):
        help_line = gro_curr_content[i_line]
        velocity_info = help_line[44:68]
        if (H3O_hop_success[i_H3O] == 'y'):
          velocity_info = H3O_OW_velocity_info[i_H3O]
        output_line = (help_line[0:20] + 
        format(H2O_OW_new_x[i_H3O], '8.3f') + 
        format(H2O_OW_new_y[i_H3O], '8.3f') + 
        format(H2O_OW_new_z[i_H3O], '8.3f') +  velocity_info)
        if (H3O_hop_success[i_H3O] == 'y'):
          output_line = output_line + (' former H3O ' + format(H3O[i_H3O], '6') + ', i_hop = ' + format(i_hop, '3'))
        output_line = output_line + '\n'
        gro_out_content[i_line] = output_line
      if ((molecule_number == H2O[i_H3O]) and (recognition_tag == 'SOL    HW1')):
        help_line = gro_curr_content[i_line]
        velocity_info = help_line[44:68]
        if (H3O_hop_success[i_H3O] == 'y'):
          velocity_info = H3O_OW_velocity_info[i_H3O] 
        output_line = (help_line[0:20] + 
        format(H2O_HW1_new_x[i_H3O], '8.3f') + 
        format(H2O_HW1_new_y[i_H3O], '8.3f') + 
        format(H2O_HW1_new_z[i_H3O], '8.3f') +  velocity_info + '\n')
        gro_out_content[i_line] = output_line
      if ((molecule_number == H2O[i_H3O]) and (recognition_tag == 'SOL    HW2')):
        help_line = gro_curr_content[i_line]
        velocity_info = help_line[44:68]
        if (H3O_hop_success[i_H3O] == 'y'):
          velocity_info = H3O_OW_velocity_info[i_H3O] 
        output_line = (help_line[0:20] + 
        format(H2O_HW2_new_x[i_H3O], '8.3f') + 
        format(H2O_HW2_new_y[i_H3O], '8.3f') + 
        format(H2O_HW2_new_z[i_H3O], '8.3f') +  velocity_info + '\n')
        gro_out_content[i_line] = output_line
      if ((molecule_number == H2O[i_H3O]) and (recognition_tag == 'SOL     MW')):
        help_line = gro_curr_content[i_line]
        velocity_info = help_line[44:68]
        if (H3O_hop_success[i_H3O] == 'y'):
          velocity_info = H3O_OW_velocity_info[i_H3O] 
        output_line = (help_line[0:20] + 
        format(H2O_MW_new_x[i_H3O], '8.3f') + 
        format(H2O_MW_new_y[i_H3O], '8.3f') + 
        format(H2O_MW_new_z[i_H3O], '8.3f') +  velocity_info + '\n')
        gro_out_content[i_line] = output_line


      # ... then write H3Os that are former H2Os to gro_out_content.
      # If hopping has taken place, the all H3O atoms are given the same velocity as the former H2O OW atom.

      recognition_line = gro_curr_content[i_line]
      molecule_number = int(recognition_line[0:5])
      recognition_tag = recognition_line[5:15]
      if ((molecule_number == H3O[i_H3O]) and (recognition_tag == 'H3O     OW')):
        help_line = gro_curr_content[i_line]
        velocity_info = help_line[44:68]
        if (H3O_hop_success[i_H3O] == 'y'):
          velocity_info = H2O_OW_velocity_info[i_H3O]
        output_line = (help_line[0:20] + 
        format(H3O_OW_new_x[i_H3O], '8.3f') + 
        format(H3O_OW_new_y[i_H3O], '8.3f') + 
        format(H3O_OW_new_z[i_H3O], '8.3f') +  velocity_info)
        if (H3O_hop_success[i_H3O] == 'y'):
          output_line = output_line + (' former H2O ' + format(H2O[i_H3O], '6') + ', i_hop = ' + format(i_hop, '3'))
        output_line = output_line + '\n'
        gro_out_content[i_line] = output_line
      if ((molecule_number == H3O[i_H3O]) and (recognition_tag == 'H3O    HW1')):
        help_line = gro_curr_content[i_line]
        velocity_info = help_line[44:68]
        if (H3O_hop_success[i_H3O] == 'y'):
          velocity_info = H2O_OW_velocity_info[i_H3O]
        output_line = (help_line[0:20] + 
        format(H3O_HW1_new_x[i_H3O], '8.3f') + 
        format(H3O_HW1_new_y[i_H3O], '8.3f') + 
        format(H3O_HW1_new_z[i_H3O], '8.3f') +  velocity_info + '\n')
        gro_out_content[i_line] = output_line
      if ((molecule_number == H3O[i_H3O]) and (recognition_tag == 'H3O    HW2')):
        help_line = gro_curr_content[i_line]
        velocity_info = help_line[44:68]
        if (H3O_hop_success[i_H3O] == 'y'):
          velocity_info = H2O_OW_velocity_info[i_H3O]
        output_line = (help_line[0:20] + 
        format(H3O_HW2_new_x[i_H3O], '8.3f') + 
        format(H3O_HW2_new_y[i_H3O], '8.3f') + 
        format(H3O_HW2_new_z[i_H3O], '8.3f') +  velocity_info + '\n')
        gro_out_content[i_line] = output_line
      if ((molecule_number == H3O[i_H3O]) and (recognition_tag == 'H3O    HW3')):
        help_line = gro_curr_content[i_line]
        velocity_info = help_line[44:68]
        if (H3O_hop_success[i_H3O] == 'y'):
          velocity_info = H2O_OW_velocity_info[i_H3O]
        output_line = (help_line[0:20] + 
        format(H3O_HW3_new_x[i_H3O], '8.3f') + 
        format(H3O_HW3_new_y[i_H3O], '8.3f') + 
        format(H3O_HW3_new_z[i_H3O], '8.3f') +  velocity_info + '\n')
        gro_out_content[i_line] = output_line
      if ((molecule_number == H3O[i_H3O]) and (recognition_tag == 'H3O     MW')):
        help_line = gro_curr_content[i_line]
        velocity_info = help_line[44:68]
        if (H3O_hop_success[i_H3O] == 'y'):
          velocity_info = H2O_OW_velocity_info[i_H3O]
        output_line = (help_line[0:20] + 
        format(H3O_MW_new_x[i_H3O], '8.3f') + 
        format(H3O_MW_new_y[i_H3O], '8.3f') + 
        format(H3O_MW_new_z[i_H3O], '8.3f') +  velocity_info + '\n')
        gro_out_content[i_line] = output_line

  # Writing output files, either for all n_hop steps or only for the last step.
  # Note that (n_hop-1 is the last frame because Python starts counting at zero.
  if (   (out_control == 'out-all') or ( (out_control == 'out-last') and (i_hop == (n_hop-1)) )   ):
    gro_out     = open(gro_out_name[i_hop],'w+')  #open output file
    gro_out.writelines(gro_out_content)           #write modified data to output file
    gro_out.close()                               #close output file
    print(" ")
    print("  output file generated:         ",gro_out_name[i_hop])
    print("  ")

  # End of function write_output():
#------------------------------------------------------------------------------------------
# Start of function H3O_to_H2O.
# This function transfers a proton from H3O to H2O, and leaves both donor and acceptor in proper orientations.
# More specifically, the H3O coordinates are turned into H2O coordinates, and
#                    the H2O coordinates are turned into H3O coordinates.
#                    The OW atoms of both H2O and H3O trade places.
# H+ hopping take place only if hop_flag = 'y'. Distances that are too large (controlled by hop_cutoff) will prevent hopping.

def  H3O_to_H2O():

  global H3O_OW,H2O_OW

  print("  ")
  print("  Checking if hopping distances are OK:")

  for i_H3O in range(n_H3O):

    hop_flag = 'y'         #proton transfer only takes place if this flag is 'y', we set it to 'y' by default
    help_line1 = gro_info_H3O_OW[i_H3O]
    help_line2 = gro_info_H2O_OW[i_H3O] 

    #assign variables for calculations
    help_line   = gro_info_H3O_OW[i_H3O]
    H3O_OW      = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])
    H3O_OW_velocity_info.append(help_line[44:68])

    help_line   = gro_info_H3O_HW1[i_H3O]
    H3O_HW1_ini = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])

    help_line   = gro_info_H3O_HW2[i_H3O]
    H3O_HW2_ini = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])

    help_line   = gro_info_H3O_HW3[i_H3O]
    H3O_HW3_ini = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])

    help_line   = gro_info_H3O_MW[i_H3O]
    H3O_MW      = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])

    help_line   = gro_info_H2O_OW[i_H3O]
    H2O_OW      = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])
    H2O_OW_velocity_info.append(help_line[44:68])

    help_line   = gro_info_H2O_HW1[i_H3O]
    H2O_HW1     = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])

    help_line   = gro_info_H2O_HW2[i_H3O]
    H2O_HW2     = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])

    help_line   = gro_info_H2O_MW[i_H3O]
    H2O_MW      = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])


    # STEP 0: Rotate H3O so that H3O_HW3 is the one that is closest to H2O_OW.
    #         For doing this, we reassign the atom positions.
    
    dist_HWi_OW_min = 100000000.

    dist_HWi_OW = np.linalg.norm(H3O_HW1_ini-H2O_OW)  # testing HW1
    if (dist_HWi_OW < dist_HWi_OW_min):
       dist_HWi_OW_min =  dist_HWi_OW
       closest_HWi = 1

    dist_HWi_OW = np.linalg.norm(H3O_HW2_ini-H2O_OW)  # testing HW2
    if (dist_HWi_OW < dist_HWi_OW_min):
       dist_HWi_OW_min =  dist_HWi_OW
       closest_HWi = 2

    dist_HWi_OW = np.linalg.norm(H3O_HW3_ini-H2O_OW)  # testing HW3
    if (dist_HWi_OW < dist_HWi_OW_min):
       dist_HWi_OW_min =  dist_HWi_OW
       closest_HWi = 3

    if (closest_HWi == 1):                 #these resassignments retain the H3O chirality
      H3O_HW1     = H3O_HW2_ini
      H3O_HW2     = H3O_HW3_ini
      H3O_HW3     = H3O_HW1_ini

    if (closest_HWi == 2):
      H3O_HW1     = H3O_HW3_ini
      H3O_HW2     = H3O_HW1_ini
      H3O_HW3     = H3O_HW2_ini

    if (closest_HWi == 3):
      H3O_HW1     = H3O_HW1_ini
      H3O_HW2     = H3O_HW2_ini
      H3O_HW3     = H3O_HW3_ini


    # By default, the proton will get transferred from H3O+ to H2O (hop_flag = 'y').
    # Now, we are testing a couple of conditions that will prevent this hopping event (hop_flag = 'n')

    # Testing if the Grotthus chain has been already ruptured at a previous step. If so, the current step will not take place.
    if (H3O_hop_success_prev[i_H3O] == 'n'):  #this is needed to prevent back-and-forth hopping
      print("  H3O",i_H3O," hopping unsuccessful because previous attempts failed")
      hop_flag = 'n'

    # Testing the distance cutoff. If distance is too large, proton transfer does not take place.
    if ( (dist_HWi_OW_min > hop_cutoff) and (H3O_hop_success_prev[i_H3O] == 'y') ):
      print("  H3O",i_H3O," hopping unsuccessful with distance : ","{:.3f}".format(dist_HWi_OW_min)," nm")
      hop_flag = 'n'
    if ( (dist_HWi_OW_min <= hop_cutoff) and (H3O_hop_success_prev[i_H3O] == 'y') ):
      print("  H3O",i_H3O," hopping   successful with distance : ","{:.3f}".format(dist_HWi_OW_min)," nm")
      #hop_flag remains as 'y', as defined previously

    # Testing if the current H3O is too close to another H3O (this can lead to Gromacs crashes).
    # We exclude O-O distances of less than 0.4 nm
      for j_H3O in range( (i_H3O+1), n_H3O):
        help_line3 = gro_info_H2O_OW[j_H3O]
        H3O_OWj    = np.array([float(help_line3[20:28]),float(help_line3[28:36]),float(help_line3[36:44])])
        OO_distance = np.sqrt((H3O_OW[0]-H3O_OWj[0])**2 + (H3O_OW[1]-H3O_OWj[1])**2 + (H3O_OW[2]-H3O_OWj[2])**2)
        if ( (OO_distance < 0.4) and (H3O_hop_success_prev[i_H3O] == 'y') ):
          print("  H3O",i_H3O," H+ hopping  unsuccessful because it would get too close to another H3O")
          hop_flag = 'n'

    H3O_hop_success.append(hop_flag)  #this list keeps track of which events were successful

    if (hop_flag == 'n'):    # Proton hopping does NOT take place if hop_flag = 'n'

      H3O_OW_new_x.append(H3O_OW[0])
      H3O_OW_new_y.append(H3O_OW[1])
      H3O_OW_new_z.append(H3O_OW[2])

      H3O_HW1_new_x.append(H3O_HW1_ini[0])
      H3O_HW1_new_y.append(H3O_HW1_ini[1])
      H3O_HW1_new_z.append(H3O_HW1_ini[2])

      H3O_HW2_new_x.append(H3O_HW2_ini[0])
      H3O_HW2_new_y.append(H3O_HW2_ini[1])
      H3O_HW2_new_z.append(H3O_HW2_ini[2])

      H3O_HW3_new_x.append(H3O_HW3_ini[0])
      H3O_HW3_new_y.append(H3O_HW3_ini[1])
      H3O_HW3_new_z.append(H3O_HW3_ini[2])

      H3O_MW_new_x.append(H3O_MW[0])
      H3O_MW_new_y.append(H3O_MW[1])
      H3O_MW_new_z.append(H3O_MW[2])

      H2O_OW_new_x.append(H2O_OW[0])
      H2O_OW_new_y.append(H2O_OW[1])
      H2O_OW_new_z.append(H2O_OW[2])

      H2O_HW1_new_x.append(H2O_HW1[0])
      H2O_HW1_new_y.append(H2O_HW1[1])
      H2O_HW1_new_z.append(H2O_HW1[2])

      H2O_HW2_new_x.append(H2O_HW2[0])
      H2O_HW2_new_y.append(H2O_HW2[1])
      H2O_HW2_new_z.append(H2O_HW2[2])

      H2O_MW_new_x.append(H2O_MW[0])
      H2O_MW_new_y.append(H2O_MW[1])
      H2O_MW_new_z.append(H2O_MW[2])

    # end of 'if (hop_flag == 'n'): '


    if (hop_flag == 'y'):    # Proton hopping DOES take place if hop_flag = 'y'

      # STEP 1: translate the H3O and the H2O such that their OW atoms are at the origin (at 0/0/0).
      # HW3 is tranlated such that it will become part of H2O, thereby making the proton "hop";
      # after this hop HW3 is located on the H3O_OW --- H2O_OW connecting line (between the two oxygens).
      H3O_HW1_cent    =   H3O_HW1 - H3O_OW
      H3O_HW2_cent    =   H3O_HW2 - H3O_OW
      H3O_HW3_cent    =   0.09686 * (H3O_OW-H2O_OW) / np.linalg.norm(H3O_OW-H2O_OW)  
      H3O_MW_cent     =   H3O_MW  - H3O_OW

      H2O_HW1_cent    =   H2O_HW1 - H2O_OW
      H2O_HW2_cent    =   H2O_HW2 - H2O_OW
      H2O_MW_cent     =   H2O_MW  - H2O_OW


      # STEP 2: Generating H3O with proper geometry, using H2O HW1 and HW2 

      # OW--HW bond length in TIP4P H3O = 0.09686 nm  [Temeslo, 2013]
      # OW--MW bond distance            = 0.015   nm
      H3O_HW1_temp    =  0.09686 * H2O_HW1_cent / np.linalg.norm(H2O_HW1_cent)
      H3O_HW2_temp    =  0.09686 * H2O_HW2_cent / np.linalg.norm(H2O_HW2_cent)
      H3O_HW3_temp    =  0.09686 * H3O_HW3_cent / np.linalg.norm(H3O_HW3_cent)
      H3O_HW_temp_sum =  H3O_HW1_temp + H3O_HW2_temp + H3O_HW3_temp
      H3O_MW_temp     =  0.01500 * H3O_HW_temp_sum / np.linalg.norm(H3O_HW_temp_sum) #H3O_MW_cent
      # All bond lenghts are now OK, and MW sits within the OW/HW1/HW2/HW3 tripod (ruling out an umbrella flip).
    
      # However, the bond angles are still off. We now tweak the geometry, leaving only H3O_HW3_temp and OW = 0/0/0 unchanged
 
      #First tweaking MW:
      #The HW3--OW--MW angle should be 74.4 deg = 1.299 rad  (average value measured by Pymol in H3Oplus_tip4p.gro)
      axis              =  np.cross(H3O_HW3_temp,H3O_MW_temp)
      axis_norm         =  axis / np.linalg.norm(axis)               #normalized vector orthogonal to HW3temp--OW--MWtemp

      # Determine current HW3--OW--MW angle:
      cos_HOMangle_temp = np.dot(H3O_HW3_temp,H3O_MW_temp) / (np.linalg.norm(H3O_HW3_temp)*np.linalg.norm(H3O_MW_temp))
      HOMangle_temp     = np.arccos(cos_HOMangle_temp)
      diff_angle        = 1.299 - HOMangle_temp

      # Rotate MW into proper position, using rotation matrix from
      # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
      ux        =   axis_norm[0]
      uy        =   axis_norm[1]
      uz        =   axis_norm[2]
      cos       =   np.cos(diff_angle)
      sin       =   np.sin(diff_angle)

      ROT       = np.array([ [  cos+ux*ux*(1-cos)    , ux*uy*(1-cos)-uz*sin , ux*uz*(1-cos)+uy*sin ],
                             [  uy*ux*(1-cos)+uz*sin , cos+uy*uy*(1-cos)    , uy*uz*(1-cos)-ux*sin ],
                             [  uz*ux*(1-cos)-uy*sin , uz*uy*(1-cos)+ux*sin , cos+uz*uz*(1-cos)    ] ])  

      H3O_MW_temp1 = np.matmul(ROT,H3O_MW_temp)

      # Rotate HW2 and HW1 into proper position (rotation of H3O_HW3 by -120 deg = -2.094 rad  around OW--H3O_MW_temp1 axis)

      axis_norm         =  H3O_MW_temp1 / np.linalg.norm(H3O_MW_temp1)  #normalized vector in H3O_MW direction

      ux        =   axis_norm[0]
      uy        =   axis_norm[1]
      uz        =   axis_norm[2]
      cos       =   np.cos(-2.094)
      sin       =   np.sin(-2.094)

      ROT       = np.array([ [  cos+ux*ux*(1-cos)    , ux*uy*(1-cos)-uz*sin , ux*uz*(1-cos)+uy*sin ],
                             [  uy*ux*(1-cos)+uz*sin , cos+uy*uy*(1-cos)    , uy*uz*(1-cos)-ux*sin ],
                             [  uz*ux*(1-cos)-uy*sin , uz*uy*(1-cos)+ux*sin , cos+uz*uz*(1-cos)    ] ])  

      H3O_HW2_temp1 = np.matmul(ROT,H3O_HW3_temp)
      H3O_HW1_temp1 = np.matmul(ROT,H3O_HW2_temp1)   #this retains the H3O chirality
 

      # STEP 3: Generating H2O with proper geometry, using H3O HW1 and HW2 

      # O--H bond length in TIP4P/2005 = 0.09572 nm  [Abascal, 2005]
      H2O_HW1_temp         =  0.09572 * H3O_HW1_cent / np.linalg.norm(H3O_HW1_cent)
      H2O_HW2_temp         =  0.09572 * H3O_HW2_cent / np.linalg.norm(H3O_HW2_cent)
      axis                    =  np.cross(H2O_HW1_temp,H2O_HW2_temp)
      axis_norm               =  axis / np.linalg.norm(axis)               #normalized vector orthogonal to H--O--H
      #calculating MW position from HW1--OW--HW2 bisector; OW--MW distance is 0.01546 nm [Abascal, 2005]
      H2O_MW_from_scratch  =  0.01546*(H3O_HW1_cent + H3O_HW2_cent) / np.linalg.norm(H3O_HW1_cent + H3O_HW2_cent)

      # All four H2O atoms are now in a plane, and all bond lengths are OK.

      # The H--O--H angle is still off. We now set it to 104.52 deg = 1.8242 rad [Abascal, 2005]
      # Determine current H--O--H angle:
      cos_HOHangle_temp = np.dot(H2O_HW1_temp,H2O_HW2_temp) / (np.linalg.norm(H2O_HW1_temp)*np.linalg.norm(H2O_HW2_temp))
      HOHangle_temp     = np.arccos(cos_HOHangle_temp)
      diff_angle        = HOHangle_temp - 1.8242

      # Rotate HW1 into proper position
      ux        =   axis_norm[0]
      uy        =   axis_norm[1]
      uz        =   axis_norm[2]
      cos       =   np.cos(diff_angle/2.)
      sin       =   np.sin(diff_angle/2.)

      ROT       = np.array([ [  cos+ux*ux*(1-cos)    , ux*uy*(1-cos)-uz*sin , ux*uz*(1-cos)+uy*sin ],
                             [  uy*ux*(1-cos)+uz*sin , cos+uy*uy*(1-cos)    , uy*uz*(1-cos)-ux*sin ],
                             [  uz*ux*(1-cos)-uy*sin , uz*uy*(1-cos)+ux*sin , cos+uz*uz*(1-cos)    ] ])  

      H2O_HW1_temp1 = np.matmul(ROT,H2O_HW1_temp)

      # Rotate HW2 into proper position
      cos       =   np.cos(-1. * diff_angle/2.)
      sin       =   np.sin(-1. * diff_angle/2.)

      ROT       = np.array([ [  cos+ux*ux*(1-cos)    , ux*uy*(1-cos)-uz*sin , ux*uz*(1-cos)+uy*sin ],
                             [  uy*ux*(1-cos)+uz*sin , cos+uy*uy*(1-cos)    , uy*uz*(1-cos)-ux*sin ],
                             [  uz*ux*(1-cos)-uy*sin , uz*uy*(1-cos)+ux*sin , cos+uz*uz*(1-cos)    ] ])  

      H2O_HW2_temp1 = np.matmul(ROT,H2O_HW2_temp)



  

      # STEP 4: Now that rotations are finished, each H3O is translocated to the H2O position, and 
      #                                          each H2O is translocated to the H3O position 
      # After this operation, the properly rotatated H3O and H2O have traded places. 


      H3O_OW_new_x.append(                          H2O_OW[0])
      H3O_OW_new_y.append(                          H2O_OW[1])
      H3O_OW_new_z.append(                          H2O_OW[2])

      H3O_HW1_new_x.append(H3O_HW1_temp1[0]       + H2O_OW[0])
      H3O_HW1_new_y.append(H3O_HW1_temp1[1]       + H2O_OW[1])
      H3O_HW1_new_z.append(H3O_HW1_temp1[2]       + H2O_OW[2])

      H3O_HW2_new_x.append(H3O_HW2_temp1[0]       + H2O_OW[0])
      H3O_HW2_new_y.append(H3O_HW2_temp1[1]       + H2O_OW[1])
      H3O_HW2_new_z.append(H3O_HW2_temp1[2]       + H2O_OW[2])

      H3O_HW3_new_x.append(H3O_HW3_temp[0]       + H2O_OW[0])
      H3O_HW3_new_y.append(H3O_HW3_temp[1]       + H2O_OW[1])
      H3O_HW3_new_z.append(H3O_HW3_temp[2]       + H2O_OW[2])

      H3O_MW_new_x.append(H3O_MW_temp1[0]        + H2O_OW[0])
      H3O_MW_new_y.append(H3O_MW_temp1[1]        + H2O_OW[1])
      H3O_MW_new_z.append(H3O_MW_temp1[2]        + H2O_OW[2])

      H2O_OW_new_x.append(                         H3O_OW[0])
      H2O_OW_new_y.append(                         H3O_OW[1])
      H2O_OW_new_z.append(                         H3O_OW[2])

      H2O_HW1_new_x.append(H2O_HW1_temp1[0]      + H3O_OW[0])
      H2O_HW1_new_y.append(H2O_HW1_temp1[1]      + H3O_OW[1])
      H2O_HW1_new_z.append(H2O_HW1_temp1[2]      + H3O_OW[2])

      H2O_HW2_new_x.append(H2O_HW2_temp1[0]      + H3O_OW[0])
      H2O_HW2_new_y.append(H2O_HW2_temp1[1]      + H3O_OW[1])
      H2O_HW2_new_z.append(H2O_HW2_temp1[2]      + H3O_OW[2])

      H2O_MW_new_x.append(H2O_MW_from_scratch[0] + H3O_OW[0])
      H2O_MW_new_y.append(H2O_MW_from_scratch[1] + H3O_OW[1])
      H2O_MW_new_z.append(H2O_MW_from_scratch[2] + H3O_OW[2])

    # end of 'if (hop_flag == 'y'):'

  # End of function orient_H3O_to_H2O
#------------------------------------------------------------------------------------------
# start of main program

global n_H3O,i_H3O
global gro_info_H3O_OW,gro_info_H3O_HW1,gro_info_H3O_HW2,gro_info_H3O_HW3,gro_info_H3O_MW
global H3O_OW_new_x, H3O_OW_new_y, H3O_OW_new_z
global H3O_HW1_new_x,H3O_HW1_new_y,H3O_HW1_new_z
global H3O_HW2_new_x,H3O_HW2_new_y,H3O_HW2_new_z
global H3O_HW3_new_x,H3O_HW3_new_y,H3O_HW3_new_z
global H3O_MW_new_x, H3O_MW_new_y, H3O_MW_new_z
global H2O_OW_new_x, H2O_OW_new_y, H2O_OW_new_z
global H2O_HW1_new_x,H2O_HW1_new_y,H2O_HW1_new_z
global H2O_HW2_new_x,H2O_HW2_new_y,H2O_HW2_new_z
global H2O_MW_new_x, H2O_MW_new_y, H2O_MW_new_z
global H3O_hop_success,H3O_hop_success_prev
global H3O #,H2O
global i_hop#,i_H3O
global H3O_OW_velocity_info, H2O_OW_velocity_info
global gro_curr_content
global water_number_1,water_number_2,water_number_3
global hop_cutoff
global out_control
global log_file_name
global H2O,gro_info_H2O_OW,gro_info_H2O_HW1,gro_info_H2O_HW2,gro_info_H2O_MW



# For measuring wall clock time
start_time = time.time()

gro_out_name  = [] # output file names

print("                         ")
print("  #######################")
print("  #                     #")
print("  #  H3O_Grotthuss.py   #")
print("  #                     #")
print("  #######################")
print("                         ")

gro_name_out  = []

n_cpu = mp.cpu_count()
n_cpu_used = int(mp.cpu_count()/2)
print("  CPUs used     : ", n_cpu_used)

# assign file names etc. from command line info
gro_name_ini  = sys.argv[1]        # input  gro file before protons have moved
n_hop         = int(sys.argv[2])   # number of times that protons are moved
hop_cutoff    = float(sys.argv[3]) # maximum distance (nm) for which H+ hopping is allowed
out_control   = sys.argv[4]
log_file_name = sys.argv[5]

if (n_hop > 5):
  sys.exit('   ERROR: nhop > 5')

print("  input gro file (prior to H+ hopping):  ",gro_name_ini)
print("  attempted hopping events for each H+:  ",n_hop)
print("  hopping distance cutoff (nm)        :  ",hop_cutoff)

# read all lines of gro
gro_ini    = open(gro_name_ini,'r')      #open input file
gro_ini_content = gro_ini.readlines()    #initial system structure is stored
gro_ini.close()                          #close input file

# determine number of atoms
n_atom = int(gro_ini_content[1])
print("  number of atoms: ",n_atom)

# count H3O
n_H3O = 0
for i_line in range(2,n_atom+2):
  help_line = gro_ini_content[i_line]
  recognition_tag = help_line[5:15]
  if(recognition_tag == 'H3O     OW'):
    n_H3O = n_H3O + 1
print("  number of H3O  : ",n_H3O)

# count H2O (SOL)
n_H2O = 0
for i_line in range(2,n_atom+2):
  help_line = gro_ini_content[i_line]
  recognition_tag = help_line[5:15]
  if(recognition_tag == 'SOL     OW'):
    n_H2O = n_H2O + 1
print("  number of H2O  : ",n_H2O)


if ( (out_control != 'out-all') and (out_control != 'out-last') ):
  sys.exit('   ERROR: out_control must be set to out-all or out-last')
if (out_control == 'out-all'):
  print("  output gro files will be generated for all ",n_hop," rounds of attempted H+ hopping")
if (out_control == 'out-last'):
  print("  output gro file will be generated only for the last round of attempted H+ hopping")

for i in range(n_hop):
  helpline = sys.argv[1]
  helpline = helpline[:-4]
  gro_out_name.append(helpline + '_grott' + str(i) + '.gro' )
print('  ')

# As a precaution for situations where there are no H3O or no H2O in the input file, the following option is included
if( (n_H3O == 0) or(n_H2O == 0) ):
  print("  ")
  print("  ### WARNING: No H3O or no H2O found, Grotthuss hopping not possible.")
  shutil.copy(gro_name_ini, gro_out_name[n_hop-1])
  print("  output file ",gro_out_name[n_hop-1]," is identical to input file ",gro_name_ini)
  print("  program terminated")
  print("  ")
  sys.exit()

H2O_prev = []   #list of integers (index number of H2O proton acceptors in gro file) from previous i_hop
for i_H3O in range(n_H3O):
  H2O_prev.append(0)

H3O_hop_success_prev = [] 
for i_H3O in range(n_H3O):  # 'y' or 'n' hopping success for each H3O in the previous iteration
  H3O_hop_success_prev.append('y')    

gro_out_content = gro_ini_content

# Now performing n_hop rounds of proton hopping
for i_hop in range (n_hop):

  print("  ")
  print(" ================ i_hop =",i_hop,"=================")

  # Initialize lists; in all of these lists the number [] refers to the index variable i_H3O.
  # These lists refer to coordinates BEFORE () proton hopping
  gro_info_H3O_OW  = []  #lists of text lines, H3O info BEFORE proton hopping
  gro_info_H3O_HW1 = []
  gro_info_H3O_HW2 = []
  gro_info_H3O_HW3 = []
  gro_info_H3O_MW  = []

  gro_info_H2O_OW  = []  #lists of text lines, H2O info BEFORE proton hopping
  gro_info_H2O_HW1 = []  #these are the closest waters (proton acceptors) for the corresponding H3O
  gro_info_H2O_HW2 = []
  gro_info_H2O_MW  = []
  for i_H3O in range(n_H3O):
    gro_info_H2O_OW.append(' ')
    gro_info_H2O_HW1.append(' ')
    gro_info_H2O_HW2.append(' ')
    gro_info_H2O_MW.append(' ')


  H3O_OW_new_x  = [] #lists of float numbers (coordinates after proton hopping)
  H3O_OW_new_y  = []
  H3O_OW_new_z  = []
  H3O_HW1_new_x = []
  H3O_HW1_new_y = []
  H3O_HW1_new_z = []
  H3O_HW2_new_x = []
  H3O_HW2_new_y = []
  H3O_HW2_new_z = []
  H3O_HW3_new_x = []
  H3O_HW3_new_y = []
  H3O_HW3_new_z = []
  H3O_MW_new_x  = []
  H3O_MW_new_y  = []
  H3O_MW_new_z  = []

  H2O_OW_new_x  = [] #lists of float numbers (coordinates after proton hopping)
  H2O_OW_new_y  = []
  H2O_OW_new_z  = []
  H2O_HW1_new_x = []
  H2O_HW1_new_y = []
  H2O_HW1_new_z = []
  H2O_HW2_new_x = []
  H2O_HW2_new_y = []
  H2O_HW2_new_z = []
  H2O_MW_new_x  = []
  H2O_MW_new_y  = []
  H2O_MW_new_z  = []

  water_number_1 = [] # lists of potential H+ acceptors (next to H3O-HW1/2/3)
  water_number_2 = []
  water_number_3 = []
  for i_H3O in range(n_H3O):
    water_number_1.append(0)
    water_number_2.append(0)
    water_number_3.append(0)

  H3O_hop_success      = []  # 'y' or 'n' hopping success for each H3O in the current iteration
  H2O      = []   #list of integers (index number of H2O proton acceptors in gro file)
  H3O      = []   #list of integers (index number of H3O proton donors in gro file)
  for i_H3O in range(n_H3O):
    H2O.append(0)         

  H3O_OW_velocity_info = []
  H2O_OW_velocity_info = []

  gro_curr_content    = gro_out_content  #this is the current system structure before proton hopping
  gro_out_content     = gro_ini_content  #this will become the system structure after proton hopping

  # determine current H3O coordinates
  for i_line in range(2,n_atom+2):
    help_line = gro_curr_content[i_line]
    recognition_tag = help_line[5:15]
    if(recognition_tag == 'H3O     OW'):
      H3O_number = int(help_line[0:5])
      help_line = gro_curr_content[i_line  ].replace("\n","")
      gro_info_H3O_OW.append(help_line)
      help_line = gro_curr_content[i_line+1].replace("\n","")
      gro_info_H3O_HW1.append(help_line)
      help_line = gro_curr_content[i_line+2].replace("\n","")
      gro_info_H3O_HW2.append(help_line)
      help_line = gro_curr_content[i_line+3].replace("\n","")
      gro_info_H3O_HW3.append(help_line)
      help_line = gro_curr_content[i_line+4].replace("\n","")
      gro_info_H3O_MW.append(help_line)
      H3O.append(H3O_number)

  print("  H3O+ [H2O]3 complexes:")

  # for each H3O, find the three closest H2O oxygens (shortest H3O  HWx ... SOL  OW distances).
  # One of these three waters in this (H3O)(H2O)3 Eigen complex will become the proton acceptor for the corresponding H3O.
  for i_H3O in range(n_H3O):
    if(H3O_hop_success_prev[i_H3O] == 'y'):  #We skip H3Os that were previously unsuccessful

      closest_distance_HW1 = 1000000.
      closest_distance_HW2 = 1000000.
      closest_distance_HW3 = 1000000.

      #testing H3O_HW1
      help_line_H3O = gro_info_H3O_HW1[i_H3O]
      H3O_HW_x = float(help_line_H3O[20:28])
      H3O_HW_y = float(help_line_H3O[28:36])
      H3O_HW_z = float(help_line_H3O[36:44])

      for i_line in range(2,n_atom+2):
        help_line_H2O   = gro_curr_content[i_line]
        recognition_tag = help_line_H2O[5:15]
        if(recognition_tag == 'SOL     OW'):
          H2O_O_x = float(help_line_H2O[20:28]) 
          H2O_O_y = float(help_line_H2O[28:36])
          H2O_O_z = float(help_line_H2O[36:44])
          distance = np.sqrt((H3O_HW_x-H2O_O_x)**2 + (H3O_HW_y-H2O_O_y)**2 + (H3O_HW_z-H2O_O_z)**2)
          water_number_test = int(help_line_H2O[0:5])
          if (distance < closest_distance_HW1) and (water_number_test != H2O_prev[i_H3O]): # H2O is excluded if it was 
            closest_distance_HW1 = distance                                                # proton acceptor in the previous step.
            water_number_1[i_H3O] = water_number_test

      #testing H3O_HW2
      help_line_H3O = gro_info_H3O_HW2[i_H3O]
      H3O_HW_x = float(help_line_H3O[20:28])
      H3O_HW_y = float(help_line_H3O[28:36])
      H3O_HW_z = float(help_line_H3O[36:44])

      for i_line in range(2,n_atom+2):
        help_line_H2O   = gro_curr_content[i_line]
        recognition_tag = help_line_H2O[5:15]
        if(recognition_tag == 'SOL     OW'):
          H2O_O_x = float(help_line_H2O[20:28]) 
          H2O_O_y = float(help_line_H2O[28:36])
          H2O_O_z = float(help_line_H2O[36:44])
          distance = np.sqrt((H3O_HW_x-H2O_O_x)**2 + (H3O_HW_y-H2O_O_y)**2 + (H3O_HW_z-H2O_O_z)**2)
          water_number_test = int(help_line_H2O[0:5])
          if (distance < closest_distance_HW2) and (water_number_test != H2O_prev[i_H3O]): # H2O is excluded if it was 
            closest_distance_HW2 = distance                                                # proton acceptor in the previous step.
            water_number_2[i_H3O] = water_number_test

      #testing H3O_HW3
      help_line_H3O = gro_info_H3O_HW3[i_H3O]
      H3O_HW_x = float(help_line_H3O[20:28])
      H3O_HW_y = float(help_line_H3O[28:36])
      H3O_HW_z = float(help_line_H3O[36:44])

      for i_line in range(2,n_atom+2):
        help_line_H2O   = gro_curr_content[i_line]
        recognition_tag = help_line_H2O[5:15]
        if(recognition_tag == 'SOL     OW'):
          H2O_O_x = float(help_line_H2O[20:28]) 
          H2O_O_y = float(help_line_H2O[28:36])
          H2O_O_z = float(help_line_H2O[36:44])
          distance = np.sqrt((H3O_HW_x-H2O_O_x)**2 + (H3O_HW_y-H2O_O_y)**2 + (H3O_HW_z-H2O_O_z)**2)
          water_number_test = int(help_line_H2O[0:5])
          if (distance < closest_distance_HW3) and (water_number_test != H2O_prev[i_H3O]): # H2O is excluded if it was 
            closest_distance_HW3 = distance                                                # proton acceptor in the previous step.
            water_number_3[i_H3O] = water_number_test
      
      print("  H3O ",i_H3O,H3O[i_H3O]," H2O_1","{:5}".format(water_number_1[i_H3O]),"(","{:.3f}".format(closest_distance_HW1),"nm),",
                                      " H2O_2","{:5}".format(water_number_2[i_H3O]),"(","{:.3f}".format(closest_distance_HW2),"nm),",
                                      " H2O_3","{:5}".format(water_number_3[i_H3O]),"(","{:.3f}".format(closest_distance_HW3),"nm)")

  #---end loop "for i_H3O in range(n_H3O):"

  print("  ")

  # Prior to performing electrostatic calculations, we have to assign charges to all atoms  
  assign_charges()

  # Now that the three potential H2O proton acceptors are identified for each H3O, we run another loop
  # that electrostatically selects the most favorable proton acceptor for each i_H3O.

##### START SEQUENTIAL PROCESSING
#  #Use the next two lines for sequential processing for all protons (slow, but perhaps better for troubleshooting)  
#  for i_H3O in range(n_H3O):
#    if(H3O_hop_success_prev[i_H3O] == 'y'):  #We skip H3Os that were previously unsuccessful
#      electrostatic_H2O_selection(i_H3O,n_atom,atom_x,atom_y,atom_z,atom_q,gro_curr_content,
#                        water_number_1,water_number_2,water_number_3,H3O,H3O_hop_success_prev)
#### END SEQUENTIAL PROCESSING


#### START PARALLEL PROCESSING _NEW_
  print("  Electrostatic acceptor selection:")
  H2O_output = []
  with  mp.Pool(n_cpu_used) as pool:
    for i_H3O in range(n_H3O):
      #H2O_acceptor_number = 123
      if(H3O_hop_success_prev[i_H3O] == 'y'):  #We skip H3Os that were previously unsuccessful
        H2O_acceptor_number = pool.apply_async( electrostatic_H2O_selection, [i_H3O,n_atom,atom_x,atom_y,atom_z,atom_q,gro_curr_content,
                                         water_number_1,water_number_2,water_number_3,H3O,H3O_hop_success_prev] )           
      H2O_output.append(H2O_acceptor_number) #m1

    H2O = [H2O_acceptor_number.get() for H2O_acceptor_number in H2O_output]

  for i_H3O in range(n_H3O):
    if(H3O_hop_success_prev[i_H3O] == 'n'):  #for H3Os that do not hop, we keep the previous numbers
      H2O[i_H3O] = H2O_prev[i_H3O]  
##### END PARALLEL PROCESSING _NEW_


  for i_H3O in range(n_H3O):
    # storing information for this H2O acceptor molecule from gro file
    for i_line in range(2,n_atom+2):
      help_line_H2O   = gro_curr_content[i_line]
      recognition_tag = help_line_H2O[5:15]
      H2O_number      = int(help_line_H2O[0:5])
      if ((H2O_number == H2O[i_H3O]) and (recognition_tag == 'SOL     OW')):
        help_line = gro_curr_content[i_line  ].replace("\n","")
        gro_info_H2O_OW[i_H3O] = help_line
        help_line = gro_curr_content[i_line+1].replace("\n","")
        gro_info_H2O_HW1[i_H3O] = help_line
        help_line = gro_curr_content[i_line+2].replace("\n","")
        gro_info_H2O_HW2[i_H3O] = help_line
        help_line = gro_curr_content[i_line+3].replace("\n","")
        gro_info_H2O_MW[i_H3O] = help_line

  # Calling function that transfers proton from H3O to H2O, and leaves both donor and acceptor in proper orientations
  H3O_to_H2O()

  # Calling function that writes output gro file(s)
  write_output()

  # Calling function that writes data to log file
  write_log_file()

  # Remembering which water molecules served as proton acceptors previously.
  # This allows us to prevent back-and-forth hopping of a proton between the same donor/acceptor for i_hop > 1.
  H2O_prev = H2O   
                   
  # This is needed as an additional measure to prevent back-and-forth hopping of a proton
  # between the same donor/acceptor for i_hop > 1
  H3O_hop_success_prev = H3O_hop_success
                                          

#--- end of loop "for i_hop in range (n_hop):"

end_time = time.time()
wall_clock_time = end_time - start_time
print("  wall clock time = ","{:.2f}".format(wall_clock_time)," seconds")
print("  ")

# End of Main Program
