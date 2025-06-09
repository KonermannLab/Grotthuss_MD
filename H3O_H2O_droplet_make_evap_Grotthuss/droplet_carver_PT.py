#!/usr/bin/python
#
# droplet_carver_PT.py
# Written by: Lars Konermann, April 2022
#
# Program to carve droplets, e.g. for removing evaporated molecules and for shaping a cube into a spherical droplet.
# The code generates a new spherical gro and a new top file, from a gro and top input structure.
#
##########################################################################################################
# When using this program, please cite:                                                                  #
# "Release of Native-Like Gaseous Proteins from Electrospray Droplets via The Charged Residue Mechanism: #
# Insights from Molecular Dynamics Simulations" R. G. McAllister, H. Metwally, Y. Sun, and L. Konermann  #
# J. Am. Chem. Soc. 137, 12667-12676 (2015).                                                             #
##########################################################################################################
#
#
# The following components are accomodated so far:
# Protein  
# NH4  	  ammonium+  
# NH3  	  ammonia
# ACH  	  acetic acid-
# ACO  	  acetate
# NA  	  sodium+  
# CL  	  chloride-
# HCL  	  HCL molecule
# H3O  	  hydronium+
# SOL     TIP4P-2005 water
# For PT to work, the components HAVE TO BE IN THIS ORDER (missing components are OK) and
# the first component should not be zero
#
# Calling program:
# $ python ~/Python/droplet_carver_PT.py <input.gro> <output.gro> <input.top> <output.top> <r_cutoff> <center_mode> <log_file.txt> <PT=no/1/2> <PT_dist_max> <PT_dist_max> <PT_prob> 
# <r_cutoff>     new system radius in nm
# <center_mode>  determines how centering is done. Available options are    OW  - TIP4P-2005 water  (this might be problematic after droplet splitting)
#                                                                           PR  - protein
#                                                                           NH4 - ammonium
#                                                                           NA  - sodium
#                                                                           H3O - hydronium  
# <log_file.txt> name of the output file that contains the molecule numbers in <output.gro/top>
# <PT_dist_max>  maximum distance for proton transfer (nm)
# <PT=1>       proton transfer from NH4 to ACO will take place
# <PT=2>       proton transfer from H3O to CL will take place
# <PT=no>                  ... will not take place
# <PT_dist_max>  maximum proton hopping distance (nm)
# <PT_prob>      probability of proton transfer (0..1); negative value means that PT_prob is concentration-dependent (only implemented for HCL)

import sys
import numpy as np
import os.path
from os import path
import random

#------------------------------------------------------------------------------------------
# Start of function execute_PT_H3O_CL()
# This function perform proton transfer from H3O to CL,
# subject to the conditions imposed by PT_dist_max and PT_prob.
# Overall reaction: H3O+  +  CL-  --->   H2O  +  HCL
# 

def execute_PT_H3O_CL():

  global n_PT
  global gro_in_content
 

  protein_residues = [
                       'ACE','GLY','ALA','VAL','LEU','ILE','SER','THR','PRO','CYS','MET','ASN','ASP','GLN','GLU','LYS',
                       'ARG','HIS','PHE','TYR','TRP'
                      ]

  #================================ STEP 1: Identify suitable CL(proton acceptor)/H3O(proton donor) pairs

  H3O_used             = []  #prevents using the same H3O as donor for two different CLs; also ensures that used donor is not included in output
  CL_used              = []  #ensures that used acceptors are not included in output

  # The following are lines from the input gro file that contain info on donor and acceptor
  PT_info_H3O          = [] #gro info for donors before PT
  PT_info_CL_accpt_ID  = [] #gro info for acceptor atom before PT
  PT_info_H3O_donor_ID = [] #gro info for donor H atom before PT
  PT_info_HCL          = [] #gro info for acceptors after PT
  PT_info_H2O          = [] #gro info for donors after PT

  n_PT                 = 0  #counts H3O/CL pairs that will be subjected to proton transfer (although the PT attempt may not be succesful) 

  for i_CL_line in range(2,n_atom_in+2):

    CL_line = gro_in_content[i_CL_line]
    recognition_tag_CL = CL_line[5:15]
    if(recognition_tag_CL == 'CL      CL'):          #looking for CL
      CL_x      = float(CL_line[20:28])
      CL_y      = float(CL_line[28:36])
      CL_z      = float(CL_line[36:44])

      dist_min = 1000.
      
      for i_H3O_line in range(2,n_atom_in+2):

        H3O_line = gro_in_content[i_H3O_line]
        recognition_tag_H3O = H3O_line[5:15]
        if(recognition_tag_H3O == 'H3O     OW'):      #looking for H3Os, to identify the one that is closest to the current CL

          H3O_line1 = gro_in_content[i_H3O_line+1]      #identifying H3O    HW1
          H3O_H1_x = float(H3O_line1[20:28])
          H3O_H1_y = float(H3O_line1[28:36])
          H3O_H1_z = float(H3O_line1[36:44])
          dist = np.sqrt((H3O_H1_x-CL_x)**2 + (H3O_H1_y-CL_y)**2 + (H3O_H1_z-CL_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_H3O_line = H3O_line1
            j_start_of_closest_H3O = i_H3O_line
            j_closest_CL = i_CL_line

          H3O_line2 = gro_in_content[i_H3O_line+2]    #identifying H3O    HW2
          H3O_H2_x = float(H3O_line2[20:28])
          H3O_H2_y = float(H3O_line2[28:36])
          H3O_H2_z = float(H3O_line2[36:44])
          dist = np.sqrt((H3O_H2_x-CL_x)**2 + (H3O_H2_y-CL_y)**2 + (H3O_H2_z-CL_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_H3O_line = H3O_line2
            j_start_of_closest_H3O = i_H3O_line
            j_closest_CL = i_CL_line

          H3O_line3 = gro_in_content[i_H3O_line+3]    #identifying H3O    HW3
          H3O_H3_x = float(H3O_line3[20:28])
          H3O_H3_y = float(H3O_line3[28:36])
          H3O_H3_z = float(H3O_line3[36:44])
          dist = np.sqrt((H3O_H3_x-CL_x)**2 + (H3O_H3_y-CL_y)**2 + (H3O_H3_z-CL_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_H3O_line = H3O_line3
            j_start_of_closest_H3O = i_H3O_line
            j_closest_CL = i_CL_line

      # Now calculating distance between CL and the closest H3O hydrogen
      if(dist_min < 1000.):   #This statement ensures that distances are only calculated if a H3O was found
        closest_H3O_closest_H_x = float(closest_H3O_line[20:28])
        closest_H3O_closest_H_y = float(closest_H3O_line[28:36])
        closest_H3O_closest_H_z = float(closest_H3O_line[36:44])

        dist_from_CL = np.sqrt( (closest_H3O_closest_H_x - CL_x)**2 + 
                                (closest_H3O_closest_H_y - CL_y)**2 + 
                                (closest_H3O_closest_H_z - CL_z)**2 )
        
      # Suitable acceptor/donor pairs have to match several criteria:
      if(dist_min < PT_dist_max): 
        print(" ")
        print("  H3O/CL with H..CL distance of",format(dist_min, "f")," nm")
        print( CL_line[0:15] )
        print( closest_H3O_line[0:15] )         
        H3O_number = closest_H3O_line[0:5]
        CL_number = CL_line[0:5]
        if(H3O_number in H3O_used):         #donor has been used before
          print("  Proton transfer for this pair is cancelled because H3O has already been used")   
        if(H3O_number not in H3O_used):         #donor has not been used before
          rand_numb = random.uniform(0,1.0)
          if(rand_numb >= PT_prob):             #acceptor/donor pair may be skipped randomly as per PT_prob
            print("  Proton transfer for this pair is cancelled because of random selection") 
          if(rand_numb < PT_prob):               
            H3O_used.append(H3O_number)
            CL_used.append(CL_number)
            n_PT = n_PT + 1
            PT_info_H3O.append(gro_in_content[j_start_of_closest_H3O  ]) #complete H3O gro info of all donors
            PT_info_H3O.append(gro_in_content[j_start_of_closest_H3O+1])
            PT_info_H3O.append(gro_in_content[j_start_of_closest_H3O+2])
            PT_info_H3O.append(gro_in_content[j_start_of_closest_H3O+3])
            PT_info_H3O.append(gro_in_content[j_start_of_closest_H3O+4])
            PT_info_CL_accpt_ID.append(CL_line)                          #acceptor CL atoms
            PT_info_H3O_donor_ID.append(closest_H3O_line)                #donor H atom     (can be 1 of 3 for each H3O)


  #================================ STEP 2: For troubleshooting, write H3O/CL coordinates to .gro file

  x_before_PT_content = []
  help_line =('H3O/CL pairs before PT' + "\n")
  x_before_PT_content.append(help_line)
  help_line =(str(n_PT * 6) + "\n")
  x_before_PT_content.append(help_line)

  for i_PT in range(0,n_PT):
    help_line = PT_info_CL_accpt_ID[i_PT]
    x_before_PT_content.append(help_line)
    for i_H3O_line in range(i_PT*5,5+i_PT*5):
      help_line = PT_info_H3O[i_H3O_line]
      x_before_PT_content.append(help_line)

  help_line =('0.0     0.0     0.0' + "\n")
  x_before_PT_content.append(help_line)

  x_before_PT = open('x_CL_H3O_before_PT.gro',"w+")
  x_before_PT.writelines(x_before_PT_content)
  x_before_PT.close()

  #================================ STEP 3: Turn CL into HCL, and turn H3O into H2O

  for i_PT in range(0,n_PT):
    
    #======================= STEP 3A: CL --> HCL

    help_line = PT_info_CL_accpt_ID[i_PT]  
    acceptor_oxygen = help_line[13:15]

    #Changing CL atom name to HCL atom name
    help_line = PT_info_CL_accpt_ID[i_PT]
    new_line = (help_line[0:5] + 'HCL    CL1' + help_line[15:68] + "\n")
    PT_info_HCL.append(new_line)

    # Generating new HCL hydrogen (HCL1) coordinates:
    # The new hydrogen is placed on the OW -  CL connecting line, with bond distance 0.127 nm
 

    help_line = PT_info_H3O[i_PT*5]
    #print(help_line) #m1
    O_xyz = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])
    help_line = PT_info_CL_accpt_ID[i_PT]
    #print(help_line) #m1
    CL_xyz = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])

    #calculating difference vector
    CL_minus_O = CL_xyz - O_xyz 

    CL_minus_O_norm = CL_minus_O / np.linalg.norm(CL_minus_O) * 0.127

    x_new_HCL1 = CL_xyz[0] - CL_minus_O_norm[0]
    y_new_HCL1 = CL_xyz[1] - CL_minus_O_norm[1]
    z_new_HCL1 = CL_xyz[2] - CL_minus_O_norm[2]

    new_line_H = (help_line[0:5] + 'HCL   HCL1' + help_line[15:20] + 
               format(x_new_HCL1, '8.3f')                               +
               format(y_new_HCL1, '8.3f')                               + 
               format(z_new_HCL1, '8.3f')                               +
               help_line[44:68] + "\n")
    PT_info_HCL.append(new_line_H)
    #print(new_line_H) #m1


  
    #======================= STEP 3B: H3O --> H2O
    # For producing H2O from H3O we remove the corresponding proton.
    # Atomic positions of the remaining two hydrogens are left unchanged.
    # We rely on Gromacs to subsequently generate a proper H2O geometry from the distorted "nascent" H2O.
    # MW is moved to its proper position in the H-O-H plane, with MW-OW distance of 0.015 nm.

    # Changing H3O names to H2O names

    help_line = PT_info_H3O[i_PT*5]       #turning H3O OW into H2O OW
    new_line = (help_line[0:5] + 'SOL     OW' + help_line[15:68] + "\n")
    PT_info_H2O.append(new_line)

    help_line = PT_info_H3O_donor_ID[i_PT]  #looking up which of the three Hs gets removed
    transfer_proton = help_line[12:15]

    if(transfer_proton == 'HW1'):
      help_line = PT_info_H3O[i_PT*5 +2]    #turning HW2 into HW1
      new_line = (help_line[0:5] + 'SOL    HW1' + help_line[15:68] + "\n")
      PT_info_H2O.append(new_line)

      help_line = PT_info_H3O[i_PT*5 +3]    #turning HW3 into HW2
      new_line = (help_line[0:5] + 'SOL    HW2' + help_line[15:68] + "\n")
      PT_info_H2O.append(new_line)

    if(transfer_proton == 'HW2'):
      help_line = PT_info_H3O[i_PT*5 +1]    #leaving HW1 as HW1
      new_line = (help_line[0:5] + 'SOL    HW1' + help_line[15:68] + "\n")
      PT_info_H2O.append(new_line)

      help_line = PT_info_H3O[i_PT*5 +3]    #turning HW3 into HW2
      new_line = (help_line[0:5] + 'SOL    HW2' + help_line[15:68] + "\n")
      PT_info_H2O.append(new_line)

    if(transfer_proton == 'HW3'):
      help_line = PT_info_H3O[i_PT*5 +1]    #leaving HW1 as HW1
      new_line = (help_line[0:5] + 'SOL    HW1' + help_line[15:68] + "\n")
      PT_info_H2O.append(new_line)

      help_line = PT_info_H3O[i_PT*5 +2]    #leaving HW2 as HW2
      new_line = (help_line[0:5] + 'SOL    HW2' + help_line[15:68] + "\n")
      PT_info_H2O.append(new_line)

    # Now calculating the new MW position. As a first step, the newly generated H2O is translated so that OW is at 0/0/0

    help_line = PT_info_H2O[ len(PT_info_H2O)-3] #OW
    H2O_OW_xyz = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])

    help_line = PT_info_H2O[ len(PT_info_H2O)-2] #HW1
    H2O_HW1_xyz_center = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])]) - H2O_OW_xyz

    help_line = PT_info_H2O[ len(PT_info_H2O)-1] #HW2
    H2O_HW2_xyz_center = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])]) - H2O_OW_xyz

    # Calculating new MW coordinates in centered coordinate system
    H2O_HW1_plus_HW2_xyz_center = H2O_HW1_xyz_center + H2O_HW2_xyz_center
    MW_center_xyz = H2O_HW1_plus_HW2_xyz_center/ np.linalg.norm(H2O_HW1_plus_HW2_xyz_center) * 0.015

    # Translating MW back to actual coordinate system, and generating a proper gro line
    MW_xyz =  MW_center_xyz + H2O_OW_xyz
    help_line = PT_info_H3O[i_PT*5+4]       #turning H3O MW into H2O MW

    new_line = (help_line[0:5] + 'SOL     MW' + help_line[15:20]       + 
               format(MW_xyz[0], '8.3f')                               +
               format(MW_xyz[1], '8.3f')                               + 
               format(MW_xyz[2], '8.3f')                               +
               help_line[44:68] + "\n")
    PT_info_H2O.append(new_line)


  #================================ STEP 4: For troubleshooting, write HCL/H2O coordinates to .gro file

  x_after_PT_content = []
  help_line =('HCL/H2O pairs after PT' + "\n")
  x_after_PT_content.append(help_line)
  help_line =(str(n_PT * 2 + n_PT * 4) + "\n")
  x_after_PT_content.append(help_line)

  for i_PT in range(0,n_PT):

    for i_HCL_line in range(i_PT*2,2+i_PT*2):
      help_line = PT_info_HCL[i_HCL_line]
      x_after_PT_content.append(help_line)
    for i_H2O_line in range(i_PT*4,4+i_PT*4):
      help_line = PT_info_H2O[i_H2O_line]
      x_after_PT_content.append(help_line)

  help_line =('0.0     0.0     0.0' + "\n")
  x_after_PT_content.append(help_line)

  x_after_PT = open('x_HCL_H2O_after_PT.gro',"w+")
  x_after_PT.writelines(x_after_PT_content)
  x_after_PT.close()


  #================================ STEP 5: filling gro_in_after_PT_no_carve_content with new information (after PT)

  gro_in_after_PT_no_carve_content = []
                                                 #header line is transferred
  help_line = ('droplet after proton transfer, prior to carving'  + "\n")                 
  gro_in_after_PT_no_carve_content.append(help_line)

  help_line = gro_in_content[1]                  #number of atoms is transferred
  gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #protein information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type in protein_residues):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #NH4 information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'NH4'):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #NH3 information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'NH3'):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #ACH information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'ACH'):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #ACO information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'ACO'):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #NA  information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'NA '):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #all CL that were *not* converted to HCL are transferred
    transfer = 'yes'
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'CL '):
      for i_PT in range (0,n_PT):
        if(CL_used[i_PT] in help_line[0:5]):
          transfer = 'no'
      if(transfer == 'yes'):
        gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #all orignal HCL are transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'HCL'):
      gro_in_after_PT_no_carve_content.append(help_line)
        
  for i_PT in range(0,n_PT):                     #all *new* HCL are transferred
    for i_HCL_line in range(i_PT*2,2+i_PT*2):
      help_line = PT_info_HCL[i_HCL_line]
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #all H3O that were *not* converted to H2O are transferred
    transfer = 'yes'
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'H3O'):
      for i_PT in range (0,n_PT):
        if(H3O_used[i_PT] in help_line[0:5]):
          transfer = 'no'
      if(transfer == 'yes'):
        gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #all original waters are transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'SOL'):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_PT in range(0,n_PT):                     #all *new* waters are transferred
    for i_H2O_line in range(i_PT*4,4+i_PT*4):
      help_line = PT_info_H2O[i_H2O_line]
      gro_in_after_PT_no_carve_content.append(help_line)

  help_line = gro_in_content[n_atom_in + 2]      #number of atoms is transferred
  gro_in_after_PT_no_carve_content.append(help_line)

  # For troubleshooting, gro_in_after_PT_no_carve_content is saved as a .gro file.
  # It contains the entire droplet after proton transfer, but prior to carving

  x_droplet_after_PT_no_carve = open('x_droplet_after_PT_no_carve.gro',"w+")
  x_droplet_after_PT_no_carve.writelines(gro_in_after_PT_no_carve_content)
  x_droplet_after_PT_no_carve.close()

  #================================ STEP 6: redefining gro_in_content. It has undergone PT, but no carving yet.

  gro_in_content =  gro_in_after_PT_no_carve_content


  # End of function execute_PT_H3O_CL()

#------------------------------------------------------------------------------------------
# Start of function execute_PT_NH4_ACO()
# This function perform proton transfer from NH4 to ACO,
# subject to the conditions imposed by PT_dist_max and PT_prob.
# Overall reaction: ACO-  +  NH4+  --->   ACH   + NH3
# 

def execute_PT_NH4_ACO():

  global n_PT
  global gro_in_content
 

  protein_residues = [
                       'ACE','GLY','ALA','VAL','LEU','ILE','SER','THR','PRO','CYS','MET','ASN','ASP','GLN','GLU','LYS',
                       'ARG','HIS','PHE','TYR','TRP'
                      ]

  #================================ STEP 1: Identify suitable ACO(proton acceptor)/NH4(proton donor) pairs

  NH4_used             = []  #prevents using the same NH4 as donor for two different ACOs; also ensures that used donor is not included in output
  ACO_used             = []  #ensures that used acceptors is not included in output

  # The following are lines from the input gro file that contain info on donor and acceptor
  PT_info_ACO          = [] #gro info for acceptors before PT
  PT_info_NH4          = [] #gro info for donors before PT
  PT_info_ACO_accpt_ID = [] #gro info for acceptor O atom before PT
  PT_info_NH4_donor_ID = [] #gro info for donor H atom before PT
  PT_info_ACH          = [] #gro info for acceptors after PT
  PT_info_NH3          = [] #gro info for donors after PT

  n_PT            = 0  #counts NH4/ACO pairs that will be subjected to proton transfer (although the PT attempt may not be succesful) 

  for i_ACO_line in range(2,n_atom_in+2):


    ACO_line = gro_in_content[i_ACO_line]
    recognition_tag_ACO = ACO_line[5:15]
    if(recognition_tag_ACO == 'ACO     C1'):          #looking for ACO
      ACO_line1     = gro_in_content[i_ACO_line+5]    #identifying ACO  O1
      ACO_O1_x      = float(ACO_line1[20:28])
      ACO_O1_y      = float(ACO_line1[28:36])
      ACO_O1_z      = float(ACO_line1[36:44])

      ACO_line2     = gro_in_content[i_ACO_line+6]    #identifying ACO  O2
      ACO_O2_x      = float(ACO_line2[20:28])
      ACO_O2_y      = float(ACO_line2[28:36])
      ACO_O2_z      = float(ACO_line2[36:44])

      ACO_Ocenter_x = 0.5*(ACO_O1_x + ACO_O2_x)       # ACO_Ocenter is the point (*) between the two O atoms, 
      ACO_Ocenter_y = 0.5*(ACO_O1_y + ACO_O2_y)       # not too far from where the H would be located after PT.
      ACO_Ocenter_z = 0.5*(ACO_O1_z + ACO_O2_z)       #    O\
                                                      #    * C-- 
                                                      #    O/
                                                      # Selecting NH4 donors based on their proximity to * reduces atom
                                                      # clashes after PT, because the NH4 will likely displace water molecules
                                                      # that could otherwise be a problem.

      dist_min = 1000.
      
      for i_NH4_line in range(2,n_atom_in+2):

        NH4_line = gro_in_content[i_NH4_line]
        recognition_tag_NH4 = NH4_line[5:15]
        if(recognition_tag_NH4 == 'NH4    HZ1'):      #looking for NH4s, to identify the one that is closest to the current ACO

          NH4_line1 = gro_in_content[i_NH4_line]      #identifying NH4  HZ1
          NH4_H1_x = float(NH4_line1[20:28])
          NH4_H1_y = float(NH4_line1[28:36])
          NH4_H1_z = float(NH4_line1[36:44])
          dist = np.sqrt((NH4_H1_x-ACO_O1_x)**2 + (NH4_H1_y-ACO_O1_y)**2 + (NH4_H1_z-ACO_O1_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_NH4_line = NH4_line1
            closest_ACO_line = ACO_line1
            j_start_of_closest_NH4 = i_NH4_line
            j_start_of_closest_ACO = i_ACO_line
          dist = np.sqrt((NH4_H1_x-ACO_O2_x)**2 + (NH4_H1_y-ACO_O2_y)**2 + (NH4_H1_z-ACO_O2_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_NH4_line = NH4_line1
            closest_ACO_line = ACO_line2
            j_start_of_closest_NH4 = i_NH4_line
            j_start_of_closest_ACO = i_ACO_line

          NH4_line2 = gro_in_content[i_NH4_line+2]    #identifying NH4  HZ2
          NH4_H2_x = float(NH4_line2[20:28])
          NH4_H2_y = float(NH4_line2[28:36])
          NH4_H2_z = float(NH4_line2[36:44])
          dist = np.sqrt((NH4_H2_x-ACO_O1_x)**2 + (NH4_H2_y-ACO_O1_y)**2 + (NH4_H2_z-ACO_O1_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_NH4_line = NH4_line2
            closest_ACO_line = ACO_line1
            j_start_of_closest_NH4 = i_NH4_line
            j_start_of_closest_ACO = i_ACO_line
          dist = np.sqrt((NH4_H2_x-ACO_O2_x)**2 + (NH4_H2_y-ACO_O2_y)**2 + (NH4_H2_z-ACO_O2_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_NH4_line = NH4_line2
            closest_ACO_line = ACO_line2
            j_start_of_closest_NH4 = i_NH4_line
            j_start_of_closest_ACO = i_ACO_line

          NH4_line3 = gro_in_content[i_NH4_line+3]    #identifying NH4  HZ3
          NH4_H3_x = float(NH4_line3[20:28])
          NH4_H3_y = float(NH4_line3[28:36])
          NH4_H3_z = float(NH4_line3[36:44])
          dist = np.sqrt((NH4_H3_x-ACO_O1_x)**2 + (NH4_H3_y-ACO_O1_y)**2 + (NH4_H3_z-ACO_O1_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_NH4_line = NH4_line3
            closest_ACO_line = ACO_line1
            j_start_of_closest_NH4 = i_NH4_line
            j_start_of_closest_ACO = i_ACO_line
          dist = np.sqrt((NH4_H3_x-ACO_O2_x)**2 + (NH4_H3_y-ACO_O2_y)**2 + (NH4_H3_z-ACO_O2_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_NH4_line = NH4_line3
            closest_ACO_line = ACO_line2
            j_start_of_closest_NH4 = i_NH4_line
            j_start_of_closest_ACO = i_ACO_line

          NH4_line4 = gro_in_content[i_NH4_line+4]    #identifying NH4  HZ41
          NH4_H4_x = float(NH4_line4[20:28])
          NH4_H4_y = float(NH4_line4[28:36])
          NH4_H4_z = float(NH4_line4[36:44])
          dist = np.sqrt((NH4_H4_x-ACO_O1_x)**2 + (NH4_H4_y-ACO_O1_y)**2 + (NH4_H4_z-ACO_O1_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_NH4_line = NH4_line4
            closest_ACO_line = ACO_line1
            j_start_of_closest_NH4 = i_NH4_line
            j_start_of_closest_ACO = i_ACO_line
          dist = np.sqrt((NH4_H4_x-ACO_O2_x)**2 + (NH4_H4_y-ACO_O2_y)**2 + (NH4_H4_z-ACO_O2_z)**2)
          if(dist < dist_min):
            dist_min = dist
            closest_NH4_line = NH4_line4
            closest_ACO_line = ACO_line2
            j_start_of_closest_NH4 = i_NH4_line
            j_start_of_closest_ACO = i_ACO_line

      # Now calculating distance between ACO_Ocenter and the closest NH4 hydrogen
      if(dist_min < 1000.):   #This statement ensures that distances are only calculated if a NH4 was found
        closest_NH4_closest_H_x = float(closest_NH4_line[20:28])
        closest_NH4_closest_H_y = float(closest_NH4_line[28:36])
        closest_NH4_closest_H_z = float(closest_NH4_line[36:44])

        dist_from_ACO_Ocenter = np.sqrt( (closest_NH4_closest_H_x - ACO_Ocenter_x)**2 + 
                                         (closest_NH4_closest_H_y - ACO_Ocenter_y)**2 + 
                                         (closest_NH4_closest_H_z - ACO_Ocenter_z)**2 )


      # Suitable acceptor/donor pairs have to match several criteria:
      #if(dist_min < PT_dist_max):               # with this line, acceptor and donor are chosen based on H..O distance
      if(dist_from_ACO_Ocenter < PT_dist_max):               # with this line, acceptor and donor are chosen based on H..* distance
        print(" ")
        print("  NH4/ACO with H..O distance of",format(dist_min, "f")," nm and H..* distance of",format(dist_from_ACO_Ocenter, "f")," nm")
        print( closest_ACO_line[0:15] )
        print( closest_NH4_line[0:15] )         
        NH4_number = closest_NH4_line[0:5]
        ACO_number = closest_ACO_line[0:5]
        if(NH4_number in NH4_used):         #donor has been used before
          print("  Proton transfer for this pair is cancelled because NH4 has already been used")   
        if(NH4_number not in NH4_used):         #donor has not been used before
          rand_numb = random.uniform(0,1.0)
          if(rand_numb >= PT_prob):             #acceptor/donor pair may be skipped randomly as per PT_prob
            print("  Proton transfer for this pair is cancelled because of random selection") 
          if(rand_numb < PT_prob):               
            NH4_used.append(NH4_number)
            ACO_used.append(ACO_number)
            n_PT = n_PT + 1
            PT_info_ACO.append(gro_in_content[j_start_of_closest_ACO  ]) #complete ACO gro info of all acceptors
            PT_info_ACO.append(gro_in_content[j_start_of_closest_ACO+1])
            PT_info_ACO.append(gro_in_content[j_start_of_closest_ACO+2])
            PT_info_ACO.append(gro_in_content[j_start_of_closest_ACO+3])
            PT_info_ACO.append(gro_in_content[j_start_of_closest_ACO+4])
            PT_info_ACO.append(gro_in_content[j_start_of_closest_ACO+5])
            PT_info_ACO.append(gro_in_content[j_start_of_closest_ACO+6])
            PT_info_NH4.append(gro_in_content[j_start_of_closest_NH4  ]) #complete NH4 gro info of all donors
            PT_info_NH4.append(gro_in_content[j_start_of_closest_NH4+1])
            PT_info_NH4.append(gro_in_content[j_start_of_closest_NH4+2])
            PT_info_NH4.append(gro_in_content[j_start_of_closest_NH4+3])
            PT_info_NH4.append(gro_in_content[j_start_of_closest_NH4+4])
            PT_info_ACO_accpt_ID.append(closest_ACO_line)                     #acceptor O atoms (can be 1 of 2 for each ACO)
            PT_info_NH4_donor_ID.append(closest_NH4_line)                     #donor H atom     (can be 1 of 4 for each NH4)


  #================================ STEP 2: For troubleshooting, write NH4/ACO coordinates to .gro file

  x_before_PT_content = []
  help_line =('NH4/ACO pairs before PT' + "\n")
  x_before_PT_content.append(help_line)
  help_line =(str(n_PT * 12) + "\n")
  x_before_PT_content.append(help_line)

  for i_PT in range(0,n_PT):

    for i_ACO_line in range(i_PT*7,7+i_PT*7):
      help_line = PT_info_ACO[i_ACO_line]
      x_before_PT_content.append(help_line)
    for i_NH4_line in range(i_PT*5,5+i_PT*5):
      help_line = PT_info_NH4[i_NH4_line]
      x_before_PT_content.append(help_line)

  help_line =('0.0     0.0     0.0' + "\n")
  x_before_PT_content.append(help_line)

  x_before_PT = open('x_ACO_NH4_before_PT.gro',"w+")
  x_before_PT.writelines(x_before_PT_content)
  x_before_PT.close()

  #================================ STEP 3: Turn ACO into ACH, and turn NH4 into NH3

  for i_PT in range(0,n_PT):
    
    #======================= STEP 3A: ACO --> ACH

    help_line = PT_info_ACO_accpt_ID[i_PT]  #looking up which of the two Os gets protonated
    acceptor_oxygen = help_line[13:15]

    #Changing ACO atom names to ACH atom names
    help_line = PT_info_ACO[i_PT*7   ]      #turning C1 into C2
    new_line = (help_line[0:5] + 'ACH     C2' + help_line[15:68] + "\n")
    PT_info_ACH.append(new_line)

    help_line = PT_info_ACO[i_PT*7 +1]      #turning C2 into C1
    new_line = (help_line[0:5] + 'ACH     C1' + help_line[15:68] + "\n")
    PT_info_ACH.append(new_line)

    help_line = PT_info_ACO[i_PT*7 +4]      #turning H3 into H21
    new_line = (help_line[0:5] + 'ACH    H21' + help_line[15:68] + "\n")
    PT_info_ACH.append(new_line)

    help_line = PT_info_ACO[i_PT*7 +3]      #turning H2 into H22
    new_line = (help_line[0:5] + 'ACH    H22' + help_line[15:68] + "\n")
    PT_info_ACH.append(new_line)

    help_line = PT_info_ACO[i_PT*7 +2]      #turning H1 into H23
    new_line = (help_line[0:5] + 'ACH    H23' + help_line[15:68] + "\n")
    PT_info_ACH.append(new_line)


    if(acceptor_oxygen == 'O1'):
      help_line = PT_info_ACO[i_PT*7 +6]      #turning O1 into O2
      new_line = (help_line[0:5] + 'ACH     O2' + help_line[15:68] + "\n")
      PT_info_ACH.append(new_line)

      help_line = PT_info_ACO[i_PT*7 +5]      #turning O2 into O1. HO1 is attached to O1 in ACH
      new_line = (help_line[0:5] + 'ACH     O1' + help_line[15:68] + "\n")
      PT_info_ACH.append(new_line)

    if(acceptor_oxygen == 'O2'):
      help_line = PT_info_ACO[i_PT*7 +5]      #turning O2 into O2
      new_line = (help_line[0:5] + 'ACH     O2' + help_line[15:68] + "\n")
      PT_info_ACH.append(new_line)

      help_line = PT_info_ACO[i_PT*7 +6]      #turning O1 into O1. HO1 is attached to O1 in ACH
      new_line = (help_line[0:5] + 'ACH     O1' + help_line[15:68] + "\n")
      PT_info_ACH.append(new_line)


    #generating new HO1, based on 'ACH     O1' coordinates
    # 
    # distance(HO1-O1) = 0.096 nm,  distance(HO1-C1) = 0.189 nm, determined from ACH.gro
    #
    # HO1-O1
    #      \
    #       C1 -- C2   (ACH atom names)
    #      /
    #     O2 
    #
 
    help_line = PT_info_ACH[i_PT*8 +1]
    C1_xyz = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])
    help_line = PT_info_ACH[i_PT*8 +5]
    O2_xyz = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])
    help_line = PT_info_ACH[i_PT*8 +6]
    O1_xyz = np.array([float(help_line[20:28]),float(help_line[28:36]),float(help_line[36:44])])
      
    #calculating difference vectors
    O1_minus_C1 = O1_xyz - C1_xyz
    O2_minus_C1 = O2_xyz - C1_xyz 
      
    #calculating new H position based on (ACH     O1) coordinates. OH bond distance = 0.096 nm
    displacement_vector = (O1_minus_C1 + 1.2 * O2_minus_C1)  #1.2 takes care of distance(HO1-C1), determined using trial and error  
    displacement_vector_norm = displacement_vector / np.linalg.norm(displacement_vector) * 0.096

    x_new_HO1 = O1_xyz[0] + displacement_vector_norm[0]
    y_new_HO1 = O1_xyz[1] + displacement_vector_norm[1]
    z_new_HO1 = O1_xyz[2] + displacement_vector_norm[2]

    new_line_H = (help_line[0:5] + 'ACH    HO1' + help_line[15:20] + 
               format(x_new_HO1, '8.3f')                               +
               format(y_new_HO1, '8.3f')                               + 
               format(z_new_HO1, '8.3f')                               +
               help_line[44:68] + "\n")
    PT_info_ACH.append(new_line_H)


  
    #======================= STEP 3B: NH4 --> NH3
    # For producing NH3 from NH4 we simply remove the corresponding proton.
    # Atomic positions of the other four atoms remain unchanged.
    # We rely on Gromacs to subsequently generate a proper NH3 geometry from the distorted "nascent" NH3.
    # Not sure if it matters, but the chirality of HZ 1/2/3/4 is retained during conversion to H1 1/2/3.

    # Changing NH4 names to NH3 names
    help_line = PT_info_NH4[i_PT*5 +1]      #turning NZ into N1
    new_N_line = (help_line[0:5] + 'NH3     N1' + help_line[15:68] + "\n")
    PT_info_NH3.append(new_N_line)

    help_line = PT_info_NH4_donor_ID[i_PT]  #looking up which of the four Hs gets removed
    transfer_proton = help_line[12:15]

    if(transfer_proton == 'HZ1'):
      help_line = PT_info_NH4[i_PT*5 +2]    #turning HZ2 into H11
      new_H_line = (help_line[0:5] + 'NH3    H11' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

      help_line = PT_info_NH4[i_PT*5 +3]    #turning HZ3 into H12
      new_H_line = (help_line[0:5] + 'NH3    H12' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

      help_line = PT_info_NH4[i_PT*5 +4]    #turning HZ4 into H13
      new_H_line = (help_line[0:5] + 'NH3    H13' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

    if(transfer_proton == 'HZ2'):
      help_line = PT_info_NH4[i_PT*5 +0]    #turning HZ1 into H11
      new_H_line = (help_line[0:5] + 'NH3    H11' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

      help_line = PT_info_NH4[i_PT*5 +4]    #turning HZ4 into H12
      new_H_line = (help_line[0:5] + 'NH3    H12' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

      help_line = PT_info_NH4[i_PT*5 +3]    #turning HZ3 into H13
      new_H_line = (help_line[0:5] + 'NH3    H13' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

    if(transfer_proton == 'HZ3'):
      help_line = PT_info_NH4[i_PT*5 +0]    #turning HZ1 into H11
      new_H_line = (help_line[0:5] + 'NH3    H11' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

      help_line = PT_info_NH4[i_PT*5 +2]    #turning HZ2 into H12
      new_H_line = (help_line[0:5] + 'NH3    H12' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

      help_line = PT_info_NH4[i_PT*5 +4]    #turning HZ4 into H13
      new_H_line = (help_line[0:5] + 'NH3    H13' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

    if(transfer_proton == 'HZ4'):
      help_line = PT_info_NH4[i_PT*5 +0]    #turning HZ1 into H11
      new_H_line = (help_line[0:5] + 'NH3    H11' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

      help_line = PT_info_NH4[i_PT*5 +3]    #turning HZ3 into H12
      new_H_line = (help_line[0:5] + 'NH3    H12' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)

      help_line = PT_info_NH4[i_PT*5 +2]    #turning HZ2 into H13
      new_H_line = (help_line[0:5] + 'NH3    H13' + help_line[15:68] + "\n")
      PT_info_NH3.append(new_H_line)


  #================================ STEP 4: For troubleshooting, write ACH/NH3 coordinates to .gro file

  x_after_PT_content = []
  help_line =('ACH/NH3 pairs after PT' + "\n")
  x_after_PT_content.append(help_line)
  help_line =(str(n_PT * 4 + n_PT * 8) + "\n")
  x_after_PT_content.append(help_line)

  for i_PT in range(0,n_PT):

    for i_ACH_line in range(i_PT*8,8+i_PT*8):
      help_line = PT_info_ACH[i_ACH_line]
      x_after_PT_content.append(help_line)
    for i_NH3_line in range(i_PT*4,4+i_PT*4):
      help_line = PT_info_NH3[i_NH3_line]
      x_after_PT_content.append(help_line)

  help_line =('0.0     0.0     0.0' + "\n")
  x_after_PT_content.append(help_line)

  x_after_PT = open('x_ACH_NH3_after_PT.gro',"w+")
  x_after_PT.writelines(x_after_PT_content)
  x_after_PT.close()


  #================================ STEP 5: filling gro_in_after_PT_no_carve_content with new information (after PT)

  gro_in_after_PT_no_carve_content = []
                                                 #header line is transferred
  help_line = ('droplet after proton transfer, prior to carving'  + "\n")                 
  gro_in_after_PT_no_carve_content.append(help_line)

  help_line = gro_in_content[1]                  #number of atoms is transferred
  gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #protein information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type in protein_residues):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #all NH4 that were *not* converted to NH3 are transferred
    transfer = 'yes'
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'NH4'):
      for i_PT in range (0,n_PT):
        if(NH4_used[i_PT] in help_line[0:5]):
          transfer = 'no'
      if(transfer == 'yes'):
        gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #all orignal NH3 are transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'NH3'):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_PT in range(0,n_PT):                     #all *new* NH3 are transferred
    for i_NH3_line in range(i_PT*4,4+i_PT*4):
      help_line = PT_info_NH3[i_NH3_line]
      gro_in_after_PT_no_carve_content.append(help_line)
 
  for i_line in range(2,n_atom_in+2):            #all orignal ACH are transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'ACH'):
      gro_in_after_PT_no_carve_content.append(help_line)
        
  for i_PT in range(0,n_PT):                     #all *new* ACH are transferred
    for i_ACH_line in range(i_PT*8,8+i_PT*8):
      help_line = PT_info_ACH[i_ACH_line]
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #all ACO that were *not* converted to ACH are transferred
    transfer = 'yes'
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'ACO'):
      for i_PT in range (0,n_PT):
        if(ACO_used[i_PT] in help_line[0:5]):
          transfer = 'no'
      if(transfer == 'yes'):
        gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #NA information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'NA '):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #CL information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'CL '):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #HCL information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'HCL'):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #H3O information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'H3O'):
      gro_in_after_PT_no_carve_content.append(help_line)

  for i_line in range(2,n_atom_in+2):            #SOL information is transferred
    help_line = gro_in_content[i_line]
    molecule_type   = help_line[5:8]
    if (molecule_type == 'SOL'):
      gro_in_after_PT_no_carve_content.append(help_line)

  help_line = gro_in_content[n_atom_in + 2]      #number of atoms is transferred
  gro_in_after_PT_no_carve_content.append(help_line)

  # For troubleshooting, gro_in_after_PT_no_carve_content is saved as a .gro file.
  # It contains the entire droplet after proton transfer, but prior to carving

  x_droplet_after_PT_no_carve = open('x_droplet_after_PT_no_carve.gro',"w+")
  x_droplet_after_PT_no_carve.writelines(gro_in_after_PT_no_carve_content)
  x_droplet_after_PT_no_carve.close()

  #================================ STEP 6: redefining gro_in_content. It has undergone PT, but no carving yet.

  gro_in_content =  gro_in_after_PT_no_carve_content


  # End of function execute_PT_NH4_ACO()

#------------------------------------------------------------------------------------------
# Start of function write_log()
# This function writes the current droplet composition to a log file.
# If the log file is not present, the program generates a new one.
# If the log file is already present, the program adds to it.

def write_log():

  operation = 0

  log_present = 'yes'
  path.exists(log_name)
  if (os.path.isfile(log_name) == False):
    log_present = 'no'

  if(log_present == 'no'):
    log = open(log_name,"w+")
    log_content = []
    help_line = ('operation    NH4     NH3     ACH     ACO      NA      CL     HCL     H3O   water    n_PT'  + "\n")
    log_content.append(help_line)
    help_line = str( format(operation, "8") + format(n_NH4_out, "8") + format(n_NH3_out, "8") 
                  +  format(n_ACH_out, "8") + format(n_ACO_out, "8") + format(n_NA_out , "8") 
                  +  format(n_CL_out , "8") + format(n_HCL_out, "8") + format(n_H3O_out, "8") 
                  +  format(n_SOL_out, "8") + format(n_PT, "8") + "\n" ) 
    log_content.append(help_line)
    log.writelines(log_content)
    log.close()

  if(log_present == 'yes'):
    log = open(log_name,"r+")
    log_content = log.readlines()
    n_log_lines = len(log_content)
    help_line   = log_content[n_log_lines-1]
    operation   = int(help_line[0:10])
    operation   = operation + 1

    help_line = str( format(operation, "8") + format(n_NH4_out, "8") + format(n_NH3_out, "8") 
                  +  format(n_ACH_out, "8") + format(n_ACO_out, "8") + format(n_NA_out , "8") 
                  +  format(n_CL_out , "8") + format(n_HCL_out, "8") + format(n_H3O_out, "8") 
                  +  format(n_SOL_out, "8") + format(n_PT, "8") + "\n" ) 
    log_content.append(help_line)
    log.close() 

    log = open(log_name,"w+")   
    log.writelines(log_content)
    log.close()


  # End of function write_log()

#------------------------------------------------------------------------------------------
# Start of function count_molecules_pre_PT():
# This function counts initial molecule numbers, prior to any PT

def count_molecules_pre_PT():

  global n_SOL_in_pre_PT  #number of TIP4P-2005 water
  global n_ACO_in_pre_PT  #number of acetate 
  global n_ACH_in_pre_PT  #number of acetatic acid
  global n_NH3_in_pre_PT  #number of NH3
  global n_NH4_in_pre_PT  #number of NH4
  global n_NA_in_pre_PT   #number of NA
  global n_CL_in_pre_PT   #number of CL
  global n_HCL_in_pre_PT  #number of HCL
  global n_H3O_in_pre_PT  #number of H3O

  n_SOL_in_pre_PT  = 0
  n_ACO_in_pre_PT  = 0 
  n_ACH_in_pre_PT  = 0 
  n_NH3_in_pre_PT  = 0 
  n_NH4_in_pre_PT  = 0 
  n_NA_in_pre_PT   = 0 
  n_CL_in_pre_PT   = 0 
  n_HCL_in_pre_PT  = 0 
  n_H3O_in_pre_PT  = 0 


  for i_line in range(2,n_atom_in+2):
    help_line = gro_in_content[i_line]
    recognition_tag = help_line[5:15]

    if(recognition_tag == 'NA      NA'):   #----------------------  NA
      n_NA_in_pre_PT = n_NA_in_pre_PT + 1

    if(recognition_tag == 'CL      CL'):   #----------------------  CL
      n_CL_in_pre_PT = n_CL_in_pre_PT + 1

    if(recognition_tag == 'HCL    CL1'):   #----------------------  HCL
      n_HCL_in_pre_PT = n_HCL_in_pre_PT + 1

    if(recognition_tag == 'H3O     OW'):   #----------------------  H3O
      n_H3O_in_pre_PT = n_H3O_in_pre_PT + 1

    if(recognition_tag == 'SOL     OW'):   #----------------------  TIP4P-2005 water (SOL)
      n_SOL_in_pre_PT = n_SOL_in_pre_PT + 1

    if(recognition_tag == 'ACO     C1'):   #----------------------  acetate (ACO)
      n_ACO_in_pre_PT = n_ACO_in_pre_PT + 1

    if(recognition_tag == 'ACH     C2'):   #----------------------  acetatic acid (ACH)
      n_ACH_in_pre_PT = n_ACH_in_pre_PT + 1

    if(recognition_tag == 'NH3     N1'):   #----------------------  ammonium (NH3)
      n_NH3_in_pre_PT = n_NH3_in_pre_PT + 1

    if(recognition_tag == 'NH4    HZ1'):   #----------------------  ammonia (NH4)
      n_NH4_in_pre_PT = n_NH4_in_pre_PT + 1



  # End of function count_molecules_pre_PT()

#------------------------------------------------------------------------------------------
# Start of function generate_top():
# This function generates a new top file by updating [ molecule ] numbers

def generate_top():

  global top_out_content     #this will turn into the output top file
  
  n_top_lines = len(top_in_content) #number of lines in top file

  top_out_content = []

  start_modifying_numbers = 'no'
  for i_line in range(0,n_top_lines):
    help_line = top_in_content[i_line]
    if('[ molecules ]' in help_line):     #modifications start once the program has read past this identifier
      start_modifying_numbers = 'yes'
    if(start_modifying_numbers == 'yes'):
      if('NH4' in help_line):
        help_line = ('NH4            ' + str(n_NH4_out) + "\n")
      if('NH3' in help_line):
        help_line = ('NH3            ' + str(n_NH3_out) + "\n")
      if('ACH' in help_line):
        help_line = ('ACH            ' + str(n_ACH_out) + "\n")
      if('ACO' in help_line):
        help_line = ('ACO            ' + str(n_ACO_out) + "\n")
      if('NA' in help_line):
        help_line = ('NA             ' + str(n_NA_out)  + "\n")
      if( ('CL' in help_line) and ('HCL' not in help_line) ):
        help_line = ('CL             ' + str(n_CL_out)  + "\n")
      if('HCL' in help_line):
        help_line = ('HCL            ' + str(n_HCL_out)  + "\n")
      if('H3O' in help_line):
        help_line = ('H3O            ' + str(n_H3O_out)  + "\n")
      if('SOL' in help_line):
        help_line = ('SOL            ' + str(n_SOL_out) + "\n")
    top_out_content.append(help_line)



  # End of function generate_top()

#------------------------------------------------------------------------------------------
# Start of function carve_gro()
# This function eliminates all molecules that are farther than r_cutoff away from center of mass (COM).
# The results are written to gro_out_content

def carve_gro():

  global gro_out_content     #this will turn into the output gro file
  global n_atom_out          #atoms in output file
  global n_SOL_in,n_SOL_out  #number of TIP4P-2005 water
  global n_ACO_in,n_ACO_out  #number of acetate 
  global n_ACH_in,n_ACH_out  #number of acetatic acid
  global n_NH3_in,n_NH3_out  #number of NH3
  global n_NH4_in,n_NH4_out  #number of NH4
  global n_NA_in ,n_NA_out   #number of NA
  global n_CL_in ,n_CL_out   #number of CL
  global n_HCL_in,n_HCL_out  #number of HCL
  global n_H3O_in,n_H3O_out  #number of HCL
  # Note that the _in numbers refer to AFTER PT


  gro_out_content = []

  allowed_molecules = [
                       'NH4','NH3','ACH','ACO','SOL','NA ','CL ','HCL','H3O',
                       'ACE','GLY','ALA','VAL','LEU','ILE','SER','THR','PRO','CYS','MET','ASN','ASP','GLN','GLU','LYS',
                       'ARG','HIS','PHE','TYR','TRP'
                      ]                   

  help_line = "gro file modified by droplet_carver_PT.py"  # header for output gro
  gro_out_content.append(help_line  + "\n")
  n_atom_out = 0                                           # number of atoms in output gro (will be updated later)
  help_line = str(n_atom_out)
  gro_out_content.append(help_line  + "\n")


  n_SOL_in  = 0
  n_SOL_out = 0
  n_ACO_in  = 0 
  n_ACO_out = 0
  n_ACH_in  = 0 
  n_ACH_out = 0
  n_NH3_in  = 0 
  n_NH3_out = 0
  n_NH4_in  = 0 
  n_NH4_out = 0
  n_NA_in   = 0 
  n_NA_out  = 0
  n_CL_in   = 0 
  n_CL_out  = 0
  n_HCL_in  = 0 
  n_HCL_out = 0
  n_H3O_in  = 0 
  n_H3O_out = 0

  do_not_include_in_output = [0,0,0,0,0,0,0,0]  #i_line numbers in this list will not be transferred to output

  for i_line in range(2,n_atom_in+2):
    help_line = gro_in_content[i_line]
    recognition_tag = help_line[5:15]
    molecule_type   = help_line[5:8]

    if(molecule_type not in allowed_molecules):
      print("  ### WARNING: molecule_type not defined: ",molecule_type)
      print("  ")

    if(recognition_tag == 'NA      NA'):   #---------------------- Carving NA
      n_NA_in = n_NA_in + 1

      line_messed_up = 'no'               #for molecules that are far out, gro lines may have faulty formatting. Those are eliminated
      if(help_line[20]  != ' '):
        line_messed_up = 'yes'
      if(help_line[28]  != ' '):
        line_messed_up = 'yes'
      if(help_line[36]  != ' '):
        line_messed_up = 'yes'

      if(line_messed_up == 'no'):         #line formatting is ok, and coordinates can be read
        x = float(help_line[20:28])
        y = float(help_line[28:36])
        z = float(help_line[36:44])
        dist_from_COM = np.sqrt((x-COM[0])**2 + (y-COM[1])**2 + (z-COM[2])**2)

      if(line_messed_up == 'yes'):        #lines with faulty formatting are assigned a distance that will cause elimination
        dist_from_COM = r_cutoff + 1.0

      if(dist_from_COM > r_cutoff):       #molecule is eliminated 
        do_not_include_in_output = [i_line,0,0,0,0,0,0,0]
      
      if(dist_from_COM <= r_cutoff):      #molecule is retained
        n_NA_out = n_NA_out + 1 

    if(recognition_tag == 'CL      CL'):   #---------------------- Carving CL
      n_CL_in = n_CL_in + 1

      line_messed_up = 'no'               #for molecules that are far out, gro lines may have faulty formatting. Those are eliminated
      if(help_line[20]  != ' '):
        line_messed_up = 'yes'
      if(help_line[28]  != ' '):
        line_messed_up = 'yes'
      if(help_line[36]  != ' '):
        line_messed_up = 'yes'

      if(line_messed_up == 'no'):         #line formatting is ok, and coordinates can be read
        x = float(help_line[20:28])
        y = float(help_line[28:36])
        z = float(help_line[36:44])
        dist_from_COM = np.sqrt((x-COM[0])**2 + (y-COM[1])**2 + (z-COM[2])**2)

      if(line_messed_up == 'yes'):        #lines with faulty formatting are assigned a distance that will cause elimination
        dist_from_COM = r_cutoff + 1.0

      if(dist_from_COM > r_cutoff):       #molecule is eliminated 
        do_not_include_in_output = [i_line,0,0,0,0,0,0,0]
      
      if(dist_from_COM <= r_cutoff):      #molecule is retained
        n_CL_out = n_CL_out + 1 

    if(recognition_tag == 'HCL    CL1'):   #---------------------- Carving HCL
      n_HCL_in = n_HCL_in + 1

      line_messed_up = 'no'               #for molecules that are far out, gro lines may have faulty formatting. Those are eliminated
      if(help_line[20]  != ' '):
        line_messed_up = 'yes'
      if(help_line[28]  != ' '):
        line_messed_up = 'yes'
      if(help_line[36]  != ' '):
        line_messed_up = 'yes'

      if(line_messed_up == 'no'):         #line formatting is ok, and coordinates can be read
        x = float(help_line[20:28])
        y = float(help_line[28:36])
        z = float(help_line[36:44])
        dist_from_COM = np.sqrt((x-COM[0])**2 + (y-COM[1])**2 + (z-COM[2])**2)

      if(line_messed_up == 'yes'):        #lines with faulty formatting are assigned a distance that will cause elimination
        dist_from_COM = r_cutoff + 1.0

      if(dist_from_COM > r_cutoff):       #molecule is eliminated 
        do_not_include_in_output = [i_line,i_line+1,0,0,0,0,0,0]
      
      if(dist_from_COM <= r_cutoff):      #molecule is retained
        n_HCL_out = n_HCL_out + 1 

    if(recognition_tag == 'H3O     OW'):   #---------------------- Carving H3O
      n_H3O_in = n_H3O_in + 1

      line_messed_up = 'no'               #for molecules that are far out, gro lines may have faulty formatting. Those are eliminated
      if(help_line[20]  != ' '):
        line_messed_up = 'yes'
      if(help_line[28]  != ' '):
        line_messed_up = 'yes'
      if(help_line[36]  != ' '):
        line_messed_up = 'yes'

      if(line_messed_up == 'no'):         #line formatting is ok, and coordinates can be read
        x = float(help_line[20:28])
        y = float(help_line[28:36])
        z = float(help_line[36:44])
        dist_from_COM = np.sqrt((x-COM[0])**2 + (y-COM[1])**2 + (z-COM[2])**2)

      if(line_messed_up == 'yes'):        #lines with faulty formatting are assigned a distance that will cause elimination
        dist_from_COM = r_cutoff + 1.0

      if(dist_from_COM > r_cutoff):       #molecule is eliminated 
        do_not_include_in_output = [i_line,i_line+1,i_line+2,i_line+3,i_line+4,0,0,0]
      
      if(dist_from_COM <= r_cutoff):      #molecule is retained
        n_H3O_out = n_H3O_out + 1 

    if(recognition_tag == 'SOL     OW'):   #---------------------- Carving TIP4P-2005 water (SOL)
      n_SOL_in = n_SOL_in + 1

      line_messed_up = 'no'               #for molecules that are far out, gro lines may have faulty formatting. Those are eliminated
      if(help_line[20]  != ' '):
        line_messed_up = 'yes'
      if(help_line[28]  != ' '):
        line_messed_up = 'yes'
      if(help_line[36]  != ' '):
        line_messed_up = 'yes'

      if(line_messed_up == 'no'):         #line formatting is ok, and coordinates can be read
        x = float(help_line[20:28])
        y = float(help_line[28:36])
        z = float(help_line[36:44])
        dist_from_COM = np.sqrt((x-COM[0])**2 + (y-COM[1])**2 + (z-COM[2])**2)

      if(line_messed_up == 'yes'):        #lines with faulty formatting are assigned a distance that will cause elimination
        dist_from_COM = r_cutoff + 1.0

      if(dist_from_COM > r_cutoff):       #molecule is eliminated 
        do_not_include_in_output = [i_line,i_line+1,i_line+2,i_line+3,0,0,0,0]
      
      if(dist_from_COM <= r_cutoff):      #molecule is retained
        n_SOL_out = n_SOL_out + 1

    if(recognition_tag == 'ACO     C1'):   #---------------------- Carving acetate (ACO)
      n_ACO_in = n_ACO_in + 1

      line_messed_up = 'no'               #for molecules that are far out, gro lines may have faulty formatting. Those are eliminated
      if(help_line[20]  != ' '):
        line_messed_up = 'yes'
      if(help_line[28]  != ' '):
        line_messed_up = 'yes'
      if(help_line[36]  != ' '):
        line_messed_up = 'yes'

      if(line_messed_up == 'no'):         #line formatting is ok, and coordinates can be read
        x = float(help_line[20:28])
        y = float(help_line[28:36])
        z = float(help_line[36:44])

      dist_from_COM = np.sqrt((x-COM[0])**2 + (y-COM[1])**2 + (z-COM[2])**2)

      if(line_messed_up == 'yes'):        #lines with faulty formatting are assigned a distance that will cause elimination
        dist_from_COM = r_cutoff + 1.0

      if(dist_from_COM > r_cutoff):
        do_not_include_in_output = [i_line,i_line+1,i_line+2,i_line+3,i_line+4,i_line+5,i_line+6,0] 
      if(dist_from_COM <= r_cutoff):
        n_ACO_out = n_ACO_out + 1

    if(recognition_tag == 'ACH     C2'):   #---------------------- Carving acetatic acid (ACH)
      n_ACH_in = n_ACH_in + 1

      line_messed_up = 'no'               #for molecules that are far out, gro lines may have faulty formatting. Those are eliminated
      if(help_line[20]  != ' '):
        line_messed_up = 'yes'
      if(help_line[28]  != ' '):
        line_messed_up = 'yes'
      if(help_line[36]  != ' '):
        line_messed_up = 'yes'

      if(line_messed_up == 'no'):         #line formatting is ok, and coordinates can be read
        x = float(help_line[20:28])
        y = float(help_line[28:36])
        z = float(help_line[36:44])

      dist_from_COM = np.sqrt((x-COM[0])**2 + (y-COM[1])**2 + (z-COM[2])**2)

      if(line_messed_up == 'yes'):        #lines with faulty formatting are assigned a distance that will cause elimination
        dist_from_COM = r_cutoff + 1.0

      if(dist_from_COM > r_cutoff):
        do_not_include_in_output = [i_line,i_line+1,i_line+2,i_line+3,i_line+4,i_line+5,i_line+6,i_line+7] 
      if(dist_from_COM <= r_cutoff):
        n_ACH_out = n_ACH_out + 1

    if(recognition_tag == 'NH3     N1'):   #---------------------- Carving ammonium (NH3)
      n_NH3_in = n_NH3_in + 1

      line_messed_up = 'no'               #for molecules that are far out, gro lines may have faulty formatting. Those are eliminated
      if(help_line[20]  != ' '):
        line_messed_up = 'yes'
      if(help_line[28]  != ' '):
        line_messed_up = 'yes'
      if(help_line[36]  != ' '):
        line_messed_up = 'yes'

      if(line_messed_up == 'no'):         #line formatting is ok, and coordinates can be read
        x = float(help_line[20:28])
        y = float(help_line[28:36])
        z = float(help_line[36:44])

      dist_from_COM = np.sqrt((x-COM[0])**2 + (y-COM[1])**2 + (z-COM[2])**2)

      if(line_messed_up == 'yes'):        #lines with faulty formatting are assigned a distance that will cause elimination
        dist_from_COM = r_cutoff + 1.0

      if(dist_from_COM > r_cutoff):
        do_not_include_in_output = [i_line,i_line+1,i_line+2,i_line+3,0,0,0,0]    
      if(dist_from_COM <= r_cutoff):
        n_NH3_out = n_NH3_out + 1

    if(recognition_tag == 'NH4    HZ1'):   #---------------------- Carving ammonia (NH4)
      n_NH4_in = n_NH4_in + 1

      line_messed_up = 'no'               #for molecules that are far out, gro lines may have faulty formatting. Those are eliminated
      if(help_line[20]  != ' '):
        line_messed_up = 'yes'
      if(help_line[28]  != ' '):
        line_messed_up = 'yes'
      if(help_line[36]  != ' '):
        line_messed_up = 'yes'

      if(line_messed_up == 'no'):         #line formatting is ok, and coordinates can be read
        x = float(help_line[20:28])
        y = float(help_line[28:36])
        z = float(help_line[36:44])

      dist_from_COM = np.sqrt((x-COM[0])**2 + (y-COM[1])**2 + (z-COM[2])**2)

      if(line_messed_up == 'yes'):        #lines with faulty formatting are assigned a distance that will cause elimination
        dist_from_COM = r_cutoff + 1.0

      if(dist_from_COM > r_cutoff):
        do_not_include_in_output = [i_line,i_line+1,i_line+2,i_line+3,i_line+4,0,0,0]
      if(dist_from_COM <= r_cutoff):
        n_NH4_out = n_NH4_out + 1
    

    if(i_line not in do_not_include_in_output):
      help_line = gro_in_content[i_line]
      gro_out_content.append(help_line)
      n_atom_out = n_atom_out + 1
 
  help_line = str(n_atom_out)         #updating number of atoms in output file
  gro_out_content[1] = (help_line  + "\n")

  help_line = gro_in_content[n_atom_in+2]
  gro_out_content.append(help_line)   #transferring box dimensions



  # End of function carve_gro()
#------------------------------------------------------------------------------------------
# Start of function calculate_COM()
# This function calculates the center of mass (COM), based on center mode.
# All atoms are treated the same (identical mass), regardless of element.

def calculate_COM():

  global COM

  COM        = [0.,0.,0.]
  COM_prelim = [0.,0.,0.]
  n_COM_atom = 0


  #----- START of PR section (centering based on protein)
  if(center_mode == "PR"):

    PR_residues = ['ACE','GLY','ALA','VAL','LEU','ILE','SER','THR','PRO','CYS','MET','ASN','ASP','GLN','GLU','LYS',
                   'ARG','HIS','PHE','TYR','TRP'] 
    #n_COM_atom = 0
    for i_line in range(2,n_atom_in+2):
      help_line = gro_in_content[i_line]
      molecule_type   = help_line[5:8]
      if(molecule_type in PR_residues):
        n_COM_atom = n_COM_atom + 1
        COM[0] = COM[0] + float(help_line[20:28]) #x of protein atom
        COM[1] = COM[1] + float(help_line[28:36]) #y of protein atom  
        COM[2] = COM[2] + float(help_line[36:44]) #z of protein atom

    COM[0] = COM[0] / n_COM_atom
    COM[1] = COM[1] / n_COM_atom 
    COM[2] = COM[2] / n_COM_atom
  #----- END of PR section


  #----- START of OW section (centering based on TIP4P-2005 water)
  if(center_mode == "OW"):

    print("  ")
    print("  ### WARNING: center_mode OW is OK for initial droplet carving.                        ")
    print("               But: This option may not be suitable for non-spherical systems, such as  ")
    print("               gro files containing a parent droplet and one or more progeny droplets.  ")
    print("  ")
    #n_COM_atom = 0
    for i_line in range(2,n_atom_in+2):
      help_line = gro_in_content[i_line]
      recognition_tag = help_line[5:15]
      if(recognition_tag == "SOL     OW"):

        line_messed_up = 'no'               # For molecules that are far out, gro lines may have faulty formatting
        if(help_line[20]  != ' '):          # and they are note considered for COM calculations
          line_messed_up = 'yes'
        if(help_line[28]  != ' '):
          line_messed_up = 'yes'
        if(help_line[36]  != ' '):
          line_messed_up = 'yes'

        if(line_messed_up == 'no'):
          n_COM_atom = n_COM_atom + 1
          COM[0] = COM[0] + float(help_line[20:28]) #x of OW atom
          COM[1] = COM[1] + float(help_line[28:36]) #y of OW atom  
          COM[2] = COM[2] + float(help_line[36:44]) #z of OW atom

    COM[0] = COM[0] / n_COM_atom
    COM[1] = COM[1] / n_COM_atom 
    COM[2] = COM[2] / n_COM_atom
  
    #print("  COM(OW) = ",COM)
  
  #----- END of OW section

  #----- START of NH4 section (centering based on NH4)
  # This code takes care of a potential problem: when one (or more) NH4 have been ejected and are far from the
  # droplet, the carving procedure might discard too many (or even all) molecules becuase the COM is somewhere
  # in empty space . This code makes sure that such outliers do not skew the COM location.
  if(center_mode == "NH4"):

    NH4_info = []                           # list of all NH4 NZ coordinates
    n_NH4    = 0                            # this counts the NH4 that have valid coordinates (no messed up lines)
    dist_min = 1000.      

    # First we calculate a "preliminary" COM, that is based on all NH4
    n_COM_atom_prelim = 0
    for i_line in range(2,n_atom_in+2):
      help_line = gro_in_content[i_line]
      recognition_tag = help_line[5:15]
      if(recognition_tag == "NH4     NZ"):

        line_messed_up = 'no'               # For molecules that are far out, gro lines may have faulty formatting
        if(help_line[20]  != ' '):          # and they are note considered for COM calculations
          line_messed_up = 'yes'
        if(help_line[28]  != ' '):
          line_messed_up = 'yes'
        if(help_line[36]  != ' '):
          line_messed_up = 'yes'

        if(line_messed_up == 'no'):
          n_COM_atom_prelim = n_COM_atom_prelim + 1
          n_NH4 = n_NH4 + 1
          NH4_info.append(help_line)
          COM_prelim[0] = COM_prelim[0] + float(help_line[20:28]) #x of NZ atom
          COM_prelim[1] = COM_prelim[1] + float(help_line[28:36]) #y of NZ atom  
          COM_prelim[2] = COM_prelim[2] + float(help_line[36:44]) #z of NZ atom


    if(n_COM_atom_prelim == 0):
      print("  ")
      print("  ### ERROR: No NH4 found. Cannot calculate COM.")
      print("  ")
      sys.exit() 

    COM_prelim[0] = COM_prelim[0] / n_COM_atom_prelim
    COM_prelim[1] = COM_prelim[1] / n_COM_atom_prelim 
    COM_prelim[2] = COM_prelim[2] / n_COM_atom_prelim

    # Now we find the coordinates of the NH4 that is closest to COM_prelim.
    # This closest NH4 will be located withint the largest droplet.

    for i_NH4 in range(0,n_NH4):
      help_line = NH4_info[i_NH4]
      x = float(help_line[20:28])
      y = float(help_line[28:36])
      z = float(help_line[36:44])
      dist_from_prelim_COM = np.sqrt((x-COM_prelim[0])**2 + (y-COM_prelim[1])**2 + (z-COM_prelim[2])**2)
      if(dist_from_prelim_COM < dist_min):
        closest_NH4_xyz = [x,y,z]
        dist_min        = dist_from_prelim_COM

    # Now calculate COM, based on only those NH4 that are in the vicinity of this closest NH4*.
    # This will eliminate any NH4 that are "far out"
    for i_NH4 in range(0,n_NH4):
      help_line = NH4_info[i_NH4]
      x = float(help_line[20:28])
      y = float(help_line[28:36])
      z = float(help_line[36:44])
      dist_from_closest_NH4 = np.sqrt((x-closest_NH4_xyz[0])**2 + (y-closest_NH4_xyz[1])**2 + (z-closest_NH4_xyz[2])**2)

      if(dist_from_closest_NH4 < r_cutoff):
        n_COM_atom = n_COM_atom + 1
        COM[0] = COM[0] + x
        COM[1] = COM[1] + y 
        COM[2] = COM[2] + z

    COM[0] = COM[0] / n_COM_atom
    COM[1] = COM[1] / n_COM_atom 
    COM[2] = COM[2] / n_COM_atom

  #----- END of NH4 section

  #----- START of NA section (centering based on NA)
  # This code takes care of a potential problem: when one (or more) NA have been ejected and are far from the
  # droplet, the carving procedure might discard too many (or even all) molecules becuase the COM is somewhere
  # in empty space . This code makes sure that such outliers do not skew the COM location.
  if(center_mode == "NA"):

    NA_info = []                           # list of all NA coordinates
    n_NA    = 0                            # this counts the NA that have valid coordinates (no messed up lines)
    dist_min = 1000.      

    # First we calculate a "preliminary" COM, that is based on all NA
    n_COM_atom_prelim = 0
    for i_line in range(2,n_atom_in+2):
      help_line = gro_in_content[i_line]
      recognition_tag = help_line[5:15]
      if(recognition_tag == "NA      NA"):

        line_messed_up = 'no'               # For molecules that are far out, gro lines may have faulty formatting
        if(help_line[20]  != ' '):          # and they are note considered for COM calculations
          line_messed_up = 'yes'
        if(help_line[28]  != ' '):
          line_messed_up = 'yes'
        if(help_line[36]  != ' '):
          line_messed_up = 'yes'

        if(line_messed_up == 'no'):
          n_COM_atom_prelim = n_COM_atom_prelim + 1
          n_NA = n_NA + 1
          NA_info.append(help_line)
          COM_prelim[0] = COM_prelim[0] + float(help_line[20:28]) #x of NA atom
          COM_prelim[1] = COM_prelim[1] + float(help_line[28:36]) #y of NA atom  
          COM_prelim[2] = COM_prelim[2] + float(help_line[36:44]) #z of NA atom


    if(n_COM_atom_prelim == 0):
      print("  ")
      print("  ### ERROR: No NA found. Cannot calculate COM.")
      print("  ")
      sys.exit() 

    COM_prelim[0] = COM_prelim[0] / n_COM_atom_prelim
    COM_prelim[1] = COM_prelim[1] / n_COM_atom_prelim 
    COM_prelim[2] = COM_prelim[2] / n_COM_atom_prelim

    # Now we find the coordinates of the NA that is closest to COM_prelim.
    # This closest NA will be located withint the largest droplet.

    for i_NA in range(0,n_NA):
      help_line = NA_info[i_NA]
      x = float(help_line[20:28])
      y = float(help_line[28:36])
      z = float(help_line[36:44])
      dist_from_prelim_COM = np.sqrt((x-COM_prelim[0])**2 + (y-COM_prelim[1])**2 + (z-COM_prelim[2])**2)
      if(dist_from_prelim_COM < dist_min):
        closest_NA_xyz = [x,y,z]
        dist_min        = dist_from_prelim_COM

    # Now calculate COM, based on only those NA that are in the vicinity of this closest NA*.
    # This will eliminate any NA that are "far out"
    for i_NA in range(0,n_NA):
      help_line = NA_info[i_NA]
      x = float(help_line[20:28])
      y = float(help_line[28:36])
      z = float(help_line[36:44])
      dist_from_closest_NA = np.sqrt((x-closest_NA_xyz[0])**2 + (y-closest_NA_xyz[1])**2 + (z-closest_NA_xyz[2])**2)

      if(dist_from_closest_NA < r_cutoff):
        n_COM_atom = n_COM_atom + 1
        COM[0] = COM[0] + x
        COM[1] = COM[1] + y 
        COM[2] = COM[2] + z

    COM[0] = COM[0] / n_COM_atom
    COM[1] = COM[1] / n_COM_atom 
    COM[2] = COM[2] / n_COM_atom

  #----- END of NA section

  #----- START of H3O section (centering based on H3O)
  # This code takes care of a potential problem: when one (or more) H3O have been ejected and are far from the
  # droplet, the carving procedure might discard too many (or even all) molecules becuase the COM is somewhere
  # in empty space . This code makes sure that such outliers do not skew the COM location.
  if(center_mode == "H3O"):

    H3O_info = []                           # list of all H3O coordinates
    n_H3O    = 0                            # this counts the H3O that have valid coordinates (no messed up lines)
    dist_min = 1000.      

    # First we calculate a "preliminary" COM, that is based on all H3O
    n_COM_atom_prelim = 0
    for i_line in range(2,n_atom_in+2):
      help_line = gro_in_content[i_line]
      recognition_tag = help_line[5:15]
      if(recognition_tag == "H3O     OW"):

        line_messed_up = 'no'               # For molecules that are far out, gro lines may have faulty formatting
        if(help_line[20]  != ' '):          # and they are note considered for COM calculations
          line_messed_up = 'yes'
        if(help_line[28]  != ' '):
          line_messed_up = 'yes'
        if(help_line[36]  != ' '):
          line_messed_up = 'yes'

        if(line_messed_up == 'no'):
          n_COM_atom_prelim = n_COM_atom_prelim + 1
          n_H3O = n_H3O + 1
          H3O_info.append(help_line)
          COM_prelim[0] = COM_prelim[0] + float(help_line[20:28]) #x of H3O
          COM_prelim[1] = COM_prelim[1] + float(help_line[28:36]) #y of H3O  
          COM_prelim[2] = COM_prelim[2] + float(help_line[36:44]) #z of H3O


    if(n_COM_atom_prelim == 0):
      print("  ")
      print("  ### ERROR: No H3O found. Cannot calculate COM.")
      print("  ")
      sys.exit() 

    COM_prelim[0] = COM_prelim[0] / n_COM_atom_prelim
    COM_prelim[1] = COM_prelim[1] / n_COM_atom_prelim 
    COM_prelim[2] = COM_prelim[2] / n_COM_atom_prelim

    # Now we find the coordinates of the H3O that is closest to COM_prelim.
    # This closest H3O will be located withint the largest droplet.

    for i_H3O in range(0,n_H3O):
      help_line = H3O_info[i_H3O]
      x = float(help_line[20:28])
      y = float(help_line[28:36])
      z = float(help_line[36:44])
      dist_from_prelim_COM = np.sqrt((x-COM_prelim[0])**2 + (y-COM_prelim[1])**2 + (z-COM_prelim[2])**2)
      if(dist_from_prelim_COM < dist_min):
        closest_H3O_xyz = [x,y,z]
        dist_min        = dist_from_prelim_COM

    # Now calculate COM, based on only those H3O that are in the vicinity of this closest H3O*.
    # This will eliminate any H3O that are "far out"
    for i_H3O in range(0,n_H3O):
      help_line = H3O_info[i_H3O]
      x = float(help_line[20:28])
      y = float(help_line[28:36])
      z = float(help_line[36:44])
      dist_from_closest_H3O = np.sqrt((x-closest_H3O_xyz[0])**2 + (y-closest_H3O_xyz[1])**2 + (z-closest_H3O_xyz[2])**2)

      if(dist_from_closest_H3O < r_cutoff):
        n_COM_atom = n_COM_atom + 1
        COM[0] = COM[0] + x
        COM[1] = COM[1] + y 
        COM[2] = COM[2] + z

    COM[0] = COM[0] / n_COM_atom
    COM[1] = COM[1] / n_COM_atom 
    COM[2] = COM[2] / n_COM_atom

  #----- END of H3O section




  # End of function calculate_COM()
#------------------------------------------------------------------------------------------
# start of main program

#global gro_in_content
global top_in_content
global center_mode
global n_atom_in
global log_name
global PT, PT_dist_max,PT_prob



print("                              ")
print("  ############################")
print("  #                          #")
print("  #   droplet_carver_PT.py   #")
print("  #                          #")
print("  ############################")
print("                         ")

n_command_line = len(sys.argv)
if(n_command_line != 11):
  print("  ")
  print("  ### ERROR: Command_line not ok.")
  print("             Expected 11 arguments, found ",n_command_line)
  print("  ")
  print(sys.argv)
  sys.exit()
  
# assign file names etc. from command line info
gro_in_name  = sys.argv[1]          # input  gro file
gro_out_name = sys.argv[2]          # output gro file
top_in_name  = sys.argv[3]          # input  top file
top_out_name = sys.argv[4]          # output top file
r_cutoff     = float(sys.argv[5])   # carving radius; all molecules beyond this radius will be eliminated
center_mode  = sys.argv[6]          # determines what molecules are used for centering
log_name     = sys.argv[7]          # name of output log file
PT           = sys.argv[8]          # determines if NH4 --> ACO proton transfer takes place

if(PT == 'PT=1'):                   # This will activate PT from NH4 to ACO
  PT = '1'
  PT_dist_max = float(sys.argv[9])  # maximum H-O distance (nm) for proton transfer
  PT_prob     = float(sys.argv[10]) # probability of PT, assuming that distance < PT_dist_max

if(PT == 'PT=2'):                   # This will activate PT from H3O to CL
  PT = '2'
  PT_dist_max = float(sys.argv[9])  # maximum H-CL distance (nm) for proton transfer
  PT_prob     = float(sys.argv[10]) # probability of PT, assuming that distance < PT_dist_max

print("  input gro file       : ",gro_in_name)
print("  input top file       : ",top_in_name)
print("  log       file       : ",log_name)
print("  cutoff radius (nm)   : ",r_cutoff)
print("  center mode          : ",center_mode)

allowed_center_mode = ['PR', 'OW', 'NH4','NA','H3O']
if( (center_mode in allowed_center_mode) == False):
  print("  ")
  print("  ### ERROR: center_mode not defined: ",center_mode)
  print("  ")
  sys.exit()
if(center_mode == "PR"):
  print("  (output will be centered around protein)")
if(center_mode == "OW"):
  print("  (output will be centered around water)")

path.exists(gro_in_name)
if (os.path.isfile(gro_in_name) == False):
  print("  ")
  print("  ### ERROR: file not found: ",gro_in_name)
  print("  ")
  sys.exit()

top_present = 'yes'
path.exists(top_in_name)  #this option allows gro carving without a top file 
if (os.path.isfile(top_in_name) == False):
  print("  ")
  print("  ### WARNING: input top file not found. Program will process gro file only.")
  print("  ")
  top_present = 'no'

# read all lines of gro_in
gro_in    = open(gro_in_name,"r")         
gro_in_content = gro_in.readlines()  
gro_in.close()

if(top_present == 'yes'):
  # read all lines of top_in
  top_in    = open(top_in_name,"r")         
  top_in_content = top_in.readlines()  
  top_in.close()                    

# determine number of atoms
n_atom_in = int(gro_in_content[1])

# The subsequent proton transfer events can redefine the composition of gro_in_content
# For this reason we now count all of the molecules (_pre_PT)
count_molecules_pre_PT()

if(PT == 'PT=no'):
  PT = 'no'
  print("  PT                   :  no proton transfer")
  
if( (PT != 'no') and (PT != '1') and (PT != '2')):
  print("  ")
  print("  ### ERROR: PT=no or PT=1 or PT=2 has to be specified in command line")
  print("  ")
  sys.exit()

if(PT == '1'):
  print("  PT option            :  1 (NH4 --> ACO proton transfer will be executed)")
  print("  PT_dist_max (nm)     : ",PT_dist_max)
  print("  PT_prob              : ",PT_prob)

if(PT == '2'):
  print("  PT option            :  2 (H3O --> CL proton transfer will be executed)")
  print("  PT_dist_max (nm)     : ",PT_dist_max)
  print("  PT_prob              : ",PT_prob)
  if(PT_prob < 0):
    print("    negative PT_prob means that H3O --> CL PT probability = |PT_prob| * n_H3O/(n_H3O + n_SOL)")
    denominator = n_H3O_in_pre_PT + n_SOL_in_pre_PT
    if (denominator > 0):
      PT_prob = abs(PT_prob) * n_H3O_in_pre_PT / denominator
      if (PT_prob > 1.0):
        PT_prob = 1.0
    if (denominator == 0):
      PT_prob = 0
    print("    current PT_prob    : ", format(PT_prob, '1.4f'))


# execute NH4 --> ACO proton transfer
n_PT = 0    #number of proton transfer events
if(PT == '1'):
  execute_PT_NH4_ACO()

# execute H3O --> CL proton transfer
n_PT = 0    #number of proton transfer events
if(PT == '2'):
  execute_PT_H3O_CL()

# calculate center of mass
calculate_COM()

# carve gro file
carve_gro()

# generate updated topology content after carving (and proton transfer, if applicable)
if(top_present == 'yes'):
  generate_top()

# write log file
write_log()

print("  ")
print("  Input file composition")
print("  Atoms ",n_atom_in)
if(n_NH4_in > 0):
  print("  NH4   ",n_NH4_in_pre_PT)
if(n_NH3_in > 0):
 print("  NH3   ",n_NH3_in_pre_PT)
if(n_ACH_in > 0):
  print("  ACH   ",n_ACH_in_pre_PT)
if(n_ACO_in > 0):
  print("  ACO   ",n_ACO_in_pre_PT)
if(n_NA_in > 0):
  print("  NA    ",n_NA_in_pre_PT)
if(n_CL_in > 0):
  print("  CL    ",n_CL_in_pre_PT)
if(n_HCL_in > 0):
  print("  HCL   ",n_HCL_in_pre_PT)
if(n_H3O_in > 0):
  print("  H3O   ",n_H3O_in_pre_PT)
if(n_SOL_in > 0):
  print("  SOL   ",n_SOL_in_pre_PT)
print("  ")

print("  Output composition (after",n_PT," PT events)")
print("  Atoms ",n_atom_out)
if(n_NH4_out > 0):
  print("  NH4   ",n_NH4_out)
if(n_NH3_out > 0):
  print("  NH3   ",n_NH3_out)
if(n_ACH_out > 0):
  print("  ACH   ",n_ACH_out)
if(n_ACO_out > 0):
  print("  ACO   ",n_ACO_out)
if(n_NA_out > 0):
  print("  NA    ",n_NA_out)
if(n_CL_out > 0):
  print("  CL    ",n_CL_out)
if(n_HCL_out > 0):
  print("  HCL   ",n_HCL_out)
if(n_H3O_out > 0):
  print("  H3O   ",n_H3O_out)
if(n_SOL_out > 0):
  print("  SOL   ",n_SOL_out)
print("  ")

# write gro output 
gro_out = open(gro_out_name,"w+")
gro_out.writelines(gro_out_content)
gro_out.close()
print("  gro output generated : ",gro_out_name)

# write top output 
if(top_present == 'yes'):
  top_out = open(top_out_name,"w+")
  top_out.writelines(top_out_content)
  top_out.close()
  print("  top output generated : ",top_out_name)

print("  log file             : ",log_name)

print("  ")
# End of Main Program
