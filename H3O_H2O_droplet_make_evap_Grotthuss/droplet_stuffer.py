#!/usr/bin/python
#
# droplet_stuffer.py
# Written by: Lars Konermann, January 2022
# konerman@uwo.ca
#
# Program to "stuff" an empty droplet with certain molecules (gro files).
# The output is a single gro files that contains all of the molecules.
# Subsequently, water can be added using gmx solvate, followed by droplet carving.
# The individual molecules are placed on grid points inside the droplet.
# The input gro files do not have to be centered.
# 
# calling program from Python directory:  $ python droplet_stuffer.py
# calling program from another directory: $ python ~/Python/droplet_stuffer.py <input_file.txt>
# 
# The program expects an input file 
# File structure:
# First line for comments
# name of output gro file
# radius:  droplet radius (nm)
# gro_mol[0]: gro file that will be placed in the center (n_mol[0] = 1)
# n_type_addtl_mol: number that specifies how many ADDITIONAL types of molecules are to be inserted
# n_mol[1]
# gro_mol[1]: number of molecules for molecule 1 and gro file
# n_mol[2]
# gro_mol[2]: number of molecules for molecule 2 and gro file
# ...

import sys
import os.path
from os import path
import numpy as np
import random

print("                           ")
print("  #########################")
print("  #                       #")
print("  #  droplet_stuffer.py   #")
print("  #                       #")
print("  #########################")
print("                           ")

input_parameter_file  = sys.argv[1]        # input  gro file

path.exists(input_parameter_file)
if (os.path.isfile(input_parameter_file) == False):
  print("  ### ERROR: input file not found: ", input_parameter_file) 
  sys.exit()

# Assign input information
# ========================
#

input_parameter         = open(input_parameter_file,"r")     #open  input file
input_parameter_content = input_parameter.readlines()        #read  input file
input_parameter.close()                                      #close input file

help_line = input_parameter_content[1]
output_file = help_line.replace("\n","")

input_parameter.close()

n_mol  = []  #list that contains the number of molecules for each type of molecule
gro_mol= []  #list that contains gro file names for each type of molecule

help_line = input_parameter_content[2]
droplet_radius = float(help_line)

n_mol.append(1) #reading information for single molecule that will be placed in center
help_line = input_parameter_content[3].replace("\n","")
gro_mol.append(help_line)
print("  molecule in center        ",gro_mol[0])

help_line = input_parameter_content[4]
n_type_addtl_mol = int(help_line)

print("  droplet radius (nm)       ",droplet_radius)
print("  types of addtl molecules  ",n_type_addtl_mol)


# reading input gro files
i_type_mol  = 1
n_mol_total = 1
for i_line in range( 5, 2*(4+n_type_addtl_mol)-3, 2 ):
  print("  ") 
  help_line = input_parameter_content[i_line]
  n_mol.append( int(help_line) )
  n_mol_total = n_mol_total + n_mol[i_type_mol]  # overall number of molecules that have to be inserted
  print("  number of molecules       ",n_mol[i_type_mol])
  help_line = input_parameter_content[i_line+1].replace("\n","")
  gro_mol.append(help_line)
  print("                            ",gro_mol[i_type_mol].replace("\n",""))
  i_type_mol = i_type_mol+1


# Generating x/y/z coordinates of grid points in droplet where molecules can be placed.
# ======================================================================================
#
# We start by arranging grid points in a cube (length = 2*droplet_radius)), but only the
# points inside the droplet get populated.
# The cube has 11 x 11 x 11 points
# The droplet center is at x = y = z = droplet_radius.
# The array grid_point_in_droplet stores these x /y /z coordinates:
#                                      i_point xi/yi/zi
grid_point_in_droplet = np.empty(shape=(1000,3))

i_point = -1

for i_x in range(0,11):
  for i_y in range(0,11):
    for i_z in range(0,11):
      x = i_x * 2. * droplet_radius/10.
      y = i_y * 2. * droplet_radius/10.
      z = i_z * 2. * droplet_radius/10.
      dist_from_center = np.sqrt((droplet_radius-x)**2 + (droplet_radius-y)**2 + (droplet_radius-z)**2)
      if (dist_from_center < droplet_radius):
        i_point = i_point + 1
        grid_point_in_droplet[i_point,0] = x
        grid_point_in_droplet[i_point,1] = y
        grid_point_in_droplet[i_point,2] = z
        
n_points_in_sphere = i_point - 1

# The very first point is now redifined and moved to the droplet center.
# This is where the first molecule will be placed.
grid_point_in_droplet[0,0] = droplet_radius
grid_point_in_droplet[0,1] = droplet_radius
grid_point_in_droplet[0,2] = droplet_radius

print("  ")
print("  number of grid points    ",n_points_in_sphere)
print("  total number of molecules",n_mol_total)

if (n_points_in_sphere < n_mol_total):
  print("  ### ERROR: not enough grid points available") 
  sys.exit()



# Starting to write to output file
# ================================
#

output         = open(output_file,"w+")     #open  input file
output_content = []
help_line = "droplet .gro file generated by droplet_stuffer.py"
output_content.append(help_line + "\n")
help_line = "0"
output_content.append(help_line + "\n")



# Reading input gro files, and writing multiples of those gro files into gro output file
# ======================================================================================
#

n_type_mol = n_type_addtl_mol + 1

n_atoms_out = 0
i_atoms_out = 1
i_mol_total = 0
collection_of_used_points = []
collection_of_used_points.append(0)


for i_type_mol in range(0,n_type_mol):
  path.exists(gro_mol[i_type_mol])
  if (os.path.isfile(gro_mol[i_type_mol]) == False):
    print("  ### ERROR: file not found: ",gro_mol[i_type_mol])
    sys.exit()


  gro_data         = open(gro_mol[i_type_mol],"r")     #open  gro file
  gro_data_content = gro_data.readlines()              #read  gro file
  gro_data.close()                                     #close gro file
  n_atom_in_gro = int(gro_data_content[1])  
  
  x_COM_curr_gro = 0.
  y_COM_curr_gro = 0.
  z_COM_curr_gro = 0.

  for i_atom_in_gro in range(0,n_atom_in_gro):         #input gro files will likely not be centered, so we calculate their COM to compensate for this.
    help_line = gro_data_content[i_atom_in_gro+2]      # (the COM calculation treats all atoms equally, regardless of their mass) 
    x = float(help_line[20:28])
    y = float(help_line[28:36])
    z = float(help_line[36:44])
    x_COM_curr_gro = x_COM_curr_gro + x
    y_COM_curr_gro = y_COM_curr_gro + y
    z_COM_curr_gro = z_COM_curr_gro + z  

  x_COM_curr_gro = x_COM_curr_gro/n_atom_in_gro
  y_COM_curr_gro = y_COM_curr_gro/n_atom_in_gro
  z_COM_curr_gro = z_COM_curr_gro/n_atom_in_gro    
   

  i_mol_success = 0                                   #number of molecules (of the current type) that have already been successfully placed in droplet
#  for i_mol in range(0,n_mol[i_type_mol]):           
  while (i_mol_success < n_mol[i_type_mol]):          #now we are placing molecules, one by one, on grid points

    if (i_type_mol == 0):                             #the first molecule gets special treatment and is placed in the droplet center
      i_point = 0
                           
    if (i_type_mol > 0):                              #all other molecules are placed on randomly selectred grid points
      while (i_point in collection_of_used_points):   #this while loop ensures that the new random point has not already been used
        i_point = random.randint(0,n_points_in_sphere)
    
    collection_of_used_points.append(i_point)         #the current grid point is added to the list, so that it will not be used again

    if (len(collection_of_used_points) > (n_points_in_sphere-1)):
      print("  ### ERROR: Program ran out of clash-free grid points during placement of molecules.")
      print("             Total number of successfully placed molecules = "),i_mol_total
      sys.exit()

    steric_clash = "no"

    clash_cutoff = 0.15 #(nm)                         #now we test if any atoms of the current gro clash with atoms that have already been placed
    for i_atom_in_gro in range(0,n_atom_in_gro):  
      help_line = gro_data_content[i_atom_in_gro+2]   
      x_test = float(help_line[20:28]) + grid_point_in_droplet[i_point,0] - x_COM_curr_gro  #coordinates of current atoms
      y_test = float(help_line[28:36]) + grid_point_in_droplet[i_point,1] - y_COM_curr_gro
      z_test = float(help_line[36:44]) + grid_point_in_droplet[i_point,2] - z_COM_curr_gro
      for j_line in range(2,n_atoms_out):
        help_line = output_content[j_line]      
        x_already_there = float(help_line[20:28])     #coordinates of atoms that have already been placed in output 
        y_already_there = float(help_line[28:36])
        z_already_there = float(help_line[36:44])
        dist = np.sqrt((x_already_there-x_test)**2 + (y_already_there-y_test)**2 + (z_already_there-z_test)**2)
        if (dist < clash_cutoff):
          steric_clash = "yes"                        # This means that there is a clash, so the current molecule will not be placed at this grid point.
                                                      # The grid point is discarded for future use by other molecules.
    if (steric_clash == "no"):                        # Only when there is no clash, will the gro be placed on the grid point and written to output.
      i_mol_success = i_mol_success + 1               # Otherwise (if there is a clash) we try another iteration of the while loop.
      i_mol_total = i_mol_total + 1
      n_atoms_out = n_atoms_out + n_atom_in_gro
      for i_atom_in_gro in range(0,n_atom_in_gro):
        help_line = gro_data_content[i_atom_in_gro+2]

        x = float(help_line[20:28]) + grid_point_in_droplet[i_point,0] - x_COM_curr_gro
        y = float(help_line[28:36]) + grid_point_in_droplet[i_point,1] - y_COM_curr_gro
        z = float(help_line[36:44]) + grid_point_in_droplet[i_point,2] - z_COM_curr_gro

        if (i_mol_total == 1):      #we need these two cases to retain protein residue numbering (this is for protein)
          output_line = (help_line[0:15] + format(i_atoms_out, "5") +
          format(x, "8.3f")                                         + 
          format(y, "8.3f")                                         + 
          format(z, "8.3f")                                         + "\n")
          output_content.append(output_line)
          i_atoms_out = i_atoms_out + 1

          last_residue_number = int(help_line[0:5])

        if (i_mol_total > 1):       #we need these two cases to retain protein residue numbering (this is for non-protein molecules)
          output_line = (format(i_mol_total+last_residue_number-1, "5") +
          help_line[5:15] + format(i_atoms_out, "5")                    + 
          format(x, "8.3f")                                             + 
          format(y, "8.3f")                                             + 
          format(z, "8.3f")                                             + "\n")
          output_content.append(output_line)
          i_atoms_out = i_atoms_out + 1

help_line = (format(n_atoms_out, "6") + "\n")  #update second line of output (number of atoms)
output_content[1] = help_line

help_line = ( format(2*droplet_radius, "12.4f") + 
              format(2*droplet_radius, "12.4f") + 
              format(2*droplet_radius, "12.4f") + "\n" )
output_content.append(help_line)              #append last line of output (box dimensions)

# Write everything to output file
output.writelines(output_content)
output.close()
print("  ")
print("  output file generated    ",output_file)
print("  ")




# End of Main Program
