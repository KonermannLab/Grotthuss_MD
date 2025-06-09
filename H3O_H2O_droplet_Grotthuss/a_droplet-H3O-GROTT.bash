#! /bin/bash
# 
# Bash script for simulating water droplets with H3O+, with Grotthuss diffusion
# written by Lars Konermann, konerman@uwo.ca
############################################################################################################
# When using this file, please cite                                                                        #
# "Grotthuss Molecular Dynamics Simulations for Modeling Proton Hopping in Electrosprayed Water Droplets"  #
# L. Konermann and S. Kim, J. Chem. Theory Comput. 18, 3781–3794 (2022)                                    #
############################################################################################################
# This bash script runs droplet simulations in a vacuum environment by using a Python program in combination 
# with Gromacs (e.g. versions 2019.3 and 2020.4).
# Some newer Gromacs version may have issues related to atom numbering in gro files 
# (all atoms have to be numbered consecutively).  
# Droplet dynamics are performed by alternating Gromacs with
# H3O_proton_hop.py (which performs H3O --> H2O Grotthuss hopping)
#
# There is no removal of evaporated molecules in this version.
# 
# For reproducing the experimental proton diffusion coefficient at 300 K, make sure that
# 1. The md_mdp... files have a 4 ps run time
# 2. n_hop=5
# 3. hop_cutoff=0.25
#
# Run this script with $ bash <filename>.bash
#


md_mdp_gen_vel_no='a_300K_4ps_gen_vel_no.mdp'          #mdp file for production run steps, with gen-vel = no
grott_log='x_LOG_grott.txt'                            #output file for H3O_Grotthuss.py (keeps track of hopping events for each H+)
n_steps_H3O_Grotthuss=10                               #number of times that H3O_Grotthuss program is called in between droplet_cleanup steps
n_hop=5                                                #number of times that each H+ can hop during each call of H3O_Grotthuss program
hop_cutoff=0.25                                        #maximum allowed O...H distance (nm) for H+ hopping

last_hop=$((n_hop-1))
Last_hop=$(printf %01d ${last_hop})  #in Python numbering the last file will have number (n_hop-1)

cp a_start.gro 00000_000.gro
cp a_start.top 00000.top

for k in $( eval echo {1..$n_steps_H3O_Grotthuss} )  #inner loop performs Grotthuss hopping 
do
  km1=$((k-1)) 
  K=$(printf %03d ${k})     #( loop step number  )    is converted to i5 text format
  Km1=$(printf %03d ${km1}) #( loop step number-1)    is converted to i5 text format
  echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@ Grotthuss step k = "${K}
  python H3O_Grotthuss.py 00000_${Km1}.gro ${n_hop} ${hop_cutoff} out-all ${grott_log} 
  gmx grompp -f ${md_mdp_gen_vel_no} -c 00000_${Km1}_grott${Last_hop}.gro -p 00000.top -o 00000_${K}.tpr
  gmx mdrun -deffnm 00000_${K} -gpu_id 0

done

echo "========================== bash DONE! ============================="







