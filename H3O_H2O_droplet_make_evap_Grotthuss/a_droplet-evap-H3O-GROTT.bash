#! /bin/bash
# 
# Bash script for simulating water droplets with H3O+, with Grotthuss diffusion, and with removal of evaporated molecules
# written by Lars Konermann, konerman@uwo.ca
############################################################################################################
# When using this file, please cite                                                                        #
# "Grotthuss Molecular Dynamics Simulations for Modeling Proton Hopping in Electrosprayed Water Droplets"  #
# L. Konermann and S. Kim, J. Chem. Theory Comput. 18, 3781–3794 (2022)                                    #
############################################################################################################
# This bash script runs droplet simulations in a vacuum environment by using 3 Python programs in combination 
# with Gromacs (e.g. versions 2019.3 and 2020.4).
# Some newer Gromacs version may have issues related to atom numbering in gro files 
# (all atoms have to be numbered consecutively).  
# Initially, droplet_stuffer.py generates the droplet, for subsequent energy minimization and equilibration.
# Droplet dynamics are then performed by alternating 
# 1. droplet_carver_PT.py (which eliminates evaporated molecules, and
# 2. H3O_proton_hop.py (which performs H3O --> H2O Grotthuss hopping)
# 
# For reproducing the experimental proton diffusion coefficient at 300 K, make sure that
# 1. The md_mdp... files have a 4 ps run time
# 2. n_hop=5
# 3. hop_cutoff=0.25
#
# Run this script with $ bash <filename>.bash
#


droplet_stuffer_input="a_droplet_stuffer_inp.txt"      #input file determines solute types and numbers in droplet 
raw_top="a_H3O_H2O_raw.top"                            #generic topology with unspecified [ molecules ] numbers
n_steps_carve=4                                        #number of carving iterations (removal of evaporated molecules)
carv_r_ini=1.2                                         #initial carving radius (has to be consistent with droplet_stuffer_input
carv_r_cleanup=5.                                      #radius for droplet cleanup operations, should be roughly 4 * carv_r_ini
em_mdp="a_em.mdp"                                      #mdp for energy minimization
equil_mdp="a_10K_300K_200ps.mdp"                       #mdp for equilibration and initial heating
md_mdp_gen_vel_no='a_300K_4ps_gen_vel_no.mdp'          #mdp file for production run steps, with gen-vel = no
md_mdp_gen_vel_yes='a_300K_4ps_gen_vel_yes.mdp'        #mdp file for production run steps, with gen-vel = yes
evap_log='x_LOG_evap.txt'                              #output file for droplet cleanup operation (number of molecules vs. time)
grott_log='x_LOG_grott.txt'                            #output file for H3O_Grotthuss.py (keeps track of hopping events for each H+)
n_steps_H3O_Grotthuss=10                               #number of times that H3O_Grotthuss program is called in between droplet_cleanup steps
n_hop=5                                                #number of times that each H+ can hop during each call of H3O_Grotthuss program
hop_cutoff=0.25                                        #maximum allowed O...H distance (nm) for H+ hopping

     
#~~~~~~~~~~~~~ PART 1 of this script generates the droplet
#Generate "dry" droplet gro file
python droplet_stuffer.py ${droplet_stuffer_input}

#Add water (water box will be cube-shaped)
gmx solvate -cp droplet_dry.gro -cs tip4p.gro -o 00000droplet_cube.gro

#Put system in large cubic box  
gmx editconf -f 00000droplet_cube.gro -o 00000droplet_cube_large.gro -bt cubic -box 999.9 -c

#Carve droplet; generate proper .gro and .top starting files
python droplet_carver_PT.py 00000droplet_cube_large.gro 00000droplet_start.gro ${raw_top} 00000droplet_start.top ${carv_r_ini} OW droplet_carver_ini_log.txt  PT=no 0.0 0.0

#Energy minimization
gmx grompp -f ${em_mdp} -c 00000droplet_start.gro -p 00000droplet_start.top -o 00000_em.tpr
gmx mdrun -v -deffnm 00000_em

#Equilibrate (raise droplet temperature)
gmx grompp -f ${equil_mdp} -c 00000_em.gro -p 00000droplet_start.top -o 00000_equil.tpr
gmx mdrun -v -deffnm 00000_equil

cp 00000_equil.gro 00000_000.gro
cp 00000droplet_start.top 00000.top


#~~~~~~~~~~~~~ PART 2 of this script performs droplet evaporation and H3O Grotthuss diffusion
n_steps_md_tot=$((n_steps_carve))
last_hop=$((n_hop-1))
Last_hop=$(printf %01d ${last_hop})  #in Python numbering the last file will have number (n_hop-1)
   
python droplet_carver_PT.py 00000_000.gro 00001_000.gro 00000.top 00001.top ${carv_r_cleanup} H3O ${evap_log} PT=no 0.0 0.0


for i in $( eval echo {1..$n_steps_md_tot} )  #outer loop performs droplet carving steps  (removal of evaporated molecules) 
do
  ip1=$((i+1)) 
  I=$(printf %05d ${i})     #(outer loop step number  ) is converted to i5 text format
  Ip1=$(printf %05d ${ip1}) #(outer loop step number+1) is converted to i5 text format
  for k in $( eval echo {1..$n_steps_H3O_Grotthuss} )  #inner loop performs Grotthuss hopping 
  do
    km1=$((k-1)) 
    K=$(printf %03d ${k})     #(inner loop step number  )    is converted to i5 text format
    Km1=$(printf %03d ${km1}) #(inner loop step number-1)    is converted to i5 text format
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@ Grotthuss step k = "${K}
    python H3O_Grotthuss.py ${I}_${Km1}.gro ${n_hop} ${hop_cutoff} out-last ${grott_log}
    

    if [ $k -eq 1 ]   # first round in each loop segment is done with gen-vel = yes to counteract evaporative cooling
    then
      gmx grompp -f ${md_mdp_gen_vel_yes} -c ${I}_${Km1}_grott${Last_hop}.gro -p ${I}.top -o ${I}_${K}.tpr
    fi
    if [ $k -gt 1 ]   # all other rounds in each loop segment are done with gen-vel = no
    then
      gmx grompp -f ${md_mdp_gen_vel_no} -c ${I}_${Km1}_grott${Last_hop}.gro -p ${I}.top -o ${I}_${K}.tpr
    fi


    gmx mdrun -deffnm ${I}_${K} -gpu_id 0

    #Keep fairly complete record of all frames by concatenating xtc files
    #For gro movie use  $ gmx trjconv -f 00001.xtc -s 00001_001.tpr  -o 00001_movie.gro
    if [ $k -eq 1 ]
    then
      cp ${I}_${K}.xtc ${I}.xtc
    fi
    if [ $k -gt 1 ]
    then
      cat ${I}_${K}.xtc >> ${I}.xtc
    fi

    
    if [ $k -ne 1 ]  #delete most files, but keep some
    then
      echo "deleting files: "${I}_${K}
      rm ${I}_${K}.tpr
      rm ${I}_${K}.log
      rm ${I}_${K}.edr
      rm ${I}_${Km1}.gro
      rm *_grott4.gro
      rm ${I}_${K}.xtc
    fi

    rm *#
    rm *.cpt
    rm *_not_centered.gro

    done
  echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ droplet cleanup step i = "${I}
  #Perform clean-up: remove evaporated molecules
  python droplet_carver_PT.py ${I}_${K}.gro ${Ip1}_000_not_centered.gro ${I}.top ${Ip1}.top ${carv_r_cleanup} H3O ${evap_log} PT=no 0.0 0.0

  #Center system after cleanup  
  gmx editconf -f ${Ip1}_000_not_centered.gro -o ${Ip1}_000.gro -bt cubic -box 999.9 -c

  rm *_001.xtc

done
rm *_000.gro



rm *_clean.gro


echo "========================== bash DONE! ============================="







