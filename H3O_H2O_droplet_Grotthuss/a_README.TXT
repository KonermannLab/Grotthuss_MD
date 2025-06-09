
This directory contains a H2O droplet with 3 H3O+ (a_start.gro, a_start.top)

The bash file a_droplet-evap-H3O-GROTT.bash performs short Gromacs MD segments,
alternating with "bursts" of Grotthuss proton hopping. These bursts take place quasi-
instantaneously, in a H2O environment that is more or less static. Individual stages
of each burst are captured in the _grott0 ... grott4 files.

Hopping information for each proton (whether a hop did occur or not) are summarized in x_LOG_grott.txt

Run bash script with 
$ bash a_droplet-evap-H3O-GROTT.bash


For more information on files etc., see the header of a_droplet-H3O-GROTT.bash

Questions?
Lars Konermann, konerman@uwo.ca

Reference:
############################################################################################################
# "Grotthuss Molecular Dynamics Simulations for Modeling Proton Hopping in Electrosprayed Water Droplets"  #
# L. Konermann and S. Kim, J. Chem. Theory Comput. 18, 3781–3794 (2022)                                    #
############################################################################################################
