# Grotthuss_MD

These files are for simulating Grotthuss diffusion of protons in water, as described in

############################################################################################################
# "Grotthuss Molecular Dynamics Simulations for Modeling Proton Hopping in Electrosprayed Water Droplets"  #
# L. Konermann and S. Kim, J. Chem. Theory Comput. 18, 3781–3794 (2022)                                    #
############################################################################################################

Simulations are performed by alternating between brief Gromacs run segments (using static protons), and a 
Python program (which performs proton hopping).

There are two directories, each with its own README file
1. In /H3O_H2O_droplet_Grotthuss we perform Grotthuss diffusion on an existing droplet with 3 H3O+ (quite simple).
2. In /H3O_H2O_droplet_make_evap_Grotthuss we generate a droplet, equilibrate, and perform Grotthuss diffusion
   with occasional removal of evaporated molecules (a bit more complicated).

Questions?
Lars Konermann, konerman@uwo.ca
