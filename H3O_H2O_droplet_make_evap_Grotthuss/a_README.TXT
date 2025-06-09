
The files in this directory
- generate a droplet containing H3O and H2O
- allow the droplet to evolve over time, with periodic removal of evaporated molecules
- allow H3O+ to undergo Grotthuss diffusion

Hopping information for each proton (whether a hop did occur or not) are summarized in x_LOG_grott.txt
Changes in droplet composition are tracked in x_LOG_evap.txt (although there is not a lot of evaporation at 300 K)

Run bash script with 
$ bash a_droplet-evap-H3O-GROTT.bash

For more information on files etc., see the header of a_droplet-evap-H3O-GROTT.bash

Questions?
Lars Konermann, konerman@uwo.ca

Reference:
############################################################################################################
# "Grotthuss Molecular Dynamics Simulations for Modeling Proton Hopping in Electrosprayed Water Droplets"  #
# L. Konermann and S. Kim, J. Chem. Theory Comput. 18, 3781–3794 (2022)                                    #
############################################################################################################
