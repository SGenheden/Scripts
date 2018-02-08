# Author: Samuel Genheden samuel.genheden@gmail.com

import sys
import math

avogrado = 6.022 * math.pow(10.0,23.0)

n_water = 5120.0
molmass_water = 18.02
mol_water = n_water / avogrado
mass_water = molmass_water * mol_water
mass_water_kg = math.pow(10.0,-3.0) * mass_water
dens_water = 997.0
vol_water = mass_water_kg / dens_water
vol_water_litre = math.pow(10.0,3.0) * vol_water

if sys.argv[1] == "eth" :
    molmass = 46.07 # g / mol
    dens = 789.0 # g / L
elif sys.argv[1] == "but" :
    molmass = 74.12 # g / mol
    dens = 810.0 # g / L

for n_arg, n_water_arg in zip(sys.argv[2::2], sys.argv[3::2]) :
    n = float(n_arg)
    mol =  n / avogrado # mol
    mass = molmass * mol # g
    vol = mass / dens # L

    n_water = float(n_water_arg)
    mol_water = n_water / avogrado # mol
    mass_water = molmass_water * mol_water # g
    vol_water = mass_water / dens_water # L

    konc = mass / (vol + vol_water) # g / L
    print "%s\t%.3f"%(sys.argv[1], konc)
