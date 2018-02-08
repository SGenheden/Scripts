# Author: Samuel Genheden samuel.genheden@gmail.com

import sys
import math

avogrado = 6.022 * math.pow(10.0,23.0) # 1 / mol
dens_water = 997.0 # kg / m3
dens = 789.0 # kg / m3
molmass_water = 18.02 # g / mol

if sys.argv[1] == "eth" :
    molmass = 46.07 # g / mol
elif sys.argv[1] == "but" :
    molmass = 74.12 # g / mol

for n_arg, n_water_arg in zip(sys.argv[2::2], sys.argv[3::2]) :
    n = float(n_arg)
    mol =  n / avogrado
    mass = molmass * mol
    mass_kg = math.pow(10.0,-3.0) * mass
    vol = mass_kg / dens
    vol_litre = math.pow(10.0,3.0) * vol

    n_water = float(n_water_arg)
    mol_water = n_water / avogrado
    mass_water = molmass_water * mol_water
    mass_water_kg = math.pow(10.0,-3.0) * mass_water
    vol_water = mass_water_kg / dens_water
    vol_water_litre = math.pow(10.0,3.0) * vol_water


    konc = mass / (vol_litre + vol_water_litre)
    konc_mol = mol / (vol_litre + vol_water_litre)
    print "%.2f %.2f"%(konc, konc_mol)
