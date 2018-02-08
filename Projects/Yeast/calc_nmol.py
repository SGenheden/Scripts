# Author: Samuel Genheden samuel.genheden@gmail.com

import sys
import math

avogrado = 6.022 * math.pow(10.0,23.0)

n_water = n_total = 5120.0
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

for konc_arg in sys.argv[2:] :
    konc = float(konc_arg)
    mass = konc * vol_water_litre
    mol = mass / molmass
    n = round(mol * avogrado)

    ratio = 1.0
    ratio_water = round(n_water / n)
    fraction = ratio / (ratio + ratio_water)
    fraction_water = ratio_water / (ratio + ratio_water)
    n_real = round(fraction * n_water)
    n_real_water = round(fraction_water * n_water)
    #print "%s\t%d\twater\t%d"%(sys.argv[1], n_real, n_real_water)

for konc_arg in sys.argv[2:] :
    konc = float(konc_arg) # g/L
    konc_mol = konc / molmass # mol/L
    konc_n   = konc_mol * avogrado # 1 / L

    n = round(n_total * molmass_water * dens / (
        dens * dens_water * avogrado / konc_n +
        molmass_water * dens - molmass * dens_water))

    print "%s\t%d\twater\t%d"%(sys.argv[1], n, n_total - n)
