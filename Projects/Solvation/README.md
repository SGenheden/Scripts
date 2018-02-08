# Solvation

This folder contains scripts created while working on this paper:

Genheden S., Predicting partition coefficients with a simple all-atom/coarse-grained hybrid model. *J. Chem. Theory Comput.*, **2015**, *12*, 297-304

and

Genheden S., Solvation free energies and partition coefficients with the coarse-grained and hybrid all-atom/coarse-grained MARTINI models *J. Comput-Aided. Mol Des.*, **2017**, *31*, 867-876

They were written to be used with python version 2.7 and requires the following packages
* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [matplotlib](http://www.matplotlib.org/)
* [MDAnalysis](http://www.mdanalysis.org/)
* [openpyxl](https://openpyxl.readthedocs.org/en/latest/)
* [croc](http://pythonhosted.org/CROC/)

In addition, you need the following programs installed and accessible through the system $PATH.
* [AmberTools 15](http://wwww.ambermd.org)
* [OpenBabel](http://openbabel.org/)
* [checkmol](http://merian.pch.univie.ac.at/~nhaider/cheminf/cmmm.html)

Many of the scripts were used to work with the data in the Minnesota solvation database,
which can be accessed from [here](http://comp.chem.umn.edu/mnsol/). In order to
use these scripts, one has to obtain a license for the database and download it.

Many of the scripts works with solvation free energy data for a list of *solutes* in
a particular *solvent*. The class `SolvDb` in `dblib.py` is used to provide an interface
to the database and useful generators to loop over solutes.  

The workflow of using the scripts are something like this:

Setup
------
* `analyse_elements.py` - is used to find solutes which contain elements not in GAFF and hence cannot be parameterised
* `analyse_weights.py` - is used to obtain a list of solutes that are has solvation free energy in two different solvents, subjected to contraints on the weight/size.
* `analyse_molradii.py` - is used to calculate a rough molecular radii of the solutes
* `param_solutes.py` - is used to parameterise the solutes. It will take the xyz coordinates from the database and send them to `Antechamber`, `ParmCheck` and `tleap` from *AmberTools* to obtain GAFF parameters and AM1-BCC charges and an Amber parameter/topology file. This can then be converted to Lammps format using public Lammps tools.
* `run_commands.py` - is used to run an external command for each solute. The commands are read from a text file and the $$ string is replaced with the solute filehandle.

This *command* was for instance be used to insert the solutes in a box of water molecules:

````
python $SCRIPTS/Lammps/insert_in_elba.py ../../Solutes/data.$$ -b ../../Boxes/Water/data.water_box -f ../../Boxes/Hexanol/forcefield.elba --center --noff -o water_$$
````

Analysis
--
* `check_results.py` - is used to check that all Lammps simulations are finished and has produced the expected result
* `collect_results.py`  - is used to read the DV/DL files produced by Lammps, perform TI and put everything in an Excel spreadsheet.
* `analyse_chemicalgroups.py` - is used to create a list of chemical groups for each solute
* `analyse_errors.py` - is used to perform error and BEDROC analysis on the results from `collect_results.py`
* `plot_correlations.py` - is used to plot predicted free energies versus experimental free energies.
