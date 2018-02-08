# Author: Samuel Genheden samuel.genheden@gmail.com

prot=3d4s
R=1

for R in {1..4}
do
$GMX/gmx insert-molecules -f ../${prot}_cg_min.gro -ci ../../../Martini/popc_single.gro -nmol 100 -box 12 12 6 -o r${R}_${prot}_popc1.gro -seed $RANDOM
python $SCRIPTS/Pdb/flip_z.py r${R}_${prot}_popc1.gro -o r${R}_${prot}_popc2.gro
$GMX/gmx insert-molecules -f r${R}_${prot}_popc2.gro  -ci ../../../Martini/popc_single.gro -nmol 100 -box 12 12 6 -o r${R}_${prot}_popc3.gro -seed $RANDOM
python $SCRIPTS/Pdb/flip_z.py r${R}_${prot}_popc3.gro -o r${R}_${prot}_popc4.gro
sed -i "" "s/W /;W /" system.top
$GMX/gmx grompp -f ../../../Martini/em_cg5.mdp -c r${R}_${prot}_popc4.gro -p system.top -o em_popc.tpr
$GMX/gmx mdrun -deffnm em_popc
mv em_popc.gro r${R}_${prot}_popc_em.gro
$GMX/gmx insert-molecules -f r${R}_${prot}_popc_em.gro -ci ../../../Martini/water_single.gro -nmol 4500 -try 500 -o r${R}_${prot}_popc_wat.gro -radius 0.21 -seed $RANDOM
sed -i "" "s/;W /W /" system.top
$GMX/gmx grompp -f ../../../Martini/em_cg5.mdp -c r${R}_${prot}_popc_wat.gro -p system.top -o em_popc_wat.tpr
$GMX/gmx mdrun -deffnm em_popc_wat
mv em_popc_wat.gro r${R}_${prot}_popc_wat_em.gro
$GMX/gmx grompp -f ../../../Martini/equil1.mdp -c r${R}_${prot}_popc_wat_em.gro -p system.top -n index.ndx -o r${R}_${prot}_equil1.tpr -maxwarn 2
done

scp -p r*_3d4s_equil1.tpr genheden@beskow.pdc.kth.se:/cfs/klemming/nobackup/g/genheden/Gpcr/Cg/B2/3d4s/Chol-0

$GMX/gmx make_ndx -f  r${R}_${prot}_popc_wat_em.gro -o index.ndx << EOF
r POPC CHOL
name 15 Membrane
quit
EOF

$GMX/gmx grompp -f ../../../Martini/equil1.mdp -c r${R}_${prot}_popc_wat_em.gro -p system.top -n index.ndx -o r${R}_${prot}_equil1.tpr -maxwarn 2
$GMX/gmx grompp -f ../../../Martini/equil2.mdp -c r${R}_${prot}_equil1.gro -p system.top -n index.ndx -o r${R}_${prot}_equil2.tpr -maxwarn 1

for R in {1..5}
do
$GMX/gmx trjconv -f r${R}_${prot}_equil2.gro -o r${R}_${prot}_equil2_cent.gro -s r${R}_${prot}_equil2.tpr -center -pbc mol << EOF
Protein
System
EOF
done

$GMX/gmx make_ndx -f r${R}_${prot}_chol_em.gro -o index.ndx << EOF
r POPC CHOL
name 16 Membrane
quit
EOF

python $SCRIPTS/Membrane/insert_chol.py -c ../../../Martini/chol_single.gro -f ../Chol-0/r${R}_${prot}_equil2_cent.gro -n 60 -o r${R}_${prot}_chol_ins.gro
$GMX/gmx grompp -f ../../../Martini/em_cg5.mdp -c r${R}_${prot}_chol_ins.gro -p system.top -o em_chol.tpr
$GMX/gmx mdrun -deffnm em_chol
mv em_chol.gro r${R}_${prot}_chol_em.gro
$GMX/gmx grompp -f ../../../Martini/equil2.mdp -c r${R}_${prot}_chol_em.gro -p system.top -n index.ndx -o r${R}_${prot}_equil3.tpr -maxwarn 2

for R in {1..5}
do
$GMX/gmx trjconv -f r${R}_${prot}_equil3.gro -o r${R}_${prot}_equil3_cent.gro -s r${R}_${prot}_equil3.tpr -center -pbc mol << EOF
Protein
System
EOF
done

$GMX/gmx grompp -f ../../../Martini/prod.mdp -c r${R}_${prot}_equil3.gro -p system.top -n index.ndx -o r${R}_${prot}_prod.tpr -maxwarn 1
