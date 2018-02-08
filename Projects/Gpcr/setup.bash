# Author: Samuel Genheden samuel.genheden@gmail.com

export PATH=$PATH:/home/sg6e12/Programs/Gromacs-4.6.2/bin

# To setup the protein
editconf -f 3qak_cg.pdb -o 3qak_cg.gro -box 12 12 12 -rotate 90 0 0
grompp -f ../../Martini/em.mdp -c 3qak_cg.gro -p a2a.top -maxwarn 10 -o em.tpr
mdrun -nt 1 -deffnm em -v -c 3qak_cg_min.gro

#
# To set up Wat
#
genbox -cp ../3d4s_cg_min.gro -ci ../../../Martini/water_single.gro -nmol 10000 -box 12 12 12 -try 500 -o r${X}_3d4s_wat.gro -seed $RANDOM
grompp -f ../../../Martini/em.mdp -c r${X}_3d4s_wat.gro -p system_en.top -maxwarn 10 -o em.tpr
mdrun -nt 1 -deffnm em -c r${X}_3d4s_wat_min.gro
grompp -f md1.mdp -c r${X}_3d4s_wat_min.gro -p system_en.top -maxwarn 10 -o r${X}_md1.tpr
mdrun -nt 4 -deffnm r1_md1.tpr

#
# To setup up Chol-0
#
make_ndx -f r1_md2_cent.gro << EOF
q
EOF

for X in {1..10}
do

#genbox -cp ../4eiy_cg_min.gro -ci ../../../Martini/popc_single.gro -nmol 200 -box 12 12 12 -try 500 -o r${X}_4eiy_popc.gro -seed $RANDOM >& log
#sed -i "s/W /;W /" system.top
#grompp -f ../../../Martini/em.mdp -c r${X}_4eiy_popc.gro -p system.top -maxwarn 10 -o em.tpr >& log
#mdrun -nt 1 -deffnm em -c r${X}_4eiy_popc_min.gro >& log

#genbox -cp r${X}_4eiy_popc_min.gro -ci ../../../Martini/water_single.gro -nmol 6000 -try 500 -o r${X}_4eiy_popc_wat.gro -vdwd 0.21 -seed $RANDOM  >& log
#sed -i "s/;W /W /" system.top
#grompp -f ../../../Martini/em.mdp -c r${X}_4eiy_popc_wat.gro -p system.top -maxwarn 10 -o em.tpr >& log
#mdrun -nt 1 -deffnm em -c r${X}_4eiy_popc_wat_min.gro >& log

grompp -f md1.mdp -c r${X}_4eiy_popc_wat_min.gro -p system.top -maxwarn 10 -o r${X}_md1.tpr >& log
#grompp -f md1.mdp -c r${X}_4eiy_popc_wat_min.gro -p system_en.top -maxwarn 10 -o r${X}_md1_en.tpr >& log

rm -rf \#* em.* mdout.mdp

done

for X in {1..10}
do
trjconv -f r${X}_md1.gro -o r${X}_md1_cent.gro -pbc mol -center -s r${X}_md1.tpr << EOF
1
0
EOF
trjconv -f r${X}_md2_en.gro -o r${X}_md2_cent_en.gro -pbc mol -center -s r${X}_md2_en.tpr << EOF
1
0
EOF
done

for X in {1..5}
do
trjconv -f r${X}_md3.gro -o r${X}_md3_cent.gro -pbc mol -center -s r${X}_md3.tpr << EOF
1
0
EOF
done

sed "s/\/Protein/\/Protein_ins/" system.top > system_ins.top
sed "s/\/Protein/\/Protein_ins/" system_en.top > system_ins_en.top

for X in {1..5}
do
grompp -f md3.mdp -c r${X}_md2_cent.gro -p system_nores.top -maxwarn 10 -o r${X}_md3.tpr  >& log
#grompp -f md3.mdp -c r${X}_md2_cent_en.gro -p system_ins_en.top -maxwarn 10 -o r${X}_md3_en.tpr  >& log
done

# For A2a
for X in {2..10}
do
python ~/Programs/a2a_solv.py 3000 r${X}_md1_cent.gro  > r${X}_md1_cent_solv.gro
grompp -f ../../../Martini/em.mdp -c r${X}_md1_cent_solv.gro -p system.top -maxwarn 10 -o em.tpr >& log
mdrun -nt 1 -deffnm em -c r${X}_md1_cent_solv_min.gro >& log
grompp -f md2.mdp -c r${X}_md1_cent_solv_min.gro -p system.top -maxwarn 10 -o r${X}_md2.tpr  >& log
done

#
# To setup Chol-10 to Chol-50
#
#


# Number of cholesterols
N=100

make_ndx -f r${X}_md1_min.gro << EOF
q
EOF

for X in {1..5}
do

python ~/Programs/insert_chol.py ../../../Martini/chol_single.gro ../Chol-0/r${X}_md2_cent.gro ${N}  > r${X}_md1.gro
grompp -f ../../../Martini/em.mdp -c r${X}_md1.gro -p system.top -maxwarn 10 -o em.tpr  >& log
mdrun -nt 1 -deffnm em -c r${X}_md1_min.gro  >& log
grompp -f md2.mdp -c r${X}_md1_min.gro -p system.top -maxwarn 10 -o r${X}_md2.tpr  >& log

#python ~/Programs/insert_chol.py ../../../Martini/chol_single.gro ../Chol-0/r${X}_md2_cent_en.gro ${N}  > r${X}_md1_en.gro
#grompp -f ../../../Martini/em.mdp -c r${X}_md1_en.gro -p system.top -maxwarn 10 -o em.tpr  >& log
#mdrun -nt 1 -deffnm em -c r${X}_md1_min_en.gro  >& log
#grompp -f md2.mdp -c r${X}_md1_min_en.gro -p system_en.top -maxwarn 10 -o r${X}_md2_en.tpr  >& log

rm -rf \#* em.* mdout.mdp step*.pdb

done


## For alternative structures
for X in {2..5}
do
python2.7 ~/Programs/gpcr_fit_alt_struct.py ../../4eiy/4eiy_cg_min.gro ../../4eiy/Chol-30/r${X}_md2_cent.gro ../3qak_cg_min.gro > r${X}_md1_fit.gro
grompp -f ../../../Martini/em.mdp -c r${X}_md1_fit.gro -p system.top -maxwarn 10 -o em.tpr >& log
mdrun -nt 1 -deffnm em -c r${X}_md1_min.gro  >& log
grompp -f md2.mdp -c r${X}_md1_min.gro -p system.top -maxwarn 10 -o r${X}_md2.tpr  >& log
rm -rf \#* em.* mdout.mdp step*.pdb
done

# A2a fitting
python2.7 ~/Programs/gpcr_overlay.py -i 4eiy 3qak -x 1 2 4 7 -p a2a -c A A

#
#
# Setup large patch
#
#
make_ndx -f r1_md2_cent.gro << EOF
q
EOF

for X in {1..10}
do

#genbox -cp 3d4s_cg_min.gro -ci ../../../Martini/popc_single.gro -nmol 400 -box 18 18 18 -try 500 -o r${X}_3d4s_popc.gro -seed $RANDOM >& log
#sed -i "s/W /;W /" system_en.top
#grompp -f ../../../Martini/em.mdp -c r${X}_3d4s_popc.gro -p system_en.top -maxwarn 10 -o em.tpr >& log
#mdrun -nt 4 -deffnm em -c r${X}_3d4s_popc_min.gro >& log

#genbox -cp r${X}_3d4s_popc_min.gro -ci ../../../Martini/water_single.gro -nmol 7000 -try 500 -o r${X}_3d4s_popc_wat.gro -vdwd 0.21 -seed $RANDOM  >& log
sed -i "s/;W /W /" system_en.top
grompp -f ../../../Martini/em.mdp -c r${X}_3d4s_popc_wat.gro -p system_en.top -maxwarn 10 -o em.tpr >& log
mdrun -nt 4 -deffnm em -c r${X}_3d4s_popc_wat_min.gro >& log

grompp -f md1.mdp -c r${X}_3d4s_popc_wat_min.gro -p system_en.top -maxwarn 10 -o r${X}_md1_en.tpr >& log

rm -rf \#* em.* mdout.mdp

done

for X in {1..10}
do
python ~/Dropbox/Scripts/insert_chol.py ../../../Martini/chol_single.gro ../Chol-0_large/r6_md2_en_cent.gro ${N}  > r${X}_md1_en.gro
grompp -f ../../../Martini/em.mdp -c r${X}_md1_en.gro -p system.top -maxwarn 10 -o em.tpr  >& log
mdrun -nt 1 -deffnm em -c r${X}_md1_min_en.gro  >& log
grompp -f md2.mdp -c r${X}_md1_min_en.gro -p system_en.top -maxwarn 10 -o r${X}_md2_en.tpr  >& log
done

for X in 6
do
grompp -f md3.mdp -c r${X}_md2_en_cent.gro -p system_ins_en.top -maxwarn 10 -o r${X}_md3_en.tpr  >& log
done
