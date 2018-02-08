echo "System" | trjconv -f /media/data3/sg6e12/Gpcr_trj/Cg/A2a/3qak/Chol-30/r1_md3_fit.xtc -o tmp.gro -sep -tu us -b 50 -e 50 -s ../../r1_md3.tpr
mv tmp0.gro r1_md3_fit_snap50.gro

echo "System" | trjconv -f /media/data3/sg6e12/Gpcr_trj/Cg/A2a/3qak/Chol-30/r1_md3_fit.xtc -o tmp.gro -sep -tu us -b 40 -e 40 -s ../../r1_md3.tpr
mv tmp0.gro r1_md3_fit_snap40.gro

echo "System" | trjconv -f /media/data3/sg6e12/Gpcr_trj/Cg/A2a/3qak/Chol-30/r1_md3_fit.xtc -o tmp.gro -sep -tu us -b 30 -e 30 -s ../../r1_md3.tpr
mv tmp0.gro r1_md3_fit_snap30.gro

echo "System" | trjconv -f /media/data3/sg6e12/Gpcr_trj/Cg/A2a/3qak/Chol-30/r1_md3_fit.xtc -o tmp.gro -sep -tu us -b 20 -e 20 -s ../../r1_md3.tpr
mv tmp0.gro r1_md3_fit_snap20.gro

echo "System" | trjconv -f /media/data3/sg6e12/Gpcr_trj/Cg/A2a/3qak/Chol-30/r1_md3_fit.xtc -o tmp.gro -sep -tu us -b 10 -e 10 -s ../../r1_md3.tpr
mv tmp0.gro r1_md3_fit_snap10.gro

python2.7 $SCRIPTS/Pdb/gro2pdb_con.py -f *gro -pb ../protein.bonds -hb chol.bonds popc.bonds
