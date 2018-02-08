# Author: Samuel Genheden samuel.genheden@gmail.com

python2.7 ~/Programs/Scripts/Gpcr/gpcr_mdanal.py -x r${X}_md3_en_fit.xtc -s r1_md2.gro -p /local/scratch/sg6e12/data/Gpcr/Gpcr_anal/B2/Chol-30/ -a contacts -c 5 --skip 2 --reslist reslist

SCRIPTS=~/Dropbox/Scripts/

# Analyse cholesterol counts
python2.7 $SCRIPTS/Gpcr/gpcr_anal_counts.py -f r1_md3_fit_chol.leaflet.dat n1 r1_md3_fit_chol.mid.dat 1an3 2an3 r1_md3_fit_chol.mstate.6.dat 6a1 6a2 6a4 6a5 6a3 -l chol/intra. chol/extra. chol/bur. chol2/intra. chol2/extra. on on/intra. on/extra. on2/intra. on2/extra. on/bur.

python2.7 $SCRIPTS/Gpcr/gpcr_anal_counts.py -f r1_md3_fit_chol.leaflet.dat n1 r1_md3_fit_chol.mid.dat 1an3 2an3 r1_md3_fit_chol-oh.mstate.6.dat 6a1 6a2 6a4 6a5 6a3 -l chol/intra. chol/extra. chol/bur. chol2/intra. chol2/extra. on on/intra. on/extra. on2/intra. on2/extra. on/bur.

# Analyse lipid counts
python $SCRIPTS/Gpcr/gpcr_anal_counts.py -f r1_md3_en_fit_lip.leaflet.dat n1 r1_md3_en_fit_short.mstate.6.dat r1_md3_en_fit_long.mstate.6.dat 3an4 4an3 3o4 3a4 3a1 4a1 5a1 6a1 3a2 4a2 5a2 6a2 -l lip/intra. lip/extra. short-on long-on short-on2 long-on2 either-on both-on short-on/intra. long-on/intra. either-on/intra. both-on/intra. short-on/extra. long-on/extra. either-on/extra. both-on/extra.

# Analyse lifetimes
python2.7 $SCRIPTS/Gpcr/gpcr_anal_lifetimes.py -f r1_md3_fit_chol.mstate.6.dat r1_md3_fit_chol.mid.dat 1a2 r1_md3_fit_short.mstate.6.dat r1_md3_fit_long.mstate.6.dat 4o5 4a5 r1_md3_fit_chol.resstate.6.dat r1_md3_fit_short.resstate.6.dat r1_md3_fit_long.resstate.6.dat 9o10 9a10 -l chol/mol chol/bur chol-on/bur short/mol long/mol either/mol both/mol chol/hlx short/hlx long/hlx either/hlx both/hlx --helical 8 9 10 11 12 --mol a2a

# Analyse cholesterol residue contacts
RESCON=~/Dropbox/Research/GPCR/Res_contacts/
MOL=b2
python $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_fit_chol.resstate.6.dat r1_md3_fit_chol-oh.resstate.6.dat -l com oh -o ${RESCON}/${MOL}_chol_6A --mol ${MOL}

python $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_en_fit_chol.resstate.6.dat r1_md3_en_fit_chol-oh.resstate.6.dat -l com oh -o ${RESCON}/${MOL}_chol_6A --mol ${MOL} --every 2

python $SCRIPTS/Plot/plot_grid_img.py -f b2_chol_6A_com_1d.png b2_a_chol_6A_com_1d.png a2a_chol_6A_com_1d.png a2a_a_chol_6A_com_1d.png  -c 4 1 -s 6.85 3.4 -o all_chol_6A_com_1d.png
python $SCRIPTS/Plot/plot_grid_img.py -f b2_chol_6A_oh_1d.png b2_a_chol_6A_oh_1d.png a2a_chol_6A_oh_1d.png a2a_a_chol_6A_oh_1d.png  -c 4 1 -s 6.85 3.4 -o all_chol_6A_oh_1d.png


python $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_fit_chol-oh.resstate.6.dat -l oh -o ${RESCON}/${MOL}_chol_6A --mol ${MOL}

python2.7 $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_fit_chol-oh.buried.rstate.6.dat -l oh -o ~/Dropbox/Research/GPCR/Res_contacts/b2_chol_buried_6A --mol b2

python2.7 $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_fit_chol-oh.buried.rstate.5.dat -l oh -o ${RESCON}/${MOL}_chol_buried_5A --mol ${MOL}

# Sites
b2 --sites 321,320,43 160,157,161,118
b2_a --sites 153,126,122,157,71 208,209
a2a --sites 93,128,44,185 133,129 182,183
--sites l 133,93,128,44,185 129,93,128,44,185
a2a_a --sites 131,132,135,180,181,184,185 93,124,128 129,133 283,286

# Analyse lipid residue contacts
python2.7 $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_fit_long.resstate.6.dat r1_md3_fit_short.resstate.6.dat -l oleoyl palmitoyl -o ~/Dropbox/Research/GPCR/Res_contacts/a2a_lip_6A --mol a2a

# Analyse reslist

python $SCRIPTS/Gpcr/gpcr_anal_reslistcontacts.py -f r1_md3_en_fit_chol-oh.resliststate.5.dat  --sites $SITES--mol ${MOL}

# Analyse joint probability

python2.7 $SCRIPTS/Gpcr/gpcr_plot_resjointcontacts.py -f r1_md3_fit_joint.5.npz -l oh -m ohburr -o ${RESCON}/${MOL}_chol_buried_5A_joint --mol ${MOL}

# Density series of cholesterol


export DENS=~/Dropbox/Research/GPCR/Density

 python $SCRIPTS/Gpcr/gpcr_plot_densityseries.py -f r{1..5}_md3_en_fit_densities2.npz -o $DENS/density_series_b2 -d chol -m b2 -r densities2. densities3. densities.20. densities100. -l  '$0.5 \mu\mathrm{s}$' '$1 \mu\mathrm{s}$' '$10 \mu\mathrm{s}$' '$50 \mu\mathrm{s}$'

# Density comparison

export DENS=???

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_b2inactive_chol -d chol2 hydroxyl2 -m b2 -l chol chol-oh -c 2:chol2

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_b2inactive_lip_head -d popc head -m b2 -l popc popc-head -c 2:popc

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_b2inactive_lip_tail -d popc tail -m b2 -l popc popc-tail -c 2:popc

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_b2inactive_lip_gly -d popc glycerol -m b2 -l popc popc-glycerol -c 2:popc

# Inactive vs active GPCRs

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -f2 Active/r{1..5}_md3_fit_densities.npz -o $DENS/density_b2_chol-oh -d hydroxyl2 hydroxyl2 -m b2 b2_a -l inactive active -c 1:chol2 2:chol2

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_fit_densities.npz -f2 Active/r{1..5}_md3_fit_densities.npz -o $DENS/density_a2a_chol-oh -d hydroxyl2 hydroxyl2 -m a2a a2a_a -l inactive active -c 1:chol2 2:chol2

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -f2 Active/r{1..5}_md3_fit_densities.npz -o $DENS/density_${MOL}_chol -d chol2 chol2 -m ${MOL} ${MOL}_a -l inactive active

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -f2 Active/r{1..5}_md3_fit_densities.npz -o $DENS/density_b2_lip -d popc popc -m b2 b2_a -l inactive active

# Density plots in paper
MOL=b2
python $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_${MOL}_ -d chol2 popc -m ${MOL} -l chol. POPC

python $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_fit_densities.npz -o $DENS/density_${MOL} -d chol2 popc -m ${MOL} -l chol. POPC

# VR
MOL=b2
python $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_${MOL}_chol2 -d chol2 -m ${MOL} -l chol.

# Misc stuff

python2.7 gpcr_plot_fractalpies.py --probes 2.0 2.5 3.0 3.5 4.0 4.5 --npies 14
gpcr_plot_gridenergy.py -p /media/data4/sg6e12/Gpcr/Cg/B2/3d4s/Chol-30/system.top -m b2

# plot_motify.py
52,155,152,156,217,310 34,38,52,49,55,155 11,7,22,61,237,248,251 12,11,14,57,239,267

### REDUX
PDB=3d4s
DENS=~/Dropbox/Research/GPCR/Redux/Density
RESCON=~/Dropbox/Research/GPCR/Redux/Res_contacts

python $SCRIPTS/Gpcr/gpcr_anal_counts.py -f r1_${PDB}_prod_fit_chol.leaflet.dat n1 r1_${PDB}_prod_fit_chol.mid.dat 1an3 2an3 r1_${PDB}_prod_fit_chol.mstate.6.dat 6a1 6a2 6a4 6a5 6a3 -l chol/intra. chol/extra. chol/bur. chol2/intra. chol2/extra. on on/intra. on/extra. on2/intra. on2/extra. on/bur.
python $SCRIPTS/Gpcr/gpcr_anal_counts.py -f r1_${PDB}_prod_fit_lip.leaflet.dat n1 r1_${PDB}_prod_fit_short.mstate.6.dat r1_${PDB}_prod_fit_long.mstate.6.dat 3an4 4an3 3o4 3a4 3a1 4a1 7a1 8a1 3a2 4a2 7a2 8a2 -l lip/intra. lip/extra. short-on long-on short-on2 long-on2 either-on both-on short-on/intra. long-on/intra. either-on/intra. both-on/intra. short-on/extra. long-on/extra. either-on/extra. both-on/extra.

python $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_${PDB}_prod_fit_chol.resstate.6.dat r1_${PDB}_prod_fit_chol-oh.resstate.6.dat -l com oh -o ${RESCON}/${mol}_chol_6A --mol ${mol}
python /Users/samuel/Dropbox/Code/Scripts/Plot/plot_grid_img.py -f b2_chol_6A_com_1d.png b2_a_chol_6A_com_1d.png a2a_chol_6A_com_1d.png a2a_a_chol_6A_com_1d.png -c 4 1 -s 6.85 3.41 -o all_chol_6A_com_1d.png
python /Users/samuel/Dropbox/Code/Scripts/Plot/plot_grid_img.py -f b2_chol_6A_oh_1d.png b2_a_chol_6A_oh_1d.png a2a_chol_6A_oh_1d.png a2a_a_chol_6A_oh_1d.png -c 4 1 -s 6.85 3.41 -o all_chol_6A_oh_1d.png

python $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_${PDB}_prod_fit_densities.npz -o $DENS/density_${mol} -d chol2 popc -m ${mol} -l chol. POPC
For A2a: --min -5.740018
For B2_a: --min -7.076167
But choose:
For A2a_a:--min -4.950703
For B2: --min -3.891713
python $SCRIPTS/Gpcr/gpcr_plot_densityseries.py -f r{1..5}_${PDB}_prod_fit_densities2.npz -o $DENS/density_series_b2 -d chol -m b2 -r densities2. densities3. densities.20. densities100. -l  '$0.5 \mu\mathrm{s}$' '$1 \mu\mathrm{s}$' '$10 \mu\mathrm{s}$' '$50 \mu\mathrm{s}$'

python $SCRIPTS/Gpcr/gpcr_anal_lifetimes.py -f r1_${PDB}_prod_fit_chol.mstate.6.dat r1_${PDB}_prod_fit_chol.mid.dat 1a2 r1_${PDB}_prod_fit_short.mstate.6.dat r1_${PDB}_prod_fit_long.mstate.6.dat 4o5 4a5 r1_${PDB}_prod_fit_chol.resstate.6.dat r1_${PDB}_prod_fit_short.resstate.6.dat r1_${PDB}_prod_fit_long.resstate.6.dat 9o10 9a10 -l chol/mol chol/bur chol-on/bur short/mol long/mol either/mol both/mol chol/hlx short/hlx long/hlx either/hlx both/hlx --helical 8 9 10 11 12 --mol ${mol}
