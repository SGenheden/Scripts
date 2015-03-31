
SCRIPTS=~/Dropbox/Scripts/

# Analyse cholesterol counts
python2.7 $SCRIPTS/Gpcr/gpcr_anal_counts.py -f r1_md3_fit_chol.leaflet.dat n1 r1_md3_fit_chol.mid.dat 1an3 2an3 r1_md3_fit_chol.mstate.6.dat 6a1 6a2 6a4 6a5 6a3 -l chol/intra. chol/extra. chol/bur. chol2/intra. chol2/extra. on on/intra. on/extra. on2/intra. on2/extra. on/bur.

python2.7 $SCRIPTS/Gpcr/gpcr_anal_counts.py -f r1_md3_fit_chol.leaflet.dat n1 r1_md3_fit_chol.mid.dat 1an3 2an3 r1_md3_fit_chol-oh.mstate.6.dat 6a1 6a2 6a4 6a5 6a3 -l chol/intra. chol/extra. chol/bur. chol2/intra. chol2/extra. on on/intra. on/extra. on2/intra. on2/extra. on/bur.

# Analyse lipid counts
python2.7 $SCRIPTS/Gpcr/gpcr_anal_counts.py -f r1_md3_fit_lip.leaflet.dat n1 r1_md3_fit_short.mstate.6.dat r1_md3_fit_long.mstate.6.dat 3o4 3a4 3a1 4a1 5a1 6a1 3a2 4a2 5a2 6a2 -l lip/intra. lip/extra. short-on long-on either-on both-on short-on/intra. long-on/intra. either-on/intra. both-on/intra. short-on/extra. long-on/extra. either-on/extra. both-on/extra.

# Analyse lifetimes
python2.7 $SCRIPTS/Gpcr/gpcr_anal_lifetimes.py -f r1_md3_fit_chol.mstate.6.dat r1_md3_fit_chol.mid.dat 1a2 r1_md3_fit_short.mstate.6.dat r1_md3_fit_long.mstate.6.dat 4o5 4a5 r1_md3_fit_chol.resstate.6.dat r1_md3_fit_short.resstate.6.dat r1_md3_fit_long.resstate.6.dat 9o10 9a10 -l chol/mol chol/bur chol-on/bur short/mol long/mol either/mol both/mol chol/hlx short/hlx long/hlx either/hlx both/hlx --helical 8 9 10 11 12 --mol a2a

# Analyse cholesterol residue contacts
RESCON=~/Dropbox/Research/GPCR/Res_contacts/
MOL=b2
python2.7 $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_fit_chol.resstate.6.dat r1_md3_fit_chol-oh.resstate.6.dat -l com oh -o ${RESCON}/${MOL}_chol_6A --mol ${MOL}

python2.7 $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_fit_chol-oh.buried.rstate.6.dat -l oh -o ~/Dropbox/Research/GPCR/Res_contacts/b2_chol_buried_6A --mol b2

python2.7 $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_fit_chol-oh.buried.rstate.5.dat -l oh -o ${RESCON}/${MOL}_chol_buried_5A --mol ${MOL}

# Analyse lipid residue contacts
python2.7 $SCRIPTS/Gpcr/gpcr_plot_rescontacts.py -f r1_md3_fit_long.resstate.6.dat r1_md3_fit_short.resstate.6.dat -l oleoyl palmitoyl -o ~/Dropbox/Research/GPCR/Res_contacts/a2a_lip_6A --mol a2a

# Analyse joint probability

python2.7 $SCRIPTS/Gpcr/gpcr_plot_resjointcontacts.py -f r1_md3_fit_joint.5.npz -l oh -m ohburr -o ${RESCON}/${MOL}_chol_buried_5A_joint --mol ${MOL}

# Density series of cholesterol


export DENS=~/Dropbox/Research/GPCR/Density

 python2.7 $SCRIPTS/Gpcr/gpcr_plot_densityseries.py -f r{1..5}_md3_en_fit_densities2.npz -o $DENS/density_series_b2 -d chol -m b2 -r densities2. densities3. densities.20. densities100. -l 500ns 1us 10us 50us

# Density comparison

export DENS=???

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_b2inactive_chol -d chol2 hydroxyl2 -m b2 -l chol chol-oh -c 2:chol2

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_b2inactive_lip_head -d popc head -m b2 -l popc popc-head -c 2:popc
 
python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_b2inactive_lip_tail -d popc tail -m b2 -l popc popc-tail -c 2:popc

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_b2inactive_lip_gly -d popc glycerol -m b2 -l popc popc-glycerol -c 2:popc

# Inactive vs active GPCRs

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -f2 Active/r{1..5}_md3_fit_densities.npz -o $DENS/density_b2_chol-oh -d hydroxyl2 hydroxyl2 -m b2 b2_a -l inactive active -c 1:chol2 2:chol2

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_fit_densities.npz -f2 Active/r{1..5}_md3_fit_densities.npz -o $DENS/density_a2a_chol-oh -d hydroxyl2 hydroxyl2 -m a2a a2a_a -l inactive active -c 1:chol2 2:chol2

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -f2 Active/r{1..5}_md3_fit_densities.npz -o $DENS/density_b2_chol -d chol2 chol2 -m b2 b2_a -l inactive active

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -f2 Active/r{1..5}_md3_fit_densities.npz -o $DENS/density_b2_lip -d popc popc -m b2 b2_a -l inactive active

# Density plots in paper
MOL=b2
python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_en_fit_densities.npz -o $DENS/density_${MOL}_ -d chol2 popc -m ${MOL} -l chol. lip.

python2.7 $SCRIPTS/Gpcr/gpcr_plot_density.py -f1 r{1..5}_md3_fit_densities.npz -o $DENS/density_${MOL} -d chol2 popc -m ${MOL} -l chol. lip.
