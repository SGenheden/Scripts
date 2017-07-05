
itp=$1

mol=`grep -A 2 moleculetype ${itp}  | tail -n1 | awk '{print $1}'`

cat ${itp}
echo "
[ system ]
; Name
Generic title

[ molecules ]
; Compound       #mols
${mol}                  1"
