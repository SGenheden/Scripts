name=$1
a=`grep "SET A " ${name}.str | awk '{print $4*0.1}'`
b=`grep "SET B " ${name}.str | awk '{print $4*0.1}'`
c=`grep "SET C " ${name}.str | awk '{print $4*0.1}'`
gmx editconf -f ${name}.pdb -o ${name}.gro -box $a $b $c
