
import sys

with open(sys.argv[1], "r") as f :
    line = f.readline()
    while line and line.find("MOD4      RE") == -1 :
        line = f.readline()

    line = f.readline()
    while line :
        col = line.strip().split()
        if len(col) < 3 : break
        atype = col[0]
        rmin = float(col[1])
        eps = float(col[2])
        print "%-9s %-8s %10.5f  %10.6f  A %13.6e %13.6e ; %6.4f %6.4f"%(
            atype, atype, 0.0, 0.0, rmin * 2**(-1.0/6.0) * 2, eps * 4.184,
            rmin, eps)
        line = f.readline()
