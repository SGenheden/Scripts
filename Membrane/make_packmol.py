
import argparse
import os

import numpy as np

from sgenlib import pdb
import build_lipid

def write_singlepdb(filename, resname, atomname):
    struct = pdb.PDBFile()
    res = pdb.Residue()
    atm = pdb.Atom()
    atm.serial = 1
    atm.name = atomname
    atm.resname = resname
    atm.residue = 1
    res.append(atm)
    struct.residues.append(res)
    struct.write(filename)

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Make packmol input to setup membranes")
    parser.add_argument('-l','--lipids',nargs="+",help="the lipids to use")
    parser.add_argument('-n','--nlipid',nargs="+",help="the number of lipids to use of each kind")
    parser.add_argument('-i','--ions',type=int,help="the number of ions to insert",default=0)
    parser.add_argument('-x','--xml',help="the definition of templates")
    parser.add_argument('-o','--out',help="the output name",default="bilayer.pdb")
    parser.add_argument('--leafletpad',type=float,help="the distance between the leaflets",default=0)
    parser.add_argument('--area',type=float,help="the area per lipid",default=60)
    parser.add_argument('--nwater',type=int,help="the number of water per lipids",default=40)
    args = parser.parse_args()

    lipidbook = build_lipid.LipidCollection()
    if args.xml is None:
        thispath = os.path.dirname(os.path.abspath(__file__))
        args.xml = os.path.join(thispath,"lipid_templates.xml")
    lipidbook.load(args.xml)

    structs = []
    lipids = []
    for l in args.lipids:
        structs.append(pdb.PDBFile())
        structs[-1].read(l)
        lipids.append(lipidbook.lipids[structs[-1].residues[0].resname])

        newfilename = os.path.splitext(l)[0]+".pdb"
        if not os.path.exists(newfilename):
            structs[-1].write(newfilename)

    zlen = np.max([l.head[0].xyz-l.tail[0].xyz for l in lipids],axis=0)[2]

    # Parse the number of lipids of each kind in each leaflet
    if len(args.nlipid) == 1 :
        args.nlipid.append(args.nlipid[0])
    nlipid = [map(int,nspec.split(":")) for nspec in args.nlipid]
    # Calculate the side of the box in x and y, assuming 60 A2 area per lipid
    xylen = np.sqrt(np.max([args.area*np.sum(n) for n in nlipid]))

    # Use 40 real water molecule / lipid as default
    nwat = np.max([args.nwater*np.sum(n) for n in nlipid])

    # Calculate the z side of the water box from ideal density
    watz = (30*nwat)/(xylen*xylen)

    if not os.path.exists("wat.pdb"):
        write_singlepdb("wat.pdb", "W", "W")

    if args.ions > 0 and not os.path.exists("ion.pdb"):
        write_singlepdb("ion.pdb", "ION", "NA+")

    # Start to produce the output

    print "tolerance 2.0"
    print "filetype pdb"
    print "seed 32535"
    print "output %s"%args.out

    for nthis, lipid, struct, filename in zip(nlipid[0], lipids, structs, args.lipids):
        print "structure %s.pdb"%os.path.splitext(filename)[0]
        print "    number %d"%nthis
        print "    inside box %.3f %.3f %.3f %.3f %.3f %.3f"%(0.0, 0.0, watz,
                xylen, xylen, watz + zlen + 10.0)
        print "    atoms %s"%" ".join("%d"%(lipid.beads.index(b)+1) for b in lipid.head)
        print "        below plane 0. 0. 1. %.2f"%(watz + 5.0)
        print "    end atoms"
        print "    atoms %s"%" ".join("%d"%(lipid.beads.index(b)+1) for b in lipid.tail)
        print "        above plane 0. 0. 1. %.2f"%(watz + 5.0 +
                (lipid.head[0].xyz[2]-lipid.tail[0].xyz[2]))
        print "    end atoms"
        print "end structure"

    for nthis, lipid, struct, filename in zip(nlipid[1], lipids, structs, args.lipids):
        print "structure %s.pdb"%os.path.splitext(filename)[0]
        print "    number %d"%nthis
        print "    inside box %.3f %.3f %.3f %.3f %.3f %.3f"%(0.0, 0.0, watz + zlen + 10.0 + args.leafletpad,
                xylen, xylen, watz + 2.0*(zlen + 10.0) + args.leafletpad)
        print "    atoms %s"%" ".join("%d"%(lipid.beads.index(b)+1) for b in lipid.head)
        print "        above plane 0. 0. 1. %.2f"%(watz + args.leafletpad + 2.0*(zlen + 10.0) - 5.0)
        print "    end atoms"
        print "    atoms %s"%" ".join("%d"%(lipid.beads.index(b)+1) for b in lipid.tail)
        print "        below plane 0. 0. 1. %.2f"%(watz + args.leafletpad + 2.0*(zlen + 10.0) - 5.0
                - (lipid.head[0].xyz[2]-lipid.tail[0].xyz[2]))
        print "    end atoms"
        print "end structure"

    print "structure wat.pdb"
    print "    number %d"%(nwat*0.25)
    print "    inside box %.3f %.3f %.3f %.3f %.3f %.3f"%(0.0, 0.0, 0.0,
                xylen, xylen, watz)
    print "end structure"

    print "structure wat.pdb"
    print "    number %d"%(nwat*0.25)
    print "    inside box %.3f %.3f %.3f %.3f %.3f %.3f"%(0.0, 0.0, watz + 2*(zlen+10) + args.leafletpad ,
                xylen, xylen, 2*watz + 2*(zlen+10) + args.leafletpad )
    print "end structure"

    if args.ions > 0 :
        print "structure ion.pdb"
        print "    number %d"%(args.ions)
        print "    inside box %.3f %.3f %.3f %.3f %.3f %.3f"%(0.0, 0.0, 0.0,
                    xylen, xylen, watz)
        print "end structure"

        print "structure ion.pdb"
        print "    number %d"%(args.ions)
        print "    inside box %.3f %.3f %.3f %.3f %.3f %.3f"%(0.0, 0.0, watz + 2*(zlen+10) + args.leafletpad,
                    xylen, xylen, 2*watz + 2*(zlen+10) + args.leafletpad)
        print "end structure"
