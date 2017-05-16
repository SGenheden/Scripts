
import argparse
import os
import warnings

import numpy as np
import parmed

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program characterise a CG topology")
    parser.add_argument('-p','--topol',nargs="+", help="the topology files")
    args = parser.parse_args()

    struct_count = []
    struct_list = []
    all_types = []
    for t in args.topol :
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", parmed.exceptions.GromacsWarning)
            topol = parmed.gromacs.GromacsTopologyFile(t, parametrize=False)
        type_dict = {}
        for atom in topol.atoms :
            if atom.type not in type_dict :
                type_dict[atom.type] = 0.0
            type_dict[atom.type] += 1.0
            all_types.append(atom.type)
        for atype in type_dict :
            type_dict[atype] /= float(len(topol.atoms))
        try :
            idx = struct_list.index(type_dict)
            struct_count[idx] += 1
        except :
            struct_list.append(type_dict)
            struct_count.append(1)

        with open(os.path.splitext(t)[0]+"_struct.txt", "w") as f:
            for atype, frac in type_dict.iteritems() :
                f.write("%s:%0.2f\n"%(atype, frac))

    for idx in np.argsort(struct_count)[::-1] :
        print struct_list[idx], struct_count[idx]
    print sum(struct_count)
    print ""

    len_stat = {}
    for struct, count in zip(struct_list, struct_count) :
        struct_len = len(struct)
        if struct_len not in len_stat :
            len_stat[struct_len] = 0
        len_stat[struct_len] += count

    for l in range(1,10) :
        if l in len_stat :
            print "%d\t%d"%(l, len_stat[l])

    with open("cg_types.txt", "w") as f :
        for atype in sorted(set(all_types)) :
            f.write("%s\n"%atype)
