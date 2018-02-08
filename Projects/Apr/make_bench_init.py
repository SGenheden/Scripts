# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Extract the host-guest structure from the benchmark prmtop+rst7 files
"""

import sys
import parmed
struct = parmed.load_file(sys.argv[1], xyz=sys.argv[2])
gh=struct[":1-2"]
gh.save(sys.argv[3],overwrite=True)
