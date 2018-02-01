# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to count the number of files in each sub-directory of
the current directory

Examples:
  count_files.py
"""

import os

paths = os.listdir(".")
total = 0
total_zip = 0
print "%-20s %5s %5s"%("Folder", "Files", "Zip-Files")
for path in paths:
    if os.path.isfile(path) : continue
    n = 0
    n_zip = 0
    for root, dirs, files in os.walk(path) :
        n += len(files)
        n_zip += len([file for file in files if os.path.splitext(file)[1] in [".gz",".tar",".zip",".bz2"]])
    print "%-20s %5d %5d"%(path, n, n_zip)
    total += n
    total_zip += n_zip
print "%-20s %5d %5d"%("total", total, total_zip)
