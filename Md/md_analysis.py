# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse a MD trajectory with one or more acions

The command file can contain all actions that are defined in sgenlib.mdactions

Examples:
    md_analysis.py -f sim.dcd -s ref.pdb -c commands
"""

import argparse
import inspect
import sys
import shlex

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    def pred(c) :
        return inspect.isclass(c) and c.__module__ == "sgenlib.mdactions" and \
                issubclass(c,mdactions.TrajectoryAction) and \
                c.__name__ != "TrajectoryAction"

    def setup_action(cls,arguments):
        action = cls(processor)
        argparser = argparse.ArgumentParser()
        action.add_arguments(argparser)
        action.setup(argparser.parse_args(arguments))

    actionclasses = {name.lower():c
        for name,c in inspect.getmembers(sys.modules["sgenlib.mdactions"],pred)}

    processor = moldyn.TrajectoryProcessor("Analyze MD trajectories",dosubsample=True)
    processor.argparser.add_argument("-c","--commands",help="a file with analysis commands")
    processor.setup()

    with open(processor.args.commands,"r") as f :
        for line in f.readlines():
            action = line.strip().split()[0].lower()
            arguments = shlex.split(line.strip())[1:]
            try :
                cls = actionclasses[action]
            except:
                try:
                    cls = actionclasses[action+"analysis"]
                except:
                    print "Skipping unknown action %s"%action
                else:
                    setup_action(cls,arguments)
            else:
                setup_action(cls,arguments)

    processor.process()
