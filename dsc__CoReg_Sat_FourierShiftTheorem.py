# -*- coding: utf-8 -*-
__author__='Daniel Scheffler'

import warnings
import os
import sys
import time

warnings.warn("\n\n"+'#'*90+"\n\nWARNING: dsc__CoReg_Sat_FourierShiftTheorem.py will be removed in future! \n\n"
              "The argument parser script 'dsc__CoReg_Sat_FourierShiftTheorem.py' has been renamed to 'coreg_cmd.py'. "
              "Please use that script instead when using CoReg_Sat from command line. \n"
              "There have also been some changes in the structure of the argument parser due to a lot of new "
              "features:\n\n"
              "The parser has now two sub-parsers: 'global' and 'local' that allow to detect and correct global X/Y "
              "shifts as well as local geometric shifts. The documentation of these sub-parsers can be accessed using "
              "the following commands:\n\n"
              "\tpython /path/to/CoReg_Sat/coreg_cmd.py global -h\n\n"
              "\tpython /path/to/CoReg_Sat/coreg_cmd.py local -h\n\n"
              "Going on in 30 seconds...\n\n"+'#'*90)

time.sleep(30)
cmd = sys.argv

cmd[0] = cmd[0].replace('dsc__CoReg_Sat_FourierShiftTheorem.py', 'coreg_cmd.py global')
os.system('python %s' %' '.join(cmd))