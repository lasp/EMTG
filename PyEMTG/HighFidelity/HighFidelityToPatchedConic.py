# Copyright (c) 2024 The Regents of the University of Colorado.
# All Other Rights Reserved.

import sys
sys.path.append("/Utilities/emtg/PyEMTG")
sys.path.append("/Utilities/emtg/PyEMTG/PEATSA")
import MissionOptions

MO = MissionOptions.MissionOptions("myOrbit1.emtgopt")

import PEATSAchef

chef = PEATSAchef.PEATSAchef()

newMO = chef.ThreeD_Flybys_to_PatchedConic(MO)

newMO.mission_name = "PatchedConicLaunchcOpen"

newMO.write_options_file("PatchedConicLaunchOpen.emtgopt")
