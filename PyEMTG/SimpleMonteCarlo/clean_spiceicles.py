# Copyright (c) 2024 The Regents of the University of Colorado.
# All Other Rights Reserved.

#command line utility to clean spicicles from an ephemeris file

def do_the_stuff(option_list):

    import ConOpsPeriod
    import SpiceyPy_Utilities as SpiceyUtil
    try:
        import spiceypy as spice
    except:
        print("spiceypy not available")
        
    working_directory = option_list[0]
    SPICE_ephem_directory = option_list[2]
    spice_handler = SpiceyUtil.SpiceHandler(SPICE_ephem_directory)
    spice_handler.loadSpiceFiles()
    
    my_ephemeris_reader = ConOpsPeriod.EphemerisFileReader()
    ephemeris_file_data = my_ephemeris_reader.parseEMTGephemerisFile(working_directory + '/' + option_list[1])

    if len(option_list) == 4:
        forward_integrated_ephemeris_minimum_timestep_kept = option_list[3]
        my_ephemeris_reader.generateCleanEphemerisFileForBSP(forward_integrated_ephemeris_minimum_timestep_kept)
    else:
        # default value is 120.0 seconds
        my_ephemeris_reader.generateCleanEphemerisFileForBSP()
    
    spice_handler.unloadSpiceFiles()

if __name__ == '__main__':
    import sys
    import os

    if len(sys.argv) != 4:
        raise Exception("Unknown number of command line options!\n\
                         Syntax: python clean_spiceicles.py working_directory input_ephemeris_file SPICE_folder forward_integrated_ephemeris_minimum_timestep_kept (optional)")
    
    thisdir = os.path.dirname(os.path.realpath(option_list[0]))
    sys.path.append(thisdir)
    sys.path.append(thisdir + '/../')
    sys.path.append(thisdir + '/../SpiceyPy_Utilities')
    
    do_the_stuff(sys.argv[1:])



