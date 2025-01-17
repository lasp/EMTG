# Copyright (c) 2024 The Regents of the University of Colorado.
# All Other Rights Reserved.

"""
make_journeyoptions_python.py
=============================
This file contains the make_PyEMTG_JourneyOptions function.

"""

def make_PyEMTG_JourneyOptions(OptionsDefinitions, now, path = '.'):  
    """
    Function for automatically creating a JourneyOptions.py file based on a csv of mission options.

    Parameters
    ----------
    OptionsDefinitions : List of dict entries
        Dict entries that define the options to be added to JourneyOptions.py
    now : string
        A string of the current time.
    path : string, optional
        Path to directory in which to write JourneyOptions.py. The default is '.'

    Returns
    -------
    None.

    """
    with open(path + "PyEMTG/JourneyOptions.py", "w") as file:
        file.write('"""\n')
        file.write('JourneyOptions.py\n')
        file.write('========================\n')
        file.write('auto-generated by make_EMTG_missionoptions_journeyoptions.py\n')
        file.write('"""')
        file.write('\n')
        file.write('class JourneyOptions(object):\n')
        file.write('\n')
        
        tab = '    '
        file.write(tab + '"""\n') # start np-style comment block
        file.write('\n')
        file.write(tab + "JourneyOptions class is used to read/write journey blocks from .emtgopt files and hold journey option variables.\n")
        
        file.write('\n')
        file.write(tab + '"""\n') # end np-style comment block

        file.write('    #************************************************************************************constructor\n')
        file.write('    def __init__(self, inputFile = None, lineNumber = None):\n')

        for option in OptionsDefinitions:
            if 'string' in option['dataType'] and option['name'] != 'user_data':
                #file.write('        self.' + option['name'] + ' = "' + str(option['defaultValue']) + '" #' + option['comment'] + '\n')
                file.write('        self.' + option['name'] + ' = "' + str(option['defaultValue']) + '"\n')
                file.write(tab + tab + '"""' + option['comment'] + '"""\n')
            else:
                #file.write('        self.' + option['name'] + ' = ' + str(option['defaultValue']) + ' #' + option['comment'] + '\n')
                file.write('        self.' + option['name'] + ' = ' + str(option['defaultValue']) + '\n')
                file.write(tab + tab + '"""' + option['comment'] + '"""\n')

        file.write('        \n')
        file.write('        #empty lists for constraint definitions and trialX\n')
        file.write('        self.ManeuverConstraintDefinitions = []\n')
        file.write('        self.BoundaryConstraintDefinitions = []\n')
        file.write('        self.PhaseDistanceConstraintDefinitions = []\n')
        file.write('        self.trialX = []\n')
        file.write('\n')
        file.write('        if inputFile != None:\n')
        file.write('            self.parse_journey(inputFile, lineNumber)\n')


        file.write('   \n')
        file.write('    #************************************************************************************parse\n')
        file.write('    def parse_journey(self, inputFile, lineNumber = 0):\n')
        file.write('        while True:\n')
        file.write('            line = inputFile.readline()\n')
        file.write('            if not line:\n')
        file.write('                break\n')
        file.write('            #strip off the newline character\n')
        file.write('            line = line.rstrip("\\n\\r ")\n')
        file.write('            \n')
        file.write('            lineNumber += 1\n')
        file.write('            if line == "END_JOURNEY":\n')
        file.write('                break\n')
        file.write('            \n')
        file.write('            #if we got this far, then this is a line worth reading\n')
        file.write('            #Note that unlike EMTG proper, PyEMTG does NOT length or bounds-check input files. But if you try to run an invalid .emtgopt, EMTG will notify you and help you fix it.\n')
        file.write('            \n')
        file.write('            if line.strip(\'\\r\') != "":\n')
        file.write('                if line[0] != "#":\n')
        file.write('                    #this is an active line, so it is space delimited\n')
        file.write('                    linecell = [entry.rstrip(" \\r\\n") for entry in line.split(" ")]\n')
        file.write('                    \n')

        ifelse = ''
        for option in OptionsDefinitions:
            length = 1
            if 'std::vector' in option['dataType']:
                length = len(option['defaultValue'])
            converter_in = ''
            converter_out = ''

            if 'double' in option['dataType']:
                converter_in = 'float('
                converter_out = ')'
            elif 'string' not in option['dataType']:#bool, size_t, int, time_t, or any enum
                converter_in = 'int('
                converter_out = ')'
            
            
            file.write('                    ' + ifelse + 'if linecell[0] == "' + option['name'] + '":\n')
            if length == 1:
                file.write('                        self.' + option['name'] + ' = ' + converter_in + 'linecell[1]' + converter_out + '\n')
            else:
                file.write('                        self.' + option['name'] + ' = [' + converter_in + 'entry' + converter_out + ' for entry in linecell[1:]]\n')
            file.write('                  \n')
            ifelse = 'el'

        file.write('                    elif linecell[0] == "BEGIN_MANEUVER_CONSTRAINT_BLOCK":\n')
        file.write('                        self.ManeuverConstraintDefinitions = []\n')
        file.write('                        while True:\n')
        file.write('                            entry = inputFile.readline()\n')
        file.write('                            lineNumber += 1\n')
        file.write('                            if "END_MANEUVER_CONSTRAINT_BLOCK" in entry:\n')
        file.write('                                break\n')
        file.write('                            self.ManeuverConstraintDefinitions.append(entry)\n')
        file.write('                    \n')

        file.write('                    elif linecell[0] == "BEGIN_BOUNDARY_CONSTRAINT_BLOCK":\n')
        file.write('                        self.BoundaryConstraintDefinitions = []\n')
        file.write('                        while True:\n')
        file.write('                            entry = inputFile.readline()\n')
        file.write('                            lineNumber += 1\n')
        file.write('                            if "END_BOUNDARY_CONSTRAINT_BLOCK" in entry:\n')
        file.write('                                break\n')
        file.write('                            self.BoundaryConstraintDefinitions.append(entry)\n')
        file.write('                    \n')

        file.write('                    elif linecell[0] == "BEGIN_PHASE_DISTANCE_CONSTRAINT_BLOCK":\n')
        file.write('                        self.PhaseDistanceConstraintDefinitions = []\n')
        file.write('                        while True:\n')
        file.write('                            entry = inputFile.readline()\n')
        file.write('                            lineNumber += 1\n')
        file.write('                            if "END_PHASE_DISTANCE_CONSTRAINT_BLOCK" in entry:\n')
        file.write('                                break\n')
        file.write('                            self.PhaseDistanceConstraintDefinitions.append(entry)\n')
        file.write('                    \n')

        file.write('                    elif linecell[0] == "BEGIN_TRIALX":\n')
        file.write('                        self.trialX = []\n')
        file.write('                        while True:\n')
        file.write('                            entry = inputFile.readline()\n')
        file.write('                            lineNumber += 1\n')
        file.write('                            if "END_TRIALX" in entry:\n')
        file.write('                                break\n')
        file.write('                            commalinecell = entry.split(\',\')\n')
        file.write('                            self.trialX.append(commalinecell)\n')
        file.write('                    \n')
        
        
        file.write('    #************************************************************************************write\n')
        file.write('    def write(self, optionsFileName, writeAll = True):\n')
        file.write('        with open(optionsFileName, "a+") as optionsFile:\n')
        file.write('            optionsFile.write("\\n")\n')
        file.write('            optionsFile.write("\\n")\n')
        file.write('            optionsFile.write("BEGIN_JOURNEY\\n")\n')
        file.write('            optionsFile.write("\\n")\n')
        file.write('            \n')

        for option in OptionsDefinitions:
            name = option['name']

            defaultValue = ''
            if option['dataType'] == 'std::string':
                defaultValue = '"' + str(option['defaultValue']) + '"'
            else:
                defaultValue = str(option['defaultValue'])

            file.write('            if (self.' + name + ' != ' + defaultValue + ' or writeAll or self.print_this_journey_options_no_matter_what):\n')
            file.write('                optionsFile.write("#' + option['description'] + '\\n")\n')
            if 'std::vector' in option['dataType']:
                elementType = option['dataType'].replace('std::vector<','').replace('>','')
                file.write('                optionsFile.write("' + name + '")\n')
                file.write('                for entry in self.' + name + ':\n')
                if 'bool' in option['dataType']:
                    file.write('                    optionsFile.write(" " + str(int(entry)))\n')
                else:
                    file.write('                    optionsFile.write(" " + str(entry))\n')
                file.write('                optionsFile.write("\\n")\n')
            else:
                if 'bool' in option['dataType']:
                    file.write('                optionsFile.write("' + name + ' " + str(int(self.' + name + ')) + "\\n")\n')
                else:
                    file.write('                optionsFile.write("' + name + ' " + str(self.' + name + ') + "\\n")\n')
            file.write('    \n')        
                
        file.write('            optionsFile.write("\\n")\n') 
        file.write('            optionsFile.write("#Maneuver constraint code\\n")\n') 
        file.write('            optionsFile.write("#Works for absolute and relative epochs and also magnitudes\\n")\n') 
        file.write('            optionsFile.write("BEGIN_MANEUVER_CONSTRAINT_BLOCK\\n")\n') 
        file.write('            for ManeuverConstraintDefinition in self.ManeuverConstraintDefinitions:\n') 
        file.write('                optionsFile.write(ManeuverConstraintDefinition + "\\n")\n') 
        file.write('            optionsFile.write("END_MANEUVER_CONSTRAINT_BLOCK\\n")\n') 
        file.write('            optionsFile.write("\\n")\n') 
        file.write('            \n') 
                
        file.write('            optionsFile.write("\\n")\n') 
        file.write('            optionsFile.write("#Boundary constraint code\\n")\n') 
        file.write('            optionsFile.write("BEGIN_BOUNDARY_CONSTRAINT_BLOCK\\n")\n') 
        file.write('            for BoundaryConstraintDefinition in self.BoundaryConstraintDefinitions:\n') 
        file.write('                optionsFile.write(BoundaryConstraintDefinition + "\\n")\n') 
        file.write('            optionsFile.write("END_BOUNDARY_CONSTRAINT_BLOCK\\n")\n') 
        file.write('            optionsFile.write("\\n")\n') 
        file.write('            \n') 
        
        file.write('            optionsFile.write("\\n")\n') 
        file.write('            optionsFile.write("#Phase distance constraint code\\n")\n') 
        file.write('            optionsFile.write("BEGIN_PHASE_DISTANCE_CONSTRAINT_BLOCK\\n")\n') 
        file.write('            for PhaseDistanceConstraintDefinition in self.PhaseDistanceConstraintDefinitions:\n') 
        file.write('                optionsFile.write(PhaseDistanceConstraintDefinition + "\\n")\n') 
        file.write('            optionsFile.write("END_PHASE_DISTANCE_CONSTRAINT_BLOCK\\n")\n') 
        file.write('            optionsFile.write("\\n")\n') 
        file.write('            \n') 

        file.write('            if len(self.trialX) > 0:\n') 
        file.write('                optionsFile.write("#trial decision vector\\n")\n') 
        file.write('                optionsFile.write("BEGIN_TRIALX\\n")\n') 
        file.write('                for entry in self.trialX:\n') 
        file.write('                    optionsFile.write(entry[0] + "," + \'%17.20f\' % float(entry[1]) + "\\n")\n') 
        file.write('                optionsFile.write("END_TRIALX\\n")\n') 
        file.write('            optionsFile.write("\\n")\n') 
        file.write('            \n') 

        file.write('            optionsFile.write("END_JOURNEY")\n')
        file.write('            \n') 

        file.write('    #************************************************************************************convert decision vector\n')
        file.write('    def ConvertDecisionVector(self, ParallelShootingStateRepresentation, PeriapseBoundaryStateRepresentation):                                                                                            \n')
        file.write('        from StateConverter import StateConverter                                                                                                                                                                                 \n')
        file.write('        myStateConverter = StateConverter()                                                                                                                                                                                 \n')
        file.write('                                                                                                                                                                                                              \n')
        file.write('        stateRepresentationNames = ["Cartesian", "SphericalRADEC", "SphericalAZFPA", "COE", "MEE", "IncomingBplane", "OutgoingBplane", "IncomingBplaneRpTA", "OutgoingBplaneRpTA"]                                                                                          \n')
        file.write('                                                                                                                                                                                                              \n')
        file.write('        mu = 1.0                                                                                                                                                                                              \n')
        file.write('        try:                                                                                                                                                                                                  \n')
        file.write('            import Universe                                                                                                                                                                                   \n')
        file.write('            myUniverse = Universe.Universe(self.universe_folder + "/" + self.journey_central_body + ".emtg_universe")                                                                                         \n')
        file.write('            mu = myUniverse.mu                                                                                                                                                                                \n')
        file.write('        except:                                                                                                                                                                                               \n')
        file.write('            print("Failed to find " + self.universe_folder + "/" + self.journey_central_body + ".emtg_universe" + "  Cannot find appropriate mu for decision vector conversion. Using 1.0. Good luck.")       \n')
        file.write('                                                                                                                                                                                                              \n')
        file.write('        #switch between MGALT/FBLT and PSBI/PSFB                                                                                                                                                              \n')
        file.write('        for entry in self.trialX:                                                                                                                                                                             \n')
        file.write('            if "MGALT" in entry[0] and (self.phase_type == 3):                                                                                                                                                \n')
        file.write('                entry[0] = entry[0].replace("MGALT","FBLT")                                                                                                                                                   \n')
        file.write('            elif "FBLT" in entry[0] and (self.phase_type == 2):                                                                                                                                               \n')  
        file.write('                entry[0] = entry[0].replace("FBLT","MGALT")                                                                                                                                                   \n')
        file.write('            elif "PSBI" in entry[0] and (self.phase_type == 5):                                                                                                                                               \n')  
        file.write('                entry[0] = entry[0].replace("PSBI","PSFB")                                                                                                                                                    \n')
        file.write('            elif "PSFB" in entry[0] and (self.phase_type == 4):                                                                                                                                               \n')  
        file.write('                entry[0] = entry[0].replace("PSFB","PSBI")                                                                                                                                                    \n')
        file.write('                                                                                                                                                                                                              \n')
        file.write('            if "xdot" in entry[0]:                                                                                                                                                                            \n')
        file.write('                entry[0] = entry[0].replace("xdot", "vx")                                                                                                                                                     \n')
        file.write('            elif "ydot" in entry[0]:                                                                                                                                                                          \n')
        file.write('                entry[0] = entry[0].replace("ydot", "vy")                                                                                                                                                     \n')
        file.write('            elif "zdot" in entry[0]:                                                                                                                                                                          \n')
        file.write('                entry[0] = entry[0].replace("zdot", "vz")                                                                                                                                                     \n')
        file.write('                                                                                                                                                                                                              \n')
        file.write('        #old launch to new launch                                                                                                          \n')
        file.write('        if self.departure_class == 3 and self.departure_type == 0: #periapse launch                                                        \n')
        file.write('            import os, sys, inspect                                                                                                        \n')
        file.write('            currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))                                         \n')
        file.write('            sys.path.append(currentdir + "/Converters")                                                                               \n')
        file.write('            import convert_old_PeriapseLaunchOrImpulsiveDeparture_to_PeriapseLaunch as lc                                                  \n')
        file.write('                                                                                                                                           \n')
        file.write('            self = lc.convert_launch(self, mu)                                                                                             \n')
        file.write('                                                                                                                                           \n')
        file.write('        #ParallelShooting                                                                                                                  \n')
        file.write('        self.trialX = myStateConverter.convertDecisionVector(self.trialX,                                                                  \n')
        file.write('                                                            stateRepresentationNames[ParallelShootingStateRepresentation],                 \n')
        file.write('                                                            ["PSBI_Step", "PSFB_Step"],                                                    \n')
        file.write('                                                            mu)                                                                            \n')
        file.write('                                                                                                                                           \n')
        file.write('        #PeriapseBoundary                                                                                                                  \n')
        file.write('        self.trialX = myStateConverter.convertDecisionVector(self.trialX,                                                                  \n')
        file.write('                                                            stateRepresentationNames[PeriapseBoundaryStateRepresentation],                 \n')
        file.write('                                                            ["Periapse"],                                                                  \n')
        file.write('                                                            mu)                                                                            \n')
        file.write('                                                                                                                                           \n')
        file.write('        #FreePointBoundary departure                                                                                                       \n')
        file.write('        self.trialX = myStateConverter.convertDecisionVector(self.trialX,                                                                  \n')
        file.write('                                                            stateRepresentationNames[self.departure_elements_state_representation],        \n')
        file.write('                                                            ["FreePointDirectInsertion", "FreePointFreeDirectDeparture"],                  \n')
        file.write('                                                            mu)                                                                            \n')
        file.write('                                                                                                                                           \n')
        file.write('        #FreePointBoundary arrival                                                                                                         \n')
        file.write('        self.trialX = myStateConverter.convertDecisionVector(self.trialX,                                                                  \n')
        file.write('                                                            stateRepresentationNames[self.arrival_elements_state_representation],          \n')
        file.write('                                                            ["FreePointChemRendezvous", "FreePointIntercept", "FreePointLTRendezvous"],    \n')
        file.write('                                                            mu,                                                                            \n')
        file.write('                                                            exceptions=["Probe"])                                                          \n')
        file.write('                                                                                                                                           \n')
        file.write('        # convert any Probe atmospheric entry states                                                                                       \n')
        file.write('        self.trialX = myStateConverter.convertDecisionVector(self.trialX,                                                                  \n')        
        file.write('                                                            stateRepresentationNames[self.Probe_AEI_elements_state_representation],        \n')          
        file.write('                                                            ["ProbeEntryPhaseProbeAEIFreePointLTRendezvous"],                              \n')
        file.write('                                                            mu)                                                                            \n')
        file.write('                                                                                                                                           \n')
        file.write('        # convert any Probe End states                                                                                                     \n')
        file.write('        self.trialX = myStateConverter.convertDecisionVector(self.trialX,                                                                  \n')        
        file.write('                                                            stateRepresentationNames[self.Probe_End_elements_state_representation],        \n')          
        file.write('                                                            ["ProbeEntryPhaseProbeEndOfMissionFreePointLTRendezvous"],                     \n')
        file.write('                                                            mu)                                                                            \n')
        file.write('                                                                                                                                           \n')
        file.write('        return                                                                                                                             \n')
        file.write('                                                                                                                                           \n')
        file.write('                                                                                                                                           \n')
        file.write('                                                                                                                                           \n')
        file.write('    #************************************************************************************getDecisionVariable()                             \n')
        file.write('    def getDecisionVariable(self, variableDefinition):                                                                                     \n')
        file.write('        for entry in self.trialX:                                                                                                          \n')
        file.write('            if entry[0].strip() == variableDefinition:                                                                                     \n')
        file.write('                return entry[1]                                                                                                            \n')
        file.write('                                                                                                                                           \n')
        file.write('        #if you made it here, something went wrong.                                                                                        \n')
        file.write('                                                                                                                                           \n')
        file.write('        raise Exception("Variable \'" + variableDefinition + "\' not found.")                                                              \n')
        file.write('                                                                                                                                           \n')
        file.write('                                                                                                                                           \n')
        file.write('    #************************************************************************************setDecisionVariable()                             \n')
        file.write('    def setDecisionVariable(self, variableDefinition, value):                                                                              \n')
        file.write('        for entryIndex in range(0, len(self.trialX)):                                                                                      \n')
        file.write('            if self.trialX[entryIndex][0].strip() == variableDefinition:                                                                   \n')
        file.write('                self.trialX[entryIndex][1] = value                                                                                         \n')
        file.write('                return                                                                                                                     \n')
        file.write('                                                                                                                                           \n')
        file.write('        #if you made it here, something went wrong.                                                                                        \n')
        file.write('                                                                                                                                           \n')
        file.write('        raise Exception("Variable \'" + variableDefinition + "\' not found.")                                                              \n')
