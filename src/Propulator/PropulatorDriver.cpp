// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Other Rights Reserved.

// Copyright (c) 2024 The Regents of the University of Colorado.
// All Other Rights Reserved.

// Licensed under the NASA Open Source License (the "License"); 
// You may not use this file except in compliance with the License. 
// You may obtain a copy of the License at:
// https://opensource.org/licenses/NASA-1.3
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
// express or implied.   See the License for the specific language
// governing permissions and limitations under the License.

#include "Propulator.h"

#include "EMTG_Matrix.h"

#include <iostream>

int main(int argc, char* argv[])
{
    try
    {
        std::cout << "Propulation starting" << std::endl;
		
        //parse the options file
        std::string options_file_name;
        std::string input_states_file_name;
        std::string output_states_file_name;
        std::string output_STM_file_name;
        bool printSTM = false;
        if (argc >= 4)
        {
            options_file_name.assign(argv[1]);
            input_states_file_name.assign(argv[2]);
            output_states_file_name.assign(argv[3]);

            if (argc == 5)
            {
                output_STM_file_name.assign(argv[4]);
                printSTM = true;
            }
        }
        else
        {
            std::cout << "Syntax is PropagatorDriver OptionsFilePath InputStatesFilePath OutputStatesFilePath [OutputSTMFilePath]" << std::endl;
            std::cout << "input columns are:" << std::endl;
            std::cout << "#Gregorian Date String, journey, x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg], propagation time [sec]" << std::endl;
            std::cout << "optional columns after these are:" << std::endl;
            std::cout << "#duty cycle, control x, control y, control z" << std::endl;
            std::cout << "#bsp file name" << std::endl;
            return 0;
        }

        //parse the input file
        std::vector< EMTG::math::Matrix<double> > InputStateLines;
        std::vector< EMTG::math::Matrix<double> > ControlVectorLines;
        std::vector<double> PropagationTimes;
        std::string inputLine;
        std::ifstream inputfile(input_states_file_name);
        std::vector<double> dutyCycle;
        std::vector<int> journeyIndices;
		std::vector<std::string> outputbsps;
		std::string outputbsp;
		std::string blank = "None";
		std::string extension = ".bsp";
		std::vector<int> journeyPtrs;
		
        while (EMTG::file_utilities::safeGetline(inputfile, inputLine))
        {
            if (inputLine.size() > 0)
            {
                if (!(inputLine.front() == *"#"))
                {
                    std::vector<std::string> linecell;
                    boost::split(linecell, inputLine, boost::is_any_of(","), boost::token_compress_on);

                    if (linecell.size() >= 10)
                    {
						bool found = false;
						int ji = std::stoi(linecell[1]);
						for (int jdx = 0; jdx < journeyIndices.size(); jdx++)
							if (journeyIndices[jdx] == ji)
							{
								found = true;
								journeyPtrs.push_back(jdx);
								break;
							}
						if (!found)
						{
	                        journeyIndices.push_back(ji);
							journeyPtrs.push_back(journeyIndices.size()-1);
						}
					}
				}
			}
		}
	    
		inputfile.clear();
	    inputfile.seekg(0, inputfile.beg);
		
        //instantiate a propulator
        EMTG::Propulator::Propulator myPropulator(options_file_name,journeyIndices);
		
        while (EMTG::file_utilities::safeGetline(inputfile, inputLine))
        {
            if (inputLine.size() > 0)
            {
                if (!(inputLine.front() == *"#"))
                {
                    std::vector<std::string> linecell;
                    boost::split(linecell, inputLine, boost::is_any_of(","), boost::token_compress_on);
					
                    if (linecell.size() >= 10)
                    {
                        //epoch, x, y, z, xdot, ydot, zdot, mass, u_x, u_y, u_z, propagationTime
                        size_t cellIndex = 0;
                        EMTG::math::Matrix<double> inputStateLine(10, 1, 0.0);
                        EMTG::math::Matrix<double> inputControlLine(3, 1, 0.0);

                        //epoch
						SpiceDouble ets;
						SpiceChar epochstring[40];
						strcpy(epochstring,linecell[cellIndex++].c_str());
						str2et_c(epochstring,&ets);
						inputStateLine(7) = ets + 51544.5 * 86400.0;

                        //journey, duty cycle, and step size
						cellIndex++; // Cycle index ptr to skip over the journey index

                        //state
                        for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                            inputStateLine(stateIndex) = std::stod(linecell[cellIndex++]);

                        InputStateLines.push_back(inputStateLine);
						
						// Time
                        PropagationTimes.push_back(std::stod(linecell[cellIndex++]));

						// control
                        if (linecell.size() == 14)
                        {
	                        dutyCycle.push_back(std::stod(linecell[cellIndex++]));
                            for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                                inputControlLine(controlIndex) = std::stod(linecell[cellIndex++]);

                            ControlVectorLines.push_back(inputControlLine);
                        }
                        else
                            ControlVectorLines.push_back(EMTG::math::Matrix<double>(1, 1, 0.0));
						
						// output bsp
						if (linecell.size() == 15 || linecell.size() == 11)
						{
							outputbsp = linecell[cellIndex++];
						
							if (!std::equal(extension.rbegin(), extension.rend(), outputbsp.rbegin()))
							{
					            std::cout << "Number of inputs suggesteded a bsp file. However entry did not end with .bsp" << std::endl;
								std::cout << "Entry: " << outputbsp << std::endl;
								std::cout << "Line: " << inputLine << std::endl;
					            return 0;
							}
							
							outputbsps.push_back(outputbsp);
						}
						else
							outputbsps.push_back(blank);
						
                    }
                    else
                    {
                        std::cout << "invalid line, '" << inputLine << "'" << std::endl;
                    }
                }
                else
                {
                    std::cout << "invalid line, '" << inputLine << "'" << std::endl;
                }
            }
            // else
            // {
            //     std::cout << "invalid line, '" << inputLine << "'" << std::endl;
            // }
        }//end loop over lines in file
		
        //run the propulator
        EMTG::math::Matrix<double> OutputState(10, 1, 0.0);
        std::ofstream outputfile(output_states_file_name, std::ios::trunc);
        outputfile << "#States output from propulator" << std::endl;

        outputfile << std::endl;
        outputfile << "Gregorian Date String [TDB], x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg]" << std::endl;

        std::ofstream STMfile;
        if (printSTM)
        {
            STMfile = std::ofstream(output_STM_file_name, std::ios::trunc);
            STMfile << "#STM output from propulator" << std::endl;
            STMfile << "row ordered, i.e. each line goes [first row, second row, ...]" << std::endl;
            STMfile << std::endl;
        }

        for (size_t caseIndex = 0; caseIndex < InputStateLines.size(); ++caseIndex)
        {
            myPropulator.setJourneyIndex(journeyPtrs[caseIndex]);
			
			myPropulator.setStoreStateHistory(true);

            if (ControlVectorLines[caseIndex].get_n() > 1) //propagate with control
            {
                if (fabs(dutyCycle[caseIndex] < 1.0e-10))
                {
                    myPropulator.propulate(InputStateLines[caseIndex],
                        PropagationTimes[caseIndex],
                        OutputState,
                        printSTM);
                }
                else
                {
                    myPropulator.propulate(InputStateLines[caseIndex],
                        ControlVectorLines[caseIndex],
                        dutyCycle[caseIndex],
                        PropagationTimes[caseIndex],
                        OutputState,
                        printSTM);
                }
            }
            else //propagate without control
            {
                myPropulator.propulate(InputStateLines[caseIndex],
                    PropagationTimes[caseIndex],
                    OutputState,
                    printSTM);
            }
            outputfile.precision(14);
			
	        SpiceChar epochstring[32];
	        timout_c(OutputState(7) _GETVALUE - (51544.5 * 86400.0), "YYYY Mon DD ::TDB HR:MN:SC.######", 32, epochstring);
	        outputfile << epochstring;
            // outputfile << OutputState(7) / 86400.0 + 2400000.5;
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                outputfile << ", " << OutputState(stateIndex);
            outputfile << std::endl;

            if (printSTM)
            {
                EMTG::math::Matrix<double> STM = myPropulator.getSTM();
                for (size_t i = 0; i < STM.get_n(); ++i)
                {
                    for (size_t j = 0; j < STM.get_m(); ++j)
                    {
                        STMfile << STM(i, j);

                        if (!(i == STM.get_n() - 1 && j == STM.get_m() - 1)) //if we're not the last entry, write a comma
                            STMfile << ", ";
                    }
                }
                STMfile << std::endl;
            }
			
			if (!std::equal(blank.rbegin(), blank.rend(), outputbsps[caseIndex].rbegin()))
				myPropulator.generateSPK(outputbsps[caseIndex]);
			
			myPropulator.clearStateHistory();
			
        }//end loop over cases


        std::cout << "Propulation complete." << std::endl;
    }
    catch (std::exception const& exception)
    {
        std::cout << "Propulator failed with error:" << std::endl;
        std::cout << exception.what() << std::endl;
    }
    return 0;
}
