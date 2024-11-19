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

//Propulator
//contains a vector of PropulatorJourney objects

#include "Propulator.h"

#include "file_utilities.h"
#include "SpacecraftOptionsFactory.h"

#include "SpiceUsr.h"

namespace EMTG
{
    namespace Propulator
    {
        Propulator::Propulator() :
            journeyIndex(0)
        {}//empty default constructor except for initialization list

        Propulator::Propulator(const std::string& options_file_name,std::vector<int> journeyIndex)
            : Propulator()
        {
            this->initialize(options_file_name,journeyIndex);
        }//end constructor

        void Propulator::initialize(const std::string& options_file_name,std::vector<int> journeyIndices)
        {
            std::cout << "Loading " << options_file_name << std::endl;
            this->myOptions = missionoptions(options_file_name);

            //set up Universe stuff
            {
                //load all ephemeris data if using SPICE
                std::vector<::boost::filesystem::path> SPICE_files_initial;
                std::vector<::boost::filesystem::path> SPICE_files_not_required;
                this->SPICE_files_required.clear();
                std::vector<int> SPICE_bodies_required;
                std::string filestring;
                if (this->myOptions.ephemeris_source >= 1)
                {
                    //load all BSP files
                    EMTG::file_utilities::get_all_files_with_extension(::boost::filesystem::path(this->myOptions.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files_initial);

                    for (size_t k = 0; k < SPICE_files_initial.size(); ++k)
                    {
                        filestring = this->myOptions.universe_folder + "/ephemeris_files/" + SPICE_files_initial[k].string();
                        furnsh_c(filestring.c_str());
                        std::cout << filestring << std::endl;
                    }

                    //disable quit-on-SPICE-error so that we can see what happens if the leap second and/or frame kernels don't load properly
                    erract_c((SpiceChar*)"SET", 100, (SpiceChar*)"RETURN");

                    //SPICE reference frame kernel
                    std::string leapsecondstring = this->myOptions.universe_folder + "/ephemeris_files/" + this->myOptions.SPICE_leap_seconds_kernel;
                    std::string referenceframestring = this->myOptions.universe_folder + "/ephemeris_files/" + this->myOptions.SPICE_reference_frame_kernel;
                    furnsh_c(leapsecondstring.c_str());
                    furnsh_c(referenceframestring.c_str());

                    //disable SPICE error printing. This is because we can, and will often, go off the edge of an ephemeris file.
                    errprt_c((SpiceChar*)"SET", 100, (SpiceChar*)"NONE");

                    SPICE_files_required = SPICE_files_initial;

                    std::cout << "Completed loading SPICE kernels." << std::endl;
                }
#ifdef SPLINE_EPHEM
                std::vector< std::tuple<int, int, int, double> > SplineUniverse_keyList;

                this->mySplineUniverse = new SplineEphem::universe(SplineUniverse_keyList);
#endif

                for (size_t j = 0; j < journeyIndices.size(); ++j)
                {
					int journeyIdx = journeyIndices[j];
#ifdef SPLINE_EPHEM
                    TheUniverse.push_back(EMTG::Astrodynamics::universe(journeyIdx, this->myOptions.universe_folder + "//" + this->myOptions.Journeys[journeyIdx].journey_central_body + ".emtg_universe", this->myOptions, this->mySplineUniverse));
#else
                    TheUniverse.push_back(EMTG::Astrodynamics::universe(journeyIdx, options.universe_folder + "//" + options.Journeys[journeyIdx].journey_central_body + ".emtg_universe", options));
#endif

                    if (TheUniverse[j].TU > this->myOptions.TU)
                        this->myOptions.TU = TheUniverse[j].TU;
                }
				
                //now that we have a Universe vector, we can use it to populate the SplineEphem::universe
                        //add every body that will we used in the mission to the SplineUniverse
#ifdef SPLINE_EPHEM
                SplineUniverse_keyList.clear();
                try
                {
                    for (size_t j = 0; j < journeyIndices.size(); ++j)
                    {
						int journeyIdx = journeyIndices[j];
						
                        std::vector<int> body_index_array;

                        //first boundary point
                        if (this->myOptions.Journeys[journeyIdx].destination_list[0] > 0)
                            body_index_array.push_back(this->myOptions.Journeys[journeyIdx].destination_list[0] - 1);

                        //last boundary point
                        if (this->myOptions.Journeys[journeyIdx].destination_list[1] > 0)
                            body_index_array.push_back(this->myOptions.Journeys[journeyIdx].destination_list[1] - 1);

                        //sequence
                        for (int body : this->myOptions.Journeys[journeyIdx].sequence)
                            if (body > 0)
                                body_index_array.push_back(body - 1);

                        //perturbation list
                        if (this->myOptions.perturb_thirdbody)
                        {
                            for (size_t b = 0; b < TheUniverse[j].perturbation_menu.size(); ++b)
                                body_index_array.push_back(TheUniverse[j].perturbation_menu[b]);
                        }

                        //distance constraint list
                        for (std::string& constraint : this->myOptions.Journeys[journeyIdx].PhaseDistanceConstraintDefinitions)
                        {
                            std::vector<std::string> ConstraintDefinitionCell;
                            boost::split(ConstraintDefinitionCell,
                                constraint,
                                boost::is_any_of("_"),
                                boost::token_compress_on);

                            if (boost::to_lower_copy(ConstraintDefinitionCell[1]) != "cb")
                            {
                                int bodyIndex = std::stoi(ConstraintDefinitionCell[1]) - 1;

                                body_index_array.push_back(bodyIndex);
                            }
                        }

                        for (size_t b = 0; b < body_index_array.size(); ++b)
                        {
                            //do we already have this body?
                            bool body_in_keylist = false;
                            for (size_t k = 0; k < SplineUniverse_keyList.size(); ++k)
                            {
                                if (std::get<0>(SplineUniverse_keyList[k]) == TheUniverse[j].bodies[body_index_array[b]].spice_ID
                                    && std::get<1>(SplineUniverse_keyList[k]) == TheUniverse[j].central_body_SPICE_ID)
                                {
                                    body_in_keylist = true;
                                    break;
                                }
                            }

                            if (!body_in_keylist && body_index_array[b] >= 0)
                            {
                                SplineUniverse_keyList.push_back(std::make_tuple(
                                    TheUniverse[j].bodies[body_index_array[b]].spice_ID,
                                    TheUniverse[j].central_body_SPICE_ID,
                                    this->myOptions.SplineEphem_points_per_period,
                                    TheUniverse[j].mu));
                            }
                        }//end loop over bodies in the universe

                        //is this universe's central body the sun? If not, let's add this body with respect to the sun. Let's add extra ephemeris points, too.
                        if (!(TheUniverse[j].central_body_SPICE_ID == 10))
                        {
                            bool body_in_keylist = false;
                            for (size_t k = 0; k < SplineUniverse_keyList.size(); ++k)
                            {
                                if (std::get<0>(SplineUniverse_keyList[k]) == TheUniverse[j].central_body_SPICE_ID
                                    && std::get<1>(SplineUniverse_keyList[k]) == 10)
                                {
                                    body_in_keylist = true;
                                    break;
                                }
                            }

                            if (!body_in_keylist)
                            {
                                SplineUniverse_keyList.push_back(std::make_tuple(
                                    TheUniverse[j].central_body_SPICE_ID,
                                    10,
                                    this->myOptions.SplineEphem_non_central_body_sun_points_per_period,
                                    1.32712440018e+11));
                            }
                        }
                    }

                    double earliestPossibleEpoch = this->myOptions.earliestPossibleEpoch * 86400.0;
                    double latestPossibleEpoch = this->myOptions.latestPossibleEpoch * 86400.0;

                    if (this->myOptions.SplineEphem_truncate_ephemeris_at_maximum_mission_epoch
                        && latestPossibleEpoch < (this->myOptions.launch_window_open_date + this->myOptions.Journeys.front().wait_time_bounds[1] + this->myOptions.total_flight_time_bounds[1]))
                        latestPossibleEpoch = this->myOptions.launch_window_open_date + this->myOptions.Journeys.front().wait_time_bounds[1] + this->myOptions.total_flight_time_bounds[1] * 86400.0;
                    /*if (earliest_possible_epoch > this->myOptions.earliestPossibleEpoch * 86400.0)
                        earliest_possible_epoch = this->myOptions.earliestPossibleEpoch * 86400.0;*/
                    this->mySplineUniverse->reinitialize(SplineUniverse_keyList,
                        earliestPossibleEpoch - 10.0 * 86400.0,
                        latestPossibleEpoch + 10.0 * 86400.0);
                }
                catch (std::exception &myError)
                {
                    std::cout << "Failure while configuring SplineEphem." << std::endl;
                    std::cout << myError.what() << std::endl;
                    std::cout << "Submit this error message to the EMTG development team, along with your .emtgopt, .emtg_universe file(s), your hardware model files, any relevant ephemeris files, and which branch you are using. This information will allow us to properly help you." << std::endl;
#ifndef BACKGROUND_MODE //macro overrides if statement
                    std::cout << "Press enter to close window." << std::endl;
                    std::cin.ignore();
#endif
                    throw;
                }
#endif
            }//end setting up universe stuff

            //spacecraft
            EMTG::HardwareModels::SpacecraftOptions mySpacecraftOptions = EMTG::HardwareModels::CreateSpacecraftOptions(this->myOptions);
            this->mySpacecraft = EMTG::HardwareModels::Spacecraft(mySpacecraftOptions);

            //set up PropulatorJourneys
            size_t stageIndex = 0;
            this->myPropulatorJourneys.clear();

            for (int j = 0; j < journeyIndices.size(); ++j)
            {
				int journeyIndex = journeyIndices[j];
                this->myPropulatorJourneys.push_back(new PropulatorJourney(journeyIndex,
                    stageIndex,
                    &this->myOptions,
                    &this->TheUniverse[j],
                    &this->mySpacecraft));
			}			
        }//end initialize()

        Propulator::~Propulator()
        {
            if (this->myOptions.ephemeris_source >= 1)
            {
                for (size_t k = 0; k < this->SPICE_files_required.size(); ++k)
                {
                    std::string filestring = this->myOptions.universe_folder + "ephemeris_files/" + this->SPICE_files_required[k].string();
                    unload_c(filestring.c_str());
                }

                unload_c((this->myOptions.universe_folder + "ephemeris_files/" + this->myOptions.SPICE_leap_seconds_kernel).c_str());
                unload_c((this->myOptions.universe_folder + "ephemeris_files/" + this->myOptions.SPICE_reference_frame_kernel).c_str());
            }
#ifdef SPLINE_EPHEM
            delete this->mySplineUniverse;
#endif
        }//end destructor

        void Propulator::propulate(const math::Matrix<double>& InputState,
            const math::Matrix<double>& ControlVector,
            const double& DutyCycle,
            const double& propagationTime,
            math::Matrix<double>& OutputState,
            const bool& needSTM)
        {
            this->myPropulatorJourneys[this->journeyIndex].propulate(InputState, 
                ControlVector,
                DutyCycle, 
                propagationTime, 
                OutputState,
                needSTM);
        }//end propulate()

        void Propulator::propulate(const math::Matrix<double>& InputState,
            const double& propagationTime,
            math::Matrix<double>& OutputState,
            const bool& needSTM)
        {
            this->myPropulatorJourneys[this->journeyIndex].propulate(InputState,
                propagationTime,
                OutputState,
                needSTM);
        }//end propulate()

        std::vector< std::vector<double> >  Propulator::getStateHistory()
        {
            return this->myPropulatorJourneys[this->journeyIndex].getStateHistory();
        }//end getStateHistory()
		
        void  Propulator::setStoreStateHistory(bool switch_in)
        {
            this->myPropulatorJourneys[this->journeyIndex].setStoreStateHistory(switch_in);
        }//end getStateHistory()
		
        bool  Propulator::getStoreStateHistory()
        {
            return this->myPropulatorJourneys[this->journeyIndex].getStoreStateHistory();
        }//end getStateHistory()
		
        void  Propulator::clearStateHistory()
        {
			this->myPropulatorJourneys[this->journeyIndex].clearStateHistory();
		}
		
		void Propulator::generateSPK(std::string filename)
		{
			std::vector< std::vector<double> > stateHistory = this->getStateHistory();
			
	        std::ofstream bspfile(filename + ".ephem", std::ios::trunc);

            bspfile.precision(14);
	        // bspfile << "#time,x,y,z,vx,vy,vz" << std::endl;

			for (int idx = 0; idx < stateHistory.size(); idx++)
			{
		        SpiceChar epochstring[32];
		        timout_c(stateHistory[idx][7] - (51544.5 * 86400.0), "YYYY Mon DD ::TDB HR:MN:SC.######", 32, epochstring);
				bspfile << epochstring << ",";
				for (int sdx = 0; sdx < 6; sdx++)
				{
					bspfile << stateHistory[idx][sdx];
					if (sdx != 5)
						bspfile << ",";
				}
				bspfile << std::endl;
			}
			
	        bspfile.close();

	        //Step 2: write .cmd file
	        std::string cmdstring = filename + ".cmd";
	        std::ofstream cmd_file(cmdstring, std::ios::out | std::ios::trunc);

	        cmd_file << "\\begindata" << std::endl;
	        cmd_file << "INPUT_DATA_TYPE = 'STATES'" << std::endl;
	        cmd_file << "OUTPUT_SPK_TYPE = 9" << std::endl;
	        cmd_file << "OBJECT_ID = " << this->myOptions.spacecraft_SPICE_ID << std::endl;
	        cmd_file << "OBJECT_NAME = '" << this->myOptions.mission_name << "'" << std::endl;
	        cmd_file << "CENTER_ID = " << this->myPropulatorJourneys[this->journeyIndex].getUniverse()->central_body_SPICE_ID << std::endl;
	        cmd_file << "CENTER_NAME = 'CENTER'" << std::endl;
	        cmd_file << "REF_FRAME_NAME = 'J2000'" << std::endl;
	        cmd_file << "PRODUCER_ID = 'EMTGv9'" << std::endl;
	        cmd_file << "DATA_ORDER = 'EPOCH X Y Z VX VY VZ'" << std::endl;
	        cmd_file << "INPUT_DATA_UNITS = ('ANGLES=DEGREES' 'DISTANCES=km')" << std::endl;
	        cmd_file << "TIME_WRAPPER = '# TDB'" << std::endl;
	        cmd_file << "DATA_DELIMITER = ','" << std::endl;
	        cmd_file << "LINES_PER_RECORD = 1" << std::endl;
	        cmd_file << "IGNORE_FIRST_LINE = 0" << std::endl;
	        cmd_file << "LEAPSECONDS_FILE = '" << boost::replace_all_copy(this->myOptions.universe_folder, "\\", "/") + "/ephemeris_files/" + this->myOptions.SPICE_leap_seconds_kernel << "'" << std::endl;
	        cmd_file << "POLYNOM_DEGREE = 3" << std::endl;
	        cmd_file << "SEGMENT_ID = 'SPK_STATES_09'" << std::endl;
	        cmd_file << "\\begintext" << std::endl;

	        cmd_file.close();
			
			std::stringstream cmd; 
			cmd << boost::replace_all_copy(this->myOptions.spice_utilities_path, "\\", "/") << "/mkspk";
			cmd << " -setup " << filename << ".cmd";
			cmd << " -input " << filename << ".ephem";
			cmd << " -output " << filename << "; rm " << filename << ".cmd " << filename << ".ephem"; 
			std::system(cmd.str().c_str());
		}
		
#ifdef PROPULATOR_PYTHON_INTERFACE
        boost::python::list Propulator::propulateThrust(const boost::python::list& InputState,
            const boost::python::list& ControlVector,
            const double& DutyCycle,
            const double& propagationTime,
            const bool& needSTM)
        {
            //Step 1: wrap inputs
            size_t nstates = boost::python::len(InputState);
            size_t ncontrols = boost::python::len(ControlVector);
            EMTG::math::Matrix<double> inputStateLine(10, 1, 0.0);
            EMTG::math::Matrix<double> inputControlLine(3, 1, 0.0);

            //position, velocity, and mass
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                inputStateLine(stateIndex) = boost::python::extract<double>(InputState[stateIndex]);

            //epoch
            SpiceChar epochstring[40];
			std::string epochstring_in = boost::python::extract<std::string>(InputState[7]);
			SpiceDouble ets;
			strcpy(epochstring,epochstring_in.c_str());
			str2et_c(epochstring,&ets);
			inputStateLine(7) = ets + 51544.5 * 86400.0;

            //leave the tanks zeroed

            //control
            for (size_t controlIndex = 0; controlIndex < ncontrols; ++controlIndex)
                inputControlLine(controlIndex) = boost::python::extract<double>(ControlVector[controlIndex]);

			
			this->clearStateHistory();
			this->setStoreStateHistory(true);
			
            //Step 2: call propulator
            EMTG::math::Matrix<double> outputStateLine(10, 1, 0.0);
            this->propulate(inputStateLine, 
                inputControlLine,
                DutyCycle, 
                propagationTime,
                outputStateLine,
                needSTM);

            //Step 3: wrap outputs
            boost::python::list outputState;
            for (size_t stateIndex = 0; stateIndex < nstates; ++stateIndex)
                outputState.append(outputStateLine(stateIndex));

            //Step 4: return the output state
            return outputState;
        }//end python propulateThrust()

        boost::python::list Propulator::propulateCoast(const boost::python::list& InputState,
            const double& propagationTime,
            const bool& needSTM)
        {
            //Step 1: wrap inputs
            size_t nstates = boost::python::len(InputState);
            EMTG::math::Matrix<double> inputStateLine(10, 1, 0.0);

            //position, velocity, and mass
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                inputStateLine(stateIndex) = boost::python::extract<double>(InputState[stateIndex]);

            //epoch
            SpiceChar epochstring[40];
			std::string epochstring_in = boost::python::extract<std::string>(InputState[7]);
			SpiceDouble ets;
			strcpy(epochstring,epochstring_in.c_str());
			str2et_c(epochstring,&ets);
			inputStateLine(7) = ets + 51544.5 * 86400.0;

            //leave the tanks zeroed
			
			this->clearStateHistory();
			this->setStoreStateHistory(true);
			
            //Step 2: call propulator
            EMTG::math::Matrix<double> outputStateLine(10, 1, 0.0);
            this->propulate(inputStateLine,
                propagationTime,
                outputStateLine,
                needSTM);

            //Step 3: wrap outputs
            boost::python::list outputState;
            for (size_t stateIndex = 0; stateIndex < nstates; ++stateIndex)
                outputState.append(outputStateLine(stateIndex));

            //Step 4: return the output state
            return outputState;
        }//end python propulateCoast()
        
        boost::python::list Propulator::getpythonSTM()
        {
            boost::python::list STMlist;

            math::Matrix<double> STM = this->getSTM();
            size_t n = STM.get_n();

            for (size_t i = 0; i < n; ++i)
            {
                boost::python::list STMrow;
                for (size_t j = 0; j < n; ++j)
                {
                    STMrow.append(STM(i, j));
                }

                STMlist.append(STMrow);
            }

            return STMlist;
        }//end getpythonSTM

        Propulator::Propulator(const std::string& options_file_name,boost::python::list journeyList)
            : Propulator()
        {
			std::vector<int> journeyIndex;
			for (int i = 0; i < len(journeyList); ++i)
				journeyIndex.push_back(boost::python::extract<int>(journeyList[i]));
            this->initialize(options_file_name,journeyIndex);
        }//end constructor

		BOOST_PYTHON_MODULE(Propulator)
		{
			boost::python::class_<Propulator>("Propulator", boost::python::init<std::string,boost::python::list>())
											  .def("propulateThrust", &Propulator::propulateThrust)
											  .def("propulateCoast", &Propulator::propulateCoast)
											  .def("getSTM", &Propulator::getpythonSTM)
											  .def("setJourneyIndex", &Propulator::setJourneyIndex)
                                              .def("setStepSize", &Propulator::setStepSize)
										      .def("generateSPK",&Propulator::generateSPK);

		}
		
#endif
    }//close namespace Propulator
}//close namespace EMTG
