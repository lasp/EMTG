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

#include "EccentricityConstraint.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

			EccentricityConstraint::EccentricityConstraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                BoundaryEventBase* myBoundaryEvent,
                const std::string& constraintDefinition) :
				OrbitElementConstraintBase::OrbitElementConstraintBase(name,
                    journeyIndex,
                    phaseIndex,
                    stageIndex,
                    Universe,
                    mySpacecraft,
                    myOptions,
                    myBoundaryEvent,
                    constraintDefinition)
            {
				// indicate to the boundary event that it should compute the orbital elements of the spacecraft
				this->myBoundaryEvent->setComputeOrbitElements(true);

                this->myReferenceFrame = ReferenceFrame::ICRF;

                this->myBoundaryEvent->add_orbit_element_reference_frame(this->myReferenceFrame);
            }

			void EccentricityConstraint::calcbounds()
			{
				// Step 1: parse the constraint definition
				// constraintDefinition is like j0p3_arrival_ecc_lb_ub
				std::vector<std::string> ConstraintDefinitionCell;
				std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
				boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

				// Step 2: create the constraint
				// figure out the lower and upper bounds - depends on units
				this->Flowerbounds->push_back(std::stod(ConstraintDefinitionCell[3]));
				this->Fupperbounds->push_back(std::stod(ConstraintDefinitionCell[4]));
				this->Fdescriptions->push_back(prefix + " " + ConstraintDefinitionCell[1] + " " + ConstraintDefinitionCell[2]);

				// Step 3: sparsity pattern
				// derivatives with respect to anything influencing the boundary event's right-hand six-state vector
				std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
				std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

				for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
				{
					// POSITION
					// non-time variables
					std::vector<size_t> state_dIndex_ecc_position_wrt_StateAfterEvent;
					std::vector<size_t> state_Gindex_ecc_position_wrt_StateAfterEvent_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex)
						{
							state_dIndex_ecc_position_wrt_StateAfterEvent.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
								state_Gindex_ecc_position_wrt_StateAfterEvent_variables);
						}
					}
					this->dIndex_ecc_position_wrt_StateAfterEvent.push_back(state_dIndex_ecc_position_wrt_StateAfterEvent);
					this->Gindex_ecc_position_wrt_StateAfterEvent_variables.push_back(state_Gindex_ecc_position_wrt_StateAfterEvent_variables);

					// time variables
					std::vector<size_t> state_dIndex_ecc_position_wrt_StateAfterEvent_wrt_Time;
					std::vector<size_t> state_Gindex_ecc_position_wrt_StateAfterEvent_time_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex)
						{
							state_dIndex_ecc_position_wrt_StateAfterEvent_wrt_Time.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
								state_Gindex_ecc_position_wrt_StateAfterEvent_time_variables);
						}
					}
					this->dIndex_ecc_position_wrt_StateAfterEvent_wrt_Time.push_back(state_dIndex_ecc_position_wrt_StateAfterEvent_wrt_Time);
					this->Gindex_ecc_position_wrt_StateAfterEvent_time_variables.push_back(state_Gindex_ecc_position_wrt_StateAfterEvent_time_variables);

					// VELOCITY
					// non-time variables
					std::vector<size_t> state_dIndex_ecc_velocity_with_respect_to_StateAfterEvent;
					std::vector<size_t> state_Gindex_ecc_velocity_wrt_StateAfterEvent_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex + 3)
						{
							state_dIndex_ecc_velocity_with_respect_to_StateAfterEvent.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
								state_Gindex_ecc_velocity_wrt_StateAfterEvent_variables);
						}
					}
					this->dIndex_ecc_velocity_wrt_StateAfterEvent.push_back(state_dIndex_ecc_velocity_with_respect_to_StateAfterEvent);
					this->Gindex_ecc_velocity_wrt_StateAfterEvent_variables.push_back(state_Gindex_ecc_velocity_wrt_StateAfterEvent_variables);

					// time variables
					std::vector<size_t> state_dIndex_ecc_velocity_wrt_StateAfterEvent_wrt_Time;
					std::vector<size_t> state_Gindex_ecc_velocity_wrt_StateAfterEvent_time_variables;
					for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
					{
						if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex + 3)
						{
							state_dIndex_ecc_velocity_wrt_StateAfterEvent_wrt_Time.push_back(dIndex);

							this->create_sparsity_entry(this->Fdescriptions->size() - 1,
								std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
								state_Gindex_ecc_velocity_wrt_StateAfterEvent_time_variables);
						}
					}
					this->dIndex_ecc_velocity_wrt_StateAfterEvent_wrt_Time.push_back(state_dIndex_ecc_velocity_wrt_StateAfterEvent_wrt_Time);
					this->Gindex_ecc_velocity_wrt_StateAfterEvent_time_variables.push_back(state_Gindex_ecc_velocity_wrt_StateAfterEvent_time_variables);
				}

			}// end calcbounds()

			void EccentricityConstraint::process_constraint(const std::vector<doubleType>& X,
				size_t& Xindex,
				std::vector<doubleType>& F,
				size_t& Findex,
				std::vector<double>& G,
				const bool& needG)
			{
				// Step 1: get the orbital elements from the boundary event
				this->boundary_orbit_elements = this->myBoundaryEvent->get_orbit_elements_after_event(ReferenceFrame::ICRF);
				this->eccentricity = this->boundary_orbit_elements(1);

				// Step 2: apply the constraint
				F[Findex++] = this->boundary_orbit_elements(1);

				// Step 3: derivatives
				if (needG)
				{

					//Step 4.1: for safety, start by zeroing out all relevant G entries
					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						// non-time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_ecc_position_wrt_StateAfterEvent_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_ecc_position_wrt_StateAfterEvent_variables[stateIndex][entryIndex]] = 0.0;
						}

						// time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_ecc_position_wrt_StateAfterEvent_time_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_ecc_position_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex]] = 0.0;
						}

						// VELOCITY
						// non-time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_ecc_velocity_wrt_StateAfterEvent_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_ecc_velocity_wrt_StateAfterEvent_variables[stateIndex][entryIndex]] = 0.0;
						}

						// time
						for (size_t entryIndex = 0; entryIndex < this->Gindex_ecc_velocity_wrt_StateAfterEvent_time_variables[stateIndex].size(); ++entryIndex)
						{
							G[this->Gindex_ecc_velocity_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex]] = 0.0;
						}
					}

					// Compute the actual Jacobian entries
					// we'll need the partials of the orbital elements w.r.t. the boundary Cartesian state
					this->boundary_orbit_elements_Jacobian = this->myBoundaryEvent->get_orbit_element_Jacobian_after_event(ReferenceFrame::ICRF);

					// first for non-time decision variables
					std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();

					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						for (size_t entryIndex = 0; entryIndex < this->dIndex_ecc_position_wrt_StateAfterEvent[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_ecc_position_wrt_StateAfterEvent[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_ecc_position_wrt_StateAfterEvent_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = this->boundary_orbit_elements_Jacobian(1, stateIndex) _GETVALUE * std::get<2>(Derivatives_of_StateAfterEvent[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry;
						}

						// VELOCITY
						for (size_t entryIndex = 0; entryIndex < this->dIndex_ecc_velocity_wrt_StateAfterEvent[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_ecc_velocity_wrt_StateAfterEvent[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_ecc_velocity_wrt_StateAfterEvent_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = this->boundary_orbit_elements_Jacobian(1, stateIndex + 3) _GETVALUE * std::get<2>(Derivatives_of_StateAfterEvent[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry;
						}
					}

					// now for the time decision variables
					std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

					for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
					{
						// POSITION
						for (size_t entryIndex = 0; entryIndex < this->dIndex_ecc_position_wrt_StateAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_ecc_position_wrt_StateAfterEvent_wrt_Time[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_ecc_position_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = this->boundary_orbit_elements_Jacobian(1, stateIndex) _GETVALUE * std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry;
						}

						// VELOCITY
						for (size_t entryIndex = 0; entryIndex < this->dIndex_ecc_velocity_wrt_StateAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
						{
							size_t dIndex = this->dIndex_ecc_velocity_wrt_StateAfterEvent_wrt_Time[stateIndex][entryIndex];
							size_t Gindex = this->Gindex_ecc_velocity_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex];
							size_t Xindex = this->jGvar->operator[](Gindex);

							double Gentry = this->boundary_orbit_elements_Jacobian(1, stateIndex + 3) _GETVALUE * std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

							G[Gindex] += this->X_scale_factors->operator[](Xindex)
								* Gentry;
						}
					}
				}
			} // end process_constraint()

			void EccentricityConstraint::output(std::ofstream& outputfile)
			{
                std::string framestring = this->myUniverse->central_body_name + " " + this->ReferenceFrameStrings[this->myReferenceFrame];

				outputfile << this->myBoundaryEvent->getName() << " eccentricity (" << framestring << "): " << this->eccentricity _GETVALUE << std::endl;
			} // end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG