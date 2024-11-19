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

//forward MGAnDSMs maneuver epoch constraint with respect to previous event
//9-26-2017

#include "MGAnDSMs_forward_maneuver_epoch_constraint_wrt_previous_event.h"

#include "EMTG_solver_utilities.h"

#include "MGAnDSMs_subphase.h"

#include "boost/algorithm/string/split.hpp"

#include <exception>

namespace EMTG
{
    namespace Phases
    {
        MGAnDSMs_forward_maneuver_epoch_constraint_wrt_previous_event::MGAnDSMs_forward_maneuver_epoch_constraint_wrt_previous_event(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& subphaseIndex,
            const size_t& stageIndex,
            MGAnDSMs_subphase* mySubPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            const std::string& ConstraintDefinition) :
            MGAnDSMs_maneuver_epoch_constraint::MGAnDSMs_maneuver_epoch_constraint(name,
                journeyIndex,
                phaseIndex,
                subphaseIndex,
                stageIndex,
                mySubPhase,
                Universe,
                mySpacecraft,
                myOptions,
                ConstraintDefinition)
        {
        }//end constructor
        //************************************calcbounds methods
        void MGAnDSMs_forward_maneuver_epoch_constraint_wrt_previous_event::calcbounds()
        {
            //Step 0: can we even create this constraint?
            if (this->mySubPhase->getBordersBoundary())
            {
                throw std::invalid_argument("You cannot constrain the epoch of an MGAnDSMs phase's first forward subphase maneuver with respect to the previous maneuver because there is no previous maneuver, you're at the beginning of the phase.\n" + this->ConstraintDefinition + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            //Step 1: parse the constraint definition
            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            //Step 2: create the constraint
            this->lowerBound = std::stod(ConstraintDefinitionCell[3]) / 100.0;
            this->upperBound = std::stod(ConstraintDefinitionCell[4]) / 100.0;
            this->Flowerbounds->push_back(this->lowerBound - this->upperBound);
            this->Fupperbounds->push_back(0.0);
            this->Fdescriptions->push_back(prefix + "maneuver epoch relative to previous event constraint");

            //sparsity pattern - derivatives with respect to phase flight time
            this->Xindex_PhaseFlightTime = this->mySubPhase->get_timeVariables().back();

            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xindex_PhaseFlightTime,
                this->Gindex_wrt_TimeVariables);

            //derivatives with respect to burn index - only my current burn index because that's the measurement between now and the previous event
            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->mySubPhase->getXindex_burnIndex().back(),
                this->Gindex_wrt_BurnIndices);
        }//end calcbounds()

         //******************************************process methods
        void MGAnDSMs_forward_maneuver_epoch_constraint_wrt_previous_event::process_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: evaluate the constraint
            doubleType SubPhaseTime = this->mySubPhase->getSubPhaseTime();
            F[Findex++] = SubPhaseTime / 86400.0 / 100.0 - this->upperBound;

            //Step 2: derivatives
            if (needG)
            {
                //with respect to time variables
                //first (n-1) time variables
                for (size_t timeIndex = 0; timeIndex < this->Gindex_wrt_TimeVariables.size() - 1; ++timeIndex)
                {
                    size_t Gindex = this->Gindex_wrt_TimeVariables[timeIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        / 86400.0 / 100.0;
                }


                double SumOfBurnIndices = 0.0;
                size_t Gindex_PhaseFlightTime = this->Gindex_wrt_TimeVariables.back();
                size_t Xindex_PhaseFlightTime = this->jGvar->operator[](Gindex_PhaseFlightTime);
                doubleType PhaseFlightTime = X[Xindex_PhaseFlightTime];

                //with respect to burn index
                for (size_t BurnIndexIndex = 0; BurnIndexIndex < this->Gindex_wrt_BurnIndices.size(); ++BurnIndexIndex)
                {
                    size_t Gindex = this->Gindex_wrt_BurnIndices[BurnIndexIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    SumOfBurnIndices += X[Xindex]_GETVALUE;

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * PhaseFlightTime _GETVALUE
                        / 86400.0 / 100.0;
                }

                //phase flight time - sum of all the burn indices
                G[Gindex_PhaseFlightTime] = this->X_scale_factors->operator[](Xindex_PhaseFlightTime)
                    * SumOfBurnIndices
                    / 86400.0 / 100.0;

            }//end derivatives
        }//end process_constraint()
    }//close namespace Phases
}//close namespace EMTG