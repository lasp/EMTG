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

#pragma once

#include "BoundaryEventBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisPeggedBoundary : public virtual BoundaryEventBase
        {
        public:
            //constructor
            EphemerisPeggedBoundary() : BoundaryEventBase::BoundaryEventBase() {};
            EphemerisPeggedBoundary(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            virtual void initialize(const std::string& name,
                                    const size_t& journeyIndex,
                                    const size_t& phaseIndex,
                                    size_t& stageIndex,
                                    Astrodynamics::universe* Universe,
                                    HardwareModels::Spacecraft* mySpacecraft,
                                    missionoptions* myOptions);


            //output
            virtual void output(std::ofstream& outputfile,
                const double& launchdate,
                size_t& eventcount) = 0;
            //calcbounds

        protected:
            //these go in the specialized event
            virtual void calcbounds_event_left_side(const std::vector<double>& MassBounds, const std::vector<double>& EpochBounds, std::vector<size_t> timeVariables = {});
            virtual void calcbounds_event_right_side() = 0;

            virtual void process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            virtual void process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //fields:
            size_t Xindex_mass;
        };
    }//close namespace events
}//close namespace EMTG