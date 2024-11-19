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

//One Step to Rule Them All (both first and last step at the same time)
//    One Step to rule them all, One step to find them; One step to bring them all
//    and in the darkness bind them.
//Jacob Englander 4-11-2018

#include "ParallelShootingOneStepToRuleThemAll.h"

namespace EMTG
{
    namespace Phases
    {
        ParallelShootingOneStepToRuleThemAll::ParallelShootingOneStepToRuleThemAll() {};

        ParallelShootingOneStepToRuleThemAll::ParallelShootingOneStepToRuleThemAll(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& stageIndex,
            ParallelShootingPhase* myPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->initialize(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                stageIndex,
                myPhase,
                myUniverse,
                mySpacecraft,
                myOptions);
        }

        void ParallelShootingOneStepToRuleThemAll::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& stageIndex,
            ParallelShootingPhase* myPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->ParallelShootingStep::initialize(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                stageIndex,
                NULL, //first step doesn't have a previous step!
                myPhase,
                myUniverse,
                mySpacecraft,
                myOptions);
        }//end initialize
    }//close namespace Phases
}//close namespace EMTG