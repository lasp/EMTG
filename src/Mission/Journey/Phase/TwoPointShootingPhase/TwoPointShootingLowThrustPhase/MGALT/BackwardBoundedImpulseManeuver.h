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

//backward bounded impulse maneuver model
//Jacob Englander 6/26/2017

#pragma once

#include "BoundedImpulseManeuver.h"




namespace EMTG
{
    namespace Phases
    {
        class BackwardBoundedImpulseManeuver : public BoundedImpulseManeuver
        {
        public:
            //constructor
            BackwardBoundedImpulseManeuver() : BoundedImpulseManeuver() {};
            BackwardBoundedImpulseManeuver(const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //methods
            void process_maneuver(const bool& needMTM);

            //get
            inline math::Matrix<double> get_dMassBeforeManeuver_dThrottleComponents() const { return this->dMassBeforeManeuver_dThrottleComponents; }

        protected:
            //fields
            math::Matrix<double> dMassBeforeManeuver_dThrottleComponents;
        };
    } //close namespace Astrodynamics
} //close namespace EMTG