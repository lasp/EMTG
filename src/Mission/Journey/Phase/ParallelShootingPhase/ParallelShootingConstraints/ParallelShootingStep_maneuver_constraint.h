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

//ParallelShootingStep maneuver constraint base class
//Jacob Englander 4-13-2018

#pragma once

#include <string>

#include "missionoptions.h"
#include "universe.h"

#include "sparsey_thing.h"

#include "ParallelShootingStep.h"

namespace EMTG
{
    namespace Phases
    {
        //forward declare parent class
        class ParallelShootingStep;

        class ParallelShootingStep_maneuver_constraint : public sparsey_thing
        {
        public:
            //constructor
            ParallelShootingStep_maneuver_constraint() {};
            ParallelShootingStep_maneuver_constraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& subStepIndex,
                const size_t& stageIndex,
                ParallelShootingStep* myStep,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                const std::string& ConstraintDefinition);

			virtual ~ParallelShootingStep_maneuver_constraint();

            //clone
            virtual ParallelShootingStep_maneuver_constraint* clone() const = 0;

            //calcbounds goes in the specialized phase
            virtual void calcbounds() = 0;

            void setup_calcbounds(std::vector<double>* Xupperbounds,
                std::vector<double>* Xlowerbounds,
                std::vector<double>* X_scale_factors,
                std::vector<double>* Fupperbounds,
                std::vector<double>* Flowerbounds,
				std::vector<double>* F_scale_factors,
                std::vector<std::string>* Xdescriptions,
                std::vector<std::string>* Fdescriptions,
                std::vector<size_t>* iGfun,
                std::vector<size_t>* jGvar,
                std::vector<std::string>* Gdescriptions,
                std::vector<size_t>* iAfun,
                std::vector<size_t>* jAvar,
                std::vector<std::string>* Adescriptions,
                std::vector<double>* A);

            //process goes in the specialized phase
            virtual void process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

        protected:

            //fields
            std::string name;
            size_t journeyIndex;
            size_t phaseIndex;
            size_t stageIndex;
            size_t stepIndex;
            size_t subStepIndex;
            ParallelShootingStep* myStep;
            Astrodynamics::universe* myUniverse;
            HardwareModels::Spacecraft* mySpacecraft;
            missionoptions* myOptions;
            JourneyOptions* myJourneyOptions;
            std::string ConstraintDefinition;

            double lowerBound, upperBound;
        };
        
        inline ParallelShootingStep_maneuver_constraint * new_clone(ParallelShootingStep_maneuver_constraint const & other)
        {
            return other.clone();
        }
    }//close namespace Phases
}//close namespace EMTG