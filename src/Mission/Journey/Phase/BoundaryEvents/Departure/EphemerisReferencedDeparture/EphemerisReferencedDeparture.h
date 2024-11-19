
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
#include "arrival.h"
#include "departure.h"
#include "EphemerisReferencedBoundary.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisReferencedDeparture : virtual public DepartureEvent,
            virtual public EphemerisReferencedBoundary
        {
        public:
            //default constructor
            EphemerisReferencedDeparture() : DepartureEvent::DepartureEvent() {};

            //specialized constructor
            EphemerisReferencedDeparture(const std::string& name,
                           const size_t& journeyIndex,
                           const size_t& phaseIndex,
                           size_t& stageIndex,
                           Astrodynamics::universe* Universe,
                           HardwareModels::Spacecraft* mySpacecraft,
                           missionoptions* myOptions,
                           ArrivalEvent* PreviousPhaseArrivalEvent);

            //initialize
            virtual void initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                ArrivalEvent* PreviousPhaseArrivalEvent);

            //destructor
            virtual ~EphemerisReferencedDeparture() {};

            //output
            virtual void output(std::ofstream& outputfile,
                const double& launchdate,
                size_t& eventcount) = 0;


            math::Matrix<doubleType> get_periapse_state() { return math::Matrix<doubleType>(6, 1, 0.0); };//TODO: eventually I might use this to derive a periapse state?

        protected:

            //method to calculate event left side

            virtual void calcbounds_event_left_side(std::vector<size_t> timeVariables);

            //virtual void calcbounds_event_right_side();

            virtual void calcbounds_specialized_constraints();

            virtual void process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_event_right_side(const std::vector<doubleType>& X, //this is just to handle the derivative with respect to wait time
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields - we need these if we are "departing" from a free point that the previous journey arrived at
            //all these do is track the derivative entries in the previous phase arrival event's "state after event" so that we can copy them
            //we only copy the position and velocity variables - time and mass are encoded directly in the departure event

            std::vector< std::vector<size_t> > dIndex_StateAfterPreviousEvent_wrt_DecisionVariables;//state, dIndex
            std::vector< std::vector<size_t> > dIndex_StateAfterPreviousEvent_wrt_DecisionVariables_wrt_Time;//state, dIndex
            std::vector< std::vector<size_t> > dIndex_StateBeforeEvent_wrt_DecisionVariables;//state, dIndex
            std::vector< std::vector<size_t> > dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time;//state, dIndex

        };

        //explicit instantiation
    }//end namespace BoundaryEvents
}//end namespace EMTG