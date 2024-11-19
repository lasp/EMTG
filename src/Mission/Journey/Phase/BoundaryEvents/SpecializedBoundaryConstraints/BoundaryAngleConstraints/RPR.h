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

#include "SpecializedBoundaryConstraintBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class BoundaryEventBase;

        namespace SpecializedConstraints
        {

            class RPRconstraint : virtual public SpecializedBoundaryConstraintBase
            {
            public:
                //constructors
                RPRconstraint() {};

                RPRconstraint(const std::string& name,
                    const size_t& journeyIndex,
                    const size_t& phaseIndex,
                    const size_t& stageIndex,
                    Astrodynamics::universe* Universe,
                    HardwareModels::Spacecraft* mySpacecraft,
                    missionoptions* myOptions,
                    BoundaryEventBase* myBoundaryEvent,
                    const std::string& constraintDefinition);

                //destructor
                virtual ~RPRconstraint() {};


                //public methods
                virtual void output(std::ofstream& outputfile);

                virtual void calcbounds();

                //process
                virtual void process_constraint(const std::vector<doubleType>& X,
                    size_t& Xindex,
                    std::vector<doubleType>& F,
                    size_t& Findex,
                    std::vector<double>& G,
                    const bool& needG);

            protected:
                //fields
                doubleType cosRPR;//centered on "ref2"
                doubleType Angle;
                math::Matrix<doubleType> R_cb_ref1;
                math::Matrix<doubleType> R_cb_ref2;
                math::Matrix<doubleType> R_probe_ref1;
                math::Matrix<doubleType> R_probe_ref2;
                Astrodynamics::body* myReference1;
                Astrodynamics::body* myReference2;
				//Astrodynamics::body* myBoundaryEventBody;

                bool ref1IsCentralBodyOfJourney;
                bool ref2IsCentralBodyOfJourney;
                std::string ref1Name;
                std::string ref2Name;
                math::Matrix<double> dcosRPR_dR_probe_ref1;
                math::Matrix<double> dcosRPR_dR_probe_ref2;
                math::Matrix<double> dR_cb_ref1_dt;
                math::Matrix<double> dR_cb_ref2_dt;
                math::Matrix<double> dR_probe_ref1_dt;
                math::Matrix<double> dR_probe_ref2_dt;

				math::Matrix<doubleType> R_cbJourney_cbEvent;
				math::Matrix<double> dR_cbJourney_cbEvent_dt;

                //this constraint has derivatives with respect to any non-time variable that affects spacecraft position and velocity
                //spacecraft position is used for "Body" vector
                //spacecraft velocity is used for "Probe" vector
                //this constraint has derivatives with respect to all time variables because they affect both spacecraft and reference
                std::vector<size_t> Gindex_constraint_wrt_time_variables;
                std::vector< std::vector<size_t> > Gindex_constraint_wrt_PositionAfterEvent_variables;//stateIndex, Gindex
                std::vector< std::vector<size_t> > Gindex_constraint_wrt_PositionAfterEvent_time_variables;//stateIndex, Gindex

                std::vector< std::vector<size_t> > dIndex_with_respect_to_PositionAfterEvent;//stateIndex, dIndex
                std::vector< std::vector<size_t> > dIndex_with_respect_to_PositionAfterEvent_wrt_Time;//stateIndex, dIndex
            };
        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG