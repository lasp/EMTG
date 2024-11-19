
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

#include "FreePointArrivalWithVinfinity.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class FreePointChemRendezvous : public FreePointArrivalWithVinfinity
        {
        public:
            //specialized constructor
            FreePointChemRendezvous(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //get
            math::Matrix<doubleType> get_Vinfinity_in() const { return this->Vinfinity_in; }

            //output
            void output(std::ofstream& outputfile,
                const double& launchdate,
                size_t& eventcount);

            //process
            void process_event(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //calcbounds
            void calcbounds(std::vector<size_t> timeVariables);

            //output
            void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

        private:
            void calcbounds_event_right_side();
            
            void calcbounds_virtual_propellant_constraints();

            void calcbounds_deltav_contribution();//capture maneuver

            void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);//capture maneuver

            //fields
            doubleType ArrivalDeltavMagnitude;

            double dChemicalFuel_ddeltavMagnitude;
            double dChemicalOxidizer_ddeltavMagnitude;

            std::vector<size_t> dIndex_mass_after_event_wrt_v_infinity;
            std::vector<size_t> Gindices_dDeltav_dVinfinity;
            std::vector<size_t> Gindices_dVirtualChemicalFuel_dVinfinity;
            std::vector<size_t> Gindices_dVirtualChemicalOxidizer_dVinfinity;
            size_t Gindices_dVirtualChemicalFuel_dLeftMass;
            size_t Gindices_dVirtualChemicalOxidizer_dLeftMass;
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG