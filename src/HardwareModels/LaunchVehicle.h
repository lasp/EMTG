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

#include "HardwareBase.h"
#include "LaunchVehicleOptions.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class LaunchVehicle : public HardwareBase
        {
        public:
            //constructor
            LaunchVehicle();
            LaunchVehicle(const LaunchVehicleOptions& launchvehicleoptions);

            //destructor

            //methods
            void initialize(const LaunchVehicleOptions& launchvehicleoptions);
            void computePerformance(const doubleType& C3, const double& LV_margin);
            
            //get
            doubleType getDeliveredMass() const { return this->DeliveredMass; }
            double getC3max() const { return this->myOptions.getC3_upperbound(); }
            double getC3min() const { return this->myOptions.getC3_lowerbound(); }
            double getdmdC3() const { return this->dmdC3; }

        private:
            //fields
            LaunchVehicleOptions myOptions;
            doubleType DeliveredMass;
            double dmdC3;
        };
    }//end namespace HardwareModels
}//end namespace EMTG