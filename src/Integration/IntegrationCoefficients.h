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

#ifndef INTEGRATIONCOEFFICIENTS_H
#define INTEGRATIONCOEFFICIENTS_H

#include "EMTG_Matrix.h"

namespace EMTG {
	namespace Integration
	{
		class IntegrationCoefficients
		{

		public:
			// constructors
			IntegrationCoefficients();

			// destructor
			~IntegrationCoefficients();

            // getters for private variables
            inline bool getHasVariableStepCoefficients() const { return this->has_variable_step_coeffs; }

            // setters for private variables
            inline void setHasVariableStepCoefficients(const bool & has_variable_step_coeffs) { this->has_variable_step_coeffs = has_variable_step_coeffs; }

		private:
			bool has_variable_step_coeffs;

		}; // end class definition

	} // end integration namespace
} // end EMTG namespace


#endif