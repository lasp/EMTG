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

//event testbed

#pragma once

#include "missionoptions.h"
#include "universe.h"

#include <random>
#include <sstream>

#include "doubleType.h"
#include "EMTG_matrix.h"
#include "StateRepresentation.h"


void StateRepresentation_testbed(EMTG::missionoptions& options,
    std::vector< EMTG::Astrodynamics::universe > TheUniverse,
    std::mt19937 RNG,
    std::uniform_real_distribution<> UniformDouble);

void testRepresentation(EMTG::Astrodynamics::StateRepresentationBase* myStateRepresentation,
    EMTG::math::Matrix<doubleType>& inputState,
    std::ofstream& outputfile);