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

#include "RungeKuttaDP87Tableau.h"

namespace EMTG {
	namespace Integration
	{
		// standard constructor
		RungeKuttaDP87Tableau::RungeKuttaDP87Tableau()
		{
			this->setHasVariableStepCoefficients(true); // variable-step method
			this->num_stages = 13;
			this->resizeArrays(); // set sizes of coefficient arrays

            this->A.assign_zeros();
            this->A(1, 0) = 1.0 / 18.0;

            this->A(2, 0) = 1.0 / 48.0;
            this->A(2, 1) = 1.0 / 16.0;

            this->A(3, 0) = 1.0 / 32.0;
            this->A(3, 1) = 0.0;
            this->A(3, 2) = 3.0 / 32.0;

            this->A(4, 0) = 5.0 / 16.0;
            this->A(4, 1) = 0.0;
            this->A(4, 2) = -75.0 / 64.0;
            this->A(4, 3) = 75.0 / 64.0;

            this->A(5, 0) = 3.0 / 80.0;
            this->A(5, 1) = 0.0;
            this->A(5, 2) = 0.0;
            this->A(5, 3) = 3.0 / 16.0;
            this->A(5, 4) = 3.0 / 20.0;

            this->A(6, 0) = 29443841.0 / 614563906.0;
            this->A(6, 1) = 0.0;
            this->A(6, 2) = 0.0;
            this->A(6, 3) = 77736538.0 / 692538347.0;
            this->A(6, 4) = -28693883.0 / 1125000000.0;
            this->A(6, 5) = 23124283.0 / 1800000000.0;

            this->A(7, 0) = 16016141.0 / 946692911.0;
            this->A(7, 1) = 0.0;
            this->A(7, 2) = 0.0;
            this->A(7, 3) = 61564180.0 / 158732637.0;
            this->A(7, 4) = 22789713.0 / 633445777.0;
            this->A(7, 5) = 545815736.0 / 2771057229.0;
            this->A(7, 6) = -180193667.0 / 1043307555.0;

            this->A(8, 0) = 39632708.0 / 573591083.0;
            this->A(8, 1) = 0.0;
            this->A(8, 2) = 0.0;
            this->A(8, 3) = -433636366.0 / 683701615.0;
            this->A(8, 4) = -421739975.0 / 2616292301.0;
            this->A(8, 5) = 100302831.0 / 723423059.0;
            this->A(8, 6) = 790204164.0 / 839813087.0;
            this->A(8, 7) = 800635310.0 / 3783071287.0;

            this->A(9, 0) = 246121993.0 / 1340847787.0;
            this->A(9, 1) = 0.0;
            this->A(9, 2) = 0.0;
            this->A(9, 3) = -37695042795.0 / 15268766246.0;
            this->A(9, 4) = -309121744.0 / 1061227803.0;
            this->A(9, 5) = -12992083.0 / 490766935.0;
            this->A(9, 6) = 6005943493.0 / 2108947869.0;
            this->A(9, 7) = 393006217.0 / 1396673457.0;
            this->A(9, 8) = 123872331.0 / 1001029789.0;

            this->A(10, 0) = -1028468189.0 / 846180014.0;
            this->A(10, 1) = 0.0;
            this->A(10, 2) = 0.0;
            this->A(10, 3) = 8478235783.0 / 508512852.0;
            this->A(10, 4) = 1311729495.0 / 1432422823.0;
            this->A(10, 5) = -10304129995.0 / 1701304382.0;
            this->A(10, 6) = -48777925059.0 / 3047939560.0;
            this->A(10, 7) = 15336726248.0 / 1032824649.0;
            this->A(10, 8) = -45442868181.0 / 3398467696.0;
            this->A(10, 9) = 3065993473.0 / 597172653.0;

            this->A(11, 0) = 185892177.0 / 718116043.0;
            this->A(11, 1) = 0.0;
            this->A(11, 2) = 0.0;
            this->A(11, 3) = -3185094517.0 / 667107341.0;
            this->A(11, 4) = -477755414.0 / 1098053517.0;
            this->A(11, 5) = -703635378.0 / 230739211.0;
            this->A(11, 6) = 5731566787.0 / 1027545527.0;
            this->A(11, 7) = 5232866602.0 / 850066563.0;
            this->A(11, 8) = -4093664535.0 / 808688257.0;
            this->A(11, 9) = 3962137247.0 / 1805957418.0;
            this->A(11, 10) = 65686358.0 / 487910083.0;

            this->A(12, 0) = 403863854.0 / 491063109.0;
            this->A(12, 1) = 0.0;
            this->A(12, 2) = 0.0;
            this->A(12, 3) = -5068492393.0 / 434740067.0;
            this->A(12, 4) = -411421997.0 / 543043805.0;
            this->A(12, 5) = 652783627.0 / 914296604.0;
            this->A(12, 6) = 11173962825.0 / 925320556.0;
            this->A(12, 7) = -13158990841.0 / 6184727034.0;
            this->A(12, 8) = 3936647629.0 / 1978049680.0;
            this->A(12, 9) = -160528059.0 / 685178525.0;
            this->A(12, 10) = 248638103.0 / 1413531060.0;
            this->A(12, 11) = 0.0;

            // populate bUpper
            this->bUpper(0) = 14005451.0 / 335480064.0;
            this->bUpper(1) = 0.0;
            this->bUpper(2) = 0.0;
            this->bUpper(3) = 0.0;
            this->bUpper(4) = 0.0;
            this->bUpper(5) = -59238493.0 / 1068277825.0;
            this->bUpper(6) = 181606767.0 / 758867731.0;
            this->bUpper(7) = 561292985.0 / 797845732.0;
            this->bUpper(8) = -1041891430.0 / 1371343529.0;
            this->bUpper(9) = 760417239.0 / 1151165299.0;
            this->bUpper(10) = 118820643.0 / 751138087.0;
            this->bUpper(11) = -528747749.0 / 2220607170.0;
            this->bUpper(12) = 1.0 / 4.0;

            // populate bLower
            this->bLower(0) = 13451932.0 / 455176623.0;
            this->bLower(1) = 0.0;
            this->bLower(2) = 0.0;
            this->bLower(3) = 0.0;
            this->bLower(4) = 0.0;
            this->bLower(5) = -808719846.0 / 976000145.0;
            this->bLower(6) = 1757004468.0 / 5645159321.0;
            this->bLower(7) = 656045339.0 / 265891186.0;
            this->bLower(8) = -3867574721.0 / 1518517206.0;
            this->bLower(9) = 465885868.0 / 322736535.0;
            this->bLower(10) = 53011238.0 / 667516719.0;
            this->bLower(11) = 2.0 / 45.0;
            this->bLower(12) = 0.0;

            // RK node constants
            // these encode the positions within the RK step of each stage
            // populate c
            this->c(0) = 0.0;
            this->c(1) = 1.0 / 18.0;
            this->c(2) = 1.0 / 12.0;
            this->c(3) = 1.0 / 8.0;
            this->c(4) = 5.0 / 16.0;
            this->c(5) = 3.0 / 8.0;
            this->c(6) = 59.0 / 400.0;
            this->c(7) = 93.0 / 200.0;
            this->c(8) = 5490023248.0 / 9719169821.0;
            this->c(9) = 13.0 / 20.0;
            this->c(10) = 1201146811.0 / 1299019798.0;
            this->c(11) = 1.0;
            this->c(12) = 1.0;

		}
	} // end integration namespace
} // end EMTG namespace
