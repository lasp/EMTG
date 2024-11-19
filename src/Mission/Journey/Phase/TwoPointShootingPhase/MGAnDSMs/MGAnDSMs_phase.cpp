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

//EMTGv9 MGAnDSMs phase
//Jacob Englander 8-25-2017

#include "MGAnDSMs_phase.h"
#include "IntegrationSchemeFactory.h"
#include "Forward_MGAnDSMs_subphase.h"
#include "Backward_MGAnDSMs_subphase.h"
#include "PhaseDistanceConstraintFactory.h"


namespace EMTG
{
    namespace Phases
    {
        MGAnDSMs_phase::MGAnDSMs_phase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions) :
            TwoPointShootingPhase::TwoPointShootingPhase(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions,
                10, //numStatesToPropagate
                9) //numMatchConstraints
        {
            //name our match point constraints
            this->stateVectorNames.push_back("virtual chemical fuel");
            this->stateVectorNames.push_back("virtual chemical oxidizer");
            this->matchPointConstraintNames.push_back("virtual chemical fuel");
            this->matchPointConstraintNames.push_back("virtual chemical oxidizer");
            this->matchPointConstraintStateIndex.push_back(8); //note we skipped the encode of boundary epoch, it does not get a match point constraint
            this->matchPointConstraintStateIndex.push_back(9);
            this->continuity_constraint_scale_factors(7) = 1.0 / (1.0 * this->myJourneyOptions->maximum_mass);
            this->continuity_constraint_scale_factors(8) = 1.0 / (1.0 * this->myJourneyOptions->maximum_mass);
            this->stateIndex_phase_propagation_variable = 10;
            this->stateVectorNames.push_back("phase flight time");

            //hey, we're a chemical phase!
            this->hasBipropManeuver = true;
            this->hasMonopropManeuver = true;

            //how many subphases do we have, and then size the subphase arrays
            this->numberOfDSMs = this->myJourneyOptions->impulses_per_phase;
            div_t divresult = div(this->numberOfDSMs, 2);
            if (divresult.rem == 0 && this->numberOfDSMs > 0)
                this->subphaseMatchPointIndex = divresult.quot - 1;
            else
                this->subphaseMatchPointIndex = divresult.quot;

            this->numberOfForwardSubphases = this->subphaseMatchPointIndex + 1;
            this->numberOfBackwardSubphases = this->numberOfDSMs == 0 ? 1 : this->numberOfDSMs - this->subphaseMatchPointIndex;

            //derivative truth table
           
            //unless SRP is enabled, only mass affects mass
            if (!this->myOptions->perturb_SRP)
            {
                //only mass variables affect mass constraints
                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[6][stateIndex] = false;
                    this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[6][stateIndex] = false;
                }
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[6][7] = false; //epoch
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[6][7] = false; //epoch

                //mass variables do not affect other constraints
                for (size_t constraintIndex = 0; constraintIndex < 6; ++constraintIndex)
                {
                    this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][6] = false; //mass variables
                    this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][6] = false; //mass variables
                }
            }
            //always true - propellant tanks are only affected by themselves and mass
            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][stateIndex] = false; //fuel mass wrt state
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[8][stateIndex] = false; //oxidizer mass wrt state
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][stateIndex] = false; //fuel mass wrt state
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[8][stateIndex] = false; //oxidizer mass wrt state
            }
            this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][7] = false; //fuel mass wrt epoch
            this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[8][7] = false; //oxidizer mass wrt epoch
            this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][9] = false; //fuel mass wrt oxidizer mass
            this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[8][8] = false; //oxidizer mass wrt fuel mass
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][7] = false; //fuel mass wrt epoch
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[8][7] = false; //oxidizer mass wrt epoch
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][9] = false; //fuel mass wrt oxidizer mass
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[8][8] = false; //oxidizer mass wrt fuel mass
            //always true - propellant tank states do not affect anything but themselves
            for (size_t constraintIndex = 0; constraintIndex < 7; ++constraintIndex)
            {
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][8] = false; //other states wrt fuel mass
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][9] = false; //other states wrt oxidizer mass
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][8] = false; //other states wrt fuel mass
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][9] = false; //other states wrt oxidizer mass
            }

            //size the vectors of state vectors

            this->SpacecraftStateForward.resize(this->numberOfForwardSubphases - 1, math::Matrix<doubleType>(this->numStatesToPropagate, 1, 0.0));
            this->SpacecraftStateBackward.resize(this->numberOfBackwardSubphases - 1, math::Matrix<doubleType>(this->numStatesToPropagate, 1, 0.0));

            //size the SPTM vectors
            this->ForwardSPTM.resize(this->numberOfForwardSubphases, math::Matrix<double>(this->numStatesToPropagate + 2, math::identity));
            this->BackwardSPTM.resize(this->numberOfBackwardSubphases, math::Matrix<double>(this->numStatesToPropagate + 2, math::identity));
            this->CumulativeForwardSPTM.resize(this->numberOfForwardSubphases, math::Matrix<double>(this->numStatesToPropagate + 2, math::identity));
            this->CumulativeBackwardSPTM.resize(this->numberOfBackwardSubphases, math::Matrix<double>(this->numStatesToPropagate + 2, math::identity));
            this->Forward_dPropagatedState_dIndependentVariable = this->myPropagatorType == PropagatorType::IntegratedPropagator ? math::Matrix<double>( 10, 2, 0.0) : math::Matrix<double>(6, 1, 0.0);
            this->Backward_dPropagatedState_dIndependentVariable = this->myPropagatorType == PropagatorType::IntegratedPropagator ? math::Matrix<double>(10, 2, 0.0) : math::Matrix<double>(6, 1, 0.0);

            //handle the initial and terminal coast
            if (this->hasInitialCoast)
            {
                this->myJourneyOptions->ManeuverConstraintDefinitions.push_back("p" + std::to_string(this->phaseIndex)
                    + "b0"
                    + "_epoch_lboundary_" + std::to_string(this->InitialCoastDuration / 86400.0) + "_100000.0");
            }
            if (this->hasTerminalCoast)
            {
                this->myJourneyOptions->ManeuverConstraintDefinitions.push_back("p" + std::to_string(this->phaseIndex)
                    + "b" + std::to_string(this->numberOfDSMs - 1)
                    + "_epoch_rboundary_" + std::to_string(this->TerminalCoastDuration / 86400.0) + "_100000.0");
            }

            //configure the forward subphases
            for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
            {
                std::string subphaseName = "forwardSubPhase" + std::to_string(forwardSubPhaseIndex);
                this->ForwardSubPhases.push_back(new Forward_MGAnDSMs_subphase(name,
                    this->journeyIndex,
                    this->phaseIndex,
                    forwardSubPhaseIndex,
                    this->stageIndex,
                    this,
                    forwardSubPhaseIndex == 0 ? NULL : &this->ForwardSubPhases[forwardSubPhaseIndex - 1],
                    this->myUniverse,
                    this->mySpacecraft,
                    this->myOptions));
                
                this->ForwardSubPhases.back().set_PhaseFlightTime(this->PhaseFlightTime);
                this->ForwardSubPhases.back().set_SPTM(this->ForwardSPTM[forwardSubPhaseIndex]);
                this->ForwardSubPhases.back().set_dPropagatedStatedIndependentVariable(this->Forward_dPropagatedState_dIndependentVariable);

                if (forwardSubPhaseIndex == 0)
                    this->ForwardSubPhases.back().set_spacecraft_state_minus(this->state_after_initial_TCM);
                else
                    this->ForwardSubPhases.back().set_spacecraft_state_minus(this->SpacecraftStateForward[forwardSubPhaseIndex - 1]);

                if (forwardSubPhaseIndex == this->numberOfForwardSubphases - 1)
                    this->ForwardSubPhases.back().set_spacecraft_state_plus(this->match_point_state_minus);
                else
                    this->ForwardSubPhases.back().set_spacecraft_state_plus(this->SpacecraftStateForward[forwardSubPhaseIndex]);
            }//end configure forward subphases

             //configure the backward subphases
            for (size_t backwardSubPhaseIndex = 0; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
            {
                std::string subphaseName = "backwardSubPhase" + std::to_string(backwardSubPhaseIndex);
                this->BackwardSubPhases.push_back(new Backward_MGAnDSMs_subphase(name,
                    this->journeyIndex,
                    this->phaseIndex,
                    backwardSubPhaseIndex,
                    this->stageIndex,
                    this,
                    backwardSubPhaseIndex == 0 ? NULL : &this->BackwardSubPhases[backwardSubPhaseIndex - 1],
                    this->myUniverse,
                    this->mySpacecraft,
                    this->myOptions));

                this->BackwardSubPhases.back().set_PhaseFlightTime(this->PhaseFlightTime);
                this->BackwardSubPhases.back().set_SPTM(this->BackwardSPTM[backwardSubPhaseIndex]);
                this->BackwardSubPhases.back().set_dPropagatedStatedIndependentVariable(this->Backward_dPropagatedState_dIndependentVariable);

                if (backwardSubPhaseIndex == 0)
                    this->BackwardSubPhases.back().set_spacecraft_state_plus(this->state_at_end_of_phase);
                else                                                    
                    this->BackwardSubPhases.back().set_spacecraft_state_plus(this->SpacecraftStateBackward[backwardSubPhaseIndex - 1]);

                if (backwardSubPhaseIndex == this->numberOfBackwardSubphases - 1)
                    this->BackwardSubPhases.back().set_spacecraft_state_minus(this->match_point_state_plus);
                else
                    this->BackwardSubPhases.back().set_spacecraft_state_minus(this->SpacecraftStateBackward[backwardSubPhaseIndex]);
            }//end configure backward subphases


            //size the TCM transition matrix
            this->TCMTM = math::Matrix<double>(this->numStatesToPropagate + 2, math::identity);

            //distance contraint
            std::string shortprefix = "p" + std::to_string(this->phaseIndex);

            for (std::string& constraint : this->myJourneyOptions->PhaseDistanceConstraintDefinitions)
            {
                if (constraint.find("#") != 0
                    && (constraint.find(shortprefix) < 1024 || (constraint.find("pEnd") < 1024 && this->isLastPhaseInJourney)))
                {
                    this->distance_constraint_definitions.push_back(CreatePhaseDistanceConstraint(constraint, this->myOptions, this->myUniverse));
                }
            }
            this->number_of_distance_constraints = this->distance_constraint_definitions.size();

            if (this->number_of_distance_constraints > 0)
            {

                this->distance_constraint_relative_position.resize(this->numberOfDSMs, std::vector< math::Matrix<doubleType> >(this->number_of_distance_constraints, math::Matrix<doubleType>(3, 1, 0.0))); //step, constraint, state variable
                this->distance_constraint_body_position_time_derivatives.resize(this->numberOfDSMs, std::vector< math::Matrix<double> >(this->number_of_distance_constraints, math::Matrix<double>(3, 1, 0.0))); //step, constraint, state 
                this->distance_from_body.resize(this->numberOfDSMs, std::vector<doubleType>(this->number_of_distance_constraints, 0.0)); //step, constraint

                //derivative indices
                this->G_indices_distance_constraints_wrt_ForwardControl.resize(this->numberOfForwardSubphases, std::vector< std::vector< std::vector<size_t> > >(this->number_of_distance_constraints)); //step, constraint, controlstep, control variable
                this->G_indices_distance_constraints_wrt_BackwardControl.resize(this->numberOfBackwardSubphases - 1, std::vector< std::vector< std::vector<size_t> > >(this->number_of_distance_constraints)); //step, constraint, controlstep, control variable
                this->G_indices_distance_constraints_wrt_LeftBoundaryState.resize(this->numberOfForwardSubphases, std::vector< std::vector<size_t> >(this->number_of_distance_constraints)); //step, constraint, state variable
                this->G_indices_distance_constraints_wrt_RightBoundaryState.resize(this->numberOfBackwardSubphases - 1, std::vector< std::vector<size_t> >(this->number_of_distance_constraints)); //step, constraint, state variable
                this->G_indices_distance_constraints_wrt_LeftBoundaryTime.resize(this->numberOfForwardSubphases, std::vector< std::vector<size_t> >(this->number_of_distance_constraints)); //step, constraint, state variable
                this->G_indices_distance_constraints_wrt_RightBoundaryTime.resize(this->numberOfBackwardSubphases - 1, std::vector< std::vector<size_t> >(this->number_of_distance_constraints)); //step, constraint, state variable
                this->G_indices_distance_constraints_wrt_PhaseFlightTime.resize(this->numberOfDSMs, std::vector<size_t>(this->number_of_distance_constraints));
                this->G_indices_distance_constraints_wrt_BurnIndex.resize(this->numberOfDSMs, std::vector< std::vector<size_t> >(this->number_of_distance_constraints));
            }
        }//end constructor

        MGAnDSMs_phase::~MGAnDSMs_phase()
        {
        }

        //************************************calcbounds methods
        void MGAnDSMs_phase::setup_calcbounds(std::vector<double>* Xupperbounds,
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
            std::vector<double>* A)
        {
            //base class
            this->phase::setup_calcbounds(Xupperbounds,
                Xlowerbounds,
                X_scale_factors,
                Fupperbounds,
                Flowerbounds,
				F_scale_factors,
                Xdescriptions,
                Fdescriptions,
                iGfun,
                jGvar,
                Gdescriptions,
                iAfun,
                jAvar,
                Adescriptions,
                A);

            //subphases
            for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
                this->ForwardSubPhases[forwardSubPhaseIndex].setup_calcbounds(Xupperbounds,
                    Xlowerbounds,
                    X_scale_factors,
                    Fupperbounds,
                    Flowerbounds,
					F_scale_factors,
                    Xdescriptions,
                    Fdescriptions,
                    iGfun,
                    jGvar,
                    Gdescriptions,
                    iAfun,
                    jAvar,
                    Adescriptions,
                    A);

            for (size_t backwardSubPhaseIndex = 0; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
                this->BackwardSubPhases[backwardSubPhaseIndex].setup_calcbounds(Xupperbounds,
                    Xlowerbounds,
                    X_scale_factors,
                    Fupperbounds,
                    Flowerbounds,
					F_scale_factors,
                    Xdescriptions,
                    Fdescriptions,
                    iGfun,
                    jGvar,
                    Gdescriptions,
                    iAfun,
                    jAvar,
                    Adescriptions,
                    A);
        }//end setup_calcbounds()

        void MGAnDSMs_phase::calcbounds()
        {
            this->calcbounds_phase_left_boundary();

            this->calcbounds_phase_flight_time();

            this->calcbounds_phase_right_boundary();

            this->calcbounds_phase_main();
        }//end calcbounds()

        void MGAnDSMs_phase::calcbounds_phase_main()
        {
            this->calcbounds_virtual_propellant_tanks();

            this->calcbounds_subphases();

            this->calcbounds_match_point_constraints();

            if (this->number_of_distance_constraints > 0)
            {
                this->calcbounds_distance_constraints();
            }

            this->calcbounds_burnindex_sum_constraint();

            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV)
                this->calcbounds_deltav_contribution();
        }//end calcbounds_phase_main()

        void MGAnDSMs_phase::calcbounds_subphases()
        {
            //subphases - give them the arrival event's time variables because they all depend on previous times plus phase flight time
            std::vector<size_t> timeVariables = this->myArrivalEvent->get_Xindices_EventRightEpoch();

            for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
                this->ForwardSubPhases[forwardSubPhaseIndex].calcbounds(timeVariables);

            for (size_t backwardSubPhaseIndex = 0; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
                this->BackwardSubPhases[backwardSubPhaseIndex].calcbounds(timeVariables);

            //maneuver constraints
            for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
                this->ForwardSubPhases[forwardSubPhaseIndex].calcbounds_maneuver_constraints();

            for (size_t backwardSubPhaseIndex = 0; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
                this->BackwardSubPhases[backwardSubPhaseIndex].calcbounds_maneuver_constraints();

        }//end calcbounds_subphases()

        void MGAnDSMs_phase::calcbounds_match_point_constraints()
        {
            //base class
            this->TwoPointShootingPhase::calcbounds_match_point_constraints();

            //forward control
            for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
            {
                
                std::vector<size_t> Xindex_DSM_components = this->ForwardSubPhases[forwardSubPhaseIndex].getXindex_DSM_components();

                std::vector< std::vector<size_t> > subphase_Gindices_match_point_constraint_ForwardControl;

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    std::vector<size_t> state_Gindices_match_point_constraint_ForwardControl;

                    for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                    {
                        if (!(forwardSubPhaseIndex == this->numberOfForwardSubphases - 1)
                            || constraintIndex == 3 + Vindex
                            || constraintIndex > 5)
                        {
                            size_t Xindex = Xindex_DSM_components[Vindex];
                            this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                                Xindex,
                                state_Gindices_match_point_constraint_ForwardControl);
                        }
                        else
                            state_Gindices_match_point_constraint_ForwardControl.push_back(0);
                    }
                    subphase_Gindices_match_point_constraint_ForwardControl.push_back(state_Gindices_match_point_constraint_ForwardControl);
                }

                this->Gindices_match_point_constraint_ForwardControl.push_back(subphase_Gindices_match_point_constraint_ForwardControl);
            }

            //backward control - skip the first one since there is no control there
            this->Gindices_match_point_constraint_BackwardControl.push_back(std::vector< std::vector<size_t> >());
            for (size_t backwardSubPhaseIndex = 1; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
            {
                std::vector<size_t> Xindex_DSM_components = this->BackwardSubPhases[backwardSubPhaseIndex].getXindex_DSM_components();

                std::vector< std::vector<size_t> > subphase_Gindices_match_point_constraint_BackwardControl;
                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    std::vector<size_t> state_Gindices_match_point_constraint_BackwardControl;

                    for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                    {
                        size_t Xindex = Xindex_DSM_components[Vindex];
                        this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                            Xindex,
                            state_Gindices_match_point_constraint_BackwardControl);
                    }
                    subphase_Gindices_match_point_constraint_BackwardControl.push_back(state_Gindices_match_point_constraint_BackwardControl);
                }

                this->Gindices_match_point_constraint_BackwardControl.push_back(subphase_Gindices_match_point_constraint_BackwardControl);
            }

            //burn indices
            for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
            {
                size_t Xindex = this->ForwardSubPhases[forwardSubPhaseIndex].getXindex_burnIndex().back();
                std::vector<size_t> subphase_Gindices_match_point_constraint_ForwardBurnIndex;

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                            Xindex,
                            subphase_Gindices_match_point_constraint_ForwardBurnIndex);
                    }
                }
                this->Gindices_match_point_constraint_ForwardBurnIndex.push_back(subphase_Gindices_match_point_constraint_ForwardBurnIndex);
            }

            for (size_t backwardSubPhaseIndex = 0; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
            {
                size_t Xindex = this->BackwardSubPhases[backwardSubPhaseIndex].getXindex_burnIndex().back();
                std::vector<size_t> subphase_Gindices_match_point_constraint_BackwardBurnIndex;

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {

                    if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                            Xindex,
                            subphase_Gindices_match_point_constraint_BackwardBurnIndex);
                    }
                }
                this->Gindices_match_point_constraint_BackwardBurnIndex.push_back(subphase_Gindices_match_point_constraint_BackwardBurnIndex);
            }


            //virtual chemical fuel
            this->create_sparsity_entry(this->Findices_match_point_constraints[7],
                this->X_index_of_first_decision_variable_in_this_phase,
                true,
                this->prefix + "virtual chemical fuel",
                this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel);

            //virtual chemical oxidizer
            this->create_sparsity_entry(this->Findices_match_point_constraints[8],
                this->X_index_of_first_decision_variable_in_this_phase,
                true,
                this->prefix + "virtual chemical oxidizer",
                this->Gindex_dVirtualChemicalOxidizerConstraint_dVirtualChemicalOxidizer);
        }//end calcbounds_match_point_constraints()

        void MGAnDSMs_phase::calcbounds_distance_constraints()
        {
            //create constraints
            for (size_t stepIndex = 0; stepIndex < this->numberOfDSMs; ++stepIndex)
            {
                for (size_t constraintIndex = 0; constraintIndex < this->number_of_distance_constraints; ++constraintIndex)
                {
                    //create the constraint
                    this->Flowerbounds->push_back((std::get<1>(this->distance_constraint_definitions[constraintIndex]) - std::get<2>(this->distance_constraint_definitions[constraintIndex])) / this->myUniverse->LU);
                    this->Fupperbounds->push_back(0.0);
                    int bodyIndex = std::get<0>(this->distance_constraint_definitions[constraintIndex]);
                    if (bodyIndex == -2)
                    {
                        this->Fdescriptions->push_back(prefix + " maneuver " + std::to_string(stepIndex) + " distance from " + this->myUniverse->central_body.name);
                    }
                    else
                    {
                        this->Fdescriptions->push_back(prefix + " maneuver " + std::to_string(stepIndex) + " distance from " + this->myUniverse->bodies[bodyIndex - 1].name);
                    }

                    //derivatives with respect to boundary states
                    for (size_t encodedStateIndex = 0; encodedStateIndex < 8; ++encodedStateIndex)
                    {
                        if (stepIndex < this->numberOfForwardSubphases) //left half-phase
                        {
                            for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingLeftBoundary.size(); ++varIndex)
                            {
                                size_t Xindex = this->ListOfVariablesAffectingLeftBoundary[varIndex];

                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->G_indices_distance_constraints_wrt_LeftBoundaryState[stepIndex][constraintIndex]);
                            }

                            for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingLeftBoundary.size(); ++varIndex)
                            {
                                size_t Xindex = this->ListOfTimeVariablesAffectingLeftBoundary[varIndex];

                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->G_indices_distance_constraints_wrt_LeftBoundaryTime[stepIndex][constraintIndex]);
                            }
                        }//end forward case
                        else //right half-phase
                        {
                            //back_burnIndex is the number of burns removed from the right-hand boundary
                            //note that the first backward subphase does not have a maneuver, so the "zeroth" backward maneuver actually corresponds to the second backward subphase
                            size_t back_burnIndex = this->numberOfDSMs - stepIndex - 1;

                            for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingRightBoundary.size(); ++varIndex)
                            {
                                size_t Xindex = this->ListOfVariablesAffectingRightBoundary[varIndex];

                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->G_indices_distance_constraints_wrt_RightBoundaryState[back_burnIndex][constraintIndex]);
                            }

                            for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingRightBoundary.size(); ++varIndex)
                            {
                                size_t Xindex = this->ListOfTimeVariablesAffectingRightBoundary[varIndex];

                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->G_indices_distance_constraints_wrt_RightBoundaryTime[back_burnIndex][constraintIndex]);
                            }
                        }//end backward case
                    }//end loop over state variables

                    //derivatives with respect to control - this is a triangle!
                    if (stepIndex < this->numberOfForwardSubphases)
                    {
                        //derivative with respect to all previous control variables
                        for (size_t controlstep = 0; controlstep < stepIndex; ++controlstep)
                        {
                            std::vector<size_t> controlstep_G_indices_distance_constraints_wrt_Control;

                            for (size_t controlIndex : {0, 1, 2})
                            {
                                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                    this->ForwardSubPhases[controlstep].getXindex_DSM_components()[controlIndex],
                                    controlstep_G_indices_distance_constraints_wrt_Control);
                            }

                            this->G_indices_distance_constraints_wrt_ForwardControl[stepIndex][constraintIndex].push_back(controlstep_G_indices_distance_constraints_wrt_Control);
                        }
                    }//end forward control
                    else
                    {
                        //back_burnIndex is the number of burns removed from the right-hand boundary
                        //note that the first backward subphase does not have a maneuver, so the "zeroth" backward maneuver actually corresponds to the second backward subphase
                        size_t back_burnIndex = this->numberOfDSMs - stepIndex - 1;

                        //derivative with respect to all previous control variables (backward subphases)
                        //note that the zeroth backward subphase does not have a maneuver
                        for (size_t controlstep = 1; controlstep <= back_burnIndex; ++controlstep)
                        {
                            std::vector<size_t> controlstep_G_indices_distance_constraints_wrt_Control;

                            for (size_t controlIndex : {0, 1, 2})
                            {
                                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                    this->BackwardSubPhases[controlstep].getXindex_DSM_components()[controlIndex],
                                    controlstep_G_indices_distance_constraints_wrt_Control);
                            }

                            this->G_indices_distance_constraints_wrt_BackwardControl[back_burnIndex][constraintIndex].push_back(controlstep_G_indices_distance_constraints_wrt_Control);
                        }
                    }//end backward control

                    //derivatives with respect to phase flight time
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        this->Xindex_PhaseFlightTime,
                        this->G_indices_distance_constraints_wrt_PhaseFlightTime[stepIndex][constraintIndex]);

                    //derivatives with respect to burn index
                    //fortunately the subphases already do the accounting for us
                    if (stepIndex < this->numberOfForwardSubphases)
                    {
                        for (size_t Xindex : this->ForwardSubPhases[stepIndex].getXindex_burnIndex())
                        {
                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                Xindex,
                                this->G_indices_distance_constraints_wrt_BurnIndex[stepIndex][constraintIndex]);
                        }
                    }//end forward case
                    else
                    {
                        //back_burnIndex is the number of burns removed from the right-hand boundary
                        //note that the first backward subphase does not have a maneuver, so the "zeroth" backward maneuver actually corresponds to the second backward subphase
                        size_t back_burnIndex = this->numberOfDSMs - stepIndex - 1;

                        for (size_t Xindex : this->BackwardSubPhases[back_burnIndex].getXindex_burnIndex())
                        {
                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                Xindex,
                                this->G_indices_distance_constraints_wrt_BurnIndex[stepIndex][constraintIndex]);
                        }
                    }//end backward case

                }//end loop over distance constraints
            }//end loop over maneuvers
        }//end calcbounds_distance_constraints

        void MGAnDSMs_phase::calcbounds_burnindex_sum_constraint()
        {
            this->Flowerbounds->push_back(-1.0e-13);
            this->Fupperbounds->push_back(1.0e-13);
            this->Fdescriptions->push_back(this->prefix + "burn index sum constraint");

            //forward dependencies
            std::vector<size_t> Xindex_burnIndex = this->ForwardSubPhases.back().getXindex_burnIndex();

            for (size_t entry = 0; entry < Xindex_burnIndex.size(); ++entry)
            {
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    Xindex_burnIndex[entry],
                    this->Gindices_BurnIndexSumConstraint);
            }

            //backward dependencies
            Xindex_burnIndex = this->BackwardSubPhases.back().getXindex_burnIndex();

            for (size_t entry = 0; entry < Xindex_burnIndex.size(); ++entry)
            {
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    Xindex_burnIndex[entry],
                    this->Gindices_BurnIndexSumConstraint);
            }

        }//end calcbounds_burnindex_sum_constraint()

        void MGAnDSMs_phase::calcbounds_virtual_propellant_tanks()
        {
            //encode virtual propellant tank variables

            //virtual chemical fuel
            this->Xlowerbounds->push_back(0.0);
            this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
            this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
            this->Xdescriptions->push_back(prefix + "virtual chemical fuel");
            size_t Xindex = this->Xdescriptions->size() - 1;

            //tell the spacecraft about this variable
            //global constraint
            this->mySpacecraft->appendGlobalChemicalFuelTank_Xindices(Xindex);
            this->mySpacecraft->appendGlobalChemicalFuelTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
            //stage constraint
            this->mySpacecraft->appendChemicalFuelTank_Xindices(this->stageIndex, Xindex);
            this->mySpacecraft->appendChemicalFuelTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());

            //virtual chemical oxidizer            
            this->Xlowerbounds->push_back(0.0);
            this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
            this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
            this->Xdescriptions->push_back(prefix + "virtual chemical oxidizer");
            Xindex = this->Xdescriptions->size() - 1;

            //tell the spacecraft about this variable
            //global constraint
            this->mySpacecraft->appendGlobalChemicalOxidizerTank_Xindices(Xindex);
            this->mySpacecraft->appendGlobalChemicalOxidizerTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
            //stage constraint
            this->mySpacecraft->appendChemicalOxidizerTank_Xindices(this->stageIndex, Xindex);
            this->mySpacecraft->appendChemicalOxidizerTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());
        }//end calcbounds_virtual_propellant_tanks()

        void MGAnDSMs_phase::calcbounds_deltav_contribution()
        {
            //this has dependencies on all of the DSM components
            //forward subphases
            for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
            {
                std::vector<size_t> Xindex_DSM_components = this->ForwardSubPhases[forwardSubPhaseIndex].getXindex_DSM_components();

                std::vector<size_t> subphase_Gindices_deltav_wrt_ForwardControl;

                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    size_t Xindex = Xindex_DSM_components[Vindex];
                    this->create_sparsity_entry(0,
                        Xindex,
                        subphase_Gindices_deltav_wrt_ForwardControl);
                }

                this->Gindices_deltav_wrt_ForwardControl.push_back(subphase_Gindices_deltav_wrt_ForwardControl);
            }


            //backward subphases - skip the first one since there is no control there
            for (size_t backwardSubPhaseIndex = 1; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
            {
                std::vector<size_t> Xindex_DSM_components = this->BackwardSubPhases[backwardSubPhaseIndex].getXindex_DSM_components();

                std::vector<size_t> subphase_Gindices_deltav_wrt_BackwardControl;

                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    size_t Xindex = Xindex_DSM_components[Vindex];
                    this->create_sparsity_entry(0,
                        Xindex,
                        subphase_Gindices_deltav_wrt_BackwardControl);
                }

                this->Gindices_deltav_wrt_BackwardControl.push_back(subphase_Gindices_deltav_wrt_BackwardControl);
            }
        }//end calcbounds_virtual_deltav_constraint()

        //************************************calcbounds methods
        void MGAnDSMs_phase::process_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_phase_left_boundary(X, Xindex, F, Findex, G, needG);

            this->process_phase_flight_time(X, Xindex, F, Findex, G, needG);

            this->process_phase_right_boundary(X, Xindex, F, Findex, G, needG);

            this->process_phase_main(X, Xindex, F, Findex, G, needG);
        }//end process_phase()

        void MGAnDSMs_phase::process_phase_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_virtual_propellant_tanks(X, Xindex, F, Findex, G, needG);

            this->process_forward_half_phase(X, Xindex, F, Findex, G, needG);

            this->process_backward_half_phase(X, Xindex, F, Findex, G, needG);

            this->process_match_point_constraints(X, Xindex, F, Findex, G, needG);

            if (this->number_of_distance_constraints > 0)
            {
                this->process_distance_constraints(X, Xindex, F, Findex, G, needG);
            }

            this->process_burnindex_sum_constraint(X, Xindex, F, Findex, G, needG);

            this->process_deltav_contribution(X, Xindex, F, Findex, G, needG);
        }//end process_phase_main()

        void MGAnDSMs_phase::process_forward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //set the tank states
            this->state_at_beginning_of_phase(8) = this->chemical_fuel_used; //TCM propellant, usually zero
            this->state_at_beginning_of_phase(9) = 0.0; //oxidizer

            //Step 2: clear dPropagatedState_dIndependentVariable - note that the left boundary NEVER has a derivative with respect to current phase flight time
            this->Forward_dPropagatedState_dIndependentVariable.assign_zeros();

            //process the forward subphases
            for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
            {
                this->ForwardSubPhases[forwardSubPhaseIndex].process_subphase(X, Xindex, F, Findex, G, needG);

                this->ForwardSubPhases[forwardSubPhaseIndex].process_maneuver_constraints(X, Xindex, F, Findex, G, needG);
            }
        }//end process_forward_half_phase()

        void MGAnDSMs_phase::process_backward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: set the tank states
            this->state_at_end_of_phase(8) = this->virtual_chemical_fuel_used;
            this->state_at_end_of_phase(9) = this->virtual_chemical_oxidizer_used;

            //Step 2: populate dPropagatedState_dIndependentVariable if we are using an integrated propagator AND the boundary had to propagate in any way
            this->Backward_dPropagatedState_dIndependentVariable.assign_zeros();

            //Step 3: process the backward subphases
            for (size_t backwardSubPhaseIndex = 0; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
            {
                this->BackwardSubPhases[backwardSubPhaseIndex].process_subphase(X, Xindex, F, Findex, G, needG);

                this->BackwardSubPhases[backwardSubPhaseIndex].process_maneuver_constraints(X, Xindex, F, Findex, G, needG);
            }
        }//end process_backward_half_phase()

        void MGAnDSMs_phase::process_match_point_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: contruct the SPTM chains and the forward and backward HPTM
            if (needG)
            {
                //Step 1.1: Forward SPTM chains
                {
                    for (int subPhaseIndex = this->numberOfForwardSubphases - 1; subPhaseIndex >= 0; --subPhaseIndex)
                    {
                        //this is done by multiplying ForwardSPTM on the right by the NEXT step's ForwardSPTM
                        if (subPhaseIndex == this->numberOfForwardSubphases - 1)
                            this->CumulativeForwardSPTM[subPhaseIndex] = this->ForwardSPTM[subPhaseIndex];
                        else
                            this->CumulativeForwardSPTM[subPhaseIndex] = this->CumulativeForwardSPTM[subPhaseIndex + 1] * this->ForwardSPTM[subPhaseIndex];

                        //IF we are using an integrator, then time derivatives must come from the LAST SPTM in the chain
						/*
                        if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
                        {
                            for (size_t i : {0, 1, 2, 3, 4, 5, 6, 8, 9})
                            {
                                this->CumulativeForwardSPTM[subPhaseIndex](i, 7) = this->ForwardSPTM[subPhaseIndex](i, 7);
                                this->CumulativeForwardSPTM[subPhaseIndex](i, 10) = this->ForwardSPTM[subPhaseIndex](i, 10);
                            }
                        }*/
                    }
                }

                //Step 1.2: Backward SPTM chains
                {
                    for (int subPhaseIndex = this->numberOfBackwardSubphases - 1; subPhaseIndex >= 0; --subPhaseIndex)
                    {
                        //this is done by multiplying BackwardSPTM on the right by the NEXT step's BackwardSPTM
                        if (subPhaseIndex == this->numberOfBackwardSubphases - 1)
                            this->CumulativeBackwardSPTM[subPhaseIndex] = this->BackwardSPTM[subPhaseIndex];
                        else
                            this->CumulativeBackwardSPTM[subPhaseIndex] = this->CumulativeBackwardSPTM[subPhaseIndex + 1] * this->BackwardSPTM[subPhaseIndex];

                        //IF we are using an integrator, then time derivatives must come from the LAST SPTM in the chain
						/*
                        if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
                        {
                            for (size_t i : {0, 1, 2, 3, 4, 5, 6, 8, 9})
                            {
                                this->CumulativeBackwardSPTM[subPhaseIndex](i, 7) = this->BackwardSPTM[subPhaseIndex](i, 7);
                                this->CumulativeBackwardSPTM[subPhaseIndex](i, 10) = this->BackwardSPTM[subPhaseIndex](i, 10);
                            }
                        }*/
                    }
                }

                //Step 1.3: form the HPTMs
                this->ForwardHPTM = this->CumulativeForwardSPTM.front() * this->TCMTM;
                this->BackwardHPTM = this->CumulativeBackwardSPTM.front();
            }

            //Step 2: call the base class
            TwoPointShootingPhase::process_match_point_constraints(X, Xindex, F, Findex, G, needG);
            

            if (needG)
            {
                //Step 4: derivatives with respect to control
                math::Matrix<double> dStateNow_dDecisionVariable(this->numStatesToPropagate + 2, 1, 0.0);
                math::Matrix<double> dMatchPointState_dDecisionVariable(this->numStatesToPropagate + 2, 1, 0.0);
                size_t constraintIndex_start = 0; //we need this to control for the fact that the last DSM before the match point does not affect position

                //Step 4.1: Forward Control
                for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
                {
                    math::Matrix<double> dMassAfterTCM_dDSMcomponents = this->ForwardSubPhases[forwardSubPhaseIndex].get_dMassAfterTCM_dDSMcomponents();
                    math::Matrix<double> dFuel_dDSMcomponents = this->ForwardSubPhases[forwardSubPhaseIndex].get_dFuel_dDSMcomponents();
                    math::Matrix<double> dOxidizer_dDSMcomponents = this->ForwardSubPhases[forwardSubPhaseIndex].get_dOxidizer_dDSMcomponents();
                    for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();
                        dStateNow_dDecisionVariable(3 + Vindex) = 1.0;
                        dStateNow_dDecisionVariable(6) = dMassAfterTCM_dDSMcomponents(Vindex);
                        dStateNow_dDecisionVariable(8) = dFuel_dDSMcomponents(Vindex);
                        dStateNow_dDecisionVariable(9) = dOxidizer_dDSMcomponents(Vindex);

                        if (forwardSubPhaseIndex == this->numberOfForwardSubphases - 1)
                        {
                            dMatchPointState_dDecisionVariable = dStateNow_dDecisionVariable;

                            //velocity
                            size_t stateIndex = this->matchPointConstraintStateIndex[3 + Vindex];
                            size_t Gindex = this->Gindices_match_point_constraint_ForwardControl[forwardSubPhaseIndex][3 + Vindex][Vindex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(3 + Vindex);

                            //masses, tanks, etc.
                            for (size_t constraintIndex = 6; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                            {
                                size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                                size_t Gindex = this->Gindices_match_point_constraint_ForwardControl[forwardSubPhaseIndex][constraintIndex][Vindex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                    * dMatchPointState_dDecisionVariable(stateIndex)
                                    * this->continuity_constraint_scale_factors(constraintIndex);
                            }//end loop over states
                        }
                        else
                        {
                            dMatchPointState_dDecisionVariable = this->CumulativeForwardSPTM[forwardSubPhaseIndex + 1] * dStateNow_dDecisionVariable;
                            
                            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                            {
                                size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                                size_t Gindex = this->Gindices_match_point_constraint_ForwardControl[forwardSubPhaseIndex][constraintIndex][Vindex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                    * dMatchPointState_dDecisionVariable(stateIndex)
                                    * this->continuity_constraint_scale_factors(constraintIndex);
                            }//end loop over states
                        }

                    }//end loop over DSM components
                }//end loop over forward substeps

                //Step 4.2: Backward Control
                for (size_t backwardSubPhaseIndex = 1; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
                {
                    math::Matrix<double> dMassBeforeDSM_dDSMcomponents = this->BackwardSubPhases[backwardSubPhaseIndex].get_dMassBeforeDSM_dDSMcomponents();
                    math::Matrix<double> dFuel_dDSMcomponents = this->BackwardSubPhases[backwardSubPhaseIndex].get_dFuel_dDSMcomponents();
                    math::Matrix<double> dOxidizer_dDSMcomponents = this->BackwardSubPhases[backwardSubPhaseIndex].get_dOxidizer_dDSMcomponents();
                    for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                    {
                        dStateNow_dDecisionVariable.assign_zeros();
                        dStateNow_dDecisionVariable(3 + Vindex) = 1.0;
                        dStateNow_dDecisionVariable(6) = dMassBeforeDSM_dDSMcomponents(Vindex);
                        dStateNow_dDecisionVariable(8) = dFuel_dDSMcomponents(Vindex);
                        dStateNow_dDecisionVariable(9) = dOxidizer_dDSMcomponents(Vindex);

                        dMatchPointState_dDecisionVariable = this->CumulativeBackwardSPTM[backwardSubPhaseIndex] * dStateNow_dDecisionVariable;

						if (backwardSubPhaseIndex == this->numberOfBackwardSubphases - 1)
						{
							dMatchPointState_dDecisionVariable(6) = dStateNow_dDecisionVariable(6);
							dMatchPointState_dDecisionVariable(8) = dStateNow_dDecisionVariable(8);
							dMatchPointState_dDecisionVariable(9) = dStateNow_dDecisionVariable(9);
						}
                        else
                        {
                            dMatchPointState_dDecisionVariable(6) /= this->BackwardSPTM[backwardSubPhaseIndex](6, 6);
                            dMatchPointState_dDecisionVariable(8) /= this->BackwardSPTM[backwardSubPhaseIndex](6, 6);
                            dMatchPointState_dDecisionVariable(9) /= this->BackwardSPTM[backwardSubPhaseIndex](6, 6);
                        }

                        for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                        {
                            size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                            size_t Gindex = this->Gindices_match_point_constraint_BackwardControl[backwardSubPhaseIndex][constraintIndex][Vindex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }//end loop over states
                    }//end loop over DSM components
                }//end loop over backward substeps

                //Step 5: Derivatives with respect to propagation time
                {
                    dStateNow_dDecisionVariable.assign_zeros();
                    dStateNow_dDecisionVariable(this->stateIndex_phase_propagation_variable) = 1.0;


                    //Forward
                    dMatchPointState_dDecisionVariable = this->ForwardHPTM * dStateNow_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                    {
                        if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                        {
                            size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                            size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime[constraintIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Backward
                    dMatchPointState_dDecisionVariable = this->BackwardHPTM * dStateNow_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                    {
                        if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                        {
                            size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                            size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime[constraintIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }

                //Step 6: Derivatives with respect to burn index
                dStateNow_dDecisionVariable.assign_zeros();
                dStateNow_dDecisionVariable(this->stateIndex_phase_propagation_variable + 1) = 1.0;

                //if we have more than one burn index, it gets more complicated...
                for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
                {
                    //multiply the first SPTM
                    dMatchPointState_dDecisionVariable = this->ForwardSPTM[forwardSubPhaseIndex] * dStateNow_dDecisionVariable;

                    //clear the burn index field because after the first multiply we are on someone else's burn index
                    dMatchPointState_dDecisionVariable(this->stateIndex_phase_propagation_variable + 1) = 0.0;

                    //finish the chain if there is more chaining to do
                    if (forwardSubPhaseIndex < this->numberOfForwardSubphases - 1)
                        dMatchPointState_dDecisionVariable = this->CumulativeForwardSPTM[forwardSubPhaseIndex + 1] * dMatchPointState_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                    {

                        if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                        {
                            size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                            size_t Gindex = this->Gindices_match_point_constraint_ForwardBurnIndex[forwardSubPhaseIndex][constraintIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }//end loop over states
                }//end loop over forward subphases

                for (size_t backwardSubPhaseIndex = 0; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
                {
                    //multiply the first SPTM
                    dMatchPointState_dDecisionVariable = this->BackwardSPTM[backwardSubPhaseIndex] * dStateNow_dDecisionVariable;
                        
                    //clear the burn index field because after the first multiply we are on someone else's burn index
                    dMatchPointState_dDecisionVariable(this->stateIndex_phase_propagation_variable + 1) = 0.0;

                    //finish the chain if there is more chaining to do
                    if (backwardSubPhaseIndex < this->numberOfBackwardSubphases - 1)
                        dMatchPointState_dDecisionVariable = this->CumulativeBackwardSPTM[backwardSubPhaseIndex + 1] * dMatchPointState_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                    {
                        if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                        {
                            size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                            size_t Gindex = this->Gindices_match_point_constraint_BackwardBurnIndex[backwardSubPhaseIndex][constraintIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }//end loop over states
                }//end loop over backward subphases
            }//end multiple DSM burn index derivatives


            //virtual chemical fuel
            {
                size_t Gindex = this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel;
                size_t Xindex = this->jGvar->operator[](Gindex);

                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * this->continuity_constraint_scale_factors(7);
            }

            //virtual chemical oxidizer
            {
                size_t Gindex = this->Gindex_dVirtualChemicalOxidizerConstraint_dVirtualChemicalOxidizer;
                size_t Xindex = this->jGvar->operator[](Gindex);

                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * this->continuity_constraint_scale_factors(8);
            }
        }//end derivatives

        void MGAnDSMs_phase::process_distance_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            std::vector< std::vector<doubleType> > body_state_check(this->numberOfDSMs, { 0.0, 0.0, 0.0 });

            //Step 1: apply the distance constraints
            for (size_t stepIndex = 0; stepIndex < this->numberOfDSMs; ++stepIndex)
            {
                //which state are we constraining?
                math::Matrix<doubleType>& state_vector = stepIndex < this->numberOfForwardSubphases
                    ? (stepIndex == this->numberOfForwardSubphases - 1 ? this->match_point_state_minus : this->SpacecraftStateForward[stepIndex])
                    : this->SpacecraftStateBackward[this->numberOfDSMs - stepIndex - 1];

                for (size_t body = 0; body < this->number_of_distance_constraints; ++body)
                {
                    //Step 1.1: get the position vector of the spacecraft relative to the body
                    if (std::get<0>(this->distance_constraint_definitions[body]) == -2)
                    {
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {
                            this->distance_constraint_relative_position[stepIndex][body](stateIndex) = state_vector(stateIndex);
                            this->distance_constraint_body_position_time_derivatives[stepIndex][body](stateIndex) = 0.0;
                        }
                    }
                    else
                    {
                        doubleType body_state[12];
                        this->myUniverse->bodies[std::get<0>(this->distance_constraint_definitions[body]) - 1].locate_body(
                            state_vector(7),
                            body_state,
                            (needG),
                            *this->myOptions);
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {
                            this->distance_constraint_relative_position[stepIndex][body](stateIndex) = state_vector(stateIndex) - body_state[stateIndex];
                            this->distance_constraint_body_position_time_derivatives[stepIndex][body](stateIndex) = body_state[6 + stateIndex] _GETVALUE;
                            body_state_check[stepIndex][stateIndex] = body_state[stateIndex];
                        }
                    }

                    //Step 1.2: apply the constraint
                    this->distance_from_body[stepIndex][body] = this->distance_constraint_relative_position[stepIndex][body].norm();
                    F[Findex] = (this->distance_from_body[stepIndex][body] - std::get<2>(this->distance_constraint_definitions[body])) / this->myUniverse->LU;
                    ++Findex;
                }
            }//end loop over steps

            //Step 2: derivatives
            if (needG)
            {
                //holder vectors
                math::Matrix<double> dStateAtDecision_dDecisionVariable(this->ForwardHPTM.get_n(), 1, 0.0);
                math::Matrix<double> dStateAtConstraint_dDecisionVariable(this->ForwardHPTM.get_n(), 1, 0.0);

                //references to boundary derivatives
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

                //Step 2.1: assemble chains of SPTMs
                for (size_t distanceStep = 0; distanceStep < this->numberOfDSMs; ++distanceStep)
                {
                    if (distanceStep < this->numberOfForwardSubphases)
                    {
                        for (int controlStep = distanceStep; controlStep >= 0; --controlStep)
                        {
                            if (controlStep == distanceStep)
                            {
                                this->STM_to_maneuver[distanceStep][controlStep] = this->ForwardSPTM[distanceStep];

                                //The maneuver is on the right-hand side of the step, so it hasn't happened yet. Remove the mass entries.
                                this->STM_to_maneuver[distanceStep][controlStep](6, 6) = 1.0;
                                this->STM_to_maneuver[distanceStep][controlStep](8, 6) = 1.0;
                                this->STM_to_maneuver[distanceStep][controlStep](9, 6) = 1.0;
                            }
                            else
                            {
                                this->STM_to_maneuver[distanceStep][controlStep] = this->STM_to_maneuver[distanceStep][controlStep + 1] * this->ForwardSPTM[controlStep];
                            }
                        }

                        //boundary STM
                        this->Boundary_STM_to_maneuver[distanceStep] = this->STM_to_maneuver[distanceStep][0] * this->TCMTM;
                    }//end forward
                    else
                    {
                        //back_burnIndex is the number of burns removed from the right-hand boundary
                        //note that the first backward subphase does not have a maneuver, so the "zeroth" backward maneuver actually corresponds to the second backward subphase
                        size_t back_burnIndex = this->numberOfDSMs - distanceStep;

                        for (int controlStep = distanceStep; controlStep <= this->numberOfDSMs; ++controlStep)
                        {
                            if (controlStep == distanceStep)
                            {
                                //backward subphases have their maneuver on the right hand side so the first "STM" is just an identify mapping
                                this->STM_to_maneuver[distanceStep][controlStep] = math::Matrix<double>(this->BackwardHPTM.get_n(), math::MatrixType::identity);
                            }
                            else
                            {
                                this->STM_to_maneuver[distanceStep][controlStep] = this->STM_to_maneuver[distanceStep][controlStep - 1] * this->BackwardSPTM[this->numberOfDSMs - controlStep];
                            }
                        }

                        //boundary STM
                        this->Boundary_STM_to_maneuver[distanceStep] = this->STM_to_maneuver[distanceStep][this->numberOfDSMs - 1] * this->BackwardSPTM[0];
                    }//end backward case
                }//end loop over distance steps

                //Step 2.2: compute and apply the derivatives themselves
                for (size_t distanceStep = 0; distanceStep < this->numberOfDSMs; ++distanceStep)
                {
                    if (distanceStep < this->numberOfForwardSubphases)
                    {
                        //control derivatives
                        for (size_t controlStep = 0; controlStep < distanceStep; ++controlStep)
                        {
                            Forward_MGAnDSMs_subphase& myManuever = this->ForwardSubPhases[controlStep];

                            math::Matrix<double> dMassAfterManeuver_dDeltaVcomponent = myManuever.get_dMassAfterTCM_dDSMcomponents();
                            math::Matrix<double> dChemicalFuel_dDeltaVcomponent = myManuever.get_dFuel_dDSMcomponents();
                            math::Matrix<double> dOxidizer_dDeltaVcomponent = myManuever.get_dOxidizer_dDSMcomponents();

                            for (size_t controlIndex : {0, 1, 2})
                            {

                                dStateAtDecision_dDecisionVariable.assign_zeros();
                                dStateAtDecision_dDecisionVariable(3 + controlIndex) = 1.0;
                                dStateAtDecision_dDecisionVariable(6) = dMassAfterManeuver_dDeltaVcomponent(controlIndex);
                                dStateAtDecision_dDecisionVariable(8) = dChemicalFuel_dDeltaVcomponent(controlIndex);
                                dStateAtDecision_dDecisionVariable(9) = dOxidizer_dDeltaVcomponent(controlIndex);

                                //forward maneuver occurs on the right-hand side of segment, so use the next step's STM
                                dStateAtConstraint_dDecisionVariable = this->STM_to_maneuver[distanceStep][controlStep + 1] * dStateAtDecision_dDecisionVariable;

                                for (size_t constraintIndex = 0; constraintIndex < this->number_of_distance_constraints; ++constraintIndex)
                                {
                                    size_t Gindex = this->G_indices_distance_constraints_wrt_ForwardControl[distanceStep][constraintIndex][controlStep][controlIndex];

                                    size_t Xindex = this->jGvar->operator[](Gindex);

                                    math::Matrix<double>& STM = this->STM_to_maneuver[distanceStep][controlStep];

                                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                        / this->distance_from_body[distanceStep][constraintIndex] _GETVALUE
                                        * (   this->distance_constraint_relative_position[distanceStep][constraintIndex](0) * dStateAtConstraint_dDecisionVariable(0)
                                            + this->distance_constraint_relative_position[distanceStep][constraintIndex](1) * dStateAtConstraint_dDecisionVariable(1)
                                            + this->distance_constraint_relative_position[distanceStep][constraintIndex](2) * dStateAtConstraint_dDecisionVariable(2)) _GETVALUE
                                        / this->myUniverse->LU;
                                }//end loop over constraints
                            }
                        }//end loop over control steps

                        //boundary state
                        for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByVariable.size(); ++varIndex)
                        {
                            dStateAtDecision_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfLeftBoundaryByVariable[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);

                                dStateAtDecision_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase[dIndex]);
                            }

                            dStateAtConstraint_dDecisionVariable = this->Boundary_STM_to_maneuver[distanceStep] * dStateAtDecision_dDecisionVariable;

                            for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                            {
                                size_t Gindex = this->G_indices_distance_constraints_wrt_LeftBoundaryState[distanceStep][bodyIndex][varIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                    / this->distance_from_body[distanceStep][bodyIndex] _GETVALUE
                                    * (   this->distance_constraint_relative_position[distanceStep][bodyIndex](0) * dStateAtConstraint_dDecisionVariable(0)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](1) * dStateAtConstraint_dDecisionVariable(1)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](2) * dStateAtConstraint_dDecisionVariable(2)) _GETVALUE
                                    / this->myUniverse->LU;
                            }
                        }

                        //boundary epoch
                        for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByTimeVariable.size(); ++varIndex)
                        {
                            dStateAtDecision_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByTimeVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfLeftBoundaryByTimeVariable[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                                dStateAtDecision_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);
                            }

                            dStateAtConstraint_dDecisionVariable = this->Boundary_STM_to_maneuver[distanceStep] * dStateAtDecision_dDecisionVariable;

                            double timeScale = dStateAtConstraint_dDecisionVariable(7);

                            for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                            {
                                size_t Gindex = this->G_indices_distance_constraints_wrt_LeftBoundaryTime[distanceStep][bodyIndex][varIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                    / this->distance_from_body[distanceStep][bodyIndex] _GETVALUE
                                    * (   this->distance_constraint_relative_position[distanceStep][bodyIndex](0) * (dStateAtConstraint_dDecisionVariable(0) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](0) * timeScale)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](1) * (dStateAtConstraint_dDecisionVariable(1) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](1) * timeScale)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](2) * (dStateAtConstraint_dDecisionVariable(2) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](2) * timeScale)) _GETVALUE
                                    / this->myUniverse->LU;
                            }
                        }

                        //phase flight time
                        {
                            dStateAtDecision_dDecisionVariable.assign_zeros();
                            dStateAtDecision_dDecisionVariable(this->stateIndex_phase_propagation_variable) = 1.0;

                            dStateAtConstraint_dDecisionVariable = this->Boundary_STM_to_maneuver[distanceStep] * dStateAtDecision_dDecisionVariable;

                            for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                            {
                                size_t Gindex = this->G_indices_distance_constraints_wrt_PhaseFlightTime[distanceStep][bodyIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                double timeScale = dStateAtConstraint_dDecisionVariable(7);

                                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                    / this->distance_from_body[distanceStep][bodyIndex] _GETVALUE
                                    * (   this->distance_constraint_relative_position[distanceStep][bodyIndex](0) * (dStateAtConstraint_dDecisionVariable(0) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](0) * timeScale)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](1) * (dStateAtConstraint_dDecisionVariable(1) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](1) * timeScale)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](2) * (dStateAtConstraint_dDecisionVariable(2) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](2) * timeScale)) _GETVALUE
                                    / this->myUniverse->LU;
                            }
                        }

                        //burn index
                        {
                            dStateAtDecision_dDecisionVariable.assign_zeros();
                            dStateAtDecision_dDecisionVariable(this->stateIndex_phase_propagation_variable + 1) = 1.0;

                            const auto& Xindex_burnIndex = this->ForwardSubPhases[distanceStep].getXindex_burnIndex();
                            for (size_t controlStep = 0; controlStep <= distanceStep; ++controlStep)
                            {
                                //multiply the first SPTM
                                dStateAtConstraint_dDecisionVariable = this->ForwardSPTM[controlStep] * dStateAtDecision_dDecisionVariable;

                                //if the burn index of interest is not associated with the step of interest, clear the burnIndex before you multiply the rest
                                if (distanceStep != controlStep)
                                {
                                    dStateAtConstraint_dDecisionVariable(this->stateIndex_phase_propagation_variable + 1) = 0.0;

                                    //multiply to the step of interest
                                    dStateAtConstraint_dDecisionVariable = this->STM_to_maneuver[distanceStep][controlStep + 1] * dStateAtConstraint_dDecisionVariable;
                                }

                                double timeScale = this->PhaseFlightTime _GETVALUE;

                                for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                                {
                                    const size_t Gindex = this->G_indices_distance_constraints_wrt_BurnIndex[distanceStep][bodyIndex][controlStep];
                                    const size_t Xindex = Xindex_burnIndex[controlStep];

                                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                        / this->distance_from_body[distanceStep][bodyIndex] _GETVALUE
                                        * (this->distance_constraint_relative_position[distanceStep][bodyIndex](0) * (dStateAtConstraint_dDecisionVariable(0) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](0) * timeScale)
                                            + this->distance_constraint_relative_position[distanceStep][bodyIndex](1) * (dStateAtConstraint_dDecisionVariable(1) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](1) * timeScale)
                                            + this->distance_constraint_relative_position[distanceStep][bodyIndex](2) * (dStateAtConstraint_dDecisionVariable(2) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](2) * timeScale)) _GETVALUE
                                        / this->myUniverse->LU;
                                }
                            }
                        }
                    }//end forward case
                    else
                    {
                        //back_burnIndex is the number of burns removed from the right-hand boundary
                        //note that the first backward subphase does not have a maneuver, so the "zeroth" backward maneuver actually corresponds to the second backward subphase
                        size_t back_maneuverIndex = this->numberOfDSMs - distanceStep - 1;

                        //control derivatives
                        for (size_t controlStep = distanceStep + 1; controlStep < this->numberOfDSMs; ++controlStep)
                        {
                            size_t controlBackStep = this->numberOfDSMs - controlStep - 1;

                            Backward_MGAnDSMs_subphase& myManuever = this->BackwardSubPhases[this->numberOfDSMs - controlStep];

                            math::Matrix<double> dMassBeforeManeuver_dDeltaVcomponent = myManuever.get_dMassBeforeDSM_dDSMcomponents();
                            math::Matrix<double> dChemicalFuel_dDeltaVcomponent = myManuever.get_dFuel_dDSMcomponents();
                            math::Matrix<double> dOxidizer_dDeltaVcomponent = myManuever.get_dOxidizer_dDSMcomponents();

                            for (size_t controlIndex : {0, 1, 2})
                            {

                                dStateAtDecision_dDecisionVariable.assign_zeros();
                                dStateAtDecision_dDecisionVariable(3 + controlIndex) = -1.0;
                                dStateAtDecision_dDecisionVariable(6) = -dMassBeforeManeuver_dDeltaVcomponent(controlIndex);
                                dStateAtDecision_dDecisionVariable(8) = -dChemicalFuel_dDeltaVcomponent(controlIndex);
                                dStateAtDecision_dDecisionVariable(9) = -dOxidizer_dDeltaVcomponent(controlIndex);

                                dStateAtConstraint_dDecisionVariable = this->STM_to_maneuver[distanceStep][controlStep] * dStateAtDecision_dDecisionVariable;

                                for (size_t constraintIndex = 0; constraintIndex < this->number_of_distance_constraints; ++constraintIndex)
                                {
                                    size_t Gindex = this->G_indices_distance_constraints_wrt_BackwardControl[back_maneuverIndex][constraintIndex][controlBackStep][controlIndex];

                                    size_t Xindex = this->jGvar->operator[](Gindex);

                                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                        / this->distance_from_body[distanceStep][constraintIndex] _GETVALUE
                                        * (   this->distance_constraint_relative_position[distanceStep][constraintIndex](0) * dStateAtConstraint_dDecisionVariable(0)
                                            + this->distance_constraint_relative_position[distanceStep][constraintIndex](1) * dStateAtConstraint_dDecisionVariable(1)
                                            + this->distance_constraint_relative_position[distanceStep][constraintIndex](2) * dStateAtConstraint_dDecisionVariable(2)) _GETVALUE
                                        / this->myUniverse->LU;
                                }//end loop over constraints
                            }
                        }//end loop over control steps

                        //boundary state
                        for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByVariable.size(); ++varIndex)
                        {
                            dStateAtDecision_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfRightBoundaryByVariable[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase[dIndex]);

                                dStateAtDecision_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterPhase[dIndex]);
                            }

                            dStateAtConstraint_dDecisionVariable = this->Boundary_STM_to_maneuver[distanceStep] * dStateAtDecision_dDecisionVariable;

                            for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                            {
                                size_t Gindex = this->G_indices_distance_constraints_wrt_RightBoundaryState[back_maneuverIndex][bodyIndex][varIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                    / this->distance_from_body[distanceStep][bodyIndex] _GETVALUE
                                    * (   this->distance_constraint_relative_position[distanceStep][bodyIndex](0) * dStateAtConstraint_dDecisionVariable(0)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](1) * dStateAtConstraint_dDecisionVariable(1)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](2) * dStateAtConstraint_dDecisionVariable(2)) _GETVALUE
                                    / this->myUniverse->LU;
                            }
                        }
                        //boundary epoch
                        for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByTimeVariable.size(); ++varIndex)
                        {
                            dStateAtDecision_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                                dStateAtDecision_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);
                            }

                            dStateAtConstraint_dDecisionVariable = this->Boundary_STM_to_maneuver[distanceStep] * dStateAtDecision_dDecisionVariable;

                            double timeScale = dStateAtConstraint_dDecisionVariable(7);

                            for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                            {
                                size_t Gindex = this->G_indices_distance_constraints_wrt_RightBoundaryTime[back_maneuverIndex][bodyIndex][varIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);
                                
                                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                    / this->distance_from_body[distanceStep][bodyIndex] _GETVALUE
                                    * (   this->distance_constraint_relative_position[distanceStep][bodyIndex](0) * (dStateAtConstraint_dDecisionVariable(0) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](0) * timeScale)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](1) * (dStateAtConstraint_dDecisionVariable(1) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](1) * timeScale)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](2) * (dStateAtConstraint_dDecisionVariable(2) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](2) * timeScale)) _GETVALUE
                                    / this->myUniverse->LU;
                            }
                        }

                        //phase flight time (this has to be a += because phase flight time is also a right boundary time variable)
                        {
                            dStateAtDecision_dDecisionVariable.assign_zeros();
                            dStateAtDecision_dDecisionVariable(this->stateIndex_phase_propagation_variable) = 1.0;

                            dStateAtConstraint_dDecisionVariable = this->Boundary_STM_to_maneuver[distanceStep] * dStateAtDecision_dDecisionVariable;

                            for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                            {
                                size_t Gindex = this->G_indices_distance_constraints_wrt_PhaseFlightTime[distanceStep][bodyIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                double timeScale = dStateAtConstraint_dDecisionVariable(7);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    / this->distance_from_body[distanceStep][bodyIndex] _GETVALUE
                                    * (   this->distance_constraint_relative_position[distanceStep][bodyIndex](0) * (dStateAtConstraint_dDecisionVariable(0) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](0) * timeScale)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](1) * (dStateAtConstraint_dDecisionVariable(1) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](1) * timeScale)
                                        + this->distance_constraint_relative_position[distanceStep][bodyIndex](2) * (dStateAtConstraint_dDecisionVariable(2) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](2) * timeScale)) _GETVALUE
                                    / this->myUniverse->LU;
                            }
                        }
                        //burn index
                        {
                            //distanceBackStep is the number of burns removed from the right-hand boundary
                            //note that the first backward subphase does not have a maneuver, so the "zeroth" backward maneuver actually corresponds to the second backward subphase
                            size_t distanceBackStep = this->numberOfDSMs - distanceStep - 1;

                            dStateAtDecision_dDecisionVariable.assign_zeros();
                            dStateAtDecision_dDecisionVariable(this->stateIndex_phase_propagation_variable + 1) = 1.0;

                            const auto& Xindex_burnIndex = this->BackwardSubPhases[back_maneuverIndex].getXindex_burnIndex();
                            for (size_t controlStep = distanceStep; controlStep < this->numberOfDSMs; ++controlStep)
                            {
                                size_t controlBackStep = this->numberOfDSMs - controlStep - 1;
                                size_t varIndex = Xindex_burnIndex.size() - controlBackStep - 1; //varIndex is index in Xindex_burnIndex. The current burnIndex is always the last one as we move backwards in the list

                                //multiply the first SPTM
                                dStateAtConstraint_dDecisionVariable = this->BackwardSPTM[controlBackStep] * dStateAtDecision_dDecisionVariable;

                                //if the burn index of interest is not associated with the step of interest, clear the burnIndex before you multiply the rest
                                if (distanceBackStep != controlBackStep)
                                {
                                    dStateAtConstraint_dDecisionVariable(this->stateIndex_phase_propagation_variable + 1) = 0.0;

                                    //multiply to the step of interest
                                    dStateAtConstraint_dDecisionVariable = this->STM_to_maneuver[distanceStep][controlStep] * dStateAtConstraint_dDecisionVariable;
                                }

                                double timeScale = this->PhaseFlightTime _GETVALUE;

                                for (size_t bodyIndex = 0; bodyIndex < this->number_of_distance_constraints; ++bodyIndex)
                                {
                                    const size_t Gindex = this->G_indices_distance_constraints_wrt_BurnIndex[distanceStep][bodyIndex][controlBackStep];
                                    const size_t Xindex = Xindex_burnIndex[controlBackStep];

                                    G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                                        / this->distance_from_body[distanceStep][bodyIndex] _GETVALUE
                                        * (   this->distance_constraint_relative_position[distanceStep][bodyIndex](0) * (dStateAtConstraint_dDecisionVariable(0) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](0) * timeScale)
                                            + this->distance_constraint_relative_position[distanceStep][bodyIndex](1) * (dStateAtConstraint_dDecisionVariable(1) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](1) * timeScale)
                                            + this->distance_constraint_relative_position[distanceStep][bodyIndex](2) * (dStateAtConstraint_dDecisionVariable(2) - this->distance_constraint_body_position_time_derivatives[distanceStep][bodyIndex](2) * timeScale)) _GETVALUE
                                        / this->myUniverse->LU;
                                }
                            }
                        }
                    }//end backward case
                }//end loop over steps
            }//end derivatives
        }//end process_distance_constraints()
        
        void MGAnDSMs_phase::process_burnindex_sum_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            doubleType SumOfBurnIndices = 0.0;
            for (size_t entry = 0; entry < this->Gindices_BurnIndexSumConstraint.size(); ++entry)
            {
                size_t Gindex = this->Gindices_BurnIndexSumConstraint[entry];
                size_t Xindex = this->jGvar->operator[](Gindex);

                SumOfBurnIndices += X[Xindex];

                G[Gindex] = this->X_scale_factors->operator[](Xindex);
            }

            F[Findex++] = SumOfBurnIndices - 1.0;
        }//end process_burnindex_sum_constraint()

        void MGAnDSMs_phase::process_virtual_propellant_tanks(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->virtual_chemical_fuel_used = X[Xindex++];
            this->virtual_chemical_oxidizer_used = X[Xindex++];
        }//end process_virtual_propellant_tanks()

        void MGAnDSMs_phase::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {

            //compute the total delta-v
            this->PhaseTotalDeterministicDeltav = this->myDepartureEvent->get_DeterministicDeltav()
                + this->myArrivalEvent->get_DeterministicDeltav();

            this->PhaseTotalStatisticalDeltav = this->myDepartureEvent->get_StatisticalDeltav()
                + this->myArrivalEvent->get_StatisticalDeltav()
                + this->initial_TCM_magnitude;

            //forward sub-phase contribution
            for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
            {
                this->PhaseTotalDeterministicDeltav += this->ForwardSubPhases[forwardSubPhaseIndex].getDSMmagnitude();
                this->PhaseTotalStatisticalDeltav += this->ForwardSubPhases[forwardSubPhaseIndex].getTCMmagnitude();
            }

            //backward sub-phase contribution
            for (size_t BackwardSubPhaseIndex = 0; BackwardSubPhaseIndex < this->numberOfBackwardSubphases; ++BackwardSubPhaseIndex)
            {
                this->PhaseTotalDeterministicDeltav += this->BackwardSubPhases[BackwardSubPhaseIndex].getDSMmagnitude();
                this->PhaseTotalStatisticalDeltav += this->BackwardSubPhases[BackwardSubPhaseIndex].getTCMmagnitude();
            }
            
            //derivatives
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV
                && needG)
            {
                //forward subphases
                for (size_t forwardSubPhaseIndex = 0; forwardSubPhaseIndex < this->numberOfForwardSubphases; ++forwardSubPhaseIndex)
                {
                    math::Matrix<doubleType> thisDSM = this->ForwardSubPhases[forwardSubPhaseIndex].getDSM();
                    doubleType DSM_magnitude = this->ForwardSubPhases[forwardSubPhaseIndex].getDSMmagnitude();

                    for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                    {
                        size_t Gindex = this->Gindices_deltav_wrt_ForwardControl[forwardSubPhaseIndex][Vindex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * (thisDSM(Vindex) / DSM_magnitude) _GETVALUE;
                    }
                }

                //backward subphases - skip the first one because it has no maneuver
                for (size_t backwardSubPhaseIndex = 1; backwardSubPhaseIndex < this->numberOfBackwardSubphases; ++backwardSubPhaseIndex)
                {
                    math::Matrix<doubleType> thisDSM = this->BackwardSubPhases[backwardSubPhaseIndex].getDSM();
                    doubleType DSM_magnitude = this->BackwardSubPhases[backwardSubPhaseIndex].getDSMmagnitude();

                    for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                    {
                        size_t Gindex = this->Gindices_deltav_wrt_BackwardControl[backwardSubPhaseIndex - 1][Vindex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * (thisDSM(Vindex) / DSM_magnitude) _GETVALUE;
                    }
                }
            }
        }

        void MGAnDSMs_phase::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            //Step 1: output the departure event
            this->myDepartureEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 2: output the initial TCM if applicable
            this->output_initial_TCM(outputfile, eventcount);

            //Step 3: output the main part of the phase
            this->output_phase_main(outputfile, eventcount);

            //Step 4: print the arrival event
            this->myArrivalEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
        }


        void MGAnDSMs_phase::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
            Astrodynamics::SpacecraftAccelerationModel mySpacecraftAccelerationModel(this->myOptions,
                this->myJourneyOptions,
                this->myUniverse,
                this->Xdescriptions,
                this->mySpacecraft,
                11); // STM size

            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 1: output the departure event
            this->myDepartureEvent->output_ephemeris(outputfile);
            
            if (this->myOptions->generate_acceleration_model_instrumentation_file && this->myPropagatorType == PropagatorType::IntegratedPropagator)
            {
                this->temp_state = this->myDepartureEvent->get_state_before_event();
                mySpacecraftAccelerationModel.populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
            }

            //Step 2: print the forward subphases
            for (size_t subPhaseIndex = 0; subPhaseIndex < this->numberOfForwardSubphases; ++subPhaseIndex)
                this->ForwardSubPhases[subPhaseIndex].output_ephemeris(outputfile, acceleration_model_file);

            //Step 3: print the backward subphases
            for (size_t subPhaseIndex = 0; subPhaseIndex < this->numberOfBackwardSubphases; ++subPhaseIndex)
            {
                size_t backSubPhaseIndex = this->numberOfBackwardSubphases - subPhaseIndex - 1;
                this->BackwardSubPhases[backSubPhaseIndex].output_ephemeris(outputfile, acceleration_model_file);
            }

            //Step 4: print the arrival event
            this->myArrivalEvent->output_ephemeris(outputfile);
            if (this->myOptions->generate_acceleration_model_instrumentation_file && this->myPropagatorType == PropagatorType::IntegratedPropagator)
            {
                this->temp_state = this->myArrivalEvent->get_state_before_event();
                mySpacecraftAccelerationModel.populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
            }
        }//end output_ephemeris()

        void MGAnDSMs_phase::output_STMs()
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::ofstream summaryfile(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + ".stm_description");
            for (size_t stateIndex = 0; stateIndex < this->stateVectorNames.size(); ++stateIndex)
                summaryfile << this->stateVectorNames[stateIndex] << std::endl;
            summaryfile.close();

            size_t STMindex = 0;
            for (size_t subPhaseIndex = 0; subPhaseIndex < this->numberOfForwardSubphases; ++subPhaseIndex)
                this->ForwardSPTM[subPhaseIndex].print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_" + std::to_string(STMindex++) + "_forward.stm");
            for (int subPhaseIndex = this->numberOfBackwardSubphases - 1; subPhaseIndex > 0; --subPhaseIndex)
                this->BackwardSPTM[subPhaseIndex].print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_" + std::to_string(STMindex++) + "_backward.stm");
        }//end output_STMs()

        void MGAnDSMs_phase::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //should write out a maneuver/target spec for the departure maneuver if there is one
            this->myDepartureEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);

            //forward phase maneuver and target spec
            for (size_t subPhaseIndex = 0; subPhaseIndex < this->numberOfForwardSubphases; ++subPhaseIndex)
                this->ForwardSubPhases[subPhaseIndex].output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);

            //backward phase maneuver and target spec
            for (size_t subPhaseIndex = 0; subPhaseIndex < this->numberOfBackwardSubphases; ++subPhaseIndex)
            {
                size_t backSubPhaseIndex = this->numberOfBackwardSubphases - subPhaseIndex - 1;
                this->BackwardSubPhases[backSubPhaseIndex].output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
            }

            //should write out a maneuver/target spec for the arrival maneuver if there is one
            this->myArrivalEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
        }//end output_maneuver_and_target_spec()

        void MGAnDSMs_phase::output_phase_main(std::ofstream& outputfile,
            size_t& eventcount)
        {
            //Step 1: print the forward subphases
            for (size_t subPhaseIndex = 0; subPhaseIndex < this->numberOfForwardSubphases; ++subPhaseIndex)
            {
                this->ForwardSubPhases[subPhaseIndex].output(outputfile, eventcount);
            }

            //Step 2: output match point
            math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
            if (this->myUniverse->central_body_SPICE_ID == 10)
            {
                R_sc_Sun = this->match_point_state_minus.getSubMatrix1D(0, 2);
            }
            else
            {
                //where is the central body relative to the sun?
                doubleType central_body_state_and_derivatives[12];
                this->myUniverse->locate_central_body(this->match_point_state_minus(7),
                    central_body_state_and_derivatives,
                    *this->myOptions,
                    false);

                math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                }

                R_sc_Sun = this->match_point_state_minus.getSubMatrix1D(0, 2) + R_CB_Sun;
            }
            doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
            this->mySpacecraft->computePowerState(r_sc_sun_AU, this->match_point_state_minus(7));

            doubleType power = this->mySpacecraft->getAvailablePower();
            doubleType active_power = this->mySpacecraft->getEPActivePower();

            phase::write_output_line(outputfile,//outputfile
                eventcount,//eventcount
                "match_point",//event_type
                "deep-space",//event_location
                0.0,// timestep_size,
                -1,//flyby_altitude,
                0,//BdotR
                0,//BdotT
                0,//angle1
                0,//angle2
                0,//C3
                this->match_point_state_minus,//state
                math::Matrix<doubleType>(3, 1, 0.0),//dV
                math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                0.0,//dVmag
                0.0,//Thrust
                0.0,//Isp
                power,//AvailPower
                0.0,//mdot
                0,//number_of_active_engines
                active_power,
                "none");//active_power)

            //Step 3: print the backward subphases
            for (size_t subPhaseIndex = 0; subPhaseIndex < this->numberOfBackwardSubphases; ++subPhaseIndex)
            {
                size_t backSubPhaseIndex = this->numberOfBackwardSubphases - subPhaseIndex - 1;
                this->BackwardSubPhases[backSubPhaseIndex].output(outputfile, eventcount);
            }
        }//end output_phase_main()
    }//end namespace Phases
}//end namespace EMTG