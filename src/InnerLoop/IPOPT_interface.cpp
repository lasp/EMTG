// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2020 United States Government as represented by the
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

//IPOPT interface
#ifndef NOIPOPT

#include "IPOPT_interface.h"
#include <time.h>
#include <cmath>

namespace EMTG
{
    namespace Solvers
    {
        IPOPT_interface::IPOPT_interface(problem* myProblem,
            const NLPoptions& myOptions)
		{
			this->initialize(myProblem,myOptions);
		}
			
        void IPOPT_interface::initialize(problem* myProblem,
            const NLPoptions& myOptions)
		{
			this->initialize_base(myProblem,myOptions);
			
			ipopt = new Ipopt::IpoptApplication(true,false);
		    ipopt->Options()->SetNumericValue("tol", this->myOptions.get_feasibility_tolerance());
			ipopt->Options()->SetStringValue("mu_strategy", "adaptive");
			ipopt->Options()->SetStringValue("output_file", myProblem->options.working_directory + "/ipopt.out");
			ipopt->Options()->SetStringValue("hessian_approximation", "limited-memory");
			if (myProblem->options.NLP_solver_type == 2)
				ipopt->Options()->SetStringValue("linear_solver", "ma27");
			else if (myProblem->options.NLP_solver_type == 3)
				ipopt->Options()->SetStringValue("linear_solver", "mumps");
			ipopt->Options()->SetNumericValue("max_wall_time", this->myOptions.get_max_run_time_seconds());
			ipopt->Options()->SetIntegerValue("max_iter", this->myOptions.get_major_iterations_limit());
			
			
			if (myProblem->options.quiet_NLP)
				ipopt->Options()->SetIntegerValue("print_level",0);
			else
				ipopt->Options()->SetIntegerValue("print_level",5);
			
			ipopt->Initialize();
			
            //set first feasibility flag
            this->first_feasibility = false;
			this->need_reinitialization = false;
		}

        void IPOPT_interface::run_NLP(const bool& X0_is_scaled)
        {	
            //compute the scaled decision vector
            if (!X0_is_scaled)
                this->scaleX0();
            else
                this->unscaleX0();

            this->X_scaled = this->X0_scaled;

            //check that the user function works, and return values of F
            this->myProblem->evaluate(this->X0_unscaled, this->F, this->myProblem->G, false);

            //set the bounds equal to the problem bounds
            for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
            {
                this->Xupperbounds[Xindex] = (this->myProblem->Xupperbounds[Xindex] - this->myProblem->Xlowerbounds[Xindex]) / this->myProblem->X_scale_factors[Xindex];
                this->Xlowerbounds[Xindex] = 0.0;
            }

            //run IPOPT
            this->NLP_start_time = time(NULL);
            this->feasibility_metric_NLP_incumbent = 1.0e+101;
            this->J_NLP_incumbent = math::LARGE;
            this->movie_frame_count = 0;

			if (this->need_reinitialization)
			{
				ipopt->Options()->SetStringValue("warm_start_same_structure","yes");
				ipopt->ReOptimizeTNLP(this);
			}
			else
			{
				ipopt->OptimizeTNLP(this);
				this->need_reinitialization = true;
			}
			
			
            //unscale the various things that might need to be uncaled
            this->unscaleX();

            //perform a feasibility check
            this->myProblem->check_feasibility(this->X_unscaled,
                this->F,
                this->worst_decision_variable,
                this->worst_constraint,
                this->feasibility_metric,
                this->normalized_feasibility_metric,
                this->distance_from_equality_filament,
                this->decision_vector_feasibility_metric);

            double worst_feasibility = fmax(normalized_feasibility_metric, decision_vector_feasibility_metric);

            //adopt the incumbent if appropriate
            if (this->myOptions.get_enable_NLP_chaperone())
            {
                this->unscaleX_NLP_incumbent();
                if (worst_feasibility < this->myOptions.get_feasibility_tolerance()
                    && this->feasibility_metric_NLP_incumbent < this->myOptions.get_feasibility_tolerance())
                {
                    // Both incumbent point and NLP exit point are feasible
                    if (this->J_NLP_incumbent < this->F.front())
                    {
                        this->X_unscaled = this->X_NLP_incumbent_unscaled;
                        this->X_scaled = this->X_NLP_incumbent_scaled;
                        this->F = this->F_NLP_incumbent;

                        std::cout << "NLP incumbent point and exit point are feasible and incumbent point is superior to exit point." << std::endl;
                    }
                    else
                    {
                        std::cout << "NLP incumbent point and exit point are feasible and exit point is superior to incumbent point." << std::endl;
                    }
                }
                //if the incumbent point is feasible but the exit point is not
                else if (this->feasibility_metric_NLP_incumbent < worst_feasibility
                    && this->feasibility_metric_NLP_incumbent < this->myOptions.get_feasibility_tolerance())
                {
                    this->X_unscaled = this->X_NLP_incumbent_unscaled;
                    this->X_scaled = this->X_NLP_incumbent_scaled;
                    this->F = this->F_NLP_incumbent;

                    std::cout << "NLP incumbent point was feasible and the exit point was not." << std::endl;
                }
                else if (this->feasibility_metric_NLP_incumbent < worst_feasibility)
                {
                    this->X_unscaled = this->X_NLP_incumbent_unscaled;
                    this->X_scaled = this->X_NLP_incumbent_scaled;
                    this->F = this->F_NLP_incumbent;

                    std::cout << "NLP incumbent point was infeasible but less infeasible than the exit point." << std::endl;
                }
                else
                {
                    std::cout << "NLP exit point was superior to incumbent point." << std::endl;
                }
            }//end NLP chaperone
        }
		
		bool IPOPT_interface::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
		                                   Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style) 
		{
			n = this->nX;
			m = this->nF;
	
			// Is this right??
			nnz_jac_g = this->nG;
	
			index_style = Ipopt::TNLP::C_STYLE;	
    
			return true;
		}

		bool IPOPT_interface::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
		                     Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) 
		{
			for (int iX = 0; iX < n; iX++)
			{
				x_l[iX] = this->Xlowerbounds[iX];
				x_u[iX] = this->Xupperbounds[iX];
			}
			for (int iF = 1; iF < m; iF++)
			{
				g_l[iF] = this->Flowerbounds[iF];
				g_u[iF] = this->Fupperbounds[iF];
			}
			return true;
		}

		bool IPOPT_interface::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
		                                bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
		                                Ipopt::Index m, bool init_lambda,
		                                Ipopt::Number* lambda)
		{
			for (int iX = 0; iX < n; iX++)
				x[iX] = this->X0_scaled[iX];
			return true;
		}


		bool IPOPT_interface::evaluate(const Ipopt::Number* x)
		{
            //Step 2: unwrap the inputs
            for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                this->X_scaled[Xindex] = x[Xindex];

            //Step 3: process the IPOPT decision vector
            this->unscaleX();

            //Step 4: Call the objective and constraint function

            if (this->myOptions.get_SolverMode() == NLPMode::FilamentFinder) //single-objective unconstrained, special Jacobian
            {

                this->myProblem->evaluate(this->X_unscaled, this->myProblem->F, this->myProblem->G, true);
                
                //then sum up the squares of the equality constraint values
                this->F.front() = 0.0;
                for (size_t Findex = 1; Findex < this->myProblem->total_number_of_constraints; ++Findex)
                {
                    if (this->myProblem->F_equality_or_inequality[Findex - 1])
                        this->F.front() += this->myProblem->F[Findex] * this->myProblem->F[Findex];
                }


                //compute dfdX
                for (size_t Gindex = 0; Gindex < this->nG; ++Gindex)
                {
                    this->G[Gindex] = 0.0;
                }

                size_t Problem_nG = this->myProblem->Gdescriptions.size();
                for (size_t Gindex = 0; Gindex < Problem_nG; ++Gindex)
                {
                    size_t Findex = this->myProblem->iGfun[Gindex];
                    size_t Xindex = this->myProblem->jGvar[Gindex];

                    if (Findex > 0)
                    {
                        if (this->myProblem->F_equality_or_inequality[Findex - 1])
                            this->G[Xindex] += 2.0 * this->myProblem->F[Findex] * this->myProblem->G[Gindex];
                    }
                }

                //inequality constraints
                for (size_t Findex = 1; Findex <= this->myProblem->F_indices_of_filament_critical_inequality_constraints.size(); ++Findex)
                {
                    this->F[Findex] = this->myProblem->F[this->myProblem->F_indices_of_filament_critical_inequality_constraints[Findex - 1]];
                }
                for (size_t Gindex = this->nX; Gindex < this->nG; ++Gindex)
                {
                    this->G[Gindex] = this->myProblem->G[this->original_G_indices_of_filament_critical_inequality_constraints[Gindex - this->nX]];
                }

            }
            else //regular problem
            {
                this->myProblem->evaluate(this->X_unscaled, this->F, this->G, true);
            }

            //Step 6: NLP chaperone
            if (true && this->myOptions.get_enable_NLP_chaperone())
            {
                //first compute the current feasibility
                double f_current;
                double f_abs_current;
                double distance_from_equality_filament;
                size_t worst_constraint;
                size_t worst_decision_variable;
                double decision_variable_feasibility_metric;

                try
                {
                    this->myProblem->check_feasibility(this->X_unscaled,
                        this->F,
                        worst_decision_variable,
                        worst_constraint,
                        f_abs_current,
                        f_current,
                        distance_from_equality_filament,
                        decision_variable_feasibility_metric,
                        true);
                }
                catch (std::runtime_error runtime_error)
                {

                    if (!this->myOptions.get_quiet_NLP())
                    {
                        std::cout << runtime_error.what() << std::endl;
                        std::cout << std::endl;
                    }
                   return false;
                }


                //stop immediately if goal met
                if (this->myOptions.get_stop_on_goal_attain())
                {
					std::cout << "snoptf0:"<<this->myProblem->getUnscaledObjective() << std::endl;
					std::cout << "goal:"<<this->myOptions.get_objective_goal() << std::endl;
                    if (f_current < this->myOptions.get_feasibility_tolerance() && this->myProblem->getUnscaledObjective() < this->myOptions.get_objective_goal())
                    {						
                        if (!this->myOptions.get_quiet_NLP())
                            std::cout << "NLP goal satisfied, exiting NLP" << std::endl;
                    return false;
                    }
                }

                //if the current point is feasible and better than the incumbent
                //raw pointer math for speed
                if (fmax(f_current, decision_variable_feasibility_metric) < this->myOptions.get_feasibility_tolerance()
                    && this->feasibility_metric_NLP_incumbent < this->myOptions.get_feasibility_tolerance()
                    && decision_variable_feasibility_metric == 0.0)
                {
                    // Both the incumbent and the current point are feasible
                    if (this->F.front() < this->J_NLP_incumbent)
                    {
                        // The current point is more optimal than the incumbent, so update the incumbent
                        this->J_NLP_incumbent = this->F.front();
                        this->feasibility_metric_NLP_incumbent = fmax(f_current, decision_variable_feasibility_metric);

                        this->X_NLP_incumbent_unscaled = this->X_unscaled;
                        this->X_NLP_incumbent_scaled = this->X_scaled;
                        this->F_NLP_incumbent = this->F;
                    }
					
					// Check if this is the first time things have been feasible. It shouldnt be, because the incumbent was already feasible
					// But do it anyway, even though this code should never be called.
					if (this->first_feasibility == false)
					{
						this->first_feasibility = true;
						this->myProblem->Xopt = this->X_NLP_incumbent_unscaled;
						this->myProblem->F = this->F_NLP_incumbent;
	                    this->myProblem->evaluate(this->X_NLP_incumbent_unscaled, this->F_NLP_incumbent, this->G, false);

                        this->myProblem->what_the_heck_am_I_called(SolutionOutputType::SUCCESS);
						this->myProblem->output(this->myProblem->options.outputfile);
					}
                }
                //if the current point is more feasible than the incumbent
                else if (fmax(f_current, decision_variable_feasibility_metric) < this->feasibility_metric_NLP_incumbent)
                {
                    this->J_NLP_incumbent = this->F.front();
                    this->feasibility_metric_NLP_incumbent = fmax(f_current, decision_variable_feasibility_metric);

                    this->X_NLP_incumbent_unscaled = this->X_unscaled;
                    this->X_NLP_incumbent_scaled = this->X_scaled;
                    this->F_NLP_incumbent = this->F;
					
					if (fmax(f_current, decision_variable_feasibility_metric) < this->myOptions.get_feasibility_tolerance() && this->first_feasibility == false)
					{
						this->first_feasibility = true;
						this->myProblem->Xopt = this->X_NLP_incumbent_unscaled;
						this->myProblem->F = this->F_NLP_incumbent;
	                    this->myProblem->evaluate(this->X_NLP_incumbent_unscaled, this->F_NLP_incumbent, this->G, false);

                        this->myProblem->what_the_heck_am_I_called(SolutionOutputType::SUCCESS);
                        this->myProblem->output(this->myProblem->options.outputfile);
					}
                }
            }//end NLP chaperone

            //Step 7: print movie frames
            if (true && this->myOptions.get_print_NLP_movie_frames())
            {
                this->myProblem->X = this->X_unscaled;
                this->myProblem->F = this->F;
                this->myProblem->output_problem_bounds_and_descriptions(this->myProblem->options.working_directory + "//" + "NLP_frame_" + std::to_string(this->movie_frame_count++) + ".csv");
            }

            // //Step 8: check the timer and kill IPOPT if necessary
            // time_t now = time(NULL);
            // if (now - this->NLP_start_time > this->myOptions.get_max_run_time_seconds())
            // {
            //     if (!this->myOptions.get_quiet_NLP())
            //         std::cout << "Exceeded NLP time limit of " << this->myOptions.get_max_run_time_seconds() << " seconds. Aborting NLP run." << std::endl;
            //     return false;
            // }

            //Step 9: do debug things if appropriate
#ifdef DEBUG_CHECK_VARIABLES_AND_CONSTRAINTS_IN_USER_FUNCTION
            if (this->myProblem->check_inputs_for_invalid_entries(x, false))
            {
                std::cout << "Invalid entry in decision vector" << endl;
                return false;
            }
            else if (this->myProblem->check_outputs_for_invalid_entries(Problem->X.data(), F, false))
            {
                std::cout << "Invalid entry in objective/constraint vector" << endl;
                return false;
            }
            else if (this->myProblem->check_derivatives_for_invalid_entries(Problem->X.data(), G, false))
            {
                std::cout << "Invalid entry in Jacobian" << endl;
                return false;
            }
#endif

			return true;
        }
		
		
		bool IPOPT_interface::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
		{
			if (new_x)
			{
				bool result = evaluate(x);		
				if (!result)
					return result;
			}
	
			obj_value = this->F[0];
			
			return true;
		}

		bool IPOPT_interface::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
		{
			if (new_x)
			{
				bool result = evaluate(x);		
				if (!result)
					return result;
			}
	

			// Put the constraint vector into a c-style array for SNOPT
			for (int i = 1; i < this->nF; i++)
				g[i] = this->F[i];

			return true;
		}
		
		bool IPOPT_interface::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
		{	
			if (new_x)
			{
				bool result = evaluate(x);		
				if (!result)
					return result;
			}
		
			for (int iX = 0; iX < n; iX++)
				grad_f[iX] = 0;
			
			for (int iG = 0; iG < this->nG; iG++)
				if (this->iGfun[iG] == 0)
					grad_f[this->jGvar[iG]] = this->G[iG];
			
			return true;
		}

		bool IPOPT_interface::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol, Ipopt::Number* values)
		{	
			if (values == NULL)
			{
				for (int iG = 0; iG < this->nG; iG++)
				{
					if (this->iGfun[iG] != 0)
					{
						iRow[iG] = this->iGfun[iG];
						jCol[iG] = this->jGvar[iG];
					}
				}
			}
			else
			{
				if (new_x)
				{
					bool result = evaluate(x);		
					if (!result)
						return result;
				}
				
				// Put the jacobian into a c-style array that IPOPT expects
				for (int iG = 0; iG < this->nG; iG++)
				{
					if (this->iGfun[iG] != 0)
					{
						values[iG] = this->G[iG];
						if (values[iG] != values[iG])
							std::cout << "NaN alert!" << std::endl;
					}
				} 
			}
			return true;
		}

		void IPOPT_interface::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
		{
			for (int iX = 0; iX < n; iX++)
				this->X_scaled[iX] = x[iX];
			return;
		}
  
    }//end namespace Solvers
}//end namespace EMTG

#endif // NOIPOPT