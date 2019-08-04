#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <limits>
#include "particleupdatesettings.h"
#include "particleswarm.h"
#include "particle.h"
#include "problem.h"
#include "topologymanager.h"
#include "particleupdatemanager.h"
#include "particleswarmsuite.h"
#include "differentialevolution.h"
#include "genome.h"
#include "crossovermanager.h"
#include "hybridalgorithm.h"
#include "hybridsuite.h"
#include "coco.h"
#include "desuite.h"

/**
 * A function type for evaluation functions, where the first argument is the vector to be evaluated and the
 * second argument the vector to which the evaluation result is stored.
 */
typedef void (*evaluate_function_t)(const double *x, double *y);

/**
 * A pointer to the problem to be optimized (needed in order to simplify the interface between the optimization
 * algorithm and the COCO platform).
 */
static coco_problem_t *PROBLEM;

/**
 * Calls coco_evaluate_function() to evaluate the objective function
 * of the problem at the point x and stores the result in the vector y
 */
static void evaluate_function(const double *x, double *y) {
  coco_evaluate_function(PROBLEM, x, y);
}

void experimentPSO(ParticleSwarm pso);

void experimentDE(DifferentialEvolution de);

void experimentHybrid(HybridAlgorithm ha);

static const long INDEPENDENT_RESTARTS = 1e5;
static const unsigned int BUDGET_MULTIPLIER = 1e4;

/* Structure and functions needed for timing the experiment */
typedef struct {
	size_t number_of_dimensions;
	size_t current_idx;
	char **output;
	size_t previous_dimension;
	size_t cumulative_evaluations;
	time_t start_time;
	time_t overall_start_time;
} 

timing_data_t;

static timing_data_t *timing_data_initialize(coco_suite_t *suite);
static void timing_data_time_problem(timing_data_t *timing_data, coco_problem_t *problem);
static void timing_data_finalize(timing_data_t *timing_data);

int main(void) {
  /* Single algorithm experiment examples */
  
  DifferentialEvolution de (RANDOM, TTB_1, BINOMIAL, false);
  experimentDE(de);
  
  ParticleSwarm ps(DECR_INERTIA_WEIGHT, VON_NEUMANN, ASYNCHRONOUS);
  experimentPSO(ps);

  ParticleSwarm ps2(INERTIA_WEIGHT, VON_NEUMANN, ASYNCHRONOUS);
  experimentPSO(ps2);

  HybridAlgorithm ha(INERTIA_WEIGHT, VON_NEUMANN, SYNCHRONOUS, BEST_1, BINOMIAL);
  experimentHybrid(ha);

  /* PSO suite experiment example */
  std::map<int,double> updateSettings;
  updateSettings[S_INER_W] = 0.8;
  /* Add mode settings if desired */

  ParticleSwarmSuite suite(updateSettings);
  suite.setUpdateManagers(std::vector<UpdateManagerType>({INERTIA_WEIGHT, DECR_INERTIA_WEIGHT, BARE_BONES, FIPS}));
   for (ParticleSwarm s : suite){
     experimentPSO(s);
  }

  return 0;
}

void experimentHybrid(HybridAlgorithm ha) {
  std::map<int, double> updateSettings;
  double const Cr = 0.7;
  double const F = 0.5;

  size_t run;
  coco_suite_t *suite;
  coco_observer_t *observer;
  timing_data_t *timing_data;
  char *observer_options =
      coco_strdupf("result_folder: %s " 
                   "algorithm_name: %s "
                   "", ha.getIdString().c_str(), ha.getIdString().c_str());

  suite = coco_suite("bbob", "instances: 1,1,1,1,1", "dimensions: 5");
  observer = coco_observer("bbob", observer_options);
  coco_free_memory(observer_options);

  timing_data = timing_data_initialize(suite);

  while ((PROBLEM = coco_suite_get_next_problem(suite, observer)) != NULL) {
    
    size_t dimension = coco_problem_get_dimension(PROBLEM);

    for (run = 1; run <= 1 + INDEPENDENT_RESTARTS; run++) {

      long evaluations_done = (long) (coco_problem_get_evaluations(PROBLEM) + 
            coco_problem_get_evaluations_constraints(PROBLEM));
      long evaluations_remaining = (long) (dimension * BUDGET_MULTIPLIER) - evaluations_done;       

      if ((coco_problem_final_target_hit(PROBLEM) && 
           coco_problem_get_number_of_constraints(PROBLEM) == 0)
           || (evaluations_remaining <= 0))
        break;

      Problem const problem(evaluate_function, PROBLEM);        
      ha.run(problem, evaluations_remaining, 6 * dimension, updateSettings, F, Cr); 

      if ((long int)coco_problem_get_evaluations(PROBLEM) == evaluations_done) {
        printf("WARNING: Budget has not been exhausted (%lu/%lu evaluations done)!\n",
            (unsigned long) evaluations_done, (unsigned long) dimension * BUDGET_MULTIPLIER);
        break;
      }
      else if ((long int)coco_problem_get_evaluations(PROBLEM) < evaluations_done)
        coco_error("Something unexpected happened - function evaluations were decreased!");
      }

      timing_data_time_problem(timing_data, PROBLEM);
  }

  timing_data_finalize(timing_data);

  coco_observer_free(observer);
  coco_suite_free(suite);  

}

void experimentDE(DifferentialEvolution de) {

  double const Cr = 0.7;
  double const F = 0.5;
  

  size_t run;
  coco_suite_t *suite;
  coco_observer_t *observer;
  timing_data_t *timing_data;
  char *observer_options =
      coco_strdupf("result_folder: %s " 
                   "algorithm_name: %s "
                   "", de.getIdString().c_str(), de.getIdString().c_str());

  suite = coco_suite("bbob", "instances: 1,1,1,1,1", "dimensions: 5");
  observer = coco_observer("bbob", observer_options);
  coco_free_memory(observer_options);

  timing_data = timing_data_initialize(suite);

  while ((PROBLEM = coco_suite_get_next_problem(suite, observer)) != NULL) {
    
    size_t dimension = coco_problem_get_dimension(PROBLEM);
    int const popSize = dimension * 8; 

    for (run = 1; run <= 1 + INDEPENDENT_RESTARTS; run++) {

      long evaluations_done = (long) (coco_problem_get_evaluations(PROBLEM) + 
            coco_problem_get_evaluations_constraints(PROBLEM));
      long evaluations_remaining = (long) (dimension * BUDGET_MULTIPLIER) - evaluations_done;       

      if ((coco_problem_final_target_hit(PROBLEM) && 
           coco_problem_get_number_of_constraints(PROBLEM) == 0)
           || (evaluations_remaining <= 0))
        break;
     
      Problem const problem(evaluate_function, PROBLEM);        
      de.run(problem, evaluations_remaining, popSize, F, Cr); 

      if ((long int)coco_problem_get_evaluations(PROBLEM) == evaluations_done) {
        printf("WARNING: Budget has not been exhausted (%lu/%lu evaluations done)!\n",
            (unsigned long) evaluations_done, (unsigned long) dimension * BUDGET_MULTIPLIER);
        break;
      }
      else if ((long int)coco_problem_get_evaluations(PROBLEM) < evaluations_done)
        coco_error("Something unexpected happened - function evaluations were decreased!");
      }

      timing_data_time_problem(timing_data, PROBLEM);
  }

  timing_data_finalize(timing_data);

  coco_observer_free(observer);
  coco_suite_free(suite);  
}

void experimentPSO(ParticleSwarm ps) {

  
  std::map<int, double> updateSettings;
  size_t run;
  coco_suite_t *suite;
  coco_observer_t *observer;
  timing_data_t *timing_data;

  char *observer_options =
      coco_strdupf("result_folder: %s " 
                   "algorithm_name: %s "
                   "", ps.getIdString().c_str(), ps.getIdString().c_str());

  suite = coco_suite("bbob", "instances: 1,1,1,1,1", "dimensions: 20");
  observer = coco_observer("bbob", observer_options);
  coco_free_memory(observer_options);

  timing_data = timing_data_initialize(suite);

  while ((PROBLEM = coco_suite_get_next_problem(suite, observer)) != NULL) {
    
    size_t dimension = coco_problem_get_dimension(PROBLEM);
    int const popSize = dimension * 6; 

    for (run = 1; run <= 1 + INDEPENDENT_RESTARTS; run++) {

      long evaluations_done = (long) (coco_problem_get_evaluations(PROBLEM) + 
            coco_problem_get_evaluations_constraints(PROBLEM));
      long evaluations_remaining = (long) (dimension * BUDGET_MULTIPLIER) - evaluations_done;       

      if ((coco_problem_final_target_hit(PROBLEM) && 
           coco_problem_get_number_of_constraints(PROBLEM) == 0)
           || (evaluations_remaining <= 0))
        break;
      
      Problem const problem(evaluate_function, PROBLEM);
      ps.run(problem, evaluations_remaining, popSize, updateSettings); 

      if ((long int)coco_problem_get_evaluations(PROBLEM) == evaluations_done) {
        printf("WARNING: Budget has not been exhausted (%lu/%lu evaluations done)!\n",
            (unsigned long) evaluations_done, (unsigned long) dimension * BUDGET_MULTIPLIER);
        break;
      }
      else if ((long int)coco_problem_get_evaluations(PROBLEM) < evaluations_done)
        coco_error("Something unexpected happened - function evaluations were decreased!");
    }

    timing_data_time_problem(timing_data, PROBLEM);
  }

  timing_data_finalize(timing_data);

  coco_observer_free(observer);
  coco_suite_free(suite);
  
}

/**
 * Allocates memory for the timing_data_t object and initializes it.
 */
static timing_data_t *timing_data_initialize(coco_suite_t *suite) {

	timing_data_t *timing_data = (timing_data_t *) coco_allocate_memory(sizeof(*timing_data));
	size_t function_idx, dimension_idx, instance_idx, i;

	/* Find out the number of all dimensions */
	coco_suite_decode_problem_index(suite, coco_suite_get_number_of_problems(suite) - 1, &function_idx,
			&dimension_idx, &instance_idx);
	timing_data->number_of_dimensions = dimension_idx + 1;
	timing_data->current_idx = 0;
	timing_data->output = (char **) coco_allocate_memory(timing_data->number_of_dimensions * sizeof(char *));
	for (i = 0; i < timing_data->number_of_dimensions; i++) {
		timing_data->output[i] = NULL;
	}
	timing_data->previous_dimension = 0;
	timing_data->cumulative_evaluations = 0;
	time(&timing_data->start_time);
	time(&timing_data->overall_start_time);

	return timing_data;
}

/**
 * Keeps track of the total number of evaluations and elapsed time. Produces an output string when the
 * current problem is of a different dimension than the previous one or when NULL.
 */
static void timing_data_time_problem(timing_data_t *timing_data, coco_problem_t *problem) {

	double elapsed_seconds = 0;

	if ((problem == NULL) || (timing_data->previous_dimension != coco_problem_get_dimension(problem))) {

		/* Output existing timing information */
		if (timing_data->cumulative_evaluations > 0) {
			time_t now;
			time(&now);
			elapsed_seconds = difftime(now, timing_data->start_time) / (double) timing_data->cumulative_evaluations;
			timing_data->output[timing_data->current_idx++] = coco_strdupf("d=%lu done in %.2e seconds/evaluation\n",
					timing_data->previous_dimension, elapsed_seconds);
		}

		if (problem != NULL) {
			/* Re-initialize the timing_data */
			timing_data->previous_dimension = coco_problem_get_dimension(problem);
			timing_data->cumulative_evaluations = coco_problem_get_evaluations(problem);
			time(&timing_data->start_time);
		}

	} else {
		timing_data->cumulative_evaluations += coco_problem_get_evaluations(problem);
	}
}

/**
 * Outputs and finalizes the given timing data.
 */
static void timing_data_finalize(timing_data_t *timing_data){

	/* Record the last problem */
	timing_data_time_problem(timing_data, NULL);

  if (timing_data) {
  	size_t i;
  	double elapsed_seconds;
		time_t now;
		int hours, minutes, seconds;

		time(&now);
		elapsed_seconds = difftime(now, timing_data->overall_start_time);

  	printf("\n");
  	for (i = 0; i < timing_data->number_of_dimensions; i++) {
    	if (timing_data->output[i]) {
				printf("%s", timing_data->output[i]);
				coco_free_memory(timing_data->output[i]);
    	}
    }

  	hours = (int) elapsed_seconds / 3600;
  	minutes = ((int) elapsed_seconds % 3600) / 60;
  	seconds = (int)elapsed_seconds - (hours * 3600) - (minutes * 60);
  	printf("Total elapsed time: %dh%02dm%02ds\n", hours, minutes, seconds);

    coco_free_memory(timing_data->output);
    coco_free_memory(timing_data);
  }
}
