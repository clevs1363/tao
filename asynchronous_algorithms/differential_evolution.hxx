/*
 * Copyright 2012, 2009 Travis Desell and the University of North Dakota.
 *
 * This file is part of the Toolkit for Asynchronous Optimization (TAO).
 *
 * TAO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TAO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TAO.  If not, see <http://www.gnu.org/licenses/>.
 * */

#ifndef TAO_DIFFERENTIAL_EVOLUTION_H
#define TAO_DIFFERENTIAL_EVOLUTION_H

#include <string>
#include <vector>

#include "asynchronous_algorithms/evolutionary_algorithm.hxx"


class DifferentialEvolution : public EvolutionaryAlgorithm {
    protected:
        uint32_t number_pairs;
        uint16_t parent_selection;
        uint16_t recombination_selection;
        bool directional;

        double parent_scaling_factor;
        double differential_scaling_factor;
        double crossover_rate;

        std::vector<double> fitnesses;
        std::vector< std::vector<double> > population;

        uint32_t initialized_individuals;

        double global_best_fitness;
        uint32_t global_best_id;

        // L-SHADE parameters
        std::vector<double> MF;   // size H
        std::vector<double> MCR;  // size H
        uint32_t memory_index;    // cycles 0..H-1
        uint32_t H = 6;           // memory size

        std::vector<double> last_Fi;
        std::vector<double> last_CRi;

        uint32_t NP_init;   // initial population size
        uint32_t NP_min = 4; // minimum active population size

        // Slots never move: BOINC results carry these stable IDs. population_size
        // is the active count; all per-slot vectors retain NP_init entries.
        std::vector<uint32_t> active_ids; // sorted by stable ID
        uint64_t reduction_interval = 10000; // reports per removal of five; 0 disables
        uint64_t completed_evaluations = 0;

        bool is_active(uint32_t id) const;
        uint32_t random_active_id();
        void reduce_population();
        void active_fitness_statistics(double &best, double &mean,
                                       double &median, double &worst) const;

        struct SuccessRecord {
          double Fi;
          double CRi;
          double df;   // fitness improvement
        };
        std::vector<SuccessRecord> success_pool;
        std::ofstream ls_log;  // log file for L-SHADE

        // end L-SHADE parameters

        DifferentialEvolution();

        void initialize();
        void parse_arguments(const std::vector<std::string> &arguments);
        void update_memory_from_successes();
        void shrink_population_to(uint32_t new_size);
        void log_generation_state(uint32_t successes_count);

    public:
        void (*print_statistics)(const std::vector<double> &);

        double get_global_best_fitness() { return global_best_fitness; }
        std::vector<double> get_global_best() { return population[global_best_id]; }

        //The following are different types parent selection
        const static uint16_t PARENT_BEST = 0;
        const static uint16_t PARENT_RANDOM = 1;
        const static uint16_t PARENT_CURRENT_TO_BEST = 2;
        const static uint16_t PARENT_CURRENT_TO_RANDOM = 3;

        //The following are different types of recombination
        const static uint16_t RECOMBINATION_BINARY = 0;
        const static uint16_t RECOMBINATION_EXPONENTIAL = 1;
        const static uint16_t RECOMBINATION_SUM = 2;
        const static uint16_t RECOMBINATION_NONE = 3;

        DifferentialEvolution( const std::vector<std::string> &arguments);

        DifferentialEvolution( const std::vector<double> &min_bound,                                    /* min bound is copied into the search */
                               const std::vector<double> &max_bound,                                    /* max bound is copied into the search */
                               const std::vector<std::string> &arguments);          /* initialize the DE from command line arguments */

        DifferentialEvolution( const std::vector<double> &min_bound,                                    /* min bound is copied into the search */
                               const std::vector<double> &max_bound,                                    /* max bound is copied into the search */
                               const uint32_t population_size,
                               const uint16_t parent_selection,                                         /* How to select the parent */
                               const uint16_t number_pairs,                                             /* How many individuals to used to calculate differntials */
                               const uint16_t recombination_selection,                                  /* How to perform recombination */
                               const double parent_scaling_factor,                                      /* weight for the parent calculation*/
                               const double differential_scaling_factor,                                /* weight for the differential calculation */
                               const double crossover_rate,                                             /* crossover rate for recombination */
                               const bool directional,                                                  /* used for directional calculation of differential (this options is not really a recombination) */
                               const uint32_t maximum_iterations                                        /* default value is 0 which means no termination */
                             );

        DifferentialEvolution( const std::vector<double> &min_bound,                                    /* min bound is copied into the search */
                               const std::vector<double> &max_bound,                                    /* max bound is copied into the search */
                               const uint32_t population_size,
                               const uint16_t parent_selection,                                         /* How to select the parent */
                               const uint16_t number_pairs,                                             /* How many individuals to used to calculate differntials */
                               const uint16_t recombination_selection,                                  /* How to perform recombination */
                               const double parent_scaling_factor,                                      /* weight for the parent calculation*/
                               const double differential_scaling_factor,                                /* weight for the differential calculation */
                               const double crossover_rate,                                             /* crossover rate for recombination */
                               const bool directional,                                                  /* used for directional calculation of differential (this options is not really a recombination) */
                               const uint32_t maximum_created,                                          /* default value is 0 which means no termination */
                               const uint32_t maximum_reported                                          /* default value is 0 which means no termination */
                             );


        virtual ~DifferentialEvolution();

        void parse_command_line(const std::vector<std::string> &arguments);
        /**
         *  The following methods are used for asynchronous optimization and are purely virtual
         */
        virtual void new_individual(uint32_t &id, std::vector<double> &parameters);
        virtual void new_individual(uint32_t &id, std::vector<double> &parameters, uint32_t &seed);
        virtual bool insert_individual(uint32_t id, const std::vector<double> &parameters, double fitness, uint32_t seed = 0);     /* Returns true if the individual is inserted. */
        virtual bool would_insert(uint32_t id, double fitness);

        /**
         *  The following method is for synchronous optimization and is purely virtual
         */
        void iterate(double (*objective_function)(const std::vector<double> &));
        void iterate(double (*objective_function)(const std::vector<double> &, const uint32_t));    //this objective function also requires a seed

        void set_print_statistics(void (*_print_statistics)(const std::vector<double> &));

        virtual void get_individuals(std::vector<Individual> &individuals);
};


#endif
