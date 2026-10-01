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

#include <string>
#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <sstream>
#include <cassert>
#include <functional>

#include "asynchronous_algorithms/evolutionary_algorithm.hxx"
#include "asynchronous_algorithms/differential_evolution.hxx"
#include "util/recombination.hxx"
#include "util/statistics.hxx"

//from undvc_common
#include "vector_io.hxx"
#include "arguments.hxx"

using namespace std;

/**
 *  Initialize a differential evolution search from command line parameters
 */

DifferentialEvolution::DifferentialEvolution() {
}

void
DifferentialEvolution::set_print_statistics(void (*_print_statistics)(const std::vector<double> &)) {
    print_statistics = _print_statistics;

}

void
DifferentialEvolution::parse_arguments(const vector<string> &arguments) {
    string parent_selection_name, recombination_selection_name;

    if (!get_argument(arguments, "--parent_scaling_factor", false, parent_scaling_factor)) {
        cerr << "Argument '--parent_scaling_factor <F>' not found, using default of 1.0." << endl;
        parent_scaling_factor = 1.0;
    }

    // with adaptive DE update, differential_scaling_factor and crossover_rate are only used as fallbacks, so they are not required to be provided
    if (!get_argument(arguments, "--differential_scaling_factor", false, differential_scaling_factor)) {
        cerr << "Argument '--differential_scaling_factor <F>' not found, using default of 1.0." << endl;
        differential_scaling_factor = 1.0;
    }

    // with adaptive DE update, differential_scaling_factor and crossover_rate are only used as fallbacks, so they are not required to be provided
    if (!get_argument(arguments, "--crossover_rate", false, crossover_rate)) {
        cerr << "Argument '--crossover_rate <F>' not found, using default of 0.5." << endl;
        crossover_rate = 0.5;
    }

    if (!get_argument(arguments, "--number_pairs", false, number_pairs)) {
        cerr << "Argument '--number_pairs <I>' not found, using default of 1." << endl;
        number_pairs = 1;
    }

    if (get_argument(arguments, "--parent_selection", false, parent_selection_name)) {
        if (parent_selection_name.compare("best") == 0) {
            parent_selection = PARENT_BEST;
        } else if (parent_selection_name.compare("random") == 0) {
            parent_selection = PARENT_RANDOM;
        } else if (parent_selection_name.compare("current-to-best") == 0) {
            parent_selection = PARENT_CURRENT_TO_BEST;
        } else if (parent_selection_name.compare("current-to-random") == 0) {
            parent_selection = PARENT_CURRENT_TO_RANDOM;
        } else {
            cerr << "Improperly specified parent selection type: '" << parent_selection_name.c_str() << "'" << endl;
            cerr << "Possibilities are:" << endl;
            cerr << "   best" << endl;
            cerr << "   random" << endl;
            cerr << "   current-to-best" << endl;
            cerr << "   current-to-random" << endl;
            exit(1);
        }   
    } else {
        cerr << "Argument '--parent_selection <S>' not found, using default of 'best'." << endl;
        parent_selection = PARENT_BEST;
    }

    if (get_argument(arguments, "--recombination_selection", false, recombination_selection_name)) {
        if (recombination_selection_name.compare("binary") == 0) {
            recombination_selection = RECOMBINATION_BINARY;
        } else if (recombination_selection_name.compare("exponential") == 0) {
            recombination_selection = RECOMBINATION_EXPONENTIAL;
        } else if (recombination_selection_name.compare("sum") == 0) {
            recombination_selection = RECOMBINATION_SUM;
        } else if (recombination_selection_name.compare("none") == 0) {
            recombination_selection = RECOMBINATION_NONE;
        } else {
            cerr << "Improperly specified recombination type: '" << recombination_selection_name.c_str() << "'" << endl;
            cerr << "Possibilities are:" << endl;
            cerr << "   binary" << endl;
            cerr << "   exponential" << endl;
            cerr << "   sum" << endl;
            cerr << "   none" << endl;
            exit(1);
        }   
    } else {
        cerr << "Argument '--recombination_selection <S>' not found, using default of 'binary'." << endl;
        recombination_selection = RECOMBINATION_BINARY;
    }

    directional = argument_exists(arguments, "directional");

    // L-SHADE-specific parameters (optional to provide)
    if (!get_argument(arguments, "--H", false, H)) {
        H = 6; // default memory size
    }
    if (!get_argument(arguments, "--NP_min", false, NP_min)) {
        NP_min = 4; // default minimum population
    }
    string interval;
    if (get_argument(arguments, "--population_reduction_interval", false, interval)) {
        if (interval.empty() || interval.find_first_not_of("0123456789") != string::npos)
            throw string("population_reduction_interval must be a nonnegative integer");
        istringstream input(interval);
        if (!(input >> reduction_interval))
            throw string("population_reduction_interval is too large");
    }
}

//
// Helper sampling functions
//
static inline double sample_uniform01(std::function<double()> &rng) {
    double u = rng();
    // clamp away from 0/1 to avoid infinities in cauchy transform
    if (u <= 0.0) u = std::numeric_limits<double>::min();
    if (u >= 1.0) u = 1.0 - std::numeric_limits<double>::epsilon();
    return u;
}

// Gaussian via Box-Muller 
static inline double sample_normal(std::function<double()> &rng, double mu, double sigma) {
    double u1 = sample_uniform01(rng);
    double u2 = sample_uniform01(rng);
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mu + sigma * z0;
}

// Cauchy using inverse CDF
static inline double sample_cauchy(std::function<double()> &rng, double location, double scale) {
    double u = sample_uniform01(rng);
    return location + scale * tan(M_PI * (u - 0.5));
}

void
DifferentialEvolution::initialize() {
    if (population_size < 4 || NP_min < 4 || NP_min > population_size)
        throw string("Differential evolution requires 4 <= NP_min <= population_size");
    this->current_individual = 0;
    this->initialized_individuals = 0;

    this->global_best_fitness = -numeric_limits<double>::max();
    this->global_best_id = 0;

    population = vector< vector<double> >(population_size, vector<double>(number_parameters, 0.0));
    fitnesses = vector<double>(population_size, -numeric_limits<double>::max());

    print_statistics = NULL;

    // --- L-SHADE additions---
    // logging 
    if (ls_log.is_open()) ls_log.close(); // prevents crashing in subsequent runs
    ls_log.open("lshade_log.csv");
    ls_log << "gen,NP,best,mean,median,worst,MF,MCR,successes\n";

    // Memory for MF and MCR (size H)
    if (H == 0) H = 6;
    MF.assign(H, 0.5);
    MCR.assign(H, 0.5);
    memory_index = 0;

    // last used Fi / CRi
    last_Fi.assign(population_size, 0.5);
    last_CRi.assign(population_size, 0.5);

    // Pool of successes collected during a generation (population_size evaluations)
    success_pool.clear();

    // population reduction parameters
    NP_init = population_size;
    active_ids.resize(NP_init);
    iota(active_ids.begin(), active_ids.end(), 0);
    completed_evaluations = 0;

    // Ensure seeds vector exists if used elsewhere
    seeds.assign(population_size, 0);

}

DifferentialEvolution::DifferentialEvolution(const vector<string> &arguments) : EvolutionaryAlgorithm(arguments) {
    parse_arguments(arguments);
    initialize();
}

DifferentialEvolution::DifferentialEvolution(const vector<double> &min_bound, const vector<double> &max_bound, const vector<string> &arguments) : EvolutionaryAlgorithm(min_bound, max_bound, arguments) {
    parse_arguments(arguments);
    initialize();
}

DifferentialEvolution::DifferentialEvolution( const std::vector<double> &min_bound,         /* min bound is copied into the search */
                                              const std::vector<double> &max_bound,         /* max bound is copied into the search */
                                              const uint32_t population_size,
                                              const uint16_t parent_selection,              /* How to select the parent */
                                              const uint16_t number_pairs,                  /* How many individuals to used to calculate differentials */
                                              const uint16_t recombination_selection,       /* How to perform recombination */
                                              const double parent_scaling_factor,           /* weight for the parent calculation*/
                                              const double differential_scaling_factor,     /* weight for the differential calculation */
                                              const double crossover_rate,                  /* crossover rate for recombination */
                                              const bool directional,                       /* used for directional calculation of differential (this options is not really a recombination) */
                                              const uint32_t maximum_iterations             /* default value is 0 which means no termination */
                                            ) : EvolutionaryAlgorithm(min_bound, max_bound, population_size, maximum_iterations) {

    if (parent_selection != PARENT_BEST && parent_selection != PARENT_RANDOM && parent_selection != PARENT_CURRENT_TO_BEST && parent_selection != PARENT_CURRENT_TO_RANDOM) {
        std::stringstream oss;
        oss << "ERROR [file: " << __FILE__ << ", line: " << __LINE__ << "]: unknown parent selection type (" << parent_selection << ")";
        throw oss.str();
    }

    if (recombination_selection != RECOMBINATION_BINARY && recombination_selection != RECOMBINATION_EXPONENTIAL && recombination_selection != RECOMBINATION_SUM && recombination_selection != RECOMBINATION_NONE) {
        std::stringstream oss;
        oss << "ERROR [file: " << __FILE__ << ", line: " << __LINE__ << "]: unknown recombination selection type (" << parent_selection << ")";
        throw oss.str();
    }

    this->parent_selection = parent_selection;
    this->number_pairs = number_pairs;
    this->recombination_selection = recombination_selection;
    this->parent_scaling_factor = parent_scaling_factor;
    this->differential_scaling_factor = differential_scaling_factor;
    this->crossover_rate = crossover_rate;
    this->directional = directional;

    maximum_created = 0;
    maximum_reported = 0;
    this->maximum_iterations = maximum_iterations;

    initialize();
}

DifferentialEvolution::DifferentialEvolution( const std::vector<double> &min_bound,         /* min bound is copied into the search */
                                              const std::vector<double> &max_bound,         /* max bound is copied into the search */
                                              const uint32_t population_size,
                                              const uint16_t parent_selection,              /* How to select the parent */
                                              const uint16_t number_pairs,                  /* How many individuals to used to calculate differentials */
                                              const uint16_t recombination_selection,       /* How to perform recombination */
                                              const double parent_scaling_factor,           /* weight for the parent calculation*/
                                              const double differential_scaling_factor,     /* weight for the differential calculation */
                                              const double crossover_rate,                  /* crossover rate for recombination */
                                              const bool directional,                       /* used for directional calculation of differential (this options is not really a recombination) */
                                              const uint32_t maximum_created,               /* default value is 0 which means no termination */
                                              const uint32_t maximum_reported               /* default value is 0 which means no termination */
                                            ) : EvolutionaryAlgorithm(min_bound, max_bound, population_size, 0) {

    if (parent_selection != PARENT_BEST && parent_selection != PARENT_RANDOM && parent_selection != PARENT_CURRENT_TO_BEST && parent_selection != PARENT_CURRENT_TO_RANDOM) {
        std::stringstream oss;
        oss << "ERROR [file: " << __FILE__ << ", line: " << __LINE__ << "]: unknown parent selection type (" << parent_selection << ")";
        throw oss.str();
    }

    if (recombination_selection != RECOMBINATION_BINARY && recombination_selection != RECOMBINATION_EXPONENTIAL && recombination_selection != RECOMBINATION_SUM && recombination_selection != RECOMBINATION_NONE) {
        std::stringstream oss;
        oss << "ERROR [file: " << __FILE__ << ", line: " << __LINE__ << "]: unknown recombination selection type (" << parent_selection << ")";
        throw oss.str();
    }

    this->parent_selection = parent_selection;
    this->number_pairs = number_pairs;
    this->recombination_selection = recombination_selection;
    this->parent_scaling_factor = parent_scaling_factor;
    this->differential_scaling_factor = differential_scaling_factor;
    this->crossover_rate = crossover_rate;
    this->directional = directional;

    this->maximum_created = maximum_created;
    this->maximum_reported = maximum_reported;
    maximum_iterations = 0;

    initialize();
}


DifferentialEvolution::~DifferentialEvolution() {
}

void
DifferentialEvolution::new_individual(uint32_t &id, std::vector<double> &parameters, uint32_t &seed) {
    DifferentialEvolution::new_individual(id, parameters);

    seeds[id] = ((*random_number_generator)() * numeric_limits<uint32_t>::max()) / 10.0;    //uint max is too large for some reason
    seed = seeds[id];
}

void
DifferentialEvolution::new_individual(uint32_t &id, std::vector<double> &parameters) {
    // A generator loaded before retirement may have an older cursor.
    current_individual %= population_size;
    id = active_ids[current_individual];
    current_individual++;
    if (current_individual >= population_size) {
        current_individual = 0;
        current_iteration++;
    }

    if (initialized_individuals < population_size) { //The search has not been fully initalized so keep generating random individuals
        Recombination::random_within(min_bound, max_bound, parameters, random_number_generator);
        population[id].assign(parameters.begin(), parameters.end());
        individuals_created++;
        return;
    }

    /**
     *  Select the parent.
     */
    vector<double> parent(number_parameters, 0);

    switch (parent_selection) {
        case PARENT_BEST:
            parent.assign(population[global_best_id].begin(), population[global_best_id].end());
            break;

        case PARENT_RANDOM:
            {   //Need a block here to avoid comiler error reusing random_individual variable
                uint32_t random_individual = random_active_id();
                parent.assign(population[random_individual].begin(), population[random_individual].end());
            }
            break;

        case PARENT_CURRENT_TO_BEST:
            for (uint32_t i = 0; i < number_parameters; i++) {
                parent[i] = parent_scaling_factor * (population[global_best_id][i] - population[id][i]);
            }
            break;

        case PARENT_CURRENT_TO_RANDOM:
            {   //Need a block here to avoid comiler error reusing random_individual variable
                uint32_t random_individual = random_active_id();
                for (uint32_t i = 0; i < number_parameters; i++) {
                    parent[i] = parent_scaling_factor * (population[random_individual][i] - population[id][i]);
                }
            }
            break;

        default:
            std::stringstream oss;
            oss << "ERROR [file: " << __FILE__ << ", line: " << __LINE__ << "]: unknown parent selection type (" << parent_selection << ")";
            throw oss.str();
            break;
    }

    /**
     *  L-SHADE adaptive parameter sampling:
     *  sample Fi and CRi using current memory MF/MCR (choose random memory index r)
     */
    uint32_t r = (uint32_t)( (*random_number_generator)() * H );
    if (r >= H) r = H - 1;

    // Wrap the RNG into a callable for the helper functions
    std::function<double()> rng = [this]() { return (*random_number_generator)(); };

    double Fi = 0.5;
    double CRi = 0.5;

    // Sample Fi from Cauchy(MF[r], 0.1) until Fi > 0 
    for (int attempts = 0; attempts < 10; ++attempts) {
        Fi = sample_cauchy(rng, MF[r], 0.1);
        if (Fi > 0.0) break;
    }
    if (Fi > 1.0) Fi = 1.0;
    if (Fi <= 0.0) Fi = 0.5; // fallback

    // Sample CRi from Normal(MCR[r], 0.1)
    CRi = sample_normal(rng, MCR[r], 0.1);
    if (CRi < 0.0) CRi = 0.0;
    if (CRi > 1.0) CRi = 1.0;

    // Store last_Fi/CRi for this individual (so insert_individual can use them if it succeeds)
    if (id >= last_Fi.size()) {
        // resize if population_size changed unexpectedly
        last_Fi.resize(id + 1, 0.5);
        last_CRi.resize(id + 1, 0.5);
    }
    last_Fi[id] = Fi;
    last_CRi[id] = CRi;

    /**
     *  Calculate the differentials.
     */
    vector<double> differential(number_parameters, 0);

    uint32_t random_individual1;
    uint32_t random_individual2; 
    for (uint32_t i = 0; i < number_pairs * 2; i++) {
        random_individual1 = random_active_id();
        random_individual2 = random_active_id();

        if (directional) { //Used for directional recombination (although that part is not the recombination step)
            if (fitnesses[random_individual2] > fitnesses[random_individual1]) {
                uint32_t temp = random_individual1;
                random_individual1 = random_individual2;
                random_individual2 = temp;
            }
        }

        for (uint32_t j = 0; j < number_parameters; j++) {
            differential[j] = population[random_individual1][j] - population[random_individual2][j];
            if (wrap_radians && (fabs(min_bound[j] - (0.0)) < 0.001) && (fabs(max_bound[j] - (M_PI)) < 0.01) ) {  
	        double tempDifferential = M_PI - differential[j];
		if(differential[j] > tempDifferential) {
                    differential[j] = tempDifferential;
                }
            }
        }
    }

    // Use Fi here instead of global differential_scaling_factor (L-SHADE uses Fi per trial)
    for (uint32_t i = 0; i < number_parameters; i++) differential[i] *= Fi / number_pairs;

    /**
     *  Perform the recombination.
     */

    for (uint32_t i = 0; i < number_parameters; i++) parent[i] += differential[i];

    switch (recombination_selection) {
        case RECOMBINATION_BINARY:
            Recombination::binary_recombination(population[id], parent, CRi, parameters, random_number_generator);
            Recombination::bound_parameters(min_bound, max_bound, parameters, wrap_radians);
            break;

        case RECOMBINATION_EXPONENTIAL:
            Recombination::exponential_recombination(population[id], parent, CRi, parameters, random_number_generator);
            Recombination::bound_parameters(min_bound, max_bound, parameters, wrap_radians);
            break;

        case RECOMBINATION_SUM:
            for (uint32_t i = 0; i < number_parameters; i++) parameters[i] = population[id][i] + parent[i];
            Recombination::bound_parameters(min_bound, max_bound, parameters, wrap_radians);
            break;

        case RECOMBINATION_NONE:
            Recombination::bound_parameters(min_bound, max_bound, parent, wrap_radians);
            parameters.assign(parent.begin(), parent.end());
            break;

        default:
            std::stringstream oss;
            oss << "ERROR [file: " << __FILE__ << ", line: " << __LINE__ << "]: unknown recombination type (" << recombination_selection << ")";
            throw oss.str();
            break;
    }

    individuals_created++;
}


bool
DifferentialEvolution::insert_individual(uint32_t id, const std::vector<double> &parameters, double fitness, uint32_t seed) {
    // Call once per canonical evaluation, not once per BOINC replica. Retired
    // slots still count as completed work but must never re-enter the population.
    if (id >= population.size() || parameters.size() != number_parameters ||
        !std::isfinite(fitness)) return false;
    ++completed_evaluations;
    ++individuals_reported;
    bool modified = false;
    if (would_insert(id, fitness)) {
        if (fitnesses[id] == -numeric_limits<double>::max()) initialized_individuals++;

        // Record success (Fi, CRi, delta f) for memory update (L-SHADE)
        {
            SuccessRecord s;
            // safety: if id out of range, use defaults
            if (id < last_Fi.size()) s.Fi = last_Fi[id]; else s.Fi = differential_scaling_factor;
            if (id < last_CRi.size()) s.CRi = last_CRi[id]; else s.CRi = crossover_rate;

            s.df = fitness - fitnesses[id];
            // Only record if positive improvement
            if (fitnesses[id] != -numeric_limits<double>::max() &&
                s.df > 0.0 && std::isfinite(s.df)) {
                ls_log << "DEBUG: Recording success Fi=" << s.Fi << " CRi=" << s.CRi << " df=" << s.df << endl; // check that successes are being recorded
                success_pool.push_back(s);
            }
        }

        fitnesses[id] = fitness;
        population[id].assign(parameters.begin(), parameters.end());
        seeds[id] = seed;

        cout.precision(10);

        if (global_best_fitness < fitness) {
            global_best_id = id;
            global_best_fitness = fitness;

            if (log_file == NULL) {
                if (!quiet) {
                    cout.precision(10);
                    cout <<  current_iteration << ":" << id << " - GLOBAL: " << global_best_fitness << " " << vector_to_string(parameters) << endl;
                }
            } else {
                double best, average, median, worst;
                active_fitness_statistics(best, average, median, worst);
                (*log_file) << individuals_reported << " -- b: " << best << ", a: " << average << ", m: " << median << ", w: " << worst << ", " << vector_to_string(parameters) << endl;
            } 
        }

        modified = true;
    }
    // Memory adaptation retains its existing report cadence. Reduction uses
    // its own evaluation schedule, independent of iterations and improvements.
    if (population_size > 0 && (individuals_reported % population_size == 0)) {
        uint32_t successes_count = (uint32_t)success_pool.size();
        update_memory_from_successes();
        log_generation_state(successes_count);
    }
    reduce_population();

    return modified;
}


bool
DifferentialEvolution::would_insert(uint32_t id, double fitness) {
    return is_active(id) && std::isfinite(fitness) && fitnesses[id] < fitness;
}
 
/**
 *  The following method is for synchronous optimization and is purely virtual
 */
void
DifferentialEvolution::iterate(double (*objective_function)(const std::vector<double> &)) {
    cout << "Initialized differential evolution. " << endl;
    cout << "   maximum_iterations:          " << maximum_iterations << endl;
    cout << "   current_iteration:           " << current_iteration << endl;
    cout << "   number_pairs:                " << number_pairs << endl;
    cout << "   parent_selection:            " << parent_selection << endl;
    cout << "   recombination_selection:     " << recombination_selection << endl;
    cout << "   parent_scaling_factor:       " << parent_scaling_factor << endl;
    cout << "   differential_scaling_factor: " << differential_scaling_factor << endl;
    cout << "   crossover_rate:              " << crossover_rate << endl;
    cout << "   directional:                 " << directional << endl;

    uint32_t id;
    vector<double> parameters(number_parameters, 0);

    while (maximum_iterations == 0 || current_iteration < maximum_iterations) {
        for (uint32_t i = 0; i < population_size; i++) {
            new_individual(id, parameters);
            double fitness = objective_function(parameters);
            insert_individual(id, parameters, fitness);
        }

        // new_individual advances current_iteration.
    }
}

void
DifferentialEvolution::iterate(double (*objective_function)(const std::vector<double> &, const uint32_t)) {
    cout << "Initialized differential evolution. " << endl;
    cout << "   maximum_iterations:          " << maximum_iterations << endl;
    cout << "   current_iteration:           " << current_iteration << endl;
    cout << "   number_pairs:                " << number_pairs << endl;
    cout << "   parent_selection:            " << parent_selection << endl;
    cout << "   recombination_selection:     " << recombination_selection << endl;
    cout << "   parent_scaling_factor:       " << parent_scaling_factor << endl;
    cout << "   differential_scaling_factor: " << differential_scaling_factor << endl;
    cout << "   crossover_rate:              " << crossover_rate << endl;
    cout << "   directional:                 " << directional << endl;

    uint32_t id;
    uint32_t seed;
    vector<double> parameters(number_parameters, 0);

    while (maximum_iterations == 0 || current_iteration < maximum_iterations) {
        for (uint32_t i = 0; i < population_size; i++) {
            new_individual(id, parameters, seed);
            double fitness = objective_function(parameters, seed);
            insert_individual(id, parameters, fitness, seed);
        }

        //This is now updated in the 'new_individual' function
//        current_iteration++;
    }
}

void
DifferentialEvolution::get_individuals(std::vector<Individual> &individuals) {
    individuals.clear();
    for (uint32_t id : active_ids) {
        individuals.push_back(Individual(id, fitnesses[id], population[id], ""));
    }
}


// ----------------- L-SHADE helper methods -----------------

// Update MF and MCR from collected successes using weighted means
void DifferentialEvolution::update_memory_from_successes() {
    if (success_pool.empty()) {
        return;
    }

    // compute sum of df
    double sum_df = 0.0;
    for (auto &s : success_pool) sum_df += s.df;
    if (sum_df <= 0.0) {
        success_pool.clear();
        return;
    }

    double weight_sum = 0.0;
    double weighted_F_lehmer_num = 0.0;
    double weighted_F_lehmer_den = 0.0;
    double weighted_CR_sum = 0.0;

    for (auto &s : success_pool) {
        double w = s.df / sum_df;
        weighted_F_lehmer_num += w * s.Fi * s.Fi;
        weighted_F_lehmer_den += w * s.Fi;
        weighted_CR_sum += w * s.CRi;
        weight_sum += w;
    }

    // Lehmer mean for MF (numerator / denominator)
    double new_MF = MF[memory_index];
    if (weighted_F_lehmer_den > 0.0) {
        new_MF = weighted_F_lehmer_num / weighted_F_lehmer_den;
    }

    // Weighted mean for MCR
    double new_MCR = MCR[memory_index];
    if (weight_sum > 0.0) {
        new_MCR = weighted_CR_sum / weight_sum;
    }

    MF[memory_index] = new_MF;
    MCR[memory_index] = new_MCR;

    memory_index = (memory_index + 1) % H;

    success_pool.clear();
}

bool DifferentialEvolution::is_active(uint32_t id) const {
    return binary_search(active_ids.begin(), active_ids.end(), id);
}

uint32_t DifferentialEvolution::random_active_id() {
    const uint32_t index = (*random_number_generator)() * active_ids.size();
    return active_ids[index];
}

void DifferentialEvolution::active_fitness_statistics(double &best, double &mean,
                                                       double &median, double &worst) const {
    vector<double> active_fitnesses;
    active_fitnesses.reserve(active_ids.size());
    for (uint32_t id : active_ids) active_fitnesses.push_back(fitnesses[id]);
    calculate_fitness_statistics(active_fitnesses, best, mean, median, worst);
}

void DifferentialEvolution::reduce_population() {
    if (reduction_interval == 0 || initialized_individuals < population_size) return;
    const uint64_t removable = NP_init - NP_min;
    // Cap before multiplying, including for very long-running searches.
    const uint64_t steps = min(completed_evaluations / reduction_interval,
                               (removable + 4) / 5);
    const uint32_t target = NP_init - min(removable, 5 * steps);
    shrink_population_to(target);
}

// Retire the worst members without moving/reusing slots referenced by workunits.
void DifferentialEvolution::shrink_population_to(uint32_t new_size) {
    new_size = max(new_size, NP_min);
    if (new_size >= population_size) return;

    vector<uint32_t> survivors(active_ids);
    sort(survivors.begin(), survivors.end(), [&](uint32_t a, uint32_t b) {
        if (fitnesses[a] != fitnesses[b]) return fitnesses[a] > fitnesses[b];
        return a < b; // deterministic ties
    });
    survivors.resize(new_size);
    global_best_id = survivors.front();
    global_best_fitness = fitnesses[global_best_id];
    sort(survivors.begin(), survivors.end());
    active_ids.swap(survivors);
    population_size = active_ids.size();
    initialized_individuals = 0;
    for (uint32_t id : active_ids)
        if (fitnesses[id] != -numeric_limits<double>::max()) ++initialized_individuals;
    current_individual %= population_size;
}

void DifferentialEvolution::log_generation_state(uint32_t successes_count) {
    double best, mean, median, worst;
    active_fitness_statistics(best, mean, median, worst);

    // memory_index points to the next slot; the most recently updated slot is (memory_index - 1)
    uint32_t last_index = 0;
    if (H > 0) last_index = (memory_index + H - 1) % H;

    ls_log << current_iteration << ","
           << population_size << ","
           << best << ","
           << mean << ","
           << median << ","
           << worst << ","
           << MF[last_index] << ","
           << MCR[last_index] << ","
           << successes_count
           << "\n";
}

// ----------------- end L-SHADE helper methods -----------------
