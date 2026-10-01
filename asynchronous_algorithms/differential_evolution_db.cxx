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

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "evolutionary_algorithm_db.hxx"
#include "differential_evolution_db.hxx"
#include "vector_io.hxx"
#include "arguments.hxx"

#include "util/statistics.hxx"

/**
 *  From MYSQL
 */
#include "mysql.h"

using namespace std;

void
DifferentialEvolutionDB::check_name(string name) throw (string) {
    if (name.substr(0,3).compare("de_") != 0) {
        ostringstream err_msg;
        err_msg << "Improper name for DifferentialEvolutionDB '" << name << "' must start with 'de_'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw err_msg.str();
    }
}

/**
 *  The following construct a DifferentialEvolution from a database entry
 */
DifferentialEvolutionDB::DifferentialEvolutionDB(MYSQL *conn, string name) throw (string) {
    check_name(name);
    this->conn = conn;

    ostringstream oss;
    oss << "SELECT * FROM differential_evolution WHERE name = '" << name << "'";
//    cout << oss.str() << endl;

    construct_from_database(oss.str());
}

DifferentialEvolutionDB::DifferentialEvolutionDB(MYSQL *conn, int id) throw (string) {
    this->conn = conn;

    ostringstream oss;
    oss << "SELECT * FROM differential_evolution WHERE id = " << id;
//    cout << oss.str() << endl;

    construct_from_database(oss.str());
}

bool
DifferentialEvolutionDB::search_exists(MYSQL *conn, string search_name) throw (string) {
    ostringstream query;
    query << "SELECT id FROM differential_evolution where name = '" << search_name << "'";

    if (mysql_query(conn, query.str().c_str())) {
        ostringstream ex_msg;
        ex_msg << "ERROR: could query database: '" << query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }

    MYSQL_RES *result = mysql_store_result(conn);
    if (result != NULL) {
        if (mysql_num_rows(result) > 0) return true;
    } else {
        ostringstream ex_msg;
        ex_msg << "ERROR: could not get result from database from query: '" << query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }

    return false;
}

void
DifferentialEvolutionDB::create_tables(MYSQL *conn) throw (string) {
    ostringstream de_query;
    de_query << "CREATE TABLE `differential_evolution` ("
                << "    `id` int(11) NOT NULL AUTO_INCREMENT,"
                << "    `name` varchar(254) NOT NULL DEFAULT '',"
                << "    `parent_selection` int(11) NOT NULL DEFAULT '0', "
                << "    `number_pairs` int(11) NOT NULL DEFAULT '0', "
                << "    `recombination_selection` int(11) NOT NULL DEFAULT '0', "
                << "    `parent_scaling_factor` double NOT NULL default  '1', "
                << "    `differential_scaling_factor` double NOT NULL default  '1', "
                << "    `crossover_rate` double NOT NULL default  '0.5', "
                << "    `directional` tinyint(1) NOT NULL default  '0', "
                << "    `current_individual` int(11) NOT NULL DEFAULT '0',"
                << "    `initialized_individuals` int(11) NOT NULL DEFAULT '0',"
                << "    `current_iteration` int(11) NOT NULL DEFAULT '0',"
                << "    `maximum_iterations` int(11) NOT NULL DEFAULT '0',"
                << "    `individuals_created` int(11) NOT NULL DEFAULT '0',"
                << "    `maximum_created` int(11) NOT NULL DEFAULT '0',"
                << "    `individuals_reported` int(11) NOT NULL DEFAULT '0',"
                << "    `maximum_reported` int(11) NOT NULL DEFAULT '0',"
                << "    `population_size` int(11) NOT NULL DEFAULT '0',"
                << "    `min_bound` varchar(2048) NOT NULL,"
                << "    `max_bound` varchar(2048) NOT NULL,"
                << "    `app_id`    int(11) NOT NULL DEFAULT '-1',"
                << "    `wrap_radians` tinyint(1) NOT NULL default '0',"
                << "    `H` int(11) NOT NULL DEFAULT '6',"
                << "    `NP_min` int(11) NOT NULL DEFAULT '4',"
                << "    `NP_init` int(11) NOT NULL DEFAULT '0',"
                << "    `memory_index` int(11) NOT NULL DEFAULT '0',"
                << "    `MF` TEXT NOT NULL,"
                << "    `MCR` TEXT NOT NULL,"
                << "    `last_Fi` TEXT NOT NULL,"
                << "    `last_CRi` TEXT NOT NULL,"
                << "    `reduction_interval` BIGINT UNSIGNED NOT NULL DEFAULT '10000',"
                << "    `active_ids` TEXT NULL,"
                << "    `completed_evaluations` BIGINT UNSIGNED NOT NULL DEFAULT '0',"
                << "PRIMARY KEY (`id`),"
                << "UNIQUE KEY `name` (`name`)"
                << ") ENGINE=InnoDB AUTO_INCREMENT=0 DEFAULT CHARSET=latin1";

    cout << "creating differential_evolution table with: " << endl << de_query.str() << endl << endl;

    if (mysql_query(conn, de_query.str().c_str())) {
        ostringstream ex_msg;
        ex_msg << "ERROR: could not create differential evolution table with query: '" << de_query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }

    ostringstream individual_query;
    individual_query  << "CREATE TABLE `de_individual` ("
                      << "    `differential_evolution_id` int(11) NOT NULL,"
                      << "    `position` int(11) NOT NULL,"
                      << "    `fitness` double NOT NULL,"
                      << "    `parameters` varchar(2048) NOT NULL,"
                      << "    `seed` int(32) UNSIGNED,"
                      << "PRIMARY KEY (`differential_evolution_id`,`position`)"
                      << ") ENGINE=InnoDB DEFAULT CHARSET=latin1";

    cout << "creating de_individual table with: " << endl << individual_query.str() << endl << endl;

    /**
     *  Create the log
     */
    //create table `particle_swarm_log` (`swarm_id` int(11) not null, `evaluation` int(11) not null, `fitness` double, `particle` int(11) not null, `r1` double, `r2` double, `global` bool not null, PRIMARY KEY(`swarm_id`, `evaluation`)) ENGINE=InnoDB;

    if (mysql_query(conn, individual_query.str().c_str())) {
        ostringstream ex_msg;
        ex_msg << "ERROR: could not create de_individual table with query: '" << individual_query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }
}

void 
DifferentialEvolutionDB::construct_from_database(string query) throw (string) {
    mysql_query(conn, query.c_str());

    if (mysql_errno(conn) != 0) {
        ostringstream ex_msg;
        ex_msg << "ERROR: could not get differential evolution from query: '" << query << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }

    MYSQL_RES *result = mysql_store_result(conn);

    if (result != NULL) {
        MYSQL_ROW row = mysql_fetch_row(result);
        
        if (row == NULL) {
            ostringstream ex_msg;
            ex_msg << "ERROR: could not construct differential evolution '" << name << "' from database, it does not exist. Thrown on " << __FILE__ << ":" << __LINE__;
            throw ex_msg.str();
        }

        if (mysql_num_fields(result) < 33) {
            mysql_free_result(result);
            throw string("Differential evolution schema requires migrations/001_de_population_reduction.sql");
        }
        try {
            construct_from_database(row);
        } catch (...) {
            mysql_free_result(result);
            throw;
        }
        mysql_free_result(result);
    } else {
        ostringstream ex_msg;
        ex_msg << "ERROR: could not get differential evolution from query: '" << query << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }
}


void 
DifferentialEvolutionDB::construct_from_database(MYSQL_ROW row) throw (string) {
    id = atoi(row[0]);
    name = row[1];

    //Inherited from DifferentialEvolution
    parent_selection = atoi(row[2]);
    number_pairs = atoi(row[3]);
    recombination_selection = atoi(row[4]);
    parent_scaling_factor = atof(row[5]);
    differential_scaling_factor = atof(row[6]);
    crossover_rate = atof(row[7]);
    directional = atoi(row[8]);

    current_individual = atoi(row[9]);
    initialized_individuals = atoi(row[10]);

    //inherited from EvolutionaryAlgorithm
    current_iteration = atoi(row[11]);
    maximum_iterations = atoi(row[12]);
    individuals_created = atoi(row[13]);
    maximum_created = atoi(row[14]);
    individuals_reported = atoi(row[15]);
    maximum_reported = atoi(row[16]);

    population_size = atoi(row[17]);
    string_to_vector<double>(row[18], min_bound);
    string_to_vector<double>(row[19], max_bound);
    app_id = atoi(row[20]);
    wrap_radians = atoi(row[21]);
    H = atoi(row[22]);
    NP_min = atoi(row[23]);
    NP_init = atoi(row[24]);
    memory_index = atoi(row[25]);
    string_to_vector<double>(row[26], MF);
    string_to_vector<double>(row[27], MCR);
    string_to_vector<double>(row[28], last_Fi);
    string_to_vector<double>(row[29], last_CRi);
    reduction_interval = strtoull(row[30], NULL, 10);
    completed_evaluations = strtoull(row[32], NULL, 10);
    if (H == 0) H = 6;
    if (NP_init == 0) NP_init = population_size;
    if (population_size < 4 || NP_min < 4 || NP_min > population_size || NP_init < population_size)
        throw string("Invalid differential evolution population sizes in database");
    active_ids.clear();
    if (row[31] == NULL) {
        // Migration of an unreduced search: original slots are all active.
        if (NP_init != population_size)
            throw string("Cannot infer active IDs for an already compacted population");
        active_ids.resize(NP_init);
        iota(active_ids.begin(), active_ids.end(), 0);
    } else {
        vector<uint64_t> stored_ids;
        string_to_vector<uint64_t>(row[31], stored_ids);
        for (uint64_t slot : stored_ids) {
            if (slot >= NP_init) throw string("Invalid active population ID in database");
            active_ids.push_back(static_cast<uint32_t>(slot));
        }
    }
    if (active_ids.size() != population_size ||
        !is_sorted(active_ids.begin(), active_ids.end()) ||
        adjacent_find(active_ids.begin(), active_ids.end()) != active_ids.end() ||
        active_ids.back() >= NP_init)
        throw string("Invalid active population IDs in database");
    current_individual %= population_size;
    if (memory_index >= H) memory_index = 0;
    if (MF.size() != H) MF.assign(H, 0.5);
    if (MCR.size() != H) MCR.assign(H, 0.5);
    if (last_Fi.size() != NP_init) last_Fi.assign(NP_init, 0.5);
    if (last_CRi.size() != NP_init) last_CRi.assign(NP_init, 0.5);
    quiet = true;
    print_statistics = NULL;
    number_parameters = min_bound.size();

    //Get the individual information from the database
    ostringstream oss;
    oss << "SELECT position, fitness, parameters, seed FROM de_individual WHERE differential_evolution_id = " << this->id << " ORDER BY position";
    mysql_query(conn, oss.str().c_str());
    MYSQL_RES *result = mysql_store_result(conn);

//    cout << oss.str() << endl;

    fitnesses.resize(NP_init, -numeric_limits<double>::max());
    population.resize(NP_init, vector<double>(number_parameters, 0.0));
    seeds.resize(NP_init, 0);

    if (random_number_generator == NULL) EvolutionaryAlgorithm::initialize_rng();    //to initialize the random number generator

    if (result != NULL) {
        uint32_t num_results = mysql_num_rows(result);
        if (num_results != NP_init) {
            ostringstream ex_msg;
            ex_msg << "ERROR: got " << num_results << " results when looking up individuals for search " << name << ", with a slot count: " << NP_init << ". Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
            throw ex_msg.str();
        }   

        MYSQL_ROW individual_row;

        while ((individual_row = mysql_fetch_row(result))) {
            int individual_id = atoi(individual_row[0]);
            if (individual_id < 0 || static_cast<uint32_t>(individual_id) >= NP_init) {
                mysql_free_result(result);
                throw string("Invalid differential evolution slot ID in database");
            }
            fitnesses[individual_id] = atof(individual_row[1]);

            if (fitnesses[individual_id] < -1.79768e+308) {
                fitnesses[individual_id] = -numeric_limits<double>::max();
            }

            string_to_vector<double>(individual_row[2], population[individual_id]);
            seeds[individual_id] = atoi(individual_row[3]);

//            cout   << "    [DEIndividual" << endl
//                   << "        position = " << individual_id << endl
//                   << "        fitness = " << fitnesses[individual_id] << endl
//                   << "        parameters = '" << vector_to_string<double>(population[individual_id]) << "'" << endl
//                   << "    ]" << endl;

         }   
        mysql_free_result(result);
    } else {
        ostringstream ex_msg;
        ex_msg << "ERROR: looking up individuals with query: '" << oss.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }

    //calculate global_best and global_best_fitness
    global_best_id = active_ids.front();
    global_best_fitness = -numeric_limits<double>::max();
    initialized_individuals = 0;
    for (uint32_t i : active_ids) {
        if (fitnesses[i] != -numeric_limits<double>::max()) ++initialized_individuals;
        if (global_best_fitness < fitnesses[i]) {
            global_best_id = i;
            global_best_fitness = fitnesses[i];
        }
    }
}


/**
 *  Insert a created differential evolution to a database.
 */
void
DifferentialEvolutionDB::insert_to_database() throw (string) {
    ostringstream query;

    query.precision(10);
    query << "INSERT INTO differential_evolution"
          << " SET "
          << "  name = '" << name << "'"
          << ", parent_selection = " << parent_selection
          << ", number_pairs = " << number_pairs
          << ", recombination_selection = " << recombination_selection
          << ", parent_scaling_factor = " << parent_scaling_factor
          << ", differential_scaling_factor = " << differential_scaling_factor
          << ", crossover_rate = " << crossover_rate
          << ", directional = " << directional
          << ", current_individual = " << current_individual
          << ", initialized_individuals = " << initialized_individuals
          << ", current_iteration = " << current_iteration  
          << ", maximum_iterations = " << maximum_iterations
          << ", individuals_created = " << individuals_created
          << ", maximum_created = " << maximum_created  
          << ", individuals_reported = " << individuals_reported
          << ", maximum_reported = " << maximum_reported 
          << ", population_size = " << population_size
          << ", min_bound = '" << vector_to_string<double>(min_bound) << "'"
          << ", max_bound = '" << vector_to_string<double>(max_bound) << "'"
          << ", app_id = " << app_id 
          << ", wrap_radians = " << wrap_radians
          << ", H = " << H
          << ", NP_min = " << NP_min
          << ", NP_init = " << NP_init
          << ", memory_index = " << memory_index
          << ", MF = '" << vector_to_string<double>(MF) << "'"
          << ", MCR = '" << vector_to_string<double>(MCR) << "'"
          << ", last_Fi = '" << vector_to_string<double>(last_Fi) << "'"
          << ", last_CRi = '" << vector_to_string<double>(last_CRi) << "'"
          << ", reduction_interval = " << reduction_interval
          << ", active_ids = '" << vector_to_string<uint32_t>(active_ids) << "'"
          << ", completed_evaluations = " << completed_evaluations;

    mysql_query(conn, query.str().c_str());

    MYSQL_RES *result;
    if ((result = mysql_store_result(conn)) == 0 && mysql_field_count(conn) == 0 && mysql_insert_id(conn) != 0) {
        id = mysql_insert_id(conn);
    } else {
        ostringstream ex_msg;
        ex_msg << "ERROR: could not get differential evolution id from insert query: '" << query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }
    mysql_free_result(result);

    for (uint32_t i = 0; i < population_size; i++) {
        ostringstream individual_query;
        individual_query.precision(10);
        individual_query << "INSERT INTO de_individual"
                         << " SET "
                         << "  differential_evolution_id = " << id
                         << ", position = " << i
                         << ", fitness = '" << setprecision(10) << fixed << fitnesses[i] << "'"
                         << ", parameters = '" << vector_to_string<double>(population[i]) << "'"
                         << ", seed = " << seeds[i];

        mysql_query(conn, individual_query.str().c_str());

        if (mysql_errno(conn) != 0) {
            ostringstream ex_msg;
            ex_msg << "ERROR: creating de_individual with insert query: '" << individual_query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
            throw ex_msg.str();
        }
        mysql_free_result(result);
    }
}

/**
 *  The following constructors create new DifferentialEvolutions and insert them into the database.
 */
//Create a differential evolution entirely from arguments
DifferentialEvolutionDB::DifferentialEvolutionDB(MYSQL *conn, const vector<string> &arguments) throw (string) : DifferentialEvolution(arguments) {
    this->conn = conn;
    this->app_id = -1;
    get_argument(arguments, "--search_name", true, name);
    check_name(name);
    insert_to_database();
}

//Create a differential evolution entirely from arguments
DifferentialEvolutionDB::DifferentialEvolutionDB(MYSQL *conn, const int32_t app_id, const vector<string> &arguments) throw (string) : DifferentialEvolution(arguments) {
    this->conn = conn;
    this->app_id = app_id;
    get_argument(arguments, "--search_name", true, name);
    check_name(name);
    insert_to_database();
}


//Create a differential evolution from arguments and a given min and max bound
DifferentialEvolutionDB::DifferentialEvolutionDB( MYSQL *conn,
                                  const vector<double> &min_bound,            /* min bound is copied into the search */
                                  const vector<double> &max_bound,            /* max bound is copied into the search */
                                  const vector<string> &arguments
                                ) throw (string) : DifferentialEvolution(min_bound, max_bound, arguments) {
    this->conn = conn;
    this->app_id = -1;
    get_argument(arguments, "--search_name", true, name);
    check_name(name);
    insert_to_database();
}


//Create a differential evolution from arguments and a given min and max bound
DifferentialEvolutionDB::DifferentialEvolutionDB( MYSQL *conn,
                                  const int32_t app_id,
                                  const vector<double> &min_bound,            /* min bound is copied into the search */
                                  const vector<double> &max_bound,            /* max bound is copied into the search */
                                  const vector<string> &arguments
                                ) throw (string) : DifferentialEvolution(min_bound, max_bound, arguments) {
    this->conn = conn;
    this->app_id = app_id;
    get_argument(arguments, "--search_name", true, name);
    check_name(name);
    insert_to_database();
}

//Create a differential evolution entirely from defined parameters.
DifferentialEvolutionDB::DifferentialEvolutionDB( MYSQL *conn,
                                                  const string name,
                                                  const vector<double> &min_bound,                                    /* min bound is copied into the search */
                                                  const vector<double> &max_bound,                                    /* max bound is copied into the search */
                                                  const uint32_t population_size,
                                                  const uint16_t parent_selection,                                         /* How to select the parent */
                                                  const uint16_t number_pairs,                                             /* How many individuals to used to calculate differntials */
                                                  const uint16_t recombination_selection,                                  /* How to perform recombination */
                                                  const double parent_scaling_factor,                                      /* weight for the parent calculation*/
                                                  const double differential_scaling_factor,                                /* weight for the differential calculation */
                                                  const double crossover_rate,                                             /* crossover rate for recombination */
                                                  const bool directional,                                                  /* used for directional calculation of differential (this options is not really a recombination) */
                                                  const uint32_t maximum_iterations                                        /* default value is 0 which means no termination */
                                                ) throw (string) : DifferentialEvolution(min_bound, max_bound, population_size, parent_selection, number_pairs, recombination_selection, parent_scaling_factor, differential_scaling_factor, crossover_rate, directional, maximum_iterations) {
    this->conn = conn;
    this->app_id = -1;
    this->name = name;
    check_name(name);
    insert_to_database();
}


//Create a differential evolution entirely from defined parameters.
DifferentialEvolutionDB::DifferentialEvolutionDB( MYSQL *conn,
                                                  const int32_t app_id,
                                                  const string name,
                                                  const vector<double> &min_bound,                                    /* min bound is copied into the search */
                                                  const vector<double> &max_bound,                                    /* max bound is copied into the search */
                                                  const uint32_t population_size,
                                                  const uint16_t parent_selection,                                         /* How to select the parent */
                                                  const uint16_t number_pairs,                                             /* How many individuals to used to calculate differntials */
                                                  const uint16_t recombination_selection,                                  /* How to perform recombination */
                                                  const double parent_scaling_factor,                                      /* weight for the parent calculation*/
                                                  const double differential_scaling_factor,                                /* weight for the differential calculation */
                                                  const double crossover_rate,                                             /* crossover rate for recombination */
                                                  const bool directional,                                                  /* used for directional calculation of differential (this options is not really a recombination) */
                                                  const uint32_t maximum_iterations                                        /* default value is 0 which means no termination */
                                                ) throw (string) : DifferentialEvolution(min_bound, max_bound, population_size, parent_selection, number_pairs, recombination_selection, parent_scaling_factor, differential_scaling_factor, crossover_rate, directional, maximum_iterations) {
    this->conn = conn;
    this->app_id = app_id;
    this->name = name;
    check_name(name);
    insert_to_database();
}

DifferentialEvolutionDB::DifferentialEvolutionDB( MYSQL *conn,
                                                 const string name,
                                                 const vector<double> &min_bound,                                    /* min bound is copied into the search */
                                                 const vector<double> &max_bound,                                    /* max bound is copied into the search */
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
                                               ) throw (string) : DifferentialEvolution(min_bound, max_bound, population_size, parent_selection, number_pairs, recombination_selection, parent_scaling_factor, differential_scaling_factor, crossover_rate, directional, maximum_created, maximum_reported) {
    this->conn = conn;
    this->app_id = -1;
    this->name = name;
    check_name(name);
    insert_to_database();
}

DifferentialEvolutionDB::DifferentialEvolutionDB( MYSQL *conn,
                                                 const int32_t app_id,
                                                 const string name,
                                                 const vector<double> &min_bound,                                    /* min bound is copied into the search */
                                                 const vector<double> &max_bound,                                    /* max bound is copied into the search */
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
                                               ) throw (string) : DifferentialEvolution(min_bound, max_bound, population_size, parent_selection, number_pairs, recombination_selection, parent_scaling_factor, differential_scaling_factor, crossover_rate, directional, maximum_created, maximum_reported) {
    this->conn = conn;
    this->app_id = app_id;
    this->name = name;
    check_name(name);
    insert_to_database();
}

DifferentialEvolutionDB::~DifferentialEvolutionDB() {
    conn = NULL;
}

void
DifferentialEvolutionDB::new_individual(uint32_t &id, vector<double> &parameters, uint32_t &seed) throw (string) {
    DifferentialEvolution::new_individual(id, parameters, seed);
    generation_dirty = true;
}

void
DifferentialEvolutionDB::new_individual(uint32_t &id, vector<double> &parameters) throw (string) {
    DifferentialEvolution::new_individual(id, parameters);
    generation_dirty = true;
}


bool
DifferentialEvolutionDB::insert_individual(uint32_t id, const vector<double> &parameters, double fitness, uint32_t seed) throw (string) {

    if (vector_to_string<double>(parameters).length() >= 2048 || !std::isfinite(fitness))
        return false;
    // Synchronous DB callers generate and insert using the same object.
    // Flush their generation cursor before reloading validator-owned state.
    if (generation_dirty) update_current_individual();
    if (mysql_query(conn, "START TRANSACTION")) throw string(mysql_error(conn));
    try {
        // Serialize validators and reload authoritative state. The generator may
        // have run since this object was cached, or another validator may shrink.
        ostringstream locked_search;
        locked_search << "SELECT * FROM differential_evolution WHERE id = "
                      << this->id << " FOR UPDATE";
        construct_from_database(locked_search.str());
        if (id >= population.size() || parameters.size() != number_parameters) {
            if (mysql_query(conn, "COMMIT")) throw string(mysql_error(conn));
            return false;
        }
        bool modified = DifferentialEvolution::insert_individual(id, parameters, fitness, seed);

        if (modified) {
            ostringstream individual_query;
            individual_query << "UPDATE de_individual"
                             << " SET "
                             << "  fitness = " << setprecision(10) << fitnesses[id]
                             << ", parameters = '" << vector_to_string<double>(population[id]) << "'"
                             << ", seed = " << seeds[id]
                             << " WHERE "
                             << "     differential_evolution_id = " << this->id
                             << " AND position = " << id;

            mysql_query(conn, individual_query.str().c_str());

            if (mysql_errno(conn) != 0) {
                ostringstream ex_msg;
                ex_msg << "ERROR: updating individual with query: '" << individual_query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__ << vector_to_string<double>(parameters) << vector_to_string<double>(parameters).length();
                throw ex_msg.str();
            }

            double best, average, median, worst;
            active_fitness_statistics(best, average, median, worst);

            ostringstream log_query;
            log_query.precision(10);
            log_query << "INSERT INTO differential_evolution_log"
                << " SET "
                << "  search_id = " << this->id
                << ", evaluation = " << this->individuals_reported
                << ", current = '" << setprecision(10) << fixed << fitnesses[id] << "'"
                << ", best = '" << best << "'"
                << ", average = '" << average << "'"
                << ", median = '" << median << "'"
                << ", worst = '" << worst << "'"
                << ", individual = " << id
                << ", seed = " << seed
                << ", global = " << setprecision(10) << (fitnesses[id] == global_best_fitness);

            mysql_query(conn, log_query.str().c_str());

            if (mysql_errno(conn) != 0) {
                ostringstream ex_msg;
                ex_msg << "ERROR: updating differential_evolution_log with query: '" << log_query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
                throw ex_msg.str();
            }


        }

        // Validator-owned fields: persist even when fitness did not improve or
        // the target slot was retired. Generation counters belong to the generator.
        ostringstream state;
        state << "UPDATE differential_evolution SET "
              << "initialized_individuals = " << initialized_individuals
              << ", individuals_reported = " << individuals_reported
              << ", completed_evaluations = " << completed_evaluations
              << ", population_size = " << population_size
              << ", active_ids = '" << vector_to_string<uint32_t>(active_ids) << "'"
              << ", memory_index = " << memory_index
              << ", MF = '" << vector_to_string<double>(MF) << "'"
              << ", MCR = '" << vector_to_string<double>(MCR) << "'"
              << " WHERE id = " << this->id;
        if (mysql_query(conn, state.str().c_str())) throw string(mysql_error(conn));
        if (mysql_query(conn, "COMMIT")) throw string(mysql_error(conn));
        return modified;
    } catch (...) {
        mysql_query(conn, "ROLLBACK");
        // Pending adaptive successes may include the rolled-back evaluation.
        success_pool.clear();
        throw;
    }
}


void
DifferentialEvolutionDB::update_current_individual() throw (string) {
    ostringstream query;
    query << " UPDATE differential_evolution"
        << " SET "
        // Only generator-owned fields. A stale generator must not undo a
        // validator's active set, report count, or adaptive memory updates.
        << "  current_individual = " << current_individual
        << ", current_iteration = GREATEST(current_iteration, " << current_iteration << ")"
        << ", individuals_created = GREATEST(individuals_created, " << individuals_created << ")"
        << ", last_Fi = '" << vector_to_string<double>(last_Fi) << "'"
        << ", last_CRi = '" << vector_to_string<double>(last_CRi) << "'"
        << " WHERE "
        << "    id = " << id << endl;

    mysql_query(conn, query.str().c_str());

    if (mysql_errno(conn) != 0) {
        ostringstream ex_msg;
        ex_msg << "ERROR: updating 'differential_evolution' with query: '" << query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }   
    generation_dirty = false;
}

void
DifferentialEvolutionDB::add_searches(MYSQL *conn, int32_t app_id, vector<EvolutionaryAlgorithmDB*> &searches) throw (string) {
    ostringstream query;
    query << "SELECT id FROM differential_evolution WHERE app_id = " << app_id;

    mysql_query(conn, query.str().c_str());
    MYSQL_RES *result = mysql_store_result(conn);

    if (mysql_errno(conn) != 0) {
        ostringstream ex_msg;
        ex_msg << "ERROR: getting unfinished searches with query: '" << query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }   

    MYSQL_ROW individual_row;

    while ((individual_row = mysql_fetch_row(result))) {
        if (mysql_errno(conn) != 0) {
            ostringstream ex_msg;
            ex_msg << "ERROR: getting row for unfinished searches with query: '" << query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
            throw ex_msg.str();
        }   

        searches.push_back(new DifferentialEvolutionDB(conn, atoi(individual_row[0])));
    }   
    mysql_free_result(result);
}


/** TODO: remove this and add an 'remote_finished_searches'/'remove_finished_searches' method to EvolutionaryAlgorithmBM **/
void
DifferentialEvolutionDB::add_unfinished_searches(MYSQL *conn, int32_t app_id, vector<EvolutionaryAlgorithmDB*> &unfinished_searches) throw (string) {
    ostringstream query;
    query << "SELECT id FROM differential_evolution WHERE app_id = " << app_id;

    mysql_query(conn, query.str().c_str());
    MYSQL_RES *result = mysql_store_result(conn);

    if (mysql_errno(conn) != 0) {
        ostringstream ex_msg;
        ex_msg << "ERROR: getting unfinished searches with query: '" << query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
        throw ex_msg.str();
    }   

    MYSQL_ROW individual_row;

    while ((individual_row = mysql_fetch_row(result))) {
        if (mysql_errno(conn) != 0) {
            ostringstream ex_msg;
            ex_msg << "ERROR: getting row for unfinished searches with query: '" << query.str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << __FILE__ << ":" << __LINE__;
            throw ex_msg.str();
        }   

        DifferentialEvolutionDB *search = new DifferentialEvolutionDB(conn, atoi(individual_row[0]));

        if (search->is_running()) {
            unfinished_searches.push_back(search);
        } else {
            delete search;
        }   
    }   
    mysql_free_result(result);
}

void
DifferentialEvolutionDB::print_to(ostream& stream) {
    stream  << "[DifferentialEvolutionDB " << endl
            << "    id = " << id << endl
            << "    name = '" << name << "'" << endl
            << "    parent_selection = " << parent_selection
            << "    number_pairs = " << number_pairs
            << "    recombination_selection = " << recombination_selection
            << "    parent_scaling_factor = " << parent_scaling_factor
            << "    differential_scaling_factor = " << differential_scaling_factor
            << "    crossover_rate = " << crossover_rate
            << "    directional = " << directional
            << "    current_individual = " << current_individual << endl
            << "    initialized_individuals = " << initialized_individuals << endl
            << "    current_iteration = " << current_iteration << endl
            << "    maximum_iterations = " << maximum_iterations << endl
            << "    individuals_created = " << individuals_created << endl
            << "    maximum_created = " << maximum_created << endl
            << "    individuals_reported = " << individuals_reported << endl
            << "    maximum_reported = " << maximum_reported << endl
            << "    population_size = " << population_size << endl
            << "    min_bound = '" << vector_to_string<double>(min_bound) << "'" << endl
            << "    max_bound = '" << vector_to_string<double>(max_bound) << "'" << endl
            << "    app_id = " << app_id << endl
            << "]" << endl;

    for (uint32_t i : active_ids) {
        stream << "    [DEIndividual" << endl
               << "        differential_evolution_id = " << id << endl
               << "        position = " << i << endl
               << "        fitness = " << setprecision(10) << fitnesses[i] << endl
               << "        parameters = '" << vector_to_string<double>(population[i]) << "'" << endl
               << "        seed = " << seeds[i] << endl
               << "    ]" << endl;
    }
}

ostream& operator<< (ostream& stream, DifferentialEvolutionDB &ps) {
    ps.print_to(stream);
    return stream;
}
