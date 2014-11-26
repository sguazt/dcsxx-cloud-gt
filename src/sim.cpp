/**
 * \file sim.cpp
 *
 * \brief Source code for the ExtremeGreen 2014 paper.
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 *
 * <hr/>
 *
 * Copyright 2013 Marco Guazzone (marco.guazzone@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
#include <cstddef>
#include <cctype>
#include <dcs/algorithm/combinatorics.hpp>
#include <dcs/assert.hpp>
#include <dcs/cli.hpp>
#include <dcs/debug.hpp>
#include <dcs/exception.hpp>
#include <dcs/macro.hpp>
#include <dcs/math/traits/float.hpp>
#include <fstream>
#include <gtpack/cooperative.hpp>
#if DCS_CLOUD_GT_HAVE_CPLEX_SOLVER
# include <ilconcert/iloalg.h>
# include <ilconcert/iloenv.h>
# include <ilconcert/iloexpression.h>
# include <ilconcert/ilomodel.h>
# include <ilcplex/ilocplex.h>
#elif DCS_CLOUD_GT_HAVE_GUROBI_SOLVER
# include <gurobi_c++.h>
#endif // DCS_CLOUD_GT_HAVE_*_SOLVER
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>


//#define CIP1 0
//#define CIP2 1
//#define CIP3 2
//#define CIP12 3
//#define CIP13 4
//#define CIP23 5
//#define CIP123 6


namespace alg = dcs::algorithm;
namespace cli = dcs::cli;
namespace math = dcs::math;


namespace detail { namespace experiment { namespace /*<unnamed>*/ {

enum coalition_formation_category
{
	merge_split_stable_coalition_formation,
	nash_stable_coalition_formation,
	pareto_optimal_coalition_formation,
	social_optimum_coalition_formation
};

enum coalition_value_division_category
{
	banzhaf_coalition_value_division,
	normalized_banzhaf_coalition_value_division,
	shapley_coalition_value_division
};

template <typename RealT>
struct options
{
	options()
	: opt_relative_gap(0),
	  opt_time_lim(-1),
	  coalition_formation(nash_stable_coalition_formation),
	  coalition_value_division(shapley_coalition_value_division),
	  rnd_gen_vms(false),
	  rnd_gen_pm_power_states(false),
	  rnd_gen_pm_on_off_costs(false),
	  rnd_gen_vm_migration_costs(false),
	  rnd_seed(5489),
	  rnd_num_iters(1),
	  csv_fname()
	{
	}

	RealT opt_relative_gap; ///< The relative gap option to set to the optimal solver
	RealT opt_time_lim; ///< The time limit (in sec) to set for each execution of the optimal solver
	coalition_formation_category coalition_formation;
	coalition_value_division_category coalition_value_division;
	bool rnd_gen_vms; ///< Tells if the number of VMs per CIP should be generated at random; if \c true, the number of VMs is randomly generated according to a integer uniform distribution in [0, scenario.cip_num_vms[cip]]
	bool rnd_gen_pm_power_states; ///< Tells if the power state of PMs per CIP should be generated at random; if \c true, the power state of PMs is randomly generated according to a Bernoulli distribution with parameter p=0.5
	bool rnd_gen_pm_on_off_costs; ///< Tells if the switch-on/off cost of PMs per CIP and PM type should be generated at random; if \c true, the switch-on/off cost of PMs is randomly generated according to a Normal distribution with parameter mu=3e-4sec and sigma=5e-5sec
	bool rnd_gen_vm_migration_costs; ///< Tells if the CIP-to-CIP migration cost of VM per CIP and VM types should be generated at random; if \c true, the CIP-to-CIP migration cost of VMs is randomly generated according to a Normal distribution with parameter mu=277sec and sigma=182sec, for the the smaller VM and double for increasing VM size.
	unsigned long rnd_seed; ///< The seed used for random generation
	std::size_t rnd_num_iters; ///< Number of iterations (used only if rnd_vms is true)
	std::string csv_fname; ///< Name of CSV file where to export coalitions enumeration.
};

template <typename RealT>
struct scenario
{
	std::size_t num_cips; ///< Number of different CIPs 
	std::size_t num_pm_types; ///< Number of different PM types
	std::size_t num_vm_types; ///< Number of different VM types
	std::vector< std::vector<std::size_t> > cip_num_pms; ///< Number of PMs per CIP and PM type
	std::vector< std::vector<std::size_t> > cip_num_vms; ///< Number of VMs per CIP and VM type
	std::vector< std::vector<bool> > cip_pm_power_states; ///< Power states of PMs per CIP and PM
	std::vector< std::vector<RealT> > cip_revenues; ///< Revenues per CIP and VM type ($/hour/VM)
	std::vector<RealT> cip_electricity_costs; ///< Energy cost per CIP (in $/kWh)
	std::vector< std::vector<RealT> > cip_pm_asleep_costs; ///< Costs to switch-off PMs, per CIP and PM type ($/hour)
	std::vector< std::vector<RealT> > cip_pm_awake_costs; ///< Costs to switch-on PMs, per CIP and PM type ($/hour)
	std::vector< std::vector< std::vector<RealT> > > cip_to_cip_vm_migration_costs; ///< Costs to migrate VMs from a CIP to another CIP, per CIP and VM type ($/hour)
	std::vector<RealT> pm_spec_min_powers; ///< Min power consumption per PM (in W)
	std::vector<RealT> pm_spec_max_powers; ///< Max power consumption per PM (in W)
	std::vector< std::vector<RealT> > vm_spec_cpus; ///< CPU share requirements per VM type and per PM type
	std::vector< std::vector<RealT> > vm_spec_rams; ///< RAM share requirements per VM type and per PM type
};

template <typename RealT>
struct optimal_allocation_info
{
	optimal_allocation_info()
	: solved(false),
	  optimal(false),
	  objective_value(std::numeric_limits<RealT>::infinity()),
	  cost(std::numeric_limits<RealT>::infinity()),
	  kwatt(std::numeric_limits<RealT>::infinity())
	{
	}


	bool solved;
	bool optimal;
	RealT objective_value;
	RealT cost;
	RealT kwatt;
	std::vector<bool> pm_power_states;
	std::vector< std::vector<bool> > pm_vm_allocations;
};

template <typename RealT>
struct cip_allocation_info
{
	std::size_t num_on_pms; // Number of powered on PMs
	std::size_t num_vms; // Number of hosted VMs
	RealT tot_watt; // Total consumed watt (in Watt)
	//RealT tot_wcost; // Cost rate due to total consumed watt (in kWh)
};

template <typename RealT>
struct coalition_info
{
	coalition_info()
	: optimal_allocation(),
	  value(::std::numeric_limits<RealT>::quiet_NaN()),
	  core_empty(true),
	  payoffs(),
	  payoffs_in_core(false),
	  cid(gtpack::empty_coalition_id)
	{
	}


	optimal_allocation_info<RealT> optimal_allocation;
	RealT value;
	bool core_empty;
	std::map<gtpack::player_type, RealT> payoffs;
	bool payoffs_in_core;
	gtpack::cid_type cid;
};

template <typename RealT>
struct partition_info
{
	RealT value;
	std::set<gtpack::cid_type> coalitions;
	std::map<gtpack::player_type, RealT> payoffs;
	std::map<gtpack::player_type, RealT> side_payments;
};

template <typename RealT>
struct coalition_formation_info
{
	std::map< gtpack::cid_type, coalition_info<RealT> > coalitions;
	std::vector< partition_info<RealT> > best_partitions;
};


template <typename RealT>
scenario<RealT> make_scenario(std::string const& fname)
{
	DCS_ASSERT(!fname.empty(),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Invalid scenario file name"));

	scenario<RealT> s;

	std::ifstream ifs(fname.c_str());

	DCS_ASSERT(ifs,
			   DCS_EXCEPTION_THROW(std::runtime_error, "Cannot open scenario file"));

	for (std::string line; std::getline(ifs, line); )
	{
		std::size_t pos(0);
		for (; pos < line.length() && std::isspace(line[pos]); ++pos)
		{
			; // empty
		}
		if (pos > 0)
		{
			line = line.substr(pos);
		}
		if (line.empty() || line.at(0) == '#')
		{
			// Skip either empty or comment lines
			continue;
		}

		boost::to_lower(line);
		if (boost::starts_with(line, "num_cips"))
		{
			std::istringstream iss(line);

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			iss >> s.num_cips;
		}
		else if (boost::starts_with(line, "num_pm_types"))
		{
			std::istringstream iss(line);

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			iss >> s.num_pm_types;
		}
		else if (boost::starts_with(line, "num_vm_types"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			iss >> s.num_vm_types;
		}
		else if (boost::starts_with(line, "cip_revenues"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

//			s.cip_revenues.resize(s.num_cips);
//			for (std::size_t c = 0; c < s.num_cips; ++c)
//			{
//				s.cip_revenues[c].resize(s.num_vm_types);
//				for (std::size_t v = 0; v < s.num_vm_types; ++v)
//				{
//					iss >> s.cip_revenues[c][v];
//				}
//			}
			s.cip_revenues.resize(s.num_cips);
			for (std::size_t c = 0; c < s.num_cips; ++c)
			{
				// Move to '['
				iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

				s.cip_revenues[c].resize(s.num_vm_types);
				for (std::size_t v = 0; v < s.num_vm_types; ++v)
				{
					iss >> s.cip_revenues[c][v];
				}

				// Move to ']'
				iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
			}
		}
		else if (boost::starts_with(line, "pm_spec_min_powers"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.pm_spec_min_powers.resize(s.num_pm_types);
			for (std::size_t p = 0; p < s.num_pm_types; ++p)
			{
				iss >> s.pm_spec_min_powers[p];
			}
		}
		else if (boost::starts_with(line, "pm_spec_max_powers"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.pm_spec_max_powers.resize(s.num_pm_types);
			for (std::size_t p = 0; p < s.num_pm_types; ++p)
			{
				iss >> s.pm_spec_max_powers[p];
			}
		}
		else if (boost::starts_with(line, "cip_num_pms"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.cip_num_pms.resize(s.num_cips);
			for (std::size_t c = 0; c < s.num_cips; ++c)
			{
				// Move to '['
				iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

				s.cip_num_pms[c].resize(s.num_pm_types);
				for (std::size_t p = 0; p < s.num_pm_types; ++p)
				{
					iss >> s.cip_num_pms[c][p];
				}

				// Move to ']'
				iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
			}
		}
		else if (boost::starts_with(line, "cip_num_vms"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.cip_num_vms.resize(s.num_cips);
			for (std::size_t c = 0; c < s.num_cips; ++c)
			{
				// Move to '['
				iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

				s.cip_num_vms[c].resize(s.num_vm_types);
				for (std::size_t v = 0; v < s.num_vm_types; ++v)
				{
					iss >> s.cip_num_vms[c][v];
				}

				// Move to ']'
				iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
			}
		}
		else if (boost::starts_with(line, "cip_pm_power_states"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.cip_pm_power_states.resize(s.num_cips);
			for (std::size_t c = 0; c < s.num_cips; ++c)
			{
				// Move to '['
				iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

				std::size_t num_pms = 0;
				for (std::size_t p = 0; p < s.num_pm_types; ++p)
				{
					num_pms += s.cip_num_pms[c][p];
				}
				s.cip_pm_power_states[c].resize(num_pms);
				for (std::size_t p = 0; p < num_pms; ++p)
				{
					bool ison = false;
					iss >> ison;
					s.cip_pm_power_states[c][p] = ison;
				}

				// Move to ']'
				iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
			}
		}
		else if (boost::starts_with(line, "cip_wcosts") || boost::starts_with(line, "cip_electricity_costs"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.cip_electricity_costs.resize(s.num_cips);
			for (std::size_t c = 0; c < s.num_cips; ++c)
			{
				iss >> s.cip_electricity_costs[c];
			}
		}
		else if (boost::starts_with(line, "cip_pm_asleep_costs"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.cip_pm_asleep_costs.resize(s.num_cips);
			for (std::size_t c = 0; c < s.num_cips; ++c)
			{
				// Move to '['
				iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

				s.cip_pm_asleep_costs[c].resize(s.num_pm_types);
				for (std::size_t p = 0; p < s.num_pm_types; ++p)
				{
					iss >> s.cip_pm_asleep_costs[c][p];
				}

				// Move to ']'
				iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
			}
		}
		else if (boost::starts_with(line, "cip_pm_awake_costs"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.cip_pm_awake_costs.resize(s.num_cips);
			for (std::size_t c = 0; c < s.num_cips; ++c)
			{
				// Move to '['
				iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

				s.cip_pm_awake_costs[c].resize(s.num_pm_types);
				for (std::size_t p = 0; p < s.num_pm_types; ++p)
				{
					iss >> s.cip_pm_awake_costs[c][p];
				}

				// Move to ']'
				iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
			}
		}
		else if (boost::starts_with(line, "vm_spec_cpus"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.vm_spec_cpus.resize(s.num_vm_types);
			for (std::size_t v = 0; v < s.num_vm_types; ++v)
			{
				// Move to '['
				iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

				s.vm_spec_cpus[v].resize(s.num_pm_types);
				for (std::size_t p = 0; p < s.num_pm_types; ++p)
				{
					iss >> s.vm_spec_cpus[v][p];
				}

				// Move to ']'
				iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
			}
		}
		else if (boost::starts_with(line, "vm_spec_rams"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.vm_spec_rams.resize(s.num_vm_types);
			for (std::size_t v = 0; v < s.num_vm_types; ++v)
			{
				// Move to '['
				iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

				s.vm_spec_rams[v].resize(s.num_pm_types);
				for (std::size_t p = 0; p < s.num_pm_types; ++p)
				{
					iss >> s.vm_spec_rams[v][p];
				}

				// Move to ']'
				iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
			}
		}
		else if (boost::starts_with(line, "cip_to_cip_vm_migration_costs"))
		{
			std::istringstream iss(line.substr(pos));

			// Move to '='
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('=' is missing)"));

			// Move to '['
			iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
			DCS_ASSERT(iss.good(),
					   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

			s.cip_to_cip_vm_migration_costs.resize(s.num_cips);
			for (std::size_t c = 0; c < s.num_cips; ++c)
			{
				// Move to '['
				iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

				s.cip_to_cip_vm_migration_costs.resize(s.num_cips);
				for (std::size_t c = 0; c < s.num_cips; ++c)
				{
					// Move to '['
					iss.ignore(std::numeric_limits<std::streamsize>::max(), '[');
					DCS_ASSERT(iss.good(),
							   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file ('[' is missing)"));

					s.cip_to_cip_vm_migration_costs[c].resize(s.num_vm_types);
					for (std::size_t p = 0; p < s.num_pm_types; ++p)
					{
						iss >> s.cip_to_cip_vm_migration_costs[c][c][p];
					}

					// Move to ']'
					iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
					DCS_ASSERT(iss.good(),
							   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
				}

				// Move to ']'
				iss.ignore(std::numeric_limits<std::streamsize>::max(), ']');
				DCS_ASSERT(iss.good(),
						   DCS_EXCEPTION_THROW(std::runtime_error, "Malformed scenario file (']' is missing)"));
			}
		}
	}

	// check: mandatory info
	DCS_ASSERT(s.num_cips > 0,
			   DCS_EXCEPTION_THROW(std::logic_error, "Number of CIP must be a positive number"));
	DCS_ASSERT(s.num_pm_types > 0,
			   DCS_EXCEPTION_THROW(std::logic_error, "Number of PM types must be a positive number"));
	DCS_ASSERT(s.num_vm_types > 0,
			   DCS_EXCEPTION_THROW(std::logic_error, "Number of VM types must be a positive number"));

	// Consistency checks
	DCS_ASSERT(s.cip_revenues.size() == 0 || s.num_cips == s.cip_revenues.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of CIP revenues"));
	DCS_ASSERT(s.cip_num_pms.size() == 0 || s.num_cips == s.cip_num_pms.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of CIP PMs"));
	DCS_ASSERT(s.cip_num_vms.size() == 0 || s.num_cips == s.cip_num_vms.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of CIP VMs"));
	DCS_ASSERT(s.cip_electricity_costs.size() == 0 || s.num_cips == s.cip_electricity_costs.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of CIP electricity costs"));
	DCS_ASSERT(s.cip_pm_power_states.size() == 0 || s.num_cips == s.cip_pm_power_states.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of CIP PM power states"));
	DCS_ASSERT(s.cip_pm_asleep_costs.size() == 0 || s.num_cips == s.cip_pm_asleep_costs.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of CIP PM switch-off costs"));
	DCS_ASSERT(s.cip_pm_awake_costs.size() == 0 || s.num_cips == s.cip_pm_awake_costs.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of CIP PM switch-on costs"));
	DCS_ASSERT(s.pm_spec_min_powers.size() == 0 || s.num_pm_types == s.pm_spec_min_powers.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of PM minimum power consumption specifications"));
	DCS_ASSERT(s.pm_spec_max_powers.size() == 0 || s.num_pm_types == s.pm_spec_max_powers.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of PM maximum power consumption specifications"));
	DCS_ASSERT(s.vm_spec_cpus.size() == 0 || s.num_vm_types == s.vm_spec_cpus.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of VM CPU share requirements"));
	DCS_ASSERT(s.vm_spec_rams.size() == 0 || s.num_vm_types == s.vm_spec_rams.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of VM RAM share requirements"));
	DCS_ASSERT(s.cip_to_cip_vm_migration_costs.size() == 0 || s.num_cips == s.cip_to_cip_vm_migration_costs.size(),
			   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of CIP-to-CIP VM migration costs"));
//	if (s.num_cips > 0)
//	{
//		DCS_ASSERT(s.cip_revenues[0].size() == 0 || s.num_vm_types == s.cip_revenues[0].size(),
//				   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of VM types in CIP revenues"));
//		DCS_ASSERT(s.cip_num_pms[0].size() == 0 || s.num_pm_types == s.cip_num_pms[0].size(),
//				   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of PM types in CIP number of PMs"));
//		DCS_ASSERT(s.cip_num_vms[0].size() == 0 || s.num_vm_types == s.cip_num_vms[0].size(),
//				   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of VM types in CIP number of VMs"));
//		DCS_ASSERT(s.cip_pm_asleep_costs[0].size() == 0 || s.num_pm_types == s.cip_pm_asleep_costs[0].size(),
//				   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of PM types in CIP switch-off costs of PMs"));
//		DCS_ASSERT(s.cip_pm_awake_costs[0].size() == 0 || s.num_pm_types == s.cip_pm_awake_costs[0].size(),
//				   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of PM types in CIP switch-on costs of PMs"));
//		DCS_ASSERT(s.cip_to_cip_vm_migration_costs[0].size() == 0 || s.num_cips == s.cip_to_cip_vm_migration_costs[0].size(),
//				   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of CIP types in CIP-to-CIP VM migration costs"));
//		DCS_ASSERT(s.cip_to_cip_vm_migration_costs[0][0].size() == 0 || s.num_vm_types == s.cip_to_cip_vm_migration_costs[0][0].size(),
//				   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of VM types in CIP-to-CIP VM migration costs"));
//	}
//	if (s.num_vm_types > 0)
//	{
//		DCS_ASSERT(s.vm_spec_cpus[0].size() == 0 || s.num_pm_types == s.vm_spec_cpus[0].size(),
//				   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of PM types in VM CPU share requirements"));
//		DCS_ASSERT(s.vm_spec_rams[0].size() == 0 || s.num_pm_types == s.vm_spec_rams[0].size(),
//				   DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected number of PM types in VM RAM share requirements"));
//	}

	//TODO: assign default values for the other info
	if (s.cip_pm_power_states.size() == 0)
	{
		// Default: all PMs are off

		s.cip_pm_power_states.resize(s.num_cips);
		for (std::size_t c = 0; c < s.num_cips; ++c)
		{
			std::size_t num_pms = 0;
			for (std::size_t p = 0; p < s.num_pm_types; ++p)
			{
				num_pms += s.cip_num_pms[c][p];
			}
			s.cip_pm_power_states[c].assign(num_pms, false);
		}
	}
	if (s.cip_pm_asleep_costs.size() == 0)
	{
		// Default: all PM switch-off costs are 0

		s.cip_pm_asleep_costs.resize(s.num_cips);
		for (std::size_t c = 0; c < s.num_cips; ++c)
		{
			s.cip_pm_asleep_costs[c].assign(s.num_pm_types, 0);
		}
	}
	if (s.cip_pm_awake_costs.size() == 0)
	{
		// Default: all PM switch-on costs are 0

		s.cip_pm_awake_costs.resize(s.num_cips);
		for (std::size_t c = 0; c < s.num_cips; ++c)
		{
			s.cip_pm_awake_costs[c].assign(s.num_pm_types, 0);
		}
	}
	if (s.cip_to_cip_vm_migration_costs.size() == 0)
	{
		// Default: all CIP-to-CIP VM migration costs are 0

		s.cip_to_cip_vm_migration_costs.resize(s.num_cips);
		for (std::size_t c1 = 0; c1 < s.num_cips; ++c1)
		{
			s.cip_to_cip_vm_migration_costs[c1].resize(s.num_cips);
			for (std::size_t c2 = 0; c2 < s.num_cips; ++c2)
			{
				s.cip_to_cip_vm_migration_costs[c1][c2].assign(s.num_vm_types, 0);
			}
		}
	}

	return s;
}

template <typename CharT, typename CharTraitsT, typename RealT>
std::basic_ostream<CharT,CharTraitsT>& operator<<(std::basic_ostream<CharT,CharTraitsT>& os, scenario<RealT> const& s)
{
	os  << "num_cips=" << s.num_cips
		<< ", " << "num_pm_types=" << s.num_pm_types
		<< ", " << "num_vm_types=" << s.num_vm_types;

	os << ", " << "cip_revenues=[";
	for (std::size_t c = 0; c < s.num_cips; ++c)
	{
		if (c > 0)
		{
			os << " ";
		}

		os << "[";
		for (std::size_t v = 0; v < s.num_vm_types; ++v)
		{
			if (v > 0)
			{
				os << ", ";
			}
			os << s.cip_revenues[c][v];
		}
		os << "]";
	}
	os << "]";
	os << ", " << "pm_spec_min_powers=[";
	for (std::size_t p = 0; p < s.num_pm_types; ++p)
	{
		if (p > 0)
		{
			os << ", ";
		}
		os << s.pm_spec_min_powers[p];
	}
	os << "]";
	os << ", " << "pm_spec_max_powers=[";
	for (std::size_t p = 0; p < s.num_pm_types; ++p)
	{
		if (p > 0)
		{
			os << ", ";
		}
		os << s.pm_spec_max_powers[p];
	}
	os << "]";
	os << ", " << "cip_num_pms=[";
	for (std::size_t c = 0; c < s.num_cips; ++c)
	{
		if (c > 0)
		{
			os << " ";
		}

		os << "[";
		for (std::size_t p = 0; p < s.num_pm_types; ++p)
		{
			if (p > 0)
			{
				os << ", ";
			}
			os << s.cip_num_pms[c][p];
		}
		os << "]";
	}
	os << "]";
	os << ", " << "cip_pm_power_states=[";
	for (std::size_t c = 0; c < s.num_cips; ++c)
	{
		if (c > 0)
		{
			os << " ";
		}

		os << "[";
		for (std::size_t p = 0; p < s.cip_pm_power_states[c].size(); ++p)
		{
			if (p > 0)
			{
				os << ", ";
			}
			os << s.cip_pm_power_states[c][p];
		}
		os << "]";
	}
	os << "]";
	os << ", " << "cip_num_vms=[";
	for (std::size_t c = 0; c < s.num_cips; ++c)
	{
		if (c > 0)
		{
			os << " ";
		}

		os << "[";
		for (std::size_t v = 0; v < s.num_vm_types; ++v)
		{
			if (v > 0)
			{
				os << ", ";
			}
			os << s.cip_num_vms[c][v];
		}
		os << "]";
	}
	os << "]";
	os << ", " << "cip_electricity_costs=[";
	for (std::size_t c = 0; c < s.num_cips; ++c)
	{
		if (c > 0)
		{
			os << ", ";
		}
		os << s.cip_electricity_costs[c];
	}
	os << "]";
	os << ", " << "cip_pm_asleep_costs=[";
	for (std::size_t c = 0; c < s.num_cips; ++c)
	{
		if (c > 0)
		{
			os << "  ";
		}

		os << "[";
		for (std::size_t p = 0; p < s.num_pm_types; ++p)
		{
			if (p > 0)
			{
				os << ", ";
			}
			os << s.cip_pm_asleep_costs[c][p];
		}
		os << "]";
	}
	os << "]";
	os << ", " << "cip_pm_awake_costs=[";
	for (std::size_t c = 0; c < s.num_cips; ++c)
	{
		if (c > 0)
		{
			os << "  ";
		}

		os << "[";
		for (std::size_t p = 0; p < s.num_pm_types; ++p)
		{
			if (p > 0)
			{
				os << ", ";
			}
			os << s.cip_pm_awake_costs[c][p];
		}
		os << "]";
	}
	os << "]";
	os << ", " << "cip_to_cip_vm_migration_costs=[";
	for (std::size_t c1 = 0; c1 < s.num_cips; ++c1)
	{
		if (c1 > 0)
		{
			os << "  ";
		}

		os << "[";
		for (std::size_t c2 = 0; c2 < s.num_cips; ++c2)
		{
			if (c2 > 0)
			{
				os << "  ";
			}

			os << "[";
			for (std::size_t p = 0; p < s.num_pm_types; ++p)
			{
				if (p > 0)
				{
					os << ", ";
				}
				os << s.cip_to_cip_vm_migration_costs[c1][c2][p];
			}
			os << "]";
		}
		os << "]";
	}
	os << "]";
	os << ", " << "vm_spec_cpus=[";
	for (std::size_t v = 0; v < s.num_vm_types; ++v)
	{
		if (v > 0)
		{
			os << " ";
		}

		os << "[";
		for (std::size_t p = 0; p < s.num_pm_types; ++p)
		{
			if (p > 0)
			{
				os << ", ";
			}
			os << s.vm_spec_cpus[v][p];
		}
		os << "]";
	}
	os << "]";
	os << ", " << "vm_spec_rams=[";
	for (std::size_t v = 0; v < s.num_vm_types; ++v)
	{
		if (v > 0)
		{
			os << " ";
		}

		os << "[";
		for (std::size_t p = 0; p < s.num_pm_types; ++p)
		{
			if (p > 0)
			{
				os << ", ";
			}
			os << s.vm_spec_rams[v][p];
		}
		os << "]";
	}
	os << "]";

	return os;
}

template <typename CharT, typename CharTraitsT, typename RealT>
std::basic_ostream<CharT,CharTraitsT>& operator<<(std::basic_ostream<CharT,CharTraitsT>& os, options<RealT> const& opt)
{
	os	<< "relative-gap: " << opt.opt_relative_gap
		<< ", time_limit: " << opt.opt_time_lim
		<< ", coalition_formation: " << opt.coalition_formation
		<< ", coalition_value_division: " << opt.coalition_value_division
		<< ", csv_file_name: " << opt.csv_fname
		<< ", random_gen_vms: " << std::boolalpha << opt.rnd_gen_vms
		<< ", random_gen_pm_power_states: " << std::boolalpha << opt.rnd_gen_pm_power_states
		<< ", random_gen_pm_on_off_costs: " << std::boolalpha << opt.rnd_gen_pm_on_off_costs
		<< ", random_gen_vm_migration_costs: " << std::boolalpha << opt.rnd_gen_vm_migration_costs
		<< ", random_seed: " << opt.rnd_seed
		<< ", random_num_iters: " << opt.rnd_num_iters;


	return os;
}

template <typename CharT, typename CharTraitsT, typename T>
std::basic_ostream<CharT,CharTraitsT>& operator<<(std::basic_ostream<CharT,CharTraitsT>& os, std::vector<T> const& v)
{
	typedef typename std::vector<T>::const_iterator iterator;

	iterator end_it(v.end());
	os << "[";
	for (iterator it = v.begin(); it != end_it; ++it)
	{
		if (it != v.begin())
		{
			os << ", ";
		}
		os << *it;
	}
	os << "]";

	return os;
}

template <typename CharT, typename CharTraitsT, typename T>
std::basic_ostream<CharT,CharTraitsT>& operator<<(std::basic_ostream<CharT,CharTraitsT>& os, std::vector< std::vector<T> > const& v)
{
	typedef typename std::vector< std::vector<T> >::const_iterator outer_iterator;
	typedef typename std::vector<T>::const_iterator inner_iterator;

	outer_iterator out_end_it(v.end());
	os << "[";
	for (outer_iterator out_it = v.begin(); out_it != out_end_it; ++out_it)
	{
		if (out_it != v.begin())
		{
			os << " ";
		}
		os << "[";
		inner_iterator in_end_it(out_it->end());
		for (inner_iterator in_it = out_it->begin(); in_it != in_end_it; ++in_it)
		{
			if (in_it != out_it->begin())
			{
				os << ", ";
			}
			os << *in_it;
		}
		os << "]";
	}
	os << "]";

	return os;
}

template <typename RealT>
inline
RealT pm_consumed_power(RealT min_power, RealT max_power, RealT u)
{
	return min_power + (max_power-min_power)*u;
}

template <typename CipsWCostVectorT,
		  typename PmsCipVectorT,
		  typename PmsCategoryVectorT,
		  typename PmSpecsMinPowerVectorT,
		  typename PmSpecsMaxPowerVectorT,
		  typename VmsCipVectorT,
		  typename VmsCategoryVectorT,
		  typename VmSpecsCpuVectorT,
		  typename VmSpecsRamVectorT,
		  typename PmPowerStatesVectorT,
		  typename CipPmAsleepCostsVectorT,
		  typename CipPmAwakeCostsVectorT,
		  typename CipToCipVmMigrationCostsCubeT,
		  typename RealT>
optimal_allocation_info<RealT> find_optimal_allocation(std::size_t ncips,
													   CipsWCostVectorT cips_electricity_cost,
													   std::size_t npms,
													   PmsCipVectorT const& pms_cip,
													   PmsCategoryVectorT const& pms_category,
													   PmSpecsMinPowerVectorT const& pm_specs_min_power,
													   PmSpecsMaxPowerVectorT const& pm_specs_max_power,
													   std::size_t nvms,
													   VmsCipVectorT const& vms_cip,
													   VmsCategoryVectorT const& vms_category,
													   VmSpecsCpuVectorT const& vm_specs_cpu,
													   VmSpecsRamVectorT const& vm_specs_ram,
													   PmPowerStatesVectorT const& pm_power_states,
													   CipPmAsleepCostsVectorT const& cip_pm_asleep_costs,
													   CipPmAwakeCostsVectorT const& cip_pm_awake_costs,
													   CipToCipVmMigrationCostsCubeT const& cip_to_cip_vm_migration_costs,
													   bool min_power,
													   RealT relative_gap = 0,
													   RealT time_lim = -1)
{
	DCS_MACRO_SUPPRESS_UNUSED_VARIABLE_WARNING( ncips );

	optimal_allocation_info<RealT> solution;

	DCS_DEBUG_TRACE("Finding optimal allocation:");
	DCS_DEBUG_TRACE("- Number of CIPs: " << ncips);
	DCS_DEBUG_TRACE("- Energy Costs per CIP: " << cips_electricity_cost);
	DCS_DEBUG_TRACE("- Number of PMs: " << npms);
	DCS_DEBUG_TRACE("- CIP per PM: " << pms_cip);
	DCS_DEBUG_TRACE("- Category per PM: " << pms_category);
	DCS_DEBUG_TRACE("- Mininimum Power Consumption per PM: " << pm_specs_min_power);
	DCS_DEBUG_TRACE("- Maximum Power Consumption per PM: " << pm_specs_max_power);
	DCS_DEBUG_TRACE("- Number of VMs: " << nvms);
	DCS_DEBUG_TRACE("- Category per VM: " << vms_category);
	DCS_DEBUG_TRACE("- CPU requirement per VM: " << vm_specs_cpu);
	DCS_DEBUG_TRACE("- RAM requirement per VM: " << vm_specs_ram);
	DCS_DEBUG_TRACE("- PM Power States: " << pm_power_states);
	DCS_DEBUG_TRACE("- PM On->Off Cost per CIP and PM Category: " << cip_pm_asleep_costs);
	DCS_DEBUG_TRACE("- PM Off->On Cost per CIP and PM Category: " << cip_pm_awake_costs);
	DCS_DEBUG_TRACE("- VM Migration Cost from CIP to CIP per VM Category: " << cip_to_cip_vm_migration_costs);
	DCS_DEBUG_TRACE("- Minimum Power: " << std::boolalpha << min_power);
	DCS_DEBUG_TRACE("- Relative Gap: " << relative_gap);


#ifdef DCS_CLOUD_GT_HAVE_CPLEX_SOLVER
	typedef IloNumVarArray var_vector_type;
	typedef IloArray<IloNumVarArray> var_matrix_type;

	//Setting up vars
	try
	{
		// Initialize the Concert Technology app
		IloEnv env;

		IloModel model(env);

		if (min_power)
		{
			model.setName("Min-Power Optimal Allocation (CPLEX)");
		}
		else
		{
			model.setName("Min-Cost Optimal Allocation (CPLEX)");
		}

		// Decision Variables

		// Variables y_{vh}: y_{vh}==1 iif VM v is on host h
		var_matrix_type y(env, nvms);
		for (std::size_t v = 0; v < nvms; ++v)
		{
			y[v] = var_vector_type(env, npms);

			for (std::size_t h = 0 ; h < npms ; ++h)
			{
				std::ostringstream oss;
				oss << "y[" << v << "][" << h << "]";
				y[v][h] = IloBoolVar(env, 0, 1, oss.str().c_str());
				model.add(y[v][h]);
			}
		}

		// Variables x_h
		var_vector_type x(env, npms);
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "x[" << h << "]";
			x[h] = IloBoolVar(env, 0, 1, oss.str().c_str());
			model.add(x[h]);
		}

		// Variables s_{h}
		var_vector_type s(env, npms);
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "s[" << h << "]";
			s[h] = IloNumVar(env, 0, 1, ILOFLOAT, oss.str().c_str());
			model.add(s[h]);
		}

		// Constraints

		std::size_t cc(0); // Constraint counter

		// C1: \forall v \in V: \sum_{h \in H} y_{vh} = 1
		++cc;
		for (std::size_t v = 0; v < nvms; ++v)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << v << "}";

			IloConstraint cons(IloSum(y[v]) == 1);
			cons.setName(oss.str().c_str());
			model.add(cons);
		}

		// C2: \forall h \in H: \sum_{v \in V} y_{vh} \le |V|*x_{h}
		++cc;
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << h << "}";

			IloExpr lhs(env);
			for (std::size_t v = 0; v < nvms; ++v)
			{
				lhs += y[v][h];
			}

			IloConstraint cons(lhs <= x[h]*IloInt(nvms));
			cons.setName(oss.str().c_str());
			model.add(cons);
		}

		// C3: \forall h \in H: \sum_{v \in V} y_{vh}M_{q(v),g(h)} \le x_{h}
		++cc;
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << h << "}";

			IloExpr lhs(env);
			for (std::size_t v = 0; v < nvms; ++v)
			{
				RealT req(vm_specs_ram[vms_category[v]][pms_category[h]]);
				lhs += y[v][h]*req;
			}

			IloConstraint cons(lhs <= x[h]);
			cons.setName(oss.str().c_str());
			model.add(cons);
		}

		// C4: \forall h \in H: \sum_{v \in V} y_{vh}S_{q(v),g(h)} == s_{h}
		++cc;
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << h << "}";

			IloExpr lhs(env);
			for (std::size_t v = 0; v < nvms; ++v)
			{
				RealT req(vm_specs_cpu[vms_category[v]][pms_category[h]]);
				lhs += y[v][h]*req;
			}

			IloConstraint cons(lhs == s[h]);
			cons.setName(oss.str().c_str());
			model.add(cons);
		}

		// C5: \forall h \in H: s_{h} \le x_{h}
		++cc;
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << h << "}";

			IloConstraint cons(s[h] <= x[h]);
			cons.setName(oss.str().c_str());
			model.add(cons);
		}

		// Set objective
		IloObjective z;
		if (min_power)
		{
			//FIXME: this does not work well when PM switch-on/off costs and VM migration costs are != zero!
			std::cerr << "(W) Power optimization does not work well when PM switch-on/off costs and VM migration costs are not zero!" << std::endl;

			// z = \min sum_{h \in H}{x_h C_{g(h)}^{min} + (C_{g(h)}^{max}-C_{g(h)}^{min})s_{h} + x_h(1-o(i)) + (1-x_h)o(i)}
			IloExpr expr(env);
			for (std::size_t h = 0; h < npms; ++h)
			{
				RealT dC = pm_specs_max_power[pms_category[h]]-pm_specs_min_power[pms_category[h]];

				expr += x[h]*pm_specs_min_power[pms_category[h]] + dC*s[h];
			}
			z = IloMinimize(env, expr);
		}
		else
		{
			IloExpr expr(env);
			for (std::size_t h = 0; h < npms; ++h)
			{
				RealT dC = pm_specs_max_power[pms_category[h]]-pm_specs_min_power[pms_category[h]];
				RealT wcost = cips_electricity_cost[pms_cip[h]]*1e-3; // Electricity cost in Wh

				expr += (x[h]*pm_specs_min_power[pms_category[h]]+dC*s[h])*wcost;
				// Add PM switch-on/off costs
				expr += x[h]*(1-pm_power_states[h])*cip_pm_awake_costs[pms_cip[h]][pms_category[h]]
					 +  (1-x[h])*pm_power_states[h]*cip_pm_asleep_costs[pms_cip[h]][pms_category[h]];
				// Add VM migration costs
				for (std::size_t v = 0; v < nvms; ++v)
				{
					expr += y[v][h]*cip_to_cip_vm_migration_costs[vms_cip[v]][pms_cip[h]][vms_category[v]];
				}
			}
			z = IloMinimize(env, expr);
		}
		model.add(z);

		// Create the CPLEX solver and make 'model' the active ("extracted") model
		IloCplex solver(model);

		//write model
#ifndef DCS_DEBUG
		solver.setOut(env.getNullStream());
		solver.setWarning(env.getNullStream());
//#else // DCS_DEBUG
//		solver.exportModel("cplex-model.lp");
#endif // DCS_DEBUG

		// Set Relative Gap to (relative_gap*100)%: CPLEX will stop as soon as it has found a feasible integer solution proved to be within (relative_gap*100)% of optimal.
		if (math::float_traits<RealT>::definitely_greater(relative_gap, 0))
		{
			solver.setParam(IloCplex::EpGap, relative_gap);
		}
		if (math::float_traits<RealT>::definitely_greater(time_lim, 0))
		{
			solver.setParam(IloCplex::TiLim, time_lim);
		}

		solution.solved = solver.solve();
		solution.optimal = false;

		IloAlgorithm::Status status = solver.getStatus();

		switch (status)
		{
			case IloAlgorithm::Optimal: // The algorithm found an optimal solution.
				solution.objective_value = static_cast<RealT>(solver.getObjValue());
				solution.optimal = true;
				break;
			case IloAlgorithm::Feasible: // The algorithm found a feasible solution, though it may not necessarily be optimal.

				solution.objective_value = static_cast<RealT>(solver.getObjValue());
				std::cerr << "(W) Optimization problem solved but non-optimal!" << std::endl;
				break;
			case IloAlgorithm::Infeasible: // The algorithm proved the model infeasible (i.e., it is not possible to find an assignment of values to variables satisfying all the constraints in the model).
			case IloAlgorithm::Unbounded: // The algorithm proved the model unbounded.
			case IloAlgorithm::InfeasibleOrUnbounded: // The model is infeasible or unbounded.
			case IloAlgorithm::Error: // An error occurred and, on platforms that support exceptions, that an exception has been thrown.
			case IloAlgorithm::Unknown: // The algorithm has no information about the solution of the model.
			{
				::std::ostringstream oss;
				std::cerr << "Optimization was stopped with status = " << status << " (CPLEX status = " << solver.getCplexStatus() << ", sub-status = " << solver.getCplexSubStatus() << ")" << std::endl;
				return solution;
			}
		}

#ifdef DCS_DEBUG
		DCS_DEBUG_TRACE( "Optimal solution: " );

		DCS_DEBUG_TRACE( "- Solved: " << std::boolalpha << solution.solved );
		DCS_DEBUG_TRACE( "- Optimal: " << std::boolalpha << solution.optimal );

		DCS_DEBUG_TRACE( "- Decision variables: " );

		// Output x_{h}
		for (std::size_t h = 0; h < npms; ++h)
		{
			DCS_DEBUG_STREAM << x[h].getName() << " = " << solver.getValue(x[h]) << " (" << static_cast<bool>(IloRound(solver.getValue(x[h]))) << ")" << ::std::endl;
		}

		// Output y_{vh}
		for (std::size_t v = 0; v < nvms; ++v)
		{
			for (std::size_t h = 0; h < npms; ++h)
			{
				DCS_DEBUG_STREAM << y[v][h].getName() << " = " << solver.getValue(y[v][h]) << " (" << static_cast<bool>(IloRound(solver.getValue(y[v][h]))) << ")" << ::std::endl;
			}
		}

		// Output s_{h}
		for (std::size_t h = 0; h < npms; ++h)
		{
			DCS_DEBUG_STREAM << s[h].getName() << " = " << solver.getValue(s[h]) << ::std::endl;
		}

		//Print z
		DCS_DEBUG_TRACE( "- Objective value: " << solution.objective_value );
#endif // DCS_DEBUG

		solution.pm_power_states.resize(npms);
		solution.pm_vm_allocations.resize(npms);
		if (min_power)
		{
			// Computed in the loop below
			solution.cost = 0;
		}
		else
		{
			solution.cost = solution.objective_value;
		}
		for (std::size_t h = 0; h < npms; ++h)
		{
			solution.pm_power_states[h] = static_cast<bool>(IloRound(solver.getValue(x[h])));

			// Compute the energy cost
			if (min_power && static_cast<bool>(IloRound(solver.getValue(x[h]))))
			{
				//RealT dC = pm_specs_max_power[pms_category[h]]-pm_specs_min_power[pms_category[h]];
				RealT wcost = cips_electricity_cost[pms_cip[h]]*1e-3; // Electricity cost in Wh

				//solution.cost += (pm_specs_min_power[pms_category[h]] + dC*static_cast<RealT>(solver.getValue(s[h])))*wcost;
				solution.cost += pm_consumed_power(pm_specs_min_power[pms_category[h]], pm_specs_max_power[pms_category[h]], static_cast<RealT>(solver.getValue(s[h])))*wcost;
			}
			solution.pm_vm_allocations[h].resize(nvms);
			for (std::size_t v = 0; v < nvms; ++v)
			{
				solution.pm_vm_allocations[h][v] = static_cast<bool>(IloRound(solver.getValue(y[v][h])));
			}
		}

		x.end();
		y.end();
		s.end();

		// Close the Concert Technology app
		env.end();
	}
	catch (IloException const& e)
	{
		::std::ostringstream oss;
		oss << "Got exception from CPLEX: " << e.getMessage();
		DCS_EXCEPTION_THROW(::std::runtime_error, oss.str());
	}
	catch (...)
	{
		DCS_EXCEPTION_THROW(::std::runtime_error,
							"Unexpected error during the optimization");
	}
#elif defined(DCS_CLOUD_GT_HAVE_GUROBI_SOLVER)
	typedef std::vector<GRBVar> var_vector_type;
	typedef std::vector< std::vector<GRBVar> > var_matrix_type;


	//Setting up vars
	try
	{
		// Initialize the Gurobi environment
		GRBEnv env;

#ifdef DCS_DEBUG
		env.set(GRB_IntParam_OutputFlag, 1);
		env.set(GRB_IntParam_LogToConsole, 1);
#else
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_LogToConsole, 0);
#endif // DCS_DEBUG

		// Set Relative Gap to (relative_gap*100)%: CPLEX will stop as soon as it has found a feasible integer solution proved to be within (relative_gap*100)% of optimal.
		if (math::float_traits<RealT>::definitely_greater(relative_gap, 0))
		{
			env.set(GRB_DoubleParam_MIPGap, relative_gap);
		}
		if (math::float_traits<RealT>::definitely_greater(time_lim, 0))
		{
			env.set(GRB_DoubleParam_TimeLimit, time_lim);
		}

		GRBModel model(env);

		if (min_power)
		{
			model.set(GRB_StringAttr_ModelName, "Min-Power Optimal Allocation (GUROBI)");
		}
		else
		{
			model.set(GRB_StringAttr_ModelName, "Min-Cost Optimal Allocation (GUROBI)");
		}

		// Decision Variables

		// Variables y_{vh}: y_{vh}==1 iif VM v is on host h
		var_matrix_type y(nvms, var_vector_type(npms));
		for (std::size_t v = 0; v < nvms; ++v)
		{
			for (std::size_t h = 0 ; h < npms ; ++h)
			{
				std::ostringstream oss;
				oss << "y[" << v << "][" << h << "]";
				y[v][h] = model.addVar(0, 1, 0, GRB_BINARY, oss.str());
			}
		}

		// Variables x_h
		var_vector_type x(npms);
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "x[" << h << "]";
			x[h] = model.addVar(0, 1, 0, GRB_BINARY, oss.str());
		}

		// Variables s_{h}
		var_vector_type s(npms);
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "s[" << h << "]";
			s[h] = model.addVar(0, 1, 0, GRB_CONTINUOUS, oss.str());
		}

		// Integrates new variables
		model.update();

		// Constraints

		std::size_t cc(0); // Constraint counter

		// C1: \forall v \in V: \sum_{h \in H} y_{vh} = 1
		++cc;
		for (std::size_t v = 0; v < nvms; ++v)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << v << "}";

			GRBLinExpr lhs(0);
			for (std::size_t h = 0; h < npms; ++h)
			{
				lhs += y[v][h];
			}

			model.addConstr(lhs, GRB_EQUAL, 1, oss.str());
		}

		// C2: \forall h \in H: \sum_{v \in V} y_{vh} \le |V|*x_{h}
		++cc;
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << h << "}";

			GRBLinExpr lhs(0);
			for (std::size_t v = 0; v < nvms; ++v)
			{
				lhs += y[v][h];
			}

			model.addConstr(lhs, GRB_LESS_EQUAL, x[h]*static_cast<RealT>(nvms), oss.str());
		}

		// C3: \forall h \in H: \sum_{v \in V} y_{vh}M_{q(v),g(h)} \le x_{h}
		++cc;
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << h << "}";

			GRBLinExpr lhs(0);
			for (std::size_t v = 0; v < nvms; ++v)
			{
				RealT req(vm_specs_ram[vms_category[v]][pms_category[h]]);
				lhs += y[v][h]*req;
			}

			model.addConstr(lhs, GRB_LESS_EQUAL, x[h], oss.str());
		}

		// C4: \forall h \in H: \sum_{v \in V} y_{vh}S_{q(v),g(h)} == s_{h}
		++cc;
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << h << "}";

			GRBLinExpr lhs(0);
			for (std::size_t v = 0; v < nvms; ++v)
			{
				RealT req(vm_specs_cpu[vms_category[v]][pms_category[h]]);
				lhs += y[v][h]*req;
			}

			model.addConstr(lhs, GRB_EQUAL, s[h], oss.str());
		}

		// C5: \forall h \in H: s_{h} \le x_{h}
		++cc;
		for (std::size_t h = 0; h < npms; ++h)
		{
			std::ostringstream oss;
			oss << "C" << cc << "_{" << h << "}";

			model.addConstr(s[h], GRB_LESS_EQUAL, x[h], oss.str());
		}

		// Set objective
		GRBLinExpr z(0);
		if (min_power)
		{
			//FIXME: this does not work well when PM switch-on/off costs and VM migration costs are != zero!
			std::cerr << "(W) Power optimization does not work well when PM switch-on/off costs and VM migration costs are not zero!" << std::endl;

			// z = \min sum_{h \in H}{x_h C_{g(h)}^{min} + (C_{g(h)}^{max}-C_{g(h)}^{min})s_{h} + x_h(1-o(i)) + (1-x_h)o(i)}
			for (std::size_t h = 0; h < npms; ++h)
			{
				RealT dC = pm_specs_max_power[pms_category[h]]-pm_specs_min_power[pms_category[h]];

				z += x[h]*pm_specs_min_power[pms_category[h]] + dC*s[h];
			}
		}
		else
		{
			for (std::size_t h = 0; h < npms; ++h)
			{
				// Add PM power costs due to computing demand
				RealT dC = pm_specs_max_power[pms_category[h]]-pm_specs_min_power[pms_category[h]];
				RealT wcost = cips_electricity_cost[pms_cip[h]]*1e-3; // Electricity cost in Wh
				z += (x[h]*pm_specs_min_power[pms_category[h]]+dC*s[h])*wcost;
				// Add PM switch-on/off costs
				z += x[h]*(1-pm_power_states[h])*cip_pm_awake_costs[pms_cip[h]][pms_category[h]]
				  +  (1-x[h])*pm_power_states[h]*cip_pm_asleep_costs[pms_cip[h]][pms_category[h]];
				// Add VM migration costs
				for (std::size_t v = 0; v < nvms; ++v)
				{
					z += y[v][h]*cip_to_cip_vm_migration_costs[vms_cip[v]][pms_cip[h]][vms_category[v]];
				}
			}
		}
		model.setObjective(z, GRB_MINIMIZE);
		model.update();

#ifdef DCS_DEBUG
		model.write("gurobi-model.lp");
#endif // DCS_DEBUG
 
		model.optimize();

		solution.solved = false;
		solution.optimal = false;

		int status = model.get(GRB_IntAttr_Status);
		switch (status)
		{
			case GRB_OPTIMAL: // The algorithm found an optimal solution.
				solution.objective_value = static_cast<RealT>(model.get(GRB_DoubleAttr_ObjVal));
				solution.solved = true;
				solution.optimal = true;
				break;
			case GRB_SUBOPTIMAL: // The algorithm found a feasible solution, though it may not necessarily be optimal.
				solution.objective_value = static_cast<RealT>(model.get(GRB_DoubleAttr_ObjVal));
				std::clog << "(W) Optimization problem solved but non-optimal" << std::endl;
				solution.solved = true;
				break;
			default:
			{
				std::ostringstream oss;
				std::clog << "Optimization was stopped with status = " << status << std::endl;
				return solution;
			}
		}

#ifdef DCS_DEBUG
		DCS_DEBUG_TRACE( "Optimal solution: " );

		DCS_DEBUG_TRACE( "- Solved: " << std::boolalpha << solution.solved );
		DCS_DEBUG_TRACE( "- Optimal: " << std::boolalpha << solution.optimal );

		DCS_DEBUG_TRACE( "- Decision variables: " );

		// Output x_{h}
		for (std::size_t h = 0; h < npms; ++h)
		{
			DCS_DEBUG_STREAM << x[h].get(GRB_StringAttr_VarName) << " = " << x[h].get(GRB_DoubleAttr_X) << " (" << static_cast<bool>(x[h].get(GRB_DoubleAttr_X)) << ")" << ::std::endl;
		}

		// Output y_{vh}
		for (std::size_t v = 0; v < nvms; ++v)
		{
			for (std::size_t h = 0; h < npms; ++h)
			{
				DCS_DEBUG_STREAM << y[v][h].get(GRB_StringAttr_VarName) << " = " << y[v][h].get(GRB_DoubleAttr_X) << " (" << static_cast<bool>(y[v][h].get(GRB_DoubleAttr_X)) << ")" << ::std::endl;
			}
		}

		// Output s_{h}
		for (std::size_t h = 0; h < npms; ++h)
		{
			DCS_DEBUG_STREAM << s[h].get(GRB_StringAttr_VarName) << " = " << s[h].get(GRB_DoubleAttr_X) << std::endl;
		}

		//Print z
		DCS_DEBUG_TRACE( "- Objective value: " << solution.objective_value );
#endif // DCS_DEBUG

		solution.pm_power_states.resize(npms);
		solution.pm_vm_allocations.resize(npms);
		if (min_power)
		{
			// Computed in the loop below
			solution.cost = 0;
		}
		else
		{
			solution.cost = solution.objective_value;
		}
		for (std::size_t h = 0; h < npms; ++h)
		{
			solution.pm_power_states[h] = static_cast<bool>(x[h].get(GRB_DoubleAttr_X));

			// Compute the energy cost
			if (min_power && static_cast<bool>(x[h].get(GRB_DoubleAttr_X)))
			{
				//RealT dC = pm_specs_max_power[pms_category[h]]-pm_specs_min_power[pms_category[h]];
				RealT wcost = cips_electricity_cost[pms_cip[h]]*1e-3; // Electricity cost in Wh

				//solution.cost += (pm_specs_min_power[pms_category[h]] + dC*static_cast<RealT>(s[h].get(GRB_DoubleAttr_X)))*wcost;
				solution.cost += pm_consumed_power(pm_specs_min_power[pms_category[h]], pm_specs_max_power[pms_category[h]], static_cast<RealT>(s[h].get(GRB_DoubleAttr_X)))*wcost;
			}
			solution.pm_vm_allocations[h].resize(nvms);
			for (std::size_t v = 0; v < nvms; ++v)
			{
				solution.pm_vm_allocations[h][v] = static_cast<bool>(y[v][h].get(GRB_DoubleAttr_X));
			}
		}
	}
	catch (GRBException const& e)
	{
		::std::ostringstream oss;
		oss << "Got exception from GUROBI: " << e.getMessage() << " (" << e.getErrorCode() << ")";
		DCS_EXCEPTION_THROW(::std::runtime_error, oss.str());
	}
	catch (...)
	{
		DCS_EXCEPTION_THROW(::std::runtime_error,
							"Unexpected error during the optimization");
	}
#else
# error Unable to find a suitable solver
#endif // DCS_CLOUD_GT_HAVE_*_SOLVER

	return solution;
}


template <typename RealT>
struct merge_split_stable_partition_selector
{
	::std::vector< partition_info<RealT> > operator()(::gtpack::cooperative_game<RealT> const& game, ::std::map< gtpack::cid_type, coalition_info<RealT> > const& visited_coalitions)
	{
		namespace alg = ::dcs::algorithm;

		// Generate all partitions and select the ones that are D_{hp}-stable, that is that are stable according to merge/split operations

		::std::vector< partition_info<RealT> > best_partitions;

		const ::std::vector<gtpack::player_type> players(game.players());
		const ::std::size_t np(players.size());

		alg::lexicographic_partition partition(np);

		while (partition.has_next())
		{
			typedef typename alg::partition_traits<gtpack::player_type>::subset_container subset_container;
			typedef typename alg::partition_traits<gtpack::player_type>::subset_const_iterator subset_iterator;

			subset_container subs;

			// Each subset is a collection of coalitions
			subs = alg::next_partition(players.begin(), players.end(), partition);

			DCS_DEBUG_TRACE("--- PARTITION: " << partition);//XXX

			partition_info<RealT> candidate_partition;

			bool Dhp_stable(true);

			std::vector<gtpack::cid_type> P;

			subset_iterator sub_end_it(subs.end());
			for (subset_iterator sub_it = subs.begin();
				 sub_it != sub_end_it && Dhp_stable;
				 ++sub_it)
			{
				const gtpack::cid_type cid = gtpack::players_coalition<RealT>::make_id(sub_it->begin(), sub_it->end());

				P.push_back(cid);

				if (visited_coalitions.count(cid) == 0)
				{
					continue;
				}

				DCS_DEBUG_TRACE("--- COALITION: " << game.coalition(cid) << " (CID=" << cid << ")");//XXX

				candidate_partition.coalitions.insert(cid);

				// Check that v(P_i) \ge \sum_{j=1}^l v(C_j), \forall partition C=\{C_1,\ldots,C_l\} of P_i

				RealT vPi = visited_coalitions.at(cid).value;

				::std::vector<gtpack::player_type> coal_players(sub_it->begin(), sub_it->end());

				alg::lexicographic_partition sub_partition(coal_players.size());
				subset_container sub_subs;
				sub_subs = alg::next_partition(coal_players.begin(), coal_players.end(), sub_partition);
				subset_iterator sub_sub_end_it(sub_subs.end());
				RealT svC = 0;
				for (subset_iterator sub_sub_it = sub_subs.begin();
					 sub_sub_it != sub_sub_end_it;
					 ++sub_sub_it)
				{
					const gtpack::cid_type sub_cid = gtpack::players_coalition<RealT>::make_id(sub_sub_it->begin(), sub_sub_it->end());

					if (visited_coalitions.count(sub_cid) == 0)
					{
						continue;
					}

					svC += visited_coalitions.at(sub_cid).value;
				}

				if (dcs::math::float_traits<RealT>::definitely_less(vPi, svC))
				{
					Dhp_stable = false;
					break;
				}

				for (::std::size_t p = 0; p < coal_players.size(); ++p)
				{
					const gtpack::player_type pid = coal_players[p];

					if (visited_coalitions.at(cid).payoffs.count(pid) > 0)
					{
						candidate_partition.payoffs[pid] = visited_coalitions.at(cid).payoffs.at(pid);
					}
					else
					{
						candidate_partition.payoffs[pid] = ::std::numeric_limits<RealT>::quiet_NaN();
					}
				}
			}

			if (!Dhp_stable)
			{
				continue;
			}

			alg::lexicographic_subset P_subset(P.size(), false);

			while (P_subset.has_next() && Dhp_stable)
			{
				typedef typename alg::subset_traits<gtpack::cid_type>::element_container cid_container;

				const cid_container sub_P = alg::next_subset(P.begin(), P.end(), P_subset);

				RealT svPi = 0;
				std::set<gtpack::player_type> UPi_players;
				for (std::size_t i = 0; i < sub_P.size(); ++i)
				{
					const gtpack::cid_type Pi_cid = sub_P[i];

					if (visited_coalitions.count(Pi_cid) == 0)
					{
						continue;
					}

					svPi += visited_coalitions.at(Pi_cid).value;

					const std::vector<gtpack::player_type> Pi_players = gtpack::players_coalition<RealT>(np, Pi_cid).players();
					UPi_players.insert(Pi_players.begin(), Pi_players.end());
				}

				const gtpack::cid_type UPi_cid = gtpack::players_coalition<RealT>(UPi_players.begin(), UPi_players.end()).id();
				if (visited_coalitions.count(UPi_cid) == 0)
				{
					continue;
				}

				const RealT vUPi = visited_coalitions.at(UPi_cid).value;

				if (dcs::math::float_traits<RealT>::definitely_less(svPi, vUPi))
				{
					Dhp_stable = false;
					break;
				}
			}

			if (Dhp_stable)
			{
				best_partitions.push_back(candidate_partition);
			}
		}

		return best_partitions;
	}
}; // nash_stable_partition_selector


template <typename RealT>
struct nash_stable_partition_selector
{
	::std::vector< partition_info<RealT> > operator()(::gtpack::cooperative_game<RealT> const& game, ::std::map< gtpack::cid_type, coalition_info<RealT> > const& visited_coalitions)
	{
		namespace alg = ::dcs::algorithm;

		// Generate all partitions and select the ones that are Nash-stable

		::std::vector< partition_info<RealT> > best_partitions;

		const ::std::vector<gtpack::player_type> players(game.players());
		const ::std::size_t np(players.size());

		alg::lexicographic_partition partition(np);

		while (partition.has_next())
		{
			typedef typename alg::partition_traits<gtpack::player_type>::subset_container subset_container;
			typedef typename alg::partition_traits<gtpack::player_type>::subset_const_iterator subset_iterator;

			subset_container subs;

			// Each subset is a collection of coalitions
			subs = alg::next_partition(players.begin(), players.end(), partition);

			DCS_DEBUG_TRACE("--- PARTITION: " << partition);//XXX

			partition_info<RealT> candidate_partition;

			subset_iterator sub_end_it(subs.end());
			for (subset_iterator sub_it = subs.begin();
				 sub_it != sub_end_it;
				 ++sub_it)
			{
				const gtpack::cid_type cid = gtpack::players_coalition<RealT>::make_id(sub_it->begin(), sub_it->end());

				if (visited_coalitions.count(cid) == 0)
				{
					continue;
				}

				DCS_DEBUG_TRACE("--- COALITION: " << game.coalition(cid) << " (CID=" << cid << ")");//XXX

				candidate_partition.coalitions.insert(cid);

				::std::vector<gtpack::player_type> coal_players(sub_it->begin(), sub_it->end());
				for (::std::size_t p = 0; p < coal_players.size(); ++p)
				{
					const gtpack::player_type pid = coal_players[p];

					if (visited_coalitions.at(cid).payoffs.count(pid) > 0)
					{
						candidate_partition.payoffs[pid] = visited_coalitions.at(cid).payoffs.at(pid);
					}
					else
					{
						candidate_partition.payoffs[pid] = ::std::numeric_limits<RealT>::quiet_NaN();
					}
				}
			}

			// Check Nash-stability

			bool nash_stable(true);

			// For all players $p$
			for (::std::size_t p = 0; p < np && nash_stable; ++p)
			{
				const gtpack::player_type pid(players[p]);

				// For all $S_k$ \in \Pi \cup \{\emptyset\}$
				bool found_singleton(false);
				subset_iterator sub_end_it(subs.end());
				for (subset_iterator sub_it = subs.begin();
					 sub_it != sub_end_it && nash_stable;
					 ++sub_it)
				{
					::std::set<gtpack::player_type> coal_players(sub_it->begin(), sub_it->end());

					if (coal_players.count(pid) == 0)
					{
						// This coalition doesn't include player pid, go on

						// Evaluate $S_k \cup {p}
						coal_players.insert(pid);

						const gtpack::cid_type cid = gtpack::players_coalition<RealT>::make_id(coal_players.begin(), coal_players.end());

						DCS_DEBUG_TRACE("--- PID: " << pid << " - AUGMENTED COALITION: " << game.coalition(cid) << " (CID=" << cid << ") - AUGMENTED PAYOFF: " << (visited_coalitions.at(cid).payoffs.count(pid) ? visited_coalitions.at(cid).payoffs.at(pid) : ::std::numeric_limits<RealT>::quiet_NaN()) << " - CANDIDATE PAYOFF: " << candidate_partition.payoffs.at(pid));///XXX

						// Check player's preference
						if (visited_coalitions.at(cid).payoffs.count(pid) == 0
							|| ::dcs::math::float_traits<RealT>::definitely_greater(visited_coalitions.at(cid).payoffs.at(pid), candidate_partition.payoffs.at(pid)))
						{
							DCS_DEBUG_TRACE("--- PID: " << pid << " - AUGMENTED COALITION: " << game.coalition(cid) << " (CID=" << cid << "): NOT NASH STABLE");//XXX
							nash_stable = false;
							break;
						}
					}
					else if (coal_players.size() == 1)
					{
						found_singleton = true;
					}
				}

				// Check singleton coalition
				if (!found_singleton)
				{
					const gtpack::cid_type cid = gtpack::players_coalition<RealT>::make_id(&(players[p]), &(players[p])+1);

					//DCS_DEBUG_TRACE("--- PID: " << pid << " - AUGMENTED COALITION: " << game.coalition(cid) << " (CID=" << cid << ")");//XXX
					DCS_DEBUG_TRACE("--- PID: " << pid << " - AUGMENTED COALITION: " << game.coalition(cid) << " (CID=" << cid << ") - AUGMENTED PAYOFF: " << (visited_coalitions.count(cid) && visited_coalitions.at(cid).payoffs.count(pid) ? visited_coalitions.at(cid).payoffs.at(pid) : ::std::numeric_limits<RealT>::quiet_NaN()) << " - CANDIDATE PAYOFF: " << candidate_partition.payoffs.at(pid));///XXX

					if (candidate_partition.coalitions.count(cid) == 0)
					{
						// This partition doesn't contain this singleton coalition
						if (visited_coalitions.at(cid).payoffs.count(pid) == 0
							|| ::dcs::math::float_traits<RealT>::definitely_greater(visited_coalitions.at(cid).payoffs.at(pid), candidate_partition.payoffs.at(pid)))
						{
							DCS_DEBUG_TRACE("--- PID: " << pid << " - AUGMENTED COALITION: " << game.coalition(cid) << " (CID=" << cid << "): NOT NASH STABLE");//XXX
							nash_stable = false;
							break;
						}
					}
				}
			}

			if (nash_stable)
			{
				best_partitions.push_back(candidate_partition);
			}
		}

		return best_partitions;
	}
}; // nash_stable_partition_selector


template <typename RealT>
struct pareto_optimal_partition_selector
{
	::std::vector< partition_info<RealT> > operator()(::gtpack::cooperative_game<RealT> const& game, ::std::map< gtpack::cid_type, coalition_info<RealT> > const& visited_coalitions)
	{
		namespace alg = ::dcs::algorithm;

		// Generate all partitions and select the ones that are Pareto optimal

		::std::vector< partition_info<RealT> > best_partitions;

		const ::std::vector<gtpack::player_type> players(game.players());
		const ::std::size_t np(players.size());

		alg::lexicographic_partition partition(np);

		::std::vector<RealT> best_payoffs(np, ::std::numeric_limits<RealT>::quiet_NaN());

		while (partition.has_next())
		{
			typedef typename alg::partition_traits<gtpack::player_type>::subset_container subset_container;
			typedef typename alg::partition_traits<gtpack::player_type>::subset_const_iterator subset_iterator;

			subset_container subs;

			// Each subset is a collection of coalitions
			subs = alg::next_partition(players.begin(), players.end(), partition);

			DCS_DEBUG_TRACE("--- PARTITION: " << partition);//XXX

			partition_info<RealT> candidate_partition;

			subset_iterator sub_end_it(subs.end());
			for (subset_iterator sub_it = subs.begin();
				 sub_it != sub_end_it;
				 ++sub_it)
			{
				const gtpack::cid_type cid = gtpack::players_coalition<RealT>::make_id(sub_it->begin(), sub_it->end());

				if (visited_coalitions.count(cid) == 0)
				{
					continue;
				}

				DCS_DEBUG_TRACE("--- COALITION: " << game.coalition(cid) << " (CID=" << cid << ")");//XXX

				candidate_partition.coalitions.insert(cid);

				::std::vector<gtpack::player_type> coal_players(sub_it->begin(), sub_it->end());
				for (::std::size_t p = 0; p < coal_players.size(); ++p)
				{
					const gtpack::player_type pid = coal_players[p];

					if (visited_coalitions.at(cid).payoffs.count(pid) > 0)
					{
						candidate_partition.payoffs[pid] = visited_coalitions.at(cid).payoffs.at(pid);
					}
					else
					{
						candidate_partition.payoffs[pid] = ::std::numeric_limits<RealT>::quiet_NaN();
					}
				}
			}

			// Check Pareto optimality

			bool pareto_optimal(true);

			// For all players $p$
			for (std::size_t p = 0; p < np && pareto_optimal; ++p)
			{
				const gtpack::player_type pid(players[p]);

				if (::std::isnan(best_payoffs[p]) || candidate_partition.payoffs.at(pid) > best_payoffs[p])
				{
					best_payoffs[p] = candidate_partition.payoffs.at(pid);
				}
				else
				{
					pareto_optimal = false;
				}
			}

			if (pareto_optimal)
			{
				best_partitions.push_back(candidate_partition);
			}
		}

		return best_partitions;
	}
}; // pareto_optimal_partition_selector


template <typename RealT>
struct social_optimum_partition_selector
{
	::std::vector< partition_info<RealT> > operator()(::gtpack::cooperative_game<RealT> const& game, ::std::map< gtpack::cid_type, coalition_info<RealT> > const& visited_coalitions)
	{
		namespace alg = ::dcs::algorithm;

		// Generate all partitions and select the ones that maximize the social welfare

		::std::vector< partition_info<RealT> > best_partitions;

		const ::std::vector<gtpack::player_type> players(game.players());
		const ::std::size_t np(players.size());
		RealT best_value(0);

		alg::lexicographic_partition partition(np);

		while (partition.has_next())
		{
			typedef typename alg::partition_traits<gtpack::player_type>::subset_container subset_container;
			typedef typename alg::partition_traits<gtpack::player_type>::subset_const_iterator subset_iterator;

			subset_container subs;

			// Each subset is a collection of coalitions
			subs = alg::next_partition(players.begin(), players.end(), partition);

			DCS_DEBUG_TRACE("--- PARTITION: " << partition);//XXX

			partition_info<RealT> candidate_partition;
			RealT candidate_partition_value(0);

			subset_iterator sub_end_it(subs.end());
			for (subset_iterator sub_it = subs.begin();
				 sub_it != sub_end_it;
				 ++sub_it)
			{
				const gtpack::cid_type cid = gtpack::players_coalition<RealT>::make_id(sub_it->begin(), sub_it->end());

				if (visited_coalitions.count(cid) == 0)
				{
					continue;
				}

				DCS_DEBUG_TRACE("--- COALITION: " << game.coalition(cid) << " (CID=" << cid << ")");//XXX

				candidate_partition.coalitions.insert(cid);

				::std::vector<gtpack::player_type> coal_players(sub_it->begin(), sub_it->end());
				for (::std::size_t p = 0; p < coal_players.size(); ++p)
				{
					const gtpack::player_type pid = coal_players[p];

					if (visited_coalitions.at(cid).payoffs.count(pid) > 0)
					{
						candidate_partition.payoffs[pid] = visited_coalitions.at(cid).payoffs.at(pid);
					}
					else
					{
						candidate_partition.payoffs[pid] = ::std::numeric_limits<RealT>::quiet_NaN();
					}
				}
				candidate_partition_value += visited_coalitions.at(cid).value;
			}

			// Check for social optimum

			if (best_partitions.size() == 0
				|| ::dcs::math::float_traits<RealT>::definitely_greater(candidate_partition_value, best_value))
			{
				best_partitions.clear();
				best_partitions.push_back(candidate_partition);
				best_value = candidate_partition_value;
			}
			else if (::dcs::math::float_traits<RealT>::essentially_equal(candidate_partition_value, best_value))
			{
				best_partitions.push_back(candidate_partition);
			}
		}

		return best_partitions;
	}
}; // social_optimum_partition_selector


template <typename RealT>
coalition_formation_info<RealT> analyze_coalitions(scenario<RealT> const& s, options<RealT> const& opts)
{
	//typedef RealT real_type;

	//const std::size_t ncoals(std::pow(2, s.num_cips)-1);

	std::vector<std::size_t> cips(s.num_cips);
	std::size_t num_pms = 0;
//	std::size_t num_vms = 0;

	gtpack::cooperative_game<RealT> game(s.num_cips, boost::make_shared< gtpack::explicit_characteristic_function<RealT> >());

	for (std::size_t k = 0; k < s.num_cips; ++k)
	{
		cips[k] = k;

		for (std::size_t p = 0; p < s.num_pm_types; ++p)
		{
			num_pms += s.cip_num_pms[k][p];
		}

//		for (std::size_t v = 0; v < s.num_vm_types; ++v)
//		{
//			num_vms += s.cip_num_vms[k][v];
//		}
	}

	std::map< gtpack::cid_type, coalition_info<RealT> > visited_coalitions;

	alg::lexicographic_subset subset(s.num_cips, false);

	while (subset.has_next())
	{
		typedef typename alg::subset_traits<std::size_t>::element_container element_container;

		DCS_DEBUG_TRACE("--- SUBSET: " << subset);//XXX

		const element_container coal_cips = alg::next_subset(cips.begin(), cips.end(), subset);
		const std::size_t coal_ncips(coal_cips.size());

		const gtpack::cid_type cid = gtpack::players_coalition<RealT>::make_id(coal_cips.begin(), coal_cips.end());

		DCS_DEBUG_TRACE("--- COALITION: " << game.coalition(cid) << " (CID=" << cid << ")");//XXX

		std::size_t coal_npms(0);
		std::size_t coal_nvms(0);

		for (std::size_t i = 0; i < coal_ncips; ++i)
		{
			std::size_t cip(coal_cips[i]);

			for (std::size_t p = 0; p < s.num_pm_types; ++p)
			{
				coal_npms += s.cip_num_pms[cip][p];
			}
			for (std::size_t v = 0; v < s.num_vm_types; ++v)
			{
				coal_nvms += s.cip_num_vms[cip][v];
			}
		}

		std::size_t pms_start(0);
		std::size_t pms_stop(0);
		std::size_t vms_start(0);
		std::size_t vms_stop(0);
		std::vector<std::size_t> coal_pms_cips(coal_npms);
		std::vector<std::size_t> coal_vms_cips(coal_nvms);
		std::vector<std::size_t> coal_pms_category(coal_npms);
		std::vector<std::size_t> coal_vms_category(coal_nvms);
		std::vector<bool> coal_pm_power_states(coal_npms);
		RealT coal_profit(0);

		for (std::size_t i = 0; i < coal_ncips; ++i)
		{
			std::size_t cip(coal_cips[i]);

			for (std::size_t p = 0; p < s.num_pm_types; ++p)
			{
				if (s.cip_num_pms[cip][p] > 0)
				{
					pms_stop += s.cip_num_pms[cip][p];
//DCS_DEBUG_TRACE("FILLING PMS -- CID: " << cid << " - CIP: " << cip << " - START: " << pms_start << " - STOP: " << pms_stop);
					std::fill(coal_pms_category.begin()+pms_start, coal_pms_category.begin()+pms_stop, p);
					std::fill(coal_pms_cips.begin()+pms_start, coal_pms_cips.begin()+pms_stop, cip);
					for (std::size_t k = 0; k < (pms_stop-pms_start); ++k)
					{
						coal_pm_power_states[k+pms_start] = s.cip_pm_power_states[cip][k];
					}
					pms_start = pms_stop;
				}
			}

			for (std::size_t v = 0; v < s.num_vm_types; ++v)
			{
				coal_profit += s.cip_revenues[cip][v]*s.cip_num_vms[cip][v];

				if (s.cip_num_vms[cip][v] > 0)
				{
					vms_stop += s.cip_num_vms[cip][v];
//DCS_DEBUG_TRACE("FILLING VMS -- CID: " << cid << " - CIP: " << cip << " - START: " << vms_start << " - STOP: " << vms_stop);
					std::fill(coal_vms_category.begin()+vms_start, coal_vms_category.begin()+vms_stop, v);
					std::fill(coal_vms_cips.begin()+vms_start, coal_vms_cips.begin()+vms_stop, cip);
					vms_start = vms_stop;
				}
			}
		}

		optimal_allocation_info<RealT> optimal_allocation;
		optimal_allocation = find_optimal_allocation(coal_ncips,
													 s.cip_electricity_costs,
													 coal_npms,
													 coal_pms_cips,
													 coal_pms_category,
													 s.pm_spec_min_powers,
													 s.pm_spec_max_powers,
													 coal_nvms,
													 coal_vms_cips,
													 coal_vms_category,
													 s.vm_spec_cpus,
													 s.vm_spec_rams,
													 coal_pm_power_states,
													 s.cip_pm_asleep_costs,
													 s.cip_pm_awake_costs,
													 s.cip_to_cip_vm_migration_costs,
													 false,
													 opts.opt_relative_gap,
													 opts.opt_time_lim);

		visited_coalitions[cid].optimal_allocation = optimal_allocation;
		if (optimal_allocation.solved)
		{
			game.value(cid, coal_profit-optimal_allocation.cost);
			visited_coalitions[cid].value = game.value(cid);

			DCS_DEBUG_TRACE( "CID: " << cid << " - Profit: " << coal_profit << " - Cost: " << optimal_allocation.cost << " => v(CID)=" << game.value(cid) );

//#ifdef DCS_DEBUG
			std::map< std::size_t, cip_allocation_info<RealT> > coal_cips_info;

			RealT tot_cost(0);
			RealT tot_watt(0);
			for (std::size_t c = 0; c < coal_ncips; ++c)
			{
				std::size_t cip(coal_cips[c]);

				coal_cips_info[cip].num_on_pms = 0;
				coal_cips_info[cip].num_vms = 0;
				coal_cips_info[cip].tot_watt = 0;
				//coal_cips_info[cip].tot_wcost = 0;
			}
			for (std::size_t p = 0; p < coal_npms; ++p)
			{
				if (optimal_allocation.pm_power_states[p])
				{
					std::size_t cip(coal_pms_cips[p]);
					std::size_t pc(coal_pms_category[p]);

					coal_cips_info[cip].num_on_pms += 1;

					RealT tot_share(0);
					for (std::size_t v = 0; v < coal_nvms; ++v)
					{
						if (optimal_allocation.pm_vm_allocations[p][v])
						{
							coal_cips_info[cip].num_vms += 1;
							tot_share += s.vm_spec_cpus[coal_vms_category[v]][pc];
						}
					}
					coal_cips_info[cip].tot_watt += pm_consumed_power(s.pm_spec_min_powers[pc], s.pm_spec_max_powers[pc], tot_share);
				}
			}
			for (std::size_t c = 0; c < coal_ncips; ++c)
			{
				std::size_t cip(coal_cips[c]);

				DCS_DEBUG_TRACE( "CID: " << cid << " - CIP: " << cip << " - # Powered-on PMs: " << coal_cips_info[cip].num_on_pms << " - # Hosted VMs: " << coal_cips_info[cip].num_vms << " - Consumed Watts: " << coal_cips_info[cip].tot_watt << " - Energy Cost: " << (coal_cips_info[cip].tot_watt*1e-3*s.cip_electricity_costs[cip]));

				tot_watt += coal_cips_info[cip].tot_watt*1e-3;
				tot_cost += coal_cips_info[cip].tot_watt*1e-3*s.cip_electricity_costs[cip];
			}
//#endif // DCS_DEBUG

			//visited_coalitions[cid].optimal_allocation.cost = tot_cost;
			visited_coalitions[cid].optimal_allocation.kwatt = tot_watt;

			// Check core existence
			gtpack::cooperative_game<RealT> subgame = game.subgame(coal_cips.begin(), coal_cips.end());
			gtpack::core<RealT> core = gtpack::find_core(subgame);

			if (core.empty())
			{
				DCS_DEBUG_TRACE( "CID: " << cid << " - The core is empty" );

				visited_coalitions[cid].core_empty = true;
				visited_coalitions[cid].payoffs_in_core = false;
				//skip_partition = true;

				if (subgame.num_players() == cips.size())
				{
					// This is the Grand coalition

					DCS_DEBUG_TRACE( "CID: " << cid << " - The Grand-Coalition has an empty core" );
				}
			}
			else
			{
				DCS_DEBUG_TRACE( "CID: " << cid << " - The core is not empty" );

				visited_coalitions[cid].core_empty = false;
			}

			// Compute the coalition value
			std::map<gtpack::player_type,RealT> coal_payoffs;
			switch (opts.coalition_value_division)
			{
				case banzhaf_coalition_value_division:
					coal_payoffs = gtpack::banzhaf_value(subgame);
					break;
				case normalized_banzhaf_coalition_value_division:
					coal_payoffs = gtpack::norm_banzhaf_value(subgame);
					break;
				case shapley_coalition_value_division:
					coal_payoffs = gtpack::shapley_value(subgame);
					break;
			}

#ifdef DCS_DEBUG
			for (std::size_t c = 0; c < coal_ncips; ++c)
			{
				std::size_t cip(coal_cips[c]);

				DCS_DEBUG_TRACE( "CID: " << cid << " - CIP: " << cip << " - Coalition payoff: " << coal_payoffs[cip] );
			}
#endif // DCS_DEBUG

			visited_coalitions[cid].payoffs = coal_payoffs;

			// Check if the value is in the core (if the core != empty)
			if (!visited_coalitions.at(cid).core_empty)
			{
				if (gtpack::belongs_to_core(game.subgame(coal_cips.begin(), coal_cips.end()), coal_payoffs.begin(), coal_payoffs.end()))
				{
					DCS_DEBUG_TRACE( "CID: " << cid << " - The Coalition value belongs to the core" );

					visited_coalitions[cid].payoffs_in_core = true;
				}
				else
				{
					DCS_DEBUG_TRACE( "CID: " << cid << " - The Coaition value does not belong to the core" );

					visited_coalitions[cid].payoffs_in_core = false;
				}
			}
		}
		else
		{
			DCS_DEBUG_TRACE( "CID: " << cid << " - The allocation problem is infeasible" );

			visited_coalitions[cid].core_empty = true;
			visited_coalitions[cid].payoffs_in_core = false;

			//skip_partition = true;
			game.value(cid, -::std::numeric_limits<RealT>::min());

			if (game.coalition(cid).num_players() == cips.size())
			{
				// This is the Grand coalition

				DCS_DEBUG_TRACE( "CID: " << cid << " - The Grand-Coalition has an infeasible solution and thus an empty core" );
			}
		}
	}

	coalition_formation_info<RealT> formed_coalitions;

	formed_coalitions.coalitions = visited_coalitions;
	switch (opts.coalition_formation)
	{
		case merge_split_stable_coalition_formation:
			formed_coalitions.best_partitions = merge_split_stable_partition_selector<RealT>()(game, visited_coalitions);
			break;
		case nash_stable_coalition_formation:
			formed_coalitions.best_partitions = nash_stable_partition_selector<RealT>()(game, visited_coalitions);
			break;
		case pareto_optimal_coalition_formation:
			formed_coalitions.best_partitions = pareto_optimal_partition_selector<RealT>()(game, visited_coalitions);
			break;
		case social_optimum_coalition_formation:
			formed_coalitions.best_partitions = social_optimum_partition_selector<RealT>()(game, visited_coalitions);
			break;
	}

/*
#ifdef DCS_DEBUG
	DCS_DEBUG_TRACE( "FORMED PARTITIONS: ");
	for (std::size_t i = 0; i < formed_coalitions.best_partitions.size(); ++i)
	{
		const partition_info<RealT> part = formed_coalitions.best_partitions[i];

		typedef typename std::set<gtpack::cid_type>::const_iterator coalition_iterator;
		coalition_iterator coal_end_it(part.coalitions.end());
		DCS_DEBUG_STREAM << "  [";
		for (coalition_iterator coal_it = part.coalitions.begin();
			 coal_it != coal_end_it;
			 ++coal_it)
		{
			const gtpack::cid_type cid(*coal_it);

			DCS_DEBUG_STREAM << cid << ",";
		}
		DCS_DEBUG_STREAM << "]" << std::endl;
	}
#endif // DCS_DEBUG
*/

	return formed_coalitions;
}

template <typename RealT>
void report(std::size_t ncips, coalition_formation_info<RealT> const& formed_coalitions)
{
	typedef typename std::map<gtpack::player_type,RealT>::const_iterator coal_val_iterator;

	std::vector<gtpack::player_type> players(ncips);
	for (std::size_t c = 0; c < ncips; ++c)
	{
		players[c] = static_cast<gtpack::player_type>(c);
	}

	// Retrieve the ID of the grand-coalition
	gtpack::cid_type gcid = gtpack::players_coalition<RealT>::make_id(players.begin(), players.end());

	std::cout << "################################################################################" << std::endl;
	std::cout << "### Report on Formed Coalitions:" << std::endl;
	std::cout << "################################################################################" << std::endl;

	std::cout << "- Best Partitions:" << std::endl;
	if (formed_coalitions.best_partitions.size() > 0)
	{
		for (std::size_t i = 0; i < formed_coalitions.best_partitions.size(); ++i)
		{
	//		const bs::partition_info<RealT> part = formed_coalitions.best_partitions[i];
	//
	//		typedef typename std::set<gtpack::cid_type>::const_iterator coalition_iterator;
	//		coalition_iterator coal_end_it(part.coalitions.end());
	//		DCS_DEBUG_STREAM << "  [";
	//		for (coalition_iterator coal_it = part.coalitions.begin();
	//			 coal_it != coal_end_it;
	//			 ++coal_it)
	//		{
	//			const gtpack::cid_type cid(*coal_it);
	//
	//			DCS_DEBUG_STREAM << cid << ",";
	//		}
	//		DCS_DEBUG_STREAM << "]" << std::endl;

			typedef typename std::set<gtpack::cid_type>::const_iterator cid_iterator;

			RealT bestpart_value(0);
			RealT grandpart_value(0);
			RealT singlepart_value(0);
			RealT bestpart_kwatt(0);
			//RealT grandpart_kwatt(0);
			RealT singlepart_kwatt(0);

			cid_iterator cid_beg_it(formed_coalitions.best_partitions[i].coalitions.begin());
			cid_iterator cid_end_it(formed_coalitions.best_partitions[i].coalitions.end());

			std::cout << " * Payoffs: {";
			for (cid_iterator cid_it = cid_beg_it; cid_it != cid_end_it; ++cid_it)
			{
				//typedef typename std::map<gtpack::player_type,RealT>::const_iterator coal_val_iterator;

				gtpack::cid_type cid(*cid_it);

				if (cid_it != cid_beg_it)
				{
					std::cout << ", ";
				}

				std::cout << "{";
				coal_val_iterator coal_val_end_it(formed_coalitions.coalitions.at(cid).payoffs.end());
				coal_val_iterator coal_val_beg_it(formed_coalitions.coalitions.at(cid).payoffs.begin());
				for (coal_val_iterator coal_val_it = coal_val_beg_it;
					 coal_val_it != coal_val_end_it;
					 ++coal_val_it)
				{
					gtpack::player_type pid(coal_val_it->first);
					RealT value(coal_val_it->second);

					if (coal_val_it != coal_val_beg_it)
					{
						std::cout << ", ";
					}

					std::cout << pid << " => " << value;
					bestpart_value += value;
				}
				std::cout << "}";

				bestpart_kwatt += formed_coalitions.coalitions.at(cid).optimal_allocation.kwatt;
			}
			std::cout << "}" << std::endl;

			std::cout << " * Value: " << bestpart_value << std::endl;
			std::cout << " * Energy Consumption: " << bestpart_kwatt << std::endl;

//			std::cout << " * Side Payments: {";
//			for (cid_iterator cid_it = cid_beg_it; cid_it != cid_end_it; ++cid_it)
//			{
//				gtpack::cid_type cid(*cid_it);
//
//				if (cid_it != cid_beg_it)
//				{
//					std::cout << ", ";
//				}
//
//				std::cout << "{";
//				coal_val_iterator coal_val_end_it(formed_coalitions.coalitions.at(cid).payoffs.end());
//				coal_val_iterator coal_val_beg_it(formed_coalitions.coalitions.at(cid).payoffs.begin());
//				for (coal_val_iterator coal_val_it = coal_val_beg_it;
//					 coal_val_it != coal_val_end_it;
//					 ++coal_val_it)
//				{
//					gtpack::player_type pid(coal_val_it->first);
//
//					if (coal_val_it != coal_val_beg_it)
//					{
//						std::cout << ", ";
//					}
//
//					std::cout << pid << " => " << formed_coalitions.best_partitions[i].side_payments.at(pid);
//				}
//				std::cout << "}";
//			}
//			std::cout << "}" << std::endl;

			std::cout << " * Core exists?: {";
			for (cid_iterator cid_it = cid_beg_it; cid_it != cid_end_it; ++cid_it)
			{
				gtpack::cid_type cid(*cid_it);

				if (cid_it != cid_beg_it)
				{
					std::cout << ", ";
				}

				std::cout << std::boolalpha << !(formed_coalitions.coalitions.at(cid).core_empty);
			}
			std::cout << "}" << std::endl;

			std::cout << " * Value inside the Core?: {";
			for (cid_iterator cid_it = cid_beg_it; cid_it != cid_end_it; ++cid_it)
			{
				gtpack::cid_type cid(*cid_it);

				if (cid_it != cid_beg_it)
				{
					std::cout << ", ";
				}

				std::cout << std::boolalpha << formed_coalitions.coalitions.at(cid).payoffs_in_core;
			}
			std::cout << "}" << std::endl;

			std::cout << " * Payoff increments wrt Grand-Coalition: {";
			for (cid_iterator cid_it = cid_beg_it; cid_it != cid_end_it; ++cid_it)
			{
				//typedef typename std::map<gtpack::player_type,RealT>::const_iterator coal_val_iterator;

				gtpack::cid_type cid(*cid_it);

				if (cid_it != cid_beg_it)
				{
					std::cout << ", ";
				}

				std::cout << "{";
				coal_val_iterator coal_val_end_it(formed_coalitions.coalitions.at(cid).payoffs.end());
				coal_val_iterator coal_val_beg_it(formed_coalitions.coalitions.at(cid).payoffs.begin());
				for (coal_val_iterator coal_val_it = coal_val_beg_it;
					 coal_val_it != coal_val_end_it;
					 ++coal_val_it)
				{
					gtpack::player_type pid(coal_val_it->first);
					RealT value(coal_val_it->second);

					if (coal_val_it != coal_val_beg_it)
					{
						std::cout << ", ";
					}

					std::cout << pid << " => " << ((value/formed_coalitions.coalitions.at(gcid).payoffs.at(pid) - 1)*100.0) << "%";
					grandpart_value += formed_coalitions.coalitions.at(gcid).payoffs.at(pid);
				}
				std::cout << "}";
			}
			std::cout << "}" << std::endl;
			std::cout << " * Value increments wrt Grand-Coalition: " << ((bestpart_value/grandpart_value-1)*100.0) << "%" << std::endl;

			std::cout << " * Payoff increments wrt Singleton Coalitions: {";
			for (cid_iterator cid_it = cid_beg_it; cid_it != cid_end_it; ++cid_it)
			{
				//typedef typename std::map<gtpack::player_type,RealT>::const_iterator coal_val_iterator;

				gtpack::cid_type cid(*cid_it);

				if (cid_it != cid_beg_it)
				{
					std::cout << ", ";
				}

				std::cout << "{";
				coal_val_iterator coal_val_end_it(formed_coalitions.coalitions.at(cid).payoffs.end());
				coal_val_iterator coal_val_beg_it(formed_coalitions.coalitions.at(cid).payoffs.begin());
				for (coal_val_iterator coal_val_it = coal_val_beg_it;
					 coal_val_it != coal_val_end_it;
					 ++coal_val_it)
				{
					gtpack::player_type pid(coal_val_it->first);
					RealT value(coal_val_it->second);

					if (coal_val_it != coal_val_beg_it)
					{
						std::cout << ", ";
					}

					gtpack::cid_type pcid = gtpack::players_coalition<RealT>::make_id(&pid, &pid+1);
					std::cout << pid << " => " << ((value/formed_coalitions.coalitions.at(pcid).payoffs.at(pid) - 1)*100.0) << "%";
					singlepart_value += formed_coalitions.coalitions.at(pcid).payoffs.at(pid);
					singlepart_kwatt += formed_coalitions.coalitions.at(pcid).optimal_allocation.kwatt;
				}
				std::cout << "}";
			}
			std::cout << "}" << std::endl;
			std::cout << " * Value increments wrt Singleton Coalitions: " << ((bestpart_value/singlepart_value-1)*100.0) << "%" << std::endl;
			std::cout << " * Energy savings wrt Singleton Coalitions: " << ((1-bestpart_kwatt/singlepart_kwatt)*100.0) << "%" << std::endl;
		}
	}
	else
	{
		std::cout << " * NOT AVAILABLE" << std::endl;
	}

	std::cout << "- Grand Coalition:" << std::endl;
	if (formed_coalitions.coalitions.count(gcid) > 0)
	{
		//typedef typename std::map<gtpack::player_type,RealT>::const_iterator coal_val_iterator;
		RealT grandpart_value(0);

		std::cout << " * Payoffs: {";
		coal_val_iterator coal_val_end_it(formed_coalitions.coalitions.at(gcid).payoffs.end());
		coal_val_iterator coal_val_beg_it(formed_coalitions.coalitions.at(gcid).payoffs.begin());
		for (coal_val_iterator coal_val_it = coal_val_beg_it;
			 coal_val_it != coal_val_end_it;
			 ++coal_val_it)
		{
			gtpack::player_type pid(coal_val_it->first);
			RealT value(coal_val_it->second);

			if (coal_val_it != coal_val_beg_it)
			{
				std::cout << ", ";
			}

			std::cout << pid << " => " << value;
			grandpart_value += value;
		}
		std::cout << "}" << std::endl;

		std::cout << " * Value: " << grandpart_value << std::endl;

		std::cout << " * Core exists?: {" << std::boolalpha << !(formed_coalitions.coalitions.at(gcid).core_empty) << "}" << std::endl;

		std::cout << " * Value inside the Core?: {" << std::boolalpha << formed_coalitions.coalitions.at(gcid).payoffs_in_core << "}" << std::endl;

#ifdef DCS_DEBUG
		if (formed_coalitions.coalitions.at(gcid).core_empty)
		{
			DCS_DEBUG_STREAM << "FOUND Grand-Coalition with empty core" << std::endl;
		}
		else
		{
			DCS_DEBUG_STREAM << "NOT FOUND Grand-Coalition with empty core" << std::endl;
		}
#endif // DCS_DEBUG
	}
	else
	{
		std::cout << " * NOT AVAILABLE" << std::endl;
	}

	std::cout << "- Singleton Coalitions:" << std::endl;
	{
		RealT singlepart_value(0);
		RealT singlepart_kwatt(0);

		std::cout << " * Payoffs: {";
		for (std::size_t c = 0; c < ncips; ++c)
		{
			gtpack::player_type pid(players[c]);

			// Retrieve the ID of the coalition for this player
			gtpack::cid_type pcid = gtpack::players_coalition<RealT>::make_id(players.begin()+c, players.begin()+c+1);

			if (c > 0)
			{
				std::cout << ", ";
			}

			std::cout << "{" << pid << " => " << formed_coalitions.coalitions.at(pcid).payoffs.at(pid) << "}";
			singlepart_value += formed_coalitions.coalitions.at(pcid).payoffs.at(pid);
			singlepart_kwatt += formed_coalitions.coalitions.at(pcid).optimal_allocation.kwatt;
		}
		std::cout << "}" << std::endl;

		std::cout << " * Value: " << singlepart_value << std::endl;
		std::cout << " * Energy Consumption: " << singlepart_kwatt << std::endl;

		std::cout << " * Core exists?: {";
		for (std::size_t c = 0; c < ncips; ++c)
		{
			gtpack::player_type pid(players[c]);

			// Retrieve the ID of the coalition for this player
			gtpack::cid_type pcid = gtpack::players_coalition<RealT>::make_id(&pid, &pid+1);

			if (c > 0)
			{
				std::cout << ", ";
			}

			std::cout << "{" << std::boolalpha << !(formed_coalitions.coalitions.at(pcid).core_empty) << "}";
		}
		std::cout << "}" << std::endl;

		std::cout << " * Value inside the Core?: {";
		for (std::size_t c = 0; c < ncips; ++c)
		{
			gtpack::player_type pid(players[c]);

			// Retrieve the ID of the coalition for this player
			gtpack::cid_type pcid = gtpack::players_coalition<RealT>::make_id(&pid, &pid+1);

			if (c > 0)
			{
				std::cout << ", ";
			}

			std::cout << "{" << std::boolalpha << formed_coalitions.coalitions.at(pcid).payoffs_in_core << "}";
		}
		std::cout << "}" << std::endl;
	}
}

template <typename RealT>
void export_csv(std::string const& fname, std::size_t ncips, coalition_formation_info<RealT> formed_coalitions, bool append = false, char field_sep=',', char line_sep='\n', char quote='"')
{
	typedef typename std::map< gtpack::cid_type, coalition_info<RealT> >::const_iterator coalition_iterator;

	std::ofstream ofs(fname.c_str(), append ? (std::ios_base::out | std::ios_base::app) : (std::ios_base::out | std::ios_base::out));
	if (!ofs)
	{
		throw std::runtime_error("Unable to open output CSV file");
	}

	// Print header or separator (i.e., an empty line)
	if (append)
	{
		for (std::size_t p = 0; p < ncips; ++p)
		{
			ofs << field_sep;
		}
	}
	else
	{
		ofs << quote << "Coalition ID" << quote;
		for (std::size_t p = 0; p < ncips; ++p)
		{
			ofs << field_sep << quote << "Payoff(CIP " << p << ")" << quote;
		}
		ofs << field_sep << quote << "Value(Coalition)" << quote;
	}
	ofs << line_sep;

	std::vector<gtpack::cid_type> cids;
	coalition_iterator coal_end_it(formed_coalitions.coalitions.end());
	for (coalition_iterator it = formed_coalitions.coalitions.begin();
		 it != coal_end_it;
		 ++it)
	{
		cids.push_back(it->first);
	}
	std::sort(cids.begin(), cids.end());

	std::size_t ncids(cids.size());
	for (std::size_t c = 0; c < ncids; ++c)
	{
		gtpack::cid_type cid(cids[c]);

		ofs << cid;

		// Output coalition value
		RealT value(0);
		for (std::size_t p = 0; p < ncips; ++p)
		{
			ofs << field_sep;
			if (formed_coalitions.coalitions.at(cid).payoffs.count(static_cast<gtpack::player_type>(p)) > 0)
			{
				ofs << formed_coalitions.coalitions.at(cid).payoffs.at(static_cast<gtpack::player_type>(p));
				value += formed_coalitions.coalitions.at(cid).payoffs.at(static_cast<gtpack::player_type>(p));
			}
		}
		ofs << field_sep << value;

		ofs << line_sep;
	}

	ofs.close();
}

template <typename RealT>
void run_experiment(scenario<RealT> const& scen, options<RealT> const& opts)
{
	const std::size_t n = opts.rnd_gen_vms ? std::max(static_cast<std::size_t>(1), opts.rnd_num_iters) : 1;

	boost::random::mt19937 rng_seed(opts.rnd_seed); // RNGs for generating seeds
	std::vector< std::vector<boost::random::mt19937> > rng_vms; // VMs RNGs
	std::vector< std::vector<boost::random::mt19937> > rng_pm_power_states; // PM power states RNGs
	std::vector< std::vector<boost::random::mt19937> > rng_pm_on_off_costs; // PM switch-on/off costs RNGs
	std::vector< std::vector< std::vector<boost::random::mt19937> > > rng_vm_migration_costs; // CIP-to-CIP VM migration costs RNGs


	if (opts.rnd_gen_vms)
	{
		boost::random::mt19937 rng(rng_seed()); // RNG used to generated random seeds for VMs RNGs
		rng_vms.resize(scen.num_cips);
		for (std::size_t c = 0; c < scen.num_cips; ++c)
		{
			rng_vms[c].resize(scen.num_vm_types);
			for (std::size_t v = 0; v < scen.num_vm_types; ++v)
			{
				rng_vms[c][v].seed(rng());
			}
		}
	}
	if (opts.rnd_gen_pm_power_states)
	{
		boost::random::mt19937 rng(rng_seed()); // RNG used to generated random seeds for PM power states RNGs
		rng_pm_power_states.resize(scen.num_cips);
		for (std::size_t c = 0; c < scen.num_cips; ++c)
		{
			rng_pm_power_states[c].resize(scen.num_pm_types);
			for (std::size_t p = 0; p < scen.num_pm_types; ++p)
			{
				rng_pm_power_states[c][p].seed(rng());
			}
		}
	}
	if (opts.rnd_gen_pm_on_off_costs)
	{
		boost::random::mt19937 rng(rng_seed()); // RNG used to generated random seeds for PM power switch-on/off costs RNGs
		rng_pm_on_off_costs.resize(scen.num_cips);
		for (std::size_t c = 0; c < scen.num_cips; ++c)
		{
			rng_pm_on_off_costs[c].resize(scen.num_pm_types);
			for (std::size_t p = 0; p < scen.num_pm_types; ++p)
			{
				rng_pm_on_off_costs[c][p].seed(rng());
			}
		}
	}
	if (opts.rnd_gen_vm_migration_costs)
	{
		boost::random::mt19937 rng(rng_seed()); // RNG used to generated random seeds for PM power switch-on/off costs RNGs
		rng_vm_migration_costs.resize(scen.num_cips);
		for (std::size_t c1 = 0; c1 < scen.num_cips; ++c1)
		{
			rng_vm_migration_costs[c1].resize(scen.num_cips);
			for (std::size_t c2 = 0; c2 < scen.num_cips; ++c2)
			{
				rng_vm_migration_costs[c1][c2].resize(scen.num_vm_types);
				for (std::size_t v = 0; v < scen.num_vm_types; ++v)
				{
					rng_vm_migration_costs[c1][c2][v].seed(rng());
				}
			}
		}
	}

	for (std::size_t i = 1; i <= n; ++i)
	{
		std::cout << "Iteration #" << i << std::endl;

		scenario<RealT> wrk_scen = scen;

		if (opts.rnd_gen_vms)
		{
			for (std::size_t c = 0; c < scen.num_cips; ++c)
			{
				for (std::size_t v = 0; v < scen.num_vm_types; ++v)
				{
					boost::random::uniform_int_distribution<std::size_t> rvg(0, scen.cip_num_vms[c][v]);
					wrk_scen.cip_num_vms[c][v] = std::max(rvg(rng_vms[c][v]), 0UL);
				}
			}
		}
		if (opts.rnd_gen_pm_power_states)
		{
			wrk_scen.cip_pm_power_states.resize(scen.num_cips);
			for (std::size_t c = 0; c < scen.num_cips; ++c)
			{
				for (std::size_t p = 0; p < scen.num_pm_types; ++p)
				{
					for (std::size_t k = 0; k < wrk_scen.cip_num_pms[c][p]; ++k)
					{
						boost::random::bernoulli_distribution<float> rvg(0.5);
						wrk_scen.cip_pm_power_states[c].push_back(rvg(rng_pm_power_states[c][p]));
					}
				}
			}
		}
		if (opts.rnd_gen_pm_on_off_costs)
		{
			// We assume the switch-on/off cost is randomly distributed as a Normal(300,50) microsec.
			// Such values are taken from [1].
			// Furthermore, we assume the switch-on cost is equal to the switch-off cost and that is
			// independent by the PM type
			//
			// - Costi switch-on degli host:
			//     <Potenza elettrica max>*<Tempo medio per passare da sleep-state a active-state>*<Costo elettricità>
			//
			// - Costi switch-off degli host:
			//     <Potenza elettrica max>*<Tempo Medio per passare da active-state a sleep-state>*<Costo elettricità>
			//
			// References:
			// 1. D. Meisner, B. Gold and T. Wenisch.
			//    "PowerNap: Eliminating Server Idle Power."
			//    In: Proc. of the ASPLOS 2009 

			const RealT norm = 3600; // normalization constant (secs in a hour)
			const RealT mu = 3e-4/norm; // Mean switch-on/off time: 300 microsec
			const RealT sigma = 5e-5/norm; // S.D. switch-on/off time: 50 microsec

			wrk_scen.cip_pm_asleep_costs.resize(scen.num_cips);
			wrk_scen.cip_pm_awake_costs.resize(scen.num_cips);
			for (std::size_t c = 0; c < scen.num_cips; ++c)
			{
				wrk_scen.cip_pm_asleep_costs[c].resize(scen.num_pm_types);
				wrk_scen.cip_pm_awake_costs[c].resize(scen.num_pm_types);
				for (std::size_t p = 0; p < scen.num_pm_types; ++p)
				{
					const RealT transition_cost_rate = scen.pm_spec_max_powers[p]*1e-3*scen.cip_electricity_costs[c]; // Transition cost in $/h

					boost::random::normal_distribution<RealT> rvg(mu, sigma);
					wrk_scen.cip_pm_asleep_costs[c][p] = wrk_scen.cip_pm_awake_costs[c][p]
													   = std::max(rvg(rng_pm_on_off_costs[c][p])*transition_cost_rate, 0.0);
				}
			}
		}
		if (opts.rnd_gen_vm_migration_costs)
		{
			// - Suppongo che il collegamento di rete fra i vari CP abbia la stessa velocità (ad es. 100Mbps).
			//
			// - Latenza tra i vari CP: suppongo sia uguale per tutti i CP
			//
			// - Tempo medio per migrare una VM da CP1 a CP2:
			//  - VM small: valore casuale Normale(277 sec, 182 sec)
			//  - VM medium: come small, ma parametri raddoppiati
			//  - VM large: come small ma parametri quadruplicati
			//  I valori di VM small li ho presi da [2], in cui sono state fatte misurazioni per collegamenti di rete a 100 Mbps, 1 Gbps e 10Gbps. Quelli per le altre due classi di VM li ho inventati, supponendo che le due classi abbiano una dimensione da migrare doppia e quadrupla di quella di una VM small.
			//  Ovviamente la migrazione tra due host di uno stesso CP ha costo zero
			//
			// - Costo upload da CP1 a CP2:
			//   Se assumo 0.01 $/GB (alla Amazon Data Transfer [3]), e quindi se trasmetto ininterroamente per T sec su una rete da 100Mbps, pagherò T*0.001 $
			//
			// - Costi migrazione VM di classe k da un host di CP1 a un host di CP2:
			//     <Tempo medio di migrazione (e downtime) di una VM  di class k da CP1 a CP2>*<Costo di upload di CP1 verso CP2>
			//
			// References:
			// 2. S. Akoush, R. Sohan, A. Rice, A.W. Moore and Andy Hopper.
			//    "Predicting the Performance of Virtual Machine Migration."
			//    In Proc. of the MASCOTS 2010.
			// 3. Amazon EC2
			//    "Amazon EC2 Data Transfer Pricing",
			//    2013, Online: https://aws.amazon.com/ec2/pricing/#DataTransfer

			const RealT norm = 3600; // normalization constant (secs in a hour)
			const RealT mu = 277/norm; // Mean migration time: 277 sec
			const RealT sigma = 61/norm; // S.D. migration time: 182 sec
			const RealT data_transfer_cost = 1e-5; // Data transfer cost per MB: 0.00001 $/MB
			const RealT activation_time = 12; // The algorithm activates every 12 hours
			const RealT data_rate = 12.5*norm; // Data rate: 12.5 MB/sec
			const RealT transfer_cost_rate = data_transfer_cost*data_rate/activation_time; // Transfer cost rate in $/hour

			wrk_scen.cip_to_cip_vm_migration_costs.resize(scen.num_cips);
			for (std::size_t c1 = 0; c1 < scen.num_cips; ++c1)
			{
				wrk_scen.cip_to_cip_vm_migration_costs[c1].resize(scen.num_cips);
				for (std::size_t c2 = 0; c2 < scen.num_cips; ++c2)
				{
					if (c1 != c2)
					{
						RealT mu2 = mu;
						RealT sigma2 = sigma;

						// We assume that VM types are ordered by increasing "size", that
						//  VMtype_1 is smaller than VMtype_2 is smaller than VMtype_3 ...
		
						wrk_scen.cip_to_cip_vm_migration_costs[c1][c2].resize(scen.num_vm_types);
						for (std::size_t v = 0; v < scen.num_vm_types; ++v)
						{
							boost::random::normal_distribution<RealT> rvg(mu2, sigma2);
							wrk_scen.cip_to_cip_vm_migration_costs[c1][c2][v] = std::max(rvg(rng_vm_migration_costs[c1][c2][v])*transfer_cost_rate, 0.0);

							mu2 *= 2;
							sigma2 *= 2;
						}
					}
					else
					{
						wrk_scen.cip_to_cip_vm_migration_costs[c1][c2].assign(scen.num_vm_types, 0);
					}
				}
			}
		}

		std::cout << "Scenario: " << wrk_scen << std::endl;
		std::cout << "Options: " << opts << std::endl;

		std::cout << "Analyzing coalitions..." << std::endl;

		detail::experiment::coalition_formation_info<RealT> formed_coalitions;
		formed_coalitions = detail::experiment::analyze_coalitions<RealT>(wrk_scen, opts);

		detail::experiment::report(wrk_scen.num_cips, formed_coalitions);

		if (!opts.csv_fname.empty())
		{
			detail::experiment::export_csv(opts.csv_fname, wrk_scen.num_cips, formed_coalitions, i > 1);
		}
	}

	std::cout << "DONE!" << std::endl;
}

}}} // Namespace detail::experiment::<unnamed>


namespace detail { namespace /*<unnamed>*/ {

void usage(char const* progname)
{
	std::cerr << "Usage: " << progname << " {options}" << std::endl
			  << "Options:" << std::endl
			  << " --csv <file>" << std::endl
			  << "   Export all the analyzed coalition onto a CSV file." << std::endl 
			  << " --formation {'nash'|'pareto'|'social'}" << std::endl
			  << "   The coalition formation strategy. Can be one of the following:" << std::endl 
			  << "   - 'merge-split': to form Merge/split-stable partitions" << std::endl
			  << "   - 'nash': to form Nash-stable partitions" << std::endl
			  << "   - 'pareto': to form Pareto-optimal partitions" << std::endl
			  << "   - 'social': to form social-optimum partitions" << std::endl
			  << " --help" << std::endl
			  << "   Show this message." << std::endl
			  << " --opt-relgap <num>" << std::endl
			  << "   A real number in [0,1] used to set the relative gap parameter of the optimal solver." << std::endl
			  << " --opt-tilim <num>" << std::endl
			  << "   A real positive number used to set the maximum number of seconds to wait for the termination of the optimal solver." << std::endl
			  << " --payoff {'banzhaf'|'norm-banzhaf'|'shapley'}" << std::endl
			  << "   The coalition value division strategy. Can be one of the following:" << std::endl 
			  << "   - 'banzhaf': the Banzhaf value" << std::endl 
			  << "   - 'norm-banzhaf': the normalized Banzhaf value" << std::endl 
			  << "   - 'shapley': the Shapley value" << std::endl 
			  << " --rnd-genvms" << std::endl
			  << "    Enable the random generation of VMs for each CIP." << std::endl
			  << " --rnd-genpmsonoff" << std::endl
			  << "    Enable the random generation of PM power states for each CIP." << std::endl
			  << " --rnd-genpmsonoffcosts" << std::endl
			  << "    Enable the random generation of switch-on/off costs of PMs for each CIP and PM type." << std::endl
			  << " --rnd-genvmsmigrcosts" << std::endl
			  << "    Enable the random generation of CIP-to-CIP migration costs of VMs for each CIP and VM type." << std::endl
			  << " --rnd-numit <number>" << std::endl
			  << "   Set the number of times that the given scenario must be run." << std::endl
			  << "   Each run will use a randomly generated number of VMs and PM power states." << std::endl 
			  << " --rnd-seed <number>" << std::endl
			  << "   Set the seed to use for random number generation." << std::endl
			  << " --scenario <file-name>" << std::endl
			  << "   The path to the file describing the scenario to use for the experiment." << std::endl;
}

}} // Namespace detail::<unnamed>


int main(int argc, char* argv[])
{
	typedef double real_type;

	if (argc < 2)
	{
		detail::usage(argv[0]);
		return -1;
	}

	std::string opt_csv_fname;
	bool opt_help;
	real_type opt_relative_gap;
	real_type opt_time_lim;
	std::string opt_scenario_file;
	detail::experiment::coalition_formation_category opt_coalition_formation;
	detail::experiment::coalition_value_division_category opt_coalition_value_division;
	bool opt_rnd_gen_vms;
	bool opt_rnd_gen_pm_power_states;
	bool opt_rnd_gen_pm_on_off_costs;
	bool opt_rnd_gen_vm_migration_costs;
	std::size_t opt_rnd_num_iters;
	unsigned long opt_rnd_seed;
	std::string opt_str;

	// Parse CLI options
	opt_csv_fname = cli::simple::get_option<std::string>(argv, argv+argc, "--csv", "");
	opt_help = cli::simple::get_option(argv, argv+argc, "--help");
	if (opt_help)
	{
		detail::usage(argv[0]);
		return 0;
	}
	opt_str = cli::simple::get_option<std::string>(argv, argv+argc, "--formation", "nash");
	if (opt_str == "nash")
	{
		opt_coalition_formation = detail::experiment::nash_stable_coalition_formation;
	}
	else if (opt_str == "merge-split")
	{
		opt_coalition_formation = detail::experiment::merge_split_stable_coalition_formation;
	}
	else if (opt_str == "pareto")
	{
		opt_coalition_formation = detail::experiment::pareto_optimal_coalition_formation;
	}
	else if (opt_str == "social")
	{
		opt_coalition_formation = detail::experiment::social_optimum_coalition_formation;
	}
	else
	{
		throw std::invalid_argument("Unknown coalition formation category");
	}
	opt_str = cli::simple::get_option<std::string>(argv, argv+argc, "--payoff", "shapley");
	if (opt_str == "banzhaf")
	{
		opt_coalition_value_division = detail::experiment::banzhaf_coalition_value_division;
	}
	else if (opt_str == "norm-banzhaf")
	{
		opt_coalition_value_division = detail::experiment::normalized_banzhaf_coalition_value_division;
	}
	else if (opt_str == "shapley")
	{
		opt_coalition_value_division = detail::experiment::shapley_coalition_value_division;
	}
	else
	{
		throw std::invalid_argument("Unknown coalition value division category");
	}
	opt_relative_gap = cli::simple::get_option<real_type>(argv, argv+argc, "--opt-relgap", 0);
	opt_time_lim = cli::simple::get_option<real_type>(argv, argv+argc, "--opt-tilim", -1);
	opt_rnd_gen_vms = cli::simple::get_option(argv, argv+argc, "--rnd-genvms");
	opt_rnd_gen_pm_power_states = cli::simple::get_option(argv, argv+argc, "--rnd-genpmsonoff");
	opt_rnd_gen_pm_on_off_costs = cli::simple::get_option(argv, argv+argc, "--rnd-genpmsonoffcosts");
	opt_rnd_gen_vm_migration_costs = cli::simple::get_option(argv, argv+argc, "--rnd-genvmsmigrcosts");
	opt_rnd_num_iters = cli::simple::get_option<std::size_t>(argv, argv+argc, "--rnd-numit", 1);
	opt_rnd_seed = cli::simple::get_option<unsigned long>(argv, argv+argc, "--rnd-seed", 5489);
	opt_scenario_file = cli::simple::get_option<std::string>(argv, argv+argc, "--scenario");

	// Check CLI options
	if (opt_scenario_file.empty())
	{
		std::cerr << "(E) Scenario file not specified" << std::endl;
		detail::usage(argv[0]);
		return -1;
	}

	// Run Experiment
	detail::experiment::options<real_type> opts;
	detail::experiment::scenario<real_type> scenario;

	scenario = detail::experiment::make_scenario<real_type>(opt_scenario_file);
	opts.opt_relative_gap = opt_relative_gap;
	opts.opt_time_lim = opt_time_lim;
	opts.coalition_formation = opt_coalition_formation;
	opts.coalition_value_division = opt_coalition_value_division;
	opts.rnd_gen_vms = opt_rnd_gen_vms;
	opts.rnd_gen_pm_power_states = opt_rnd_gen_pm_power_states;
	opts.rnd_gen_pm_on_off_costs = opt_rnd_gen_pm_on_off_costs;
	opts.rnd_gen_vm_migration_costs = opt_rnd_gen_vm_migration_costs;
	opts.rnd_num_iters = opt_rnd_num_iters;
	opts.rnd_seed = opt_rnd_seed;
	opts.csv_fname = opt_csv_fname;

	detail::experiment::run_experiment<real_type>(scenario, opts);
}
