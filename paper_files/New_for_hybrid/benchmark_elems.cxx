/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2025 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#include <array>

#include <t8.h>

#include <sc_options.h>
#include <sc_statistics.h>
#include <sc_functions.h>

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_examples.h>

#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_io.h>

#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid.hxx>


#include <t8_types/t8_vec.hxx>


/**
 * Create a partitioned cmesh. If no msh_file is given, a new hybrid cmesh is created.
 * 
 * \param[in] msh_file 
 * \param[in] mesh_dim 
 * \param[in] comm 
 * \param[in] init_level 
 * \return t8_cmesh_t 
 */
t8_cmesh_t
t8_benchmark_forest_create_cmesh ( sc_MPI_Comm comm, const int init_level, const t8_eclass_t eclass, const int num_elems)
{
  T8_ASSERT (eclass != T8_ECLASS_INVALID);
  t8_cmesh_t cmesh = t8_cmesh_new_bigmesh ( eclass, num_elems, comm);
  t8_cmesh_t cmesh_partition;
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_set_partition_uniform (cmesh_partition, init_level, t8_scheme_new_default ());
  t8_cmesh_set_profiling (cmesh_partition, 1);
  t8_cmesh_commit (cmesh_partition, comm);
  t8_cmesh_destroy (&cmesh);
  return cmesh_partition;
}


/* refine the forest in a band, given by a plane E and two constants
 * c_min, c_max. We refine the cells in the band c_min*E, c_max*E */
static int
t8_adapt_pyramid ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, [[maybe_unused]] t8_locidx_t which_tree, 
                [[maybe_unused]]t8_eclass_t tree_class,
               [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]]const t8_scheme *scheme, [[maybe_unused]]const int is_family,
               [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
    const t8_dpyramid_t *pyra = (const t8_dpyramid_t *) elements[0];
    const int type = pyra->pyramid.type;
    if (type == 6 || type == 0 || type == 2 || type == 4 ){
        return 1;
    }
    else {
        return 0;
    }
}

static int
t8_adapt_second ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, [[maybe_unused]] t8_locidx_t which_tree, 
                [[maybe_unused]]t8_eclass_t tree_class,
               [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]]const t8_scheme *scheme, [[maybe_unused]]const int is_family,
               [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
    const int child_id = scheme->element_get_child_id (tree_class, elements[0]);
    return child_id % 2;
}

static void
benchmark_band_adapt(t8_cmesh_t cmesh, sc_MPI_Comm comm, const int init_level, const t8_eclass_t eclass)
{
  double adapt_time = 0;
  double partition_time = 0;
  double new_time = 0;
  double total_time = 0;
  double ghost_time = 0;
  double balance_time = 0;
  const int num_stats = 6;
  std::array<sc_statinfo_t, num_stats> times;
  sc_stats_init (&times[0], "new");
  sc_stats_init (&times[1], "adapt");
  sc_stats_init (&times[2], "partition");
  sc_stats_init (&times[3], "ghost");
  sc_stats_init (&times[4], "balance");
  sc_stats_init (&times[5], "total");

  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_set_cmesh(forest, cmesh, comm);
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  t8_forest_set_level (forest, init_level); 

  total_time -= sc_MPI_Wtime ();

  new_time -= sc_MPI_Wtime ();
  t8_forest_commit (forest);
  new_time += sc_MPI_Wtime ();

  sc_stats_set1 (&times[0], new_time, "new");


  t8_forest_t forest_adapt, forest_partition;
  t8_forest_init (&forest_adapt);
  if (eclass == T8_ECLASS_PYRAMID) {
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_pyramid, 0);
  }
  else {
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_second, 0);
  }
  t8_forest_set_profiling (forest_adapt, 1);


  t8_forest_commit (forest_adapt);
  adapt_time += t8_forest_profile_get_adapt_time(forest_adapt);

  t8_forest_ref (forest_adapt);

  t8_forest_init (&forest_partition);
  t8_forest_set_partition(forest_partition, forest_adapt, 0);
  t8_forest_set_profiling (forest_partition, 1);
  t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
  t8_forest_set_balance (forest_partition, NULL, 0);

  t8_forest_commit (forest_partition);
  forest = forest_partition;
  int ghost_sent = 0;
  int procs_sent = 0;
  int balance_rounds = 0;
  partition_time += t8_forest_profile_get_partition_time (forest_partition, &procs_sent);
  ghost_time += t8_forest_profile_get_ghost_time (forest_partition, &ghost_sent);
  balance_time += t8_forest_profile_get_balance_time (forest_partition, &balance_rounds);


  t8_forest_unref (&forest_adapt);
  
  total_time += sc_MPI_Wtime ();

  sc_stats_accumulate (&times[0], new_time);
  sc_stats_accumulate (&times[1], adapt_time);
  sc_stats_accumulate (&times[2], partition_time);
  sc_stats_accumulate (&times[3], ghost_time);
  sc_stats_accumulate (&times[4], balance_time);
  sc_stats_accumulate (&times[5], total_time);
  sc_stats_compute (comm, num_stats, times.data ());
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, num_stats, times.data (), 1, 1);
  t8_forest_unref (&forest_partition);
}

int
main (int argc, char **argv)
{

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);
  int help = 0;
  int initial_level;
  int eclass_int;
  int num_runs;
  int num_elems;

  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  sc_options_t *options = sc_options_new (argv[0]);

  sc_options_add_switch (options, 'h', "help", &help, "Print this help message and exit");
  sc_options_add_int (options, 'e', "eclass", &eclass_int, 0,
                      "0: Tetrahedron, 1: Hexahedron, 2: Prism, 3: Pyramid");
  sc_options_add_int (options, 'l', "level", &initial_level, 0, "The initial uniform refinement level of the forest.");
  sc_options_add_int (options, 'n', "num-runs", &num_runs, 1,
                          "The number of runs to perform. If specified, the program will run num_runs times with the same parameters. ");
  sc_options_add_int (options, 'N', "num-elems", &num_elems, 1,
                          "The number of elements in the forest. If specified, the program will create a forest with num_elems elements. ");
  const int options_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, options, argc, argv);

  if( options_argc <= 0 || options_argc != argc || help )
  {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, options, NULL);
    return 1;
  }
  t8_eclass_t eclass = T8_ECLASS_INVALID;

  switch (eclass_int)
  {
  case 0:
    eclass = T8_ECLASS_TET;
    break;
  case 1:
    eclass = T8_ECLASS_HEX;
    break;
  case 2:
    eclass = T8_ECLASS_PRISM;
    break;
  case 3:
    eclass = T8_ECLASS_PYRAMID;
    break;
  default:
    break;
  }

  for (int irun = 0; irun < num_runs; ++irun) {
    t8_global_essentialf ("#################### Run %d of %d ####################\n", irun + 1, num_runs);
    t8_cmesh_t cmesh = t8_benchmark_forest_create_cmesh (sc_MPI_COMM_WORLD, initial_level, eclass, num_elems);


    benchmark_band_adapt (cmesh, sc_MPI_COMM_WORLD, initial_level, eclass);
  }

  sc_options_destroy (options);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
