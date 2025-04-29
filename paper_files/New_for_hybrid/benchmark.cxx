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
#include <t8_forest/t8_forest_types.h>

#include <t8_schemes/t8_default/t8_default.hxx>

#include <t8_types/t8_vec.hxx>

typedef struct
{
  double c_min, c_max; /* constants that define the thickness of the refinement region */
  t8_3D_vec normal;    /* normal vector to the plane E */
  int base_level;      /* A given level that is not coarsend further, see -l argument */
  int max_level;       /* A max level that is not refined further, see -L argument */
} adapt_data_t;

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
t8_benchmark_forest_create_cmesh (const char *msh_file, const int mesh_dim, sc_MPI_Comm comm, const int init_level, [[ maybe_unused ]]const t8_eclass_t eclass)
{
  t8_cmesh_t cmesh;
  if (msh_file != NULL){
    cmesh = t8_cmesh_from_msh_file ((char *) msh_file, 1, comm, mesh_dim, 0, false);
  }
  else {
    T8_ASSERT (eclass != T8_ECLASS_INVALID);
    cmesh = t8_cmesh_new_hypercube ( eclass, comm, 0, 0, 0);
  }
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
t8_band_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
               [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family,
               [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  t8_3D_vec elem_midpoint;
  T8_ASSERT (!is_family || num_elements == scheme->element_get_num_children (tree_class, elements[0]));
  const int level = scheme->element_get_level (tree_class, elements[0]);
  /* Get the minimum and maximum x-coordinate from the user data pointer of forest */
  const adapt_data_t * adapt_data = (adapt_data_t *) t8_forest_get_user_data (forest);
  const t8_3D_vec normal = adapt_data->normal;
  const int base_level = adapt_data->base_level;
  const int max_level = adapt_data->max_level;
  /* Compute the coordinates of the anchor node. */
  t8_forest_element_centroid (forest_from, which_tree, elements[0], elem_midpoint.data ());

  /* Calculate elem_midpoint - c_min n */
  t8_axpy (normal, elem_midpoint, -adapt_data->c_min);

  /* The purpose of the factor C*h is that the levels get smaller, the
   * closer we get to the interface. We refine a cell if it is at most
   * C times its own height away from the interface */
  if (t8_dot (elem_midpoint, normal) >= 0) {
    /* if the anchor node is to the right of c_min*E,
     * check if it is to the left of c_max*E */

    /* set elem_midpoint to the original anchor - c_max*normal */
    t8_axpy (normal, elem_midpoint, -(adapt_data->c_max - adapt_data->c_min));
    if (t8_dot (elem_midpoint, normal) <= 0) {
      if (level < max_level) {
        /* We do refine if level smaller 1+base level and the anchor is
         * to the left of c_max*E */
        return 1;
      }
    }
    else if (is_family && level > base_level) {
      /* Otherwise, we coarse if we have a family and level is greater
       * than the base level. */
      return -1;
    }
  }
  else if (is_family && level > base_level) {
    /* If element lies out of the refinement region and a family was given
     * as argument, we coarsen to level base level */
    /* set elem_midpoint to the original midpoint - c_max*normal */
    return -1;
  }
  return 0;
}

static void
benchmark_band_adapt(t8_cmesh_t cmesh, sc_MPI_Comm comm, const int init_level, const int max_level, const std::array<double, 2> &x_min_max, 
                      const double delta_t, const double max_time, const int do_ghost)
{
  double adapt_time = 0;
  double partition_time = 0;
  double new_time = 0;
  double total_time = 0;
  double ghost_time = 0;
  const int num_stats = 5;
  std::array<sc_statinfo_t, num_stats> times;
  sc_stats_init (&times[0], "new");
  sc_stats_init (&times[1], "adapt");
  sc_stats_init (&times[2], "partition");
  sc_stats_init (&times[3], "ghost");
  sc_stats_init (&times[4], "total");

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

  t8_3D_vec normal({0.8, 0.3, 0.0});
  adapt_data_t adapt_data = {x_min_max[0], x_min_max[1], normal, init_level, max_level};
  t8_normalize (adapt_data.normal);
  int num_steps = 0;
  t8_forest_t forest_adapt, forest_partition;
  for (double time = 0; time < max_time; time += delta_t, ++num_steps) {
    t8_forest_init (&forest_adapt);
    t8_forest_set_adapt (forest_adapt, forest, t8_band_adapt, 1);
    t8_forest_set_profiling (forest_adapt, 1);

    adapt_data.c_min = x_min_max[0] + time ;
    adapt_data.c_max = x_min_max[1] + time ;

    t8_forest_set_user_data (forest_adapt, (void *)&adapt_data);
    adapt_time -= sc_MPI_Wtime ();
    t8_forest_commit (forest_adapt);
    adapt_time += sc_MPI_Wtime ();

    t8_forest_ref (forest_adapt);

    t8_forest_init (&forest_partition);
    t8_forest_set_partition(forest_partition, forest_adapt, 0);
    t8_forest_set_profiling (forest_partition, 1);
    if (do_ghost) {
      t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
    }
 
    t8_forest_commit (forest_partition);
    forest = forest_partition;

    t8_profile_t *profile = forest_partition->profile;

    partition_time = profile->partition_runtime;
    ghost_time = profile->ghost_runtime;


    t8_forest_unref (&forest_adapt);
  }
  
  total_time += sc_MPI_Wtime ();

  t8_global_productionf ("Num steps: %d\n", num_steps);

  sc_stats_accumulate (&times[0], new_time);
  sc_stats_accumulate (&times[1], adapt_time);
  sc_stats_accumulate (&times[2], partition_time);
  sc_stats_accumulate (&times[3], ghost_time);
  sc_stats_accumulate (&times[4], total_time);
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
  const char *mshfileprefix = NULL;
  int dim;
  int initial_level;
  int level_diff;
  std::array<double, 2> x_min_max;
  double T;
  double cfl = 0;
  int eclass_int;
  int num_runs;
  int do_ghost;

  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  sc_options_t *options = sc_options_new (argv[0]);

  sc_options_add_switch (options, 'h', "help", &help, "Print this help message and exit");
  sc_options_add_string (options, 'f', "mshfile", &mshfileprefix, NULL,
                         "If specified, the cmesh is constructed from a .msh file with the given prefix. "
                         "The files must end in .msh and be created with gmsh.");
  sc_options_add_int (options, 'e', "eclass", &eclass_int, 0,
                      "If no mshfile is given, the cmesh is created with the given element class. "
                      "0: Tetrahedron, 1: Hexahedron, 2: Prism, 3: Pyramid");
  sc_options_add_int (options, 'd', "dim", &dim, 2, "Together with -f: The dimension of the coarse mesh. 2 or 3.");
  sc_options_add_int (options, 'l', "level", &initial_level, 0, "The initial uniform refinement level of the forest.");
  sc_options_add_int (options, 'r', "rlevel", &level_diff, 1,
                      "The number of levels that the forest is refined from the initial level.");
  sc_options_add_switch (options, 'g', "ghost", &do_ghost, "If specified, the forest is created with ghost cells.");
  sc_options_add_double (options, 'x', "xmin", &x_min_max[0], 0, "The minimum x coordinate in the mesh.");
  sc_options_add_double (options, 'X', "xmax", &x_min_max[1], 1, "The maximum x coordinate in the mesh.");
  sc_options_add_double (options, 'T', "time", &T, 1,
                         "The simulated time span. We simulate the time from 0 to T. T has to be > 0.");
  /* CFL number. delta_t = CFL * 0.64 / 2^level */
  sc_options_add_double (options, 'C', "cfl", &cfl, 0,
                         "The CFL number. If specified, then delta_t is set to CFL * 0.64 / 2^level. ");
  sc_options_add_int (options, 'n', "num-runs", &num_runs, 1,
                          "The number of runs to perform. If specified, the program will run num_runs times with the same parameters. ");
  
  const int options_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, options, argc, argv);

  if( options_argc <= 0 || options_argc != argc || help || initial_level < 0 || level_diff <= 0 || cfl == 0)
  {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, options, NULL);
    return 1;
  }
  const double delta_t = cfl * 0.64 / (1 << initial_level);
  t8_global_productionf ("Using CFL %f, delta_t = %f\n", cfl, delta_t);

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

  T8_ASSERT (mshfileprefix != NULL || eclass != T8_ECLASS_INVALID);

  t8_global_productionf ("Using mshfileprefix %s with dim %d\n", mshfileprefix, dim);
  const int max_level = initial_level + level_diff;
  for (int irun = 0; irun < num_runs; ++irun) {
    t8_global_essentialf ("#################### Run %d of %d ####################\n", irun + 1, num_runs);
    t8_cmesh_t cmesh = t8_benchmark_forest_create_cmesh (mshfileprefix, dim, sc_MPI_COMM_WORLD, initial_level, eclass);


    benchmark_band_adapt (cmesh,  sc_MPI_COMM_WORLD, initial_level, max_level, x_min_max, delta_t, T, do_ghost);
  }

  sc_options_destroy (options);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
