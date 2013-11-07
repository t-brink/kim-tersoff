/*
  Copyright (c) 2012,2013 Tobias Brink

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstdlib>
#include <string>
#include <map>
#include <set>
#include <stdexcept>

#include "KIM_API.h"
#include "KIM_API_status.h"

#include "pair_tersoff.hpp"
#include "ndarray.hpp"

using namespace std;
using namespace model_driver_Tersoff;


extern "C" {
  // Can't be both static and extern.  But this is not needed: extern
  // just tells the compiler to use C naming conventions for the
  // function instead of C++ (name mangling).  The functions destroy
  // and compute are not linked at compile time but passed at runtime
  // by a pointer, so they do not need to be extern (and should not
  // be, there could be namespace clashes).
  /*
  static int destroy(KIM_API_model**);
  static int compute(KIM_API_model**);
  */
  int model_driver_init(void*, char*, int*, int*);
}


// WRAPPERS AND INTERFACE TO KIM ///////////////////////////////////////


static int destroy(KIM_API_model** kimmdl) {
  KIM_API_model& kim_model = **kimmdl;
  int error;
  PairTersoff* tersoff =
    static_cast<PairTersoff*>(kim_model.get_model_buffer(&error));
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer",
                           error);
    return error;
  }
  delete tersoff;

  return KIM_STATUS_OK;
}


static int reinit(KIM_API_model** kimmdl) {
  KIM_API_model& kim_model = **kimmdl;

  int error;
  PairTersoff* tersoff =
    static_cast<PairTersoff*>(kim_model.get_model_buffer(&error));
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer",
                           error);
    return error;
  }

  // Copy KIM-published parameters to internal storage format and
  // re-compute some pre-calculated values from new parameters.
  try {
    tersoff->update_params();
  } catch (const exception& e) {
    kim_model.report_error(__LINE__, __FILE__, e.what(), KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }

  // Re-publish cutoff.
  double& cutoff =
    *static_cast<double*>(kim_model.get_data("cutoff", &error));
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_data(\"cutoff\")",
                           error);
    return error;
  }
  cutoff = tersoff->cutoff();

  return KIM_STATUS_OK;
}


static int compute(KIM_API_model** kimmdl) {
  KIM_API_model& kim_model = **kimmdl;
  int error;
  PairTersoff* tersoff =
    static_cast<PairTersoff*>(kim_model.get_model_buffer(&error));
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer",
                           error);
    return error;
  }

  const KimIndices& ki = tersoff->kim_indices;
  int compute_energy, compute_particleEnergy, compute_forces,
    compute_virial, compute_particleVirial;
  kim_model.getm_compute_by_index(&error, 5*3,
                                  ki.energy, &compute_energy, 1,
                                  ki.particleEnergy, &compute_particleEnergy, 1,
                                  ki.forces, &compute_forces, 1,
                                  ki.virial, &compute_virial, 1,
                                  ki.particleVirial, &compute_particleVirial, 1
                                  );
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_getm_compute",
                           error);
    return error;
  }

  // Get model input and requested outputs.
  int* n_atoms;
  int* atom_types; // 1D array
  double* c;
  double* energy;
  double* particle_energy;
  double* f;
  double* virial;
  double* particleVirial_;
  double* boxSideLengths;
  kim_model.getm_data_by_index(&error, 9*3,
                               ki.numberOfParticles, &n_atoms, 1,
                               ki.particleTypes, &atom_types, 1,
                               ki.coordinates, &c, 1,
                               ki.boxSideLengths, &boxSideLengths,
                                                  (ki.nbc == KIM_MI_OPBC_F),
                               ki.energy, &energy, compute_energy,
                               ki.particleEnergy, &particle_energy,
                                                  compute_particleEnergy,
                               ki.forces, &f, compute_forces,
                               ki.virial, &virial, compute_virial,
                               ki.particleVirial, &particleVirial_,
                                                  compute_particleVirial
                               );
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index",
                           error);
    return error;
  }
  Array2D<double> coord(c, *n_atoms, 3);
  if (!compute_energy)
    energy = NULL;
  if (ki.nbc != KIM_MI_OPBC_F)
    boxSideLengths = NULL;
  if (!compute_particleEnergy)
    particle_energy = NULL;
  Array2D<double> forces(f, *n_atoms, 3);
  Array2D<double>* forces_ptr = compute_forces ? &forces : NULL;
  if (!compute_virial)
    virial = NULL;
  Array2D<double> particleVirial(particleVirial_, *n_atoms, 6);
  Array2D<double>* particleVirial_ptr =
    compute_particleVirial ? &particleVirial : NULL;

  // Calculate values.
  try {
    bool use_neighbor_list, use_distvec;
    switch (ki.nbc) {
    case KIM_CLUSTER:
    case KIM_MI_OPBC_F:
      use_neighbor_list = false;
      use_distvec = false;
      break;
    case KIM_NEIGH_PURE_F:
      use_neighbor_list = true;
      use_distvec = false;
      break;
    case KIM_NEIGH_RVEC_F:
      use_neighbor_list = true;
      use_distvec = true;
      break;
    default:
      kim_model.report_error(__LINE__, __FILE__,
                             "Unknown NBC, probably bug in model driver!",
                             KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
    tersoff->compute(kim_model, use_neighbor_list, use_distvec,
                     ki.neigh_access_mode,
                     *n_atoms, atom_types, coord, boxSideLengths,
                     energy, particle_energy, forces_ptr,
                     virial, particleVirial_ptr);
  } catch (const exception& e) {
    kim_model.report_error(__LINE__, __FILE__, e.what(), KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }

  return KIM_STATUS_OK;
}



int model_driver_init(void* km, // The KIM model object
                      // Array[numparamfiles][nmstrlen]
                      char* paramfile_names,
                      int* nmstrlen,
                      int* numparamfiles) {
  /* KIM variables */
  KIM_API_model& kim_model = **static_cast<KIM_API_model**>(km);

  int error = KIM_STATUS_OK;

  // Get number and names of species.
  map<string,int> spec_idx;
  set<int> indices;
  int n_spec;
  char* partcl_types = kim_model.get_model_partcl_typs(&n_spec, &error);
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_model_partcl_typs",
                           error);
    return error;
  } else if (partcl_types == NULL) {
    kim_model.report_error(__LINE__, __FILE__,
                           "No species names specified, "
                           "but this model driver needs at least one "
                           "species name.", KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }
  for (int i = 0; i != n_spec; ++i) {
    string spec_name(partcl_types + i*KIM_KEY_STRING_LENGTH);
    int type_idx =
      kim_model.get_partcl_type_code(spec_name.c_str(), &error);
    if (error != KIM_STATUS_OK) {
      kim_model.report_error(__LINE__, __FILE__,
                             "KIM_API_get_partcl_type_code", error);
      return error;
    }
    spec_idx[spec_name] = type_idx;
    indices.insert(type_idx);
  }
  free(partcl_types);
  // Check if the indices are continuous (as set members are unique,
  // the minimum must be 0, the maximum n_spec-1).  The set is
  // guaranteed to contain at least one entry, as we check for empty
  // species list above.
  int min_idx = *indices.begin();  // First key is the smallest.
  int max_idx = *indices.rbegin(); // Last key is the biggest.
  if (min_idx != 0 || max_idx != n_spec - 1) {
    kim_model.report_error(__LINE__, __FILE__,
                           "Species indices are not continuous.  "
                           "Models using this driver MUST use continuous "
                           "indices for the atom species!",
                           KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }

  // Prepare parameter file.
  if (*numparamfiles != 1) {
    error = KIM_STATUS_FAIL;
    kim_model.report_error(__LINE__, __FILE__,
                           "Incorrect number of parameter files.",
                           error);
    return error;
  }
  string param_filename(paramfile_names); // NUL-terminated, can treat
                                          // as simple string here
                                          // (only one file)

  // Init unit conversion factors.  We are using LAMMPS "metal" units,
  // i.e. eV and Ångstroms.
  const double energy_conv =
    kim_model.convert_to_act_unit("A", "eV", "e", "K", "ps",
                                  0.0,  1.0, 0.0, 0.0,  0.0,
                                  &error);
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                           error);
    return error;
  }
  const double length_conv =
    kim_model.convert_to_act_unit("A", "eV", "e", "K", "ps",
                                  1.0,  0.0, 0.0, 0.0,  0.0,
                                  &error);
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                           error);
    return error;
  }
  const double inv_length_conv =
    kim_model.convert_to_act_unit("A", "eV", "e", "K", "ps",
                                  -1.0, 0.0, 0.0, 0.0,  0.0,
                                  &error);
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                           error);
    return error;
  }

  // Get KIM indices for efficiency.
  KimIndices kim_indices;
  kim_model.getm_index(&error, 9*3,
                       "numberOfParticles", &kim_indices.numberOfParticles, 1,
                       "particleTypes", &kim_indices.particleTypes, 1,
                       "coordinates", &kim_indices.coordinates, 1,
                       "boxSideLengths", &kim_indices.boxSideLengths, 1,
                       "energy", &kim_indices.energy, 1,
                       "particleEnergy", &kim_indices.particleEnergy, 1,
                       "forces", &kim_indices.forces, 1,
                       "virial", &kim_indices.virial, 1,
                       "particleVirial", &kim_indices.particleVirial, 1
                       );
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_getm_index",
                           error);
    return error;
  }

  // Get neighbor list style/boundary conditions.
  char* nbc_str = kim_model.get_NBC_method(&error);
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method",
                           error);
    return error;
  }
  string nbc(nbc_str);
  free(nbc_str);

  if (nbc == "NEIGH_RVEC_F")
    kim_indices.nbc = KIM_NEIGH_RVEC_F;
  else if (nbc == "NEIGH_PURE_F")
    kim_indices.nbc = KIM_NEIGH_PURE_F;
  else if (nbc == "MI_OPBC_F")
    kim_indices.nbc = KIM_MI_OPBC_F;
  else if (nbc == "CLUSTER")
    kim_indices.nbc = KIM_CLUSTER;
  else {
    kim_model.report_error(__LINE__, __FILE__,
                           ("Unknown NBC: " + nbc +
                            ". What are you doing?").c_str(),
                           KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }

  // Iterator or locator mode?
  int iter_loca = kim_model.get_neigh_mode(&error);
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_neigh_mode",
                           error);
    return error;
  }
  switch (iter_loca) {
  case 1:
    kim_indices.neigh_access_mode = KIM_ITERATOR_MODE;
    break;
  case 2:
    kim_indices.neigh_access_mode = KIM_LOCATOR_MODE;
    break;
  /*
  case 3:
    ????
    break;
  */
  default:
    kim_model.report_error(__LINE__, __FILE__,
                           "Unknown/unsupported neighbor mode.",
                           KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }

  // Init the model object (reads parameters from file).
  PairTersoff* tersoff;
  try {
    tersoff = new PairTersoff(param_filename, n_spec, spec_idx,
                              energy_conv, length_conv, inv_length_conv,
                              kim_indices);
  } catch (const exception& e) {
    kim_model.report_error(__LINE__, __FILE__, e.what(), KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }

  // Pass functions and parameters to KIM.
  kim_model.setm_data(&error, 17*4,
                      "compute", 1, &compute, 1,
                      "reinit",  1, &reinit,  1,
                      "destroy", 1, &destroy, 1,
                      "PARAM_FREE_A", 1, &tersoff->kim_params.A(0,0,0), 1,
                      "PARAM_FREE_B", 1, &tersoff->kim_params.B(0,0,0), 1,
                      "PARAM_FREE_lambda1", 1, &tersoff->kim_params.lam1(0,0,0), 1,
                      "PARAM_FREE_lambda2", 1, &tersoff->kim_params.lam2(0,0,0), 1,
                      "PARAM_FREE_lambda3", 1, &tersoff->kim_params.lam3(0,0,0), 1,
                      "PARAM_FREE_beta", 1, &tersoff->kim_params.beta(0,0,0), 1,
                      "PARAM_FREE_n", 1, &tersoff->kim_params.n(0,0,0), 1,
                      "PARAM_FREE_m", 1, &tersoff->kim_params.m(0,0,0), 1,
                      "PARAM_FREE_gamma", 1, &tersoff->kim_params.gamma(0,0,0), 1,
                      "PARAM_FREE_c", 1, &tersoff->kim_params.c(0,0,0), 1,
                      "PARAM_FREE_d", 1, &tersoff->kim_params.d(0,0,0), 1,
                      "PARAM_FREE_h", 1, &tersoff->kim_params.h(0,0,0), 1,
                      "PARAM_FREE_Rc", 1, &tersoff->kim_params.R(0,0,0), 1,
                      "PARAM_FREE_Dc", 1, &tersoff->kim_params.D(0,0,0), 1
                      );
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_setm_data", error);
    delete tersoff;
    return error;
  }

  // Set up cutoff for KIM.
  double& cutoff =
    *static_cast<double*>(kim_model.get_data("cutoff", &error));
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_data(\"cutoff\")",
                           error);
    delete tersoff;
    return error;
  }
  cutoff = tersoff->cutoff();

  // Store model object.
  kim_model.set_model_buffer(static_cast<void*>(tersoff), &error);
  if (error != KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer",
                           error);
    delete tersoff;
    return error;
  }

  return KIM_STATUS_OK;
}
