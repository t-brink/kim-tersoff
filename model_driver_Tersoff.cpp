/*
  Copyright (c) 2012,2013,2014,2018,2019,2020 Tobias Brink

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
#include <utility>
#include <fstream>
#include <sstream>

#include "KIM_LogMacros.hpp"
#include "KIM_ModelDriverHeaders.hpp"

#include "pair_tersoff.hpp"
#include "ndarray.hpp"

using namespace std;
using namespace model_driver_Tersoff;


extern "C" {
  // Can't be both static and extern.  But this is not needed: extern
  // just tells the compiler to use C naming conventions for the
  // function instead of C++ (name mangling).  All functions except
  // for the dirver creation are not linked at compile time but passed
  // at runtime by a pointer, so they do not need to be extern (and
  // should not be, there could be namespace clashes).
  int model_driver_create(KIM::ModelDriverCreate * const,
                          const KIM::LengthUnit,
                          const KIM::EnergyUnit,
                          const KIM::ChargeUnit,
                          const KIM::TemperatureUnit,
                          const KIM::TimeUnit);
}

// For some reason, KIM wants a pointer to that instead of just
// copying a bool. So be it, I'll define a constant, which is false
// since we do use ghost particles' neighbors.
static const int doesnt_use_ghost_neighbors = 0;

// LOCAL DEFINITIONS ///////////////////////////////////////////////////

enum PotentialVariant { standard, zbl };

// Helper to trim a string. For some reason C++ doesn't provide this.
string trim(const string &s)
{
    string::const_iterator it = s.begin();
    while (it != s.end() && isspace(*it))
        it++;

    string::const_reverse_iterator rit = s.rbegin();
    while (rit.base() != it && isspace(*rit))
        rit++;

    return string(it, rit.base());
}


// WRAPPERS AND INTERFACE TO KIM ///////////////////////////////////////

static int
compute_arguments_create(const KIM::ModelCompute * const, // unused
                         KIM::ModelComputeArgumentsCreate * const
                                              model_compute_arguments_create) {
  int error =

    // Tell KIM what inputs/outputs are supported and how.
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialForces,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
      KIM::SUPPORT_STATUS::optional)

    // Register callback support.
    ||
    model_compute_arguments_create->SetCallbackSupportStatus(
      KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetCallbackSupportStatus(
      KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,
      KIM::SUPPORT_STATUS::notSupported)
    ;

  return error;
}
#undef KIM_LOGGER_OBJECT_NAME

#define KIM_LOGGER_OBJECT_NAME model_compute

static int
compute(const KIM::ModelCompute * const model_compute,
        const KIM::ModelComputeArguments * const model_compute_arguments) {
  PairTersoff* tersoff;
  model_compute->GetModelBufferPointer(reinterpret_cast<void **>(&tersoff));

  // Unpack data.
  const int * n_atoms;
  const int * atom_types;
  const int * contributing;
  const double * atom_coords_ptr;

  double * energy;
  double * atom_energy;
  double * forces_ptr;
  double * virial;
  double * particle_virial_ptr;

  int error =

    // Input
    model_compute_arguments->GetArgumentPointer(
      KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles, &n_atoms)
    ||
    model_compute_arguments->GetArgumentPointer(
      KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes, &atom_types)
    ||
    model_compute_arguments->GetArgumentPointer(
      KIM::COMPUTE_ARGUMENT_NAME::particleContributing, &contributing)
    ||
    model_compute_arguments->GetArgumentPointer(
      KIM::COMPUTE_ARGUMENT_NAME::coordinates, &atom_coords_ptr)

    // Output
    ||
    model_compute_arguments->GetArgumentPointer(
      KIM::COMPUTE_ARGUMENT_NAME::partialEnergy, &energy)
    ||
    model_compute_arguments->GetArgumentPointer(
      KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy, &atom_energy)
    ||
    model_compute_arguments->GetArgumentPointer(
      KIM::COMPUTE_ARGUMENT_NAME::partialForces, &forces_ptr)
    ||
    model_compute_arguments->GetArgumentPointer(
      KIM::COMPUTE_ARGUMENT_NAME::partialVirial, &virial)
    ||
    model_compute_arguments->GetArgumentPointer(
      KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial, &particle_virial_ptr);

  if (error) return error;

  int compute_process_dEdr;
  error =
    model_compute_arguments->IsCallbackPresent(
      KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm, &compute_process_dEdr);

  if (error) return error;

  // Wrap some stuff for convenience.
  Array2D<const double> atom_coords(atom_coords_ptr, *n_atoms, 3);
  Array2D<double> f(forces_ptr, *n_atoms, 3);
  Array2D<double>* forces = forces_ptr ? &f : NULL;
  Array2D<double> v(particle_virial_ptr, *n_atoms, 6);
  Array2D<double>* particle_virial = particle_virial_ptr ? &v : NULL;

  // Do the compute.
  try {
    tersoff->compute(*model_compute_arguments,
                     *n_atoms,
                     atom_types,
                     contributing,
                     atom_coords,
                     energy,
                     atom_energy,
                     forces,
                     virial,
                     particle_virial,
                     compute_process_dEdr);
  } catch (const exception& e) {
    LOG_ERROR(string("compute: ") + e.what());
    return 1;
  }

  return 0;
}
#undef KIM_LOGGER_OBJECT_NAME


#define KIM_LOGGER_OBJECT_NAME model_refresh

static int
refresh(KIM::ModelRefresh * const model_refresh) {
  PairTersoff* tersoff;
  model_refresh->GetModelBufferPointer(reinterpret_cast<void **>(&tersoff));

  // Recalculate internal derived values.
  try {
    tersoff->update_params();
  } catch (const exception& e) {
    LOG_ERROR(string("refresh: ") + e.what());
    return 1;
  }

  // Republish cutoff.
  model_refresh->SetInfluenceDistancePointer(tersoff->cutoff_ptr());
  // Pass info about neighbor lists.
  model_refresh->SetNeighborListPointers(1, tersoff->cutoff_ptr(),
                                         &doesnt_use_ghost_neighbors);

  return 0;
}
#undef KIM_LOGGER_OBJECT_NAME


/*
static int
write_parameterized_model(.....) {     
  TODO: would be super-nice, but have to find example    
}
*/

static int
compute_arguments_destroy(const KIM::ModelCompute * const, // all ununsed
                          KIM::ModelComputeArgumentsDestroy * const) {
  // We did not allocate anything for compute_arguments_create(), thus
  // no cleanup is needed.

  return 0;
}


#define KIM_LOGGER_OBJECT_NAME model_destroy

static int destroy(KIM::ModelDestroy * const model_destroy) {
  PairTersoff* tersoff;
  model_destroy->GetModelBufferPointer(reinterpret_cast<void **>(&tersoff));

  if (tersoff != NULL) {
    delete tersoff;
  } else {
    LOG_ERROR("destroy: tried to destroy a model driver that is already null");
  }

  return 0;
}
#undef KIM_LOGGER_OBJECT_NAME


// Init stuff must be last, so that the other functions are already defined.

#define KIM_LOGGER_OBJECT_NAME model_driver_create

static int
read_settings(KIM::ModelDriverCreate * const model_driver_create,
              const string& settings_filename,
              int& n_spec, map<string,int>& type_map,
              PotentialVariant& potential_variant) {
  ifstream settings_file(settings_filename.c_str()); // passing the std::string
                                                     // is C++11
  bool got_line;

  // Get the list of species. //////////////////////////////////////////
  string species_line;
  got_line = getline(settings_file, species_line);
  if (!got_line) {
    LOG_ERROR("The settings file ("
              + settings_filename
              + ") does not contain a line with supported particle types.");
    return 1;
  }

  // Parse the list.
  istringstream iss(species_line);
  string species_name;
  int species_id = 0;
  while (iss >> species_name) {
    // Collect in map.
    pair<map<string,int>::iterator, bool> insertion_result =
      type_map.insert(pair<string,int>(species_name, species_id));
    if (!insertion_result.second) {
      LOG_ERROR("Particle type \"" + species_name + "\" occurs twice in file "
                + settings_filename);
      return 1;
    }
    // Register to KIM.
    const KIM::SpeciesName kim_spec(species_name);
    const int error = model_driver_create->SetSpeciesCode(kim_spec, species_id);
    if (error) {
      LOG_ERROR("Error returned by KIM's SetSpeciesCode().");
      return error;
    }
    //
    ++species_id;
  }
  n_spec = type_map.size();

  // See if there is a variant (e.g. ZBL) requested ////////////////////
  potential_variant = standard;
  string variant_line;
  got_line = getline(settings_file, variant_line);
  if (got_line) {
    variant_line = trim(variant_line);
    if (!variant_line.empty()) {
      if (variant_line == "ZBL") {
        potential_variant = zbl;
      } else {
        LOG_ERROR("Illegal potential variant ("
                  + variant_line
                  + ") specified in second line of settings file ("
                  + settings_filename
                  + ")");
        return 1;
      }
    }
  }

  return 0;
}
#undef KIM_LOGGER_OBJECT_NAME

#define KIM_LOGGER_OBJECT_NAME model_driver_create

static int
init_unit_conv(KIM::ModelDriverCreate * const model_driver_create,
               const KIM::LengthUnit length_unit,
               const KIM::EnergyUnit energy_unit,
               const KIM::ChargeUnit charge_unit,
               const KIM::TemperatureUnit temperature_unit,
               const KIM::TimeUnit time_unit,
               double& length_conv,
               double& inv_length_conv,
               double& energy_conv) {
  int error;

  // Length ////////////////////////////////////////////////////////////
  error = model_driver_create->ConvertUnit(KIM::LENGTH_UNIT::A,
                                           KIM::ENERGY_UNIT::eV,
                                           KIM::CHARGE_UNIT::e,
                                           KIM::TEMPERATURE_UNIT::K,
                                           KIM::TIME_UNIT::ps,
                                           length_unit, energy_unit,
                                           charge_unit,
                                           temperature_unit, time_unit,
                                           1.0, 0.0, 0.0, 0.0, 0.0,
                                           &length_conv);
  if (error) {
    LOG_ERROR("Error returned by KIM's ConvertUnit() when trying to "
              "get length units.");
    return error;
  }

  // Inverse length ////////////////////////////////////////////////////
  error = model_driver_create->ConvertUnit(KIM::LENGTH_UNIT::A,
                                           KIM::ENERGY_UNIT::eV,
                                           KIM::CHARGE_UNIT::e,
                                           KIM::TEMPERATURE_UNIT::K,
                                           KIM::TIME_UNIT::ps,
                                           length_unit, energy_unit,
                                           charge_unit,
                                           temperature_unit, time_unit,
                                           -1.0, 0.0, 0.0, 0.0, 0.0,
                                           &inv_length_conv);
  if (error) {
    LOG_ERROR("Error returned by KIM's ConvertUnit() when trying to "
              "get inverse length units.");
    return error;
  }

  // Energy ////////////////////////////////////////////////////////////
  error = model_driver_create->ConvertUnit(KIM::LENGTH_UNIT::A,
                                           KIM::ENERGY_UNIT::eV,
                                           KIM::CHARGE_UNIT::e,
                                           KIM::TEMPERATURE_UNIT::K,
                                           KIM::TIME_UNIT::ps,
                                           length_unit, energy_unit,
                                           charge_unit,
                                           temperature_unit, time_unit,
                                           0.0, 1.0, 0.0, 0.0, 0.0,
                                           &energy_conv);
  if (error) {
    LOG_ERROR("Error returned by KIM's ConvertUnit() when trying to "
              "get energy units.");
    return error;
  }

  // KIM wants to know what happened, I guess. /////////////////////////
  error = model_driver_create->SetUnits(length_unit, energy_unit,
                                        KIM::CHARGE_UNIT::unused,
                                        KIM::TEMPERATURE_UNIT::unused,
                                        KIM::TIME_UNIT::unused);
  if (error) {
    LOG_ERROR("Error returned by KIM's SetUnits().");
    return error;
  }

  return 0;
}
#undef KIM_LOGGER_OBJECT_NAME


#define REG2BODY(memb, name, expl)                                           \
  error =                                                                    \
    model_driver_create->SetParameterPointer(tersoff->kim_params.size2,      \
                                             &tersoff->kim_params.memb(0,0), \
                                             #name,                          \
                                             "The two-body parameter "       \
                                             #name " " #expl ". "            \
                                             "Size N*N, where N is the "     \
                                             "number of species supported "  \
                                             "by the model. Storage in "     \
                                             "row-major order by ascending " \
                                             "species code.");               \
  if (error) {                                                               \
    return 1;                                                                \
  }

#define REG3BODY(memb, name, expl)                                           \
  error =                                                                    \
    model_driver_create->SetParameterPointer(tersoff->kim_params.size3,      \
                                             &tersoff->kim_params.memb(0,0,0),\
                                             #name,                          \
                                             "The three-body parameter "     \
                                             #name " " #expl ". "            \
                                             "Size N*N*N, where N is the "   \
                                             "number of species supported "  \
                                             "by the model. Storage in "     \
                                             "row-major order by ascending " \
                                             "species code.");               \
  if (error) {                                                               \
    return 1;                                                                \
  }

static int
reg_params(KIM::ModelDriverCreate * const model_driver_create,
           PairTersoff * const tersoff) {
  int error;

  // Two-body parameters.
  REG2BODY(A, A, in units of energy);
  REG2BODY(B, B, in units of energy);
  REG2BODY(lam1, lambda1, in units of inverse length);
  REG2BODY(lam2, lambda2, in units of inverse length);
  REG2BODY(beta, beta, (unitless));
  REG2BODY(n, n, (unitless));

  // Three-body parameters.
  REG3BODY(lam3, lambda3, in units of inverse length);
  REG3BODY(m, m, (unitless). This parameter is an integer exponent of
           value 1 or 3 that is used to implement slightly different
           variants of the Tersoff potential);
  REG3BODY(gamma, gamma, (unitless));
  REG3BODY(c, c, (unitless));
  REG3BODY(d, d, (unitless));
  REG3BODY(h, h, (unitless));
  REG3BODY(R, Rc, in units of length. This is a cutoff parameter);
  REG3BODY(D, Dc, in units of length. This is a cutoff parameter);

  return 0;
}

#undef REG2BODY
#undef REG3BODY

#define KIM_LOGGER_OBJECT_NAME model_driver_create

int
model_driver_create(KIM::ModelDriverCreate * const model_driver_create,
                    const KIM::LengthUnit length_unit,
                    const KIM::EnergyUnit energy_unit,
                    const KIM::ChargeUnit charge_unit,
                    const KIM::TemperatureUnit temperature_unit,
                    const KIM::TimeUnit time_unit) {
  int error;

  // Get parameter files. //////////////////////////////////////////////
  int n_param_files;
  model_driver_create->GetNumberOfParameterFiles(&n_param_files);
  if (n_param_files != 2) {
    // Since we cannot use C++11's to_string (which would be dead
    // easy), we have to this dance instead:
    ostringstream s;
    s << n_param_files;
    LOG_ERROR("This model driver requires exactly two parameter files, but "
              + s.str() + " were provided.");
    return 1;
  }

  const string * settings_filename;
  error = model_driver_create->GetParameterFileName(0, &settings_filename);
  if (error) {
    LOG_ERROR("Error returned by KIM's GetParameterFileName() "
              "for the first parameter file.");
    return 1;
  }

  const string * param_filename;
  error = model_driver_create->GetParameterFileName(1, &param_filename);
  if (error) {
    LOG_ERROR("Error returned by KIM's GetParameterFileName() "
              "for the second parameter file.");
    return 1;
  }

  // Get number and name of species. ///////////////////////////////////
  int n_spec = 0;
  map<string,int> type_map;
  PotentialVariant potential_variant;      
  error =
    read_settings(model_driver_create, *settings_filename,
                  n_spec, type_map, potential_variant);
  if (error) {
    return error; // already logged.
  }

  // Get variant of potential (e.g. ZBL). //////////////////////////////
  //TODO (will go into read_settings, too)  

  // Init unit conversion factors. /////////////////////////////////////
  // We are using LAMMPS "metal" units, i.e., eV and Ångstroms.
  double length_conv;
  double inv_length_conv;
  double energy_conv;
  error = init_unit_conv(model_driver_create,
                         length_unit, energy_unit, charge_unit,
                         temperature_unit, time_unit,
                         length_conv,
                         inv_length_conv,
                         energy_conv);
  if (error) {
    return error; // already logged.
  }


  // Init the core class. //////////////////////////////////////////////
  PairTersoff* tersoff;
  try {
    tersoff = new PairTersoff(*param_filename, n_spec, type_map,
                              energy_conv, length_conv, inv_length_conv);
  } catch (const exception& e) {
    LOG_ERROR(string("model_driver_create: ") + e.what());
    return 1; // error
  }

  // Pass stuff to KIM. ////////////////////////////////////////////////
  model_driver_create->SetModelBufferPointer(static_cast<void *>(tersoff));
  model_driver_create->SetInfluenceDistancePointer(tersoff->cutoff_ptr());
  model_driver_create->SetNeighborListPointers(1, tersoff->cutoff_ptr(),
                                               &doesnt_use_ghost_neighbors);
  error = model_driver_create->SetModelNumbering(KIM::NUMBERING::zeroBased);
  if (error) {
    LOG_ERROR("Error returned by KIM's SetModelNumbering().");
    delete tersoff;
    return 1;
  }

  // Register parameters.
  error = reg_params(model_driver_create, tersoff);
  if (error) {
    delete tersoff;
    return error; // logging already done in reg_params()
  }

  // Use function pointer definitions to statically verify correct
  // prototypes.
  KIM::ModelComputeArgumentsCreateFunction * kim_ca_create
    = &compute_arguments_create;
  KIM::ModelComputeFunction * kim_compute = &compute;
  KIM::ModelRefreshFunction * kim_refresh = &refresh;
  //KIM::ModelWriteParameterizedModelFunction * kim_write_params
  //  = &write_parameterized_model;
  KIM::ModelComputeArgumentsDestroyFunction * kim_ca_destroy
    = &compute_arguments_destroy;
  KIM::ModelDestroyFunction * kim_destroy = &destroy;

  // Register the function pointers.
  error =
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::ComputeArgumentsCreate,
      KIM::LANGUAGE_NAME::cpp, true,
      reinterpret_cast<KIM::Function *>(kim_ca_create))
    ||
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::Compute,
      KIM::LANGUAGE_NAME::cpp, true,
      reinterpret_cast<KIM::Function *>(kim_compute))
    ||
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::Refresh,
      KIM::LANGUAGE_NAME::cpp, false,
      reinterpret_cast<KIM::Function *>(kim_refresh))
    ||
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::ComputeArgumentsDestroy,
      KIM::LANGUAGE_NAME::cpp, true,
      reinterpret_cast<KIM::Function *>(kim_ca_destroy))
    ||
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::Destroy,
      KIM::LANGUAGE_NAME::cpp, true,
      reinterpret_cast<KIM::Function *>(kim_destroy));
  if (error) {
    LOG_ERROR("Error returned by KIM's SetRoutinePointer().");
    delete tersoff;
    return 1;
  }


  return 0;
}

#undef KIM_LOGGER_OBJECT_NAME
