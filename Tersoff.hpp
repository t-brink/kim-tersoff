/*
  Copyright (c) 2012 Tobias Brink

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

#ifndef TERSOFF_TERSOFF_HPP
#define TERSOFF_TERSOFF_HPP

#include <string>
#include <map>
#include <fstream>
#include <istream>

#include <KIM_API.h>

#include "ndarray.hpp"

namespace model_driver_Tersoff {

  /*!
    Internal use enum to quickly check for neighbor list/boundary
    condition style.
  */
  enum KIM_NBC { // Only for internal use no mapping to any integer
    // returned from the KIM API.
    KIM_CLUSTER,
    KIM_NEIGH_PURE_H,
    KIM_NEIGH_PURE_F,
    KIM_NEIGH_RVEC_F,
    KIM_MI_OPBC_H,
    KIM_MI_OPBC_F
  };

  /*!
    Constants for integers used to tell the get_neigh() KIM API
    function if we want an iterator or random access.
  */
  enum KIM_IterLoca { // constants defined by the KIM API!
    KIM_ITERATOR_MODE = 0, // Sequential access.
    KIM_LOCATOR_MODE = 1   // This is random access.
  };

  /*!
    KIM related data that can be stored at least until reinit is called.

    These are mostly KIM indices for fast acces with getm/setm.
  */
  struct KimIndices {
    // Required inputs.
    int numberOfParticles;
    int particleTypes;
    int coordinates;
    // Optional outputs.
    int energy;
    //int compute_energy; // was energy requested?
    int particleEnergy;
    //int compute_particleEnergy; // was particleEnergy requested?
    int forces;
    //int compute_forces; // was forces requested?

    // Non-index stuff.
    KIM_NBC nbc;
    KIM_IterLoca neigh_access_mode;
  };

  //////////////////////////////////////////////////////////////////////

  /*!
    This is the actual implementation of the model driver.
  */
  class Tersoff {
  public:
    const KimIndices kim_indices;
    /*!
      This constructor reads the parameters from an input file, the
      name of which is given as a parameter.
    */
    explicit Tersoff(std::string parameter_file,
                     int n_spec,
                     std::map<std::string,int> type_map,
                     // Conversion factors.
                     double energy_conv,
                     double length_conv,
                     double inv_length_conv,
                     // KIM indices.
                     const KimIndices& ki
                     )
      : kim_indices(ki), n_spec(n_spec), params(n_spec, n_spec, n_spec)
    {
      std::fstream f(parameter_file.c_str(), std::ios_base::in);
      read_params(f, type_map, energy_conv, length_conv, inv_length_conv);
    }

    /*!
      Destructor.
    */
    ~Tersoff() {
      // Nothing to do for now.
    }

    /*!
      Return maximum cutoff for creating neighbor lists.
    */
    double cutoff() const {
      return max_cutoff;
    }

    /*!
      Do the actual computation for a cluster.
    */
    void compute_cluster(// Input.
                         int n_atoms,
                         const int* atom_types,
                         const Array2D<double>& atom_coords,
                         // Output. All optional.
                         double* energy,
                         double* atom_energy,
                         Array2D<double>* forces
                         ) const;
    /*!
      Do the actual computation with methods that do neighbor lists.
    */
    void compute_neigh_list(// Input.
                            KIM_API_model& kim_model,
                            KIM_IterLoca access_mode,
                            int n_atoms,
                            const int* atom_types,
                            const Array2D<double>& atom_coords,
                            // Output. All optional.
                            double* energy,
                            double* atom_energy
                            ) const;

  private:
    /*!
      Storage for parameters.

      Holds parameters for one interaction (e.g. Si-C-C, with Si being
      the central atom).
    */
    struct Params {
      double A, B;
      double lam1, lam2, lam3;
      double c, d, h;
      double gamma;
      int m;
      double n, beta;
      // Cutoff related.
      double D, R;
      double cut, cutmin, cutsq;
    };

    //! Number of atom types.
    int n_spec;
    //! Holds all parameters.
    Array3D<Params> params;
    //! The biggest cutoff distance.  Used by the test.
    double max_cutoff;

    /*!
      Read and init parameters from file.
    */
    void read_params(std::istream& infile,
                     std::map<std::string,int> type_map,
                     double energy_conv,
                     double length_conv,
                     double inv_length_conv);
  };

}

#endif /* TERSOFF_TERSOFF_HPP */
