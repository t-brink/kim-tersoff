/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Modified for use with KIM by Tobias Brink (2012).
------------------------------------------------------------------------- */


#ifndef LMP_PAIR_TERSOFF_H
#define LMP_PAIR_TERSOFF_H

#include <string>
#include <map>
#include <istream>

#include <KIM_API.h>

#include "ndarray.hpp"

namespace model_driver_Tersoff {

// The M_PI* macros are not part of the C++ standard, so to be safe we
// use these constants instead.
const double pi   = 3.14159265358979323846;  /* pi */
const double pi_2 = 1.57079632679489661923;  /* pi/2 */
const double pi_4 = 0.78539816339744830962;  /* pi/4 */

  enum KIM_NBC { // Only for internal use no mapping to any integer
                 // returned from the KIM API.
    KIM_CLUSTER,
    KIM_NEIGH_PURE_H,
    KIM_NEIGH_PURE_F,
    KIM_NEIGH_RVEC_F,
    KIM_MI_OPBC_H,
    KIM_MI_OPBC_F
  };

  enum KIM_IterLoca { // constants defined by the KIM API!
    KIM_ITERATOR_MODE = 0, // Sequential access.
    KIM_LOCATOR_MODE = 1   // This is random access.
  };

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


class PairTersoff /*: public Pair*/ {
 public:
  const KimIndices kim_indices;
  PairTersoff(std::string parameter_file,
              int n_spec,
              std::map<std::string,int> type_map,
              // Conversion factors.
              double energy_conv,
              double length_conv,
              double inv_length_conv,
              // KIM indices.
              const KimIndices& ki
              );
  virtual ~PairTersoff();
  void compute(KIM_API_model&, bool, bool, KIM_IterLoca,
               int, const int*, const Array2D<double>&,
               double*, double*, Array2D<double>*) const;
  void prepare_params();
  double cutoff() const {
    return max_cutoff;
  }


 private:
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
    // Pre-computed.
    double n_precomp[4];
    double c2;    // c^2
    double d2;    // d^2
    double c2_d2; // c^2 / d^2
  };

  int n_spec;                   // number of species
  Array3D<Params> params;       // n_spec*n_spec*n_spec array of parameters
  double max_cutoff;            // max cutoff for all elements
  std::map<int,std::string> to_spec;  // map element index to element
                                      // name, needed for user-
                                      // friendly error messages

  void read_params(std::istream&, std::map<std::string,int>,
                   double, double, double);
  double repulsive(double, double, double, double, double,
                   bool, double&) const;
  double zeta(double, double,
              int, double, double, double,
              double, double, double, double,
              double,
              double*, double*) const;
  double force_zeta(double, double, double, double,
                    double, double,
                    double, double,
                    const double[4],
                    double&, bool, double&) const;
  void attractive(double, double, double,
                  double, double, int, double,
                  double, double, double, double,
                  double,
                  double *, double *,
                  double *, double *, double *) const;

  double ters_fc(double, double, double) const;
  double ters_fc_d(double, double, double) const;
  double ters_fa(double, double, double, double) const;
  double ters_fa_d(double, double, double, double, double) const;
  double ters_bij(double, double, double, const double[4]) const;
  double ters_bij_d(double, double, double, const double[4]) const;

  void ters_zetaterm_d(double,
                       double, double, int, double,
                       double, double, double,
                       double, double,
                       double *, double, double *, double,
                       double *, double *, double *) const;
  void costheta_d(double *, double, double *, double,
                  double *, double *, double *) const;

  // inlined functions for efficiency

  inline double ters_gijk(double costheta,
                          double gamma, double c2, double d2, double c2_d2,
                          double h) const {
    const double hcth = h - costheta;
    return gamma*(1.0 + c2_d2 - c2 / (d2 + hcth*hcth));
  }

  inline double ters_gijk_d(const double costheta,
                            double gamma, double c2, double d2,
                            double h) const {
    const double hcth = h - costheta;
    const double numerator = -2.0 * c2 * hcth;
    const double denominator = 1.0/(d2 + hcth*hcth);
    return gamma*numerator*denominator*denominator;
  }

  inline double vec3_dot(const double x[3], const double y[3]) const {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

  inline void vec3_add(const double x[3], const double y[3],
                       double * const z) const {
    z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
  }

  inline void vec3_scale(const double k, const double x[3],
                         double y[3]) const {
    y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
  }

  inline void vec3_scaleadd(const double k, const double x[3],
                            const double y[3], double * const z) const {
    z[0] = k*x[0]+y[0];
    z[1] = k*x[1]+y[1];
    z[2] = k*x[2]+y[2];
  }
};

}

#endif
