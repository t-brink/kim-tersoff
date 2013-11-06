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
   Modified for use with KIM by Tobias Brink (2012,2013).
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
    // Optional inputs.
    int boxSideLengths;
    // Optional outputs.
    int energy;
    //int compute_energy; // was energy requested?
    int particleEnergy;
    //int compute_particleEnergy; // was particleEnergy requested?
    int forces;
    //int compute_forces; // was forces requested?
    int virial;
    int particleVirial;

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
               double*, double*, double*, Array2D<double>*,
               double*, Array2D<double>*) const;
  void prepare_params();
  double cutoff() const {
    return max_cutoff;
  }


 private:
  // Parameters are stored together in a struct because I found that
  // using the layout that conforms to KIM (one array per parameter)
  // can lead to a slowdown of up to 5% because the data is not in
  // adjacent memory.  The parameter data is not so much, so we will
  // copy between the KIM-published parameters and this internal data.
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
  struct KIMParams {
    explicit KIMParams(int N) // Number of particle types
      : A(N,N,N), B(N,N,N),
        lam1(N,N,N), lam2(N,N,N), lam3(N,N,N),
        c(N,N,N), d(N,N,N), h(N,N,N),
        gamma(N,N,N),
        m(N,N,N),
        n(N,N,N), beta(N,N,N),
        D(N,N,N), R(N,N,N) {}
    // Copy data from a Params array.
    void from_params(const Array3D<Params>& p) {
      for (int i = 0; i < A.extent(0); ++i)
        for (int j = 0; j < A.extent(1); ++j)
          for (int k = 0; k < A.extent(2); ++k) {
            A(i,j,k) = p(i,j,k).A;
            B(i,j,k) = p(i,j,k).B;
            lam1(i,j,k) = p(i,j,k).lam1;
            lam2(i,j,k) = p(i,j,k).lam2;
            lam3(i,j,k) = p(i,j,k).lam3;
            c(i,j,k) = p(i,j,k).c;
            d(i,j,k) = p(i,j,k).d;
            h(i,j,k) = p(i,j,k).h;
            gamma(i,j,k) = p(i,j,k).gamma;
            m(i,j,k) = p(i,j,k).m;
            n(i,j,k) = p(i,j,k).n;
            beta(i,j,k) = p(i,j,k).beta;
            D(i,j,k) = p(i,j,k).D;
            R(i,j,k) = p(i,j,k).R;
          }
    };
    // Copy data to a Params array.
    void to_params(Array3D<Params>& p) const {
      for (int i = 0; i < A.extent(0); ++i)
        for (int j = 0; j < A.extent(1); ++j)
          for (int k = 0; k < A.extent(2); ++k) {
            p(i,j,k).A = A(i,j,k);
            p(i,j,k).B = B(i,j,k);
            p(i,j,k).lam1 = lam1(i,j,k);
            p(i,j,k).lam2 = lam2(i,j,k);
            p(i,j,k).lam3 = lam3(i,j,k);
            p(i,j,k).c = c(i,j,k);
            p(i,j,k).d = d(i,j,k);
            p(i,j,k).h = h(i,j,k);
            p(i,j,k).gamma = gamma(i,j,k);
            p(i,j,k).m = m(i,j,k);
            p(i,j,k).n = n(i,j,k);
            p(i,j,k).beta = beta(i,j,k);
            p(i,j,k).D = D(i,j,k);
            p(i,j,k).R = R(i,j,k);
          }
    }
    Array3D<double> A, B;
    Array3D<double> lam1, lam2, lam3;
    Array3D<double> c, d, h;
    Array3D<double> gamma;
    Array3D<int> m;
    Array3D<double> n, beta;
    // Cutoff related.
    Array3D<double> D, R;
  };

  int n_spec;                   // number of species
  Array3D<Params> params;       // n_spec*n_spec*n_spec array of parameters
  KIMParams kim_params;         // Parameters published to KIM, see above
                                // why we keep two copies.
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
