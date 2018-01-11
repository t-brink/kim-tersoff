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
   Modified for use with KIM by Tobias Brink (2012,2013,2017).
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

enum KIM_NBC { // Only for internal use; no mapping to any integer
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
  // Parameters are stored together in a struct because I found that
  // using the layout that conforms to KIM (one array per parameter)
  // can lead to a slowdown of up to 5% because the data is not in
  // adjacent memory.  The parameter data is not so much, so we will
  // copy between the KIM-published parameters and this internal data.
  struct Params2 {
    // Cutoff related.
    double cutsq;
    double R, D;
    // Two-body parameters.
    double lam1;
    double A;
    double B;
    double lam2;
    double beta;
    double n;
    // Pre-computed.
    double n_precomp[4];
  };
  struct Params3 {
    // Cutoff related.
    double cutsq;
    double R, D;
    // Three-body parameters.
    int m;
    double lam3;
    double gamma;
    double h;
    // Pre-computed.
    double c2;    // c^2
    double d2;    // d^2
    double c2_d2; // c^2 / d^2
    // Not used in computation.
    //double c, d;
  };
  struct KIMParams {
    explicit KIMParams(int N) // Number of particle types
      : A(N,N), B(N,N),
        lam1(N,N), lam2(N,N), lam3(N,N,N),
        c(N,N,N), d(N,N,N), h(N,N,N),
        gamma(N,N,N),
        m(N,N,N),
        n(N,N), beta(N,N),
        D(N,N,N), R(N,N,N)
    {
      shape2[0] = N; shape2[1] = N;
      shape3[0] = N; shape3[1] = N; shape3[2] = N;
    }
    // Copy data from a Params array.
    void from_params(const Array2D<Params2>& p2, const Array3D<Params3>& p3) {
      for (int i = 0; i < A.extent(0); ++i)
        for (int j = 0; j < A.extent(1); ++j) {
          A(i,j) = p2(i,j).A;
          B(i,j) = p2(i,j).B;
          lam1(i,j) = p2(i,j).lam1;
          lam2(i,j) = p2(i,j).lam2;
          n(i,j) = p2(i,j).n;
          beta(i,j) = p2(i,j).beta;
          for (int k = 0; k < lam3.extent(2); ++k) {
            lam3(i,j,k) = p3(i,j,k).lam3;
            //c(i,j,k) = p3(i,j,k).c; // those are not kept there,
            //d(i,j,k) = p3(i,j,k).d; // but only in derived form c², d²
            h(i,j,k) = p3(i,j,k).h;
            gamma(i,j,k) = p3(i,j,k).gamma;
            m(i,j,k) = p3(i,j,k).m;
            D(i,j,k) = p3(i,j,k).D;
            R(i,j,k) = p3(i,j,k).R;
          }
        }
    };
    // Copy data to a Params array.
    void to_params(Array2D<Params2>& p2, Array3D<Params3>& p3) const {
      for (int i = 0; i < lam3.extent(0); ++i)
        for (int j = 0; j < lam3.extent(1); ++j)
          for (int k = 0; k < lam3.extent(2); ++k) {
            p2(i,j).A = A(i,j);
            p2(i,j).B = B(i,j);
            p2(i,j).lam1 = lam1(i,j);
            p2(i,j).lam2 = lam2(i,j);
            p3(i,j,k).lam3 = lam3(i,j,k);
            //p3(i,j,k).c = c(i,j,k); // those are not kept there,
            //p3(i,j,k).d = d(i,j,k); // but only in derived form c², d²
            p3(i,j,k).h = h(i,j,k);
            p3(i,j,k).gamma = gamma(i,j,k);
            p3(i,j,k).m = m(i,j,k);
            p2(i,j).n = n(i,j);
            p2(i,j).beta = beta(i,j);
            p3(i,j,k).D = D(i,j,k);
            p3(i,j,k).R = R(i,j,k);
          }
    }
    Array2D<double> A, B;
    Array2D<double> lam1, lam2;
    Array3D<double> lam3;
    Array3D<double> c, d, h;
    Array3D<double> gamma;
    Array3D<int> m;
    Array2D<double> n, beta;
    // Cutoff related.
    Array3D<double> D, R;
    // The shape of all parameter arrays. Needed for calls to
    // KIM_API_model::set_shape().
    int shape2[2];
    int shape3[3];
  };

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
  void update_params(); // Copy from KIM-published parameters to internal.
  double cutoff() const {
    return max_cutoff;
  }
  const KimIndices kim_indices;
  KIMParams kim_params; // Parameters published to KIM, see above why
                        // we keep two copies.

 private:
  int n_spec;                   // number of species
  Array2D<Params2> params2;     // n_spec*n_spec array of parameters
  Array3D<Params3> params3;     // n_spec*n_spec*n_spec array of parameters
  double max_cutoff;            // max cutoff for all elements
  std::map<int,std::string> to_spec;  // map element index to element
                                      // name, needed for user-
                                      // friendly error messages

  void read_params(std::istream&, std::map<std::string,int>,
                   double, double, double);
  void prepare_params();
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
