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
   Modified for use with KIM by Tobias Brink (2020).
------------------------------------------------------------------------- */


#ifndef LMP_PAIR_TERSOFF_ZBL_H
#define LMP_PAIR_TERSOFF_ZBL_H

#include <string>
#include <map>
#include <istream>

#include "KIM_ModelDriverHeaders.hpp"

#include "pair_tersoff.hpp"
#include "ndarray.hpp"

namespace model_driver_Tersoff {


class PairTersoffZBL : public PairTersoff {
 public:
  // Parameters are stored together in a struct because I found that
  // using the layout that conforms to KIM (one array per parameter)
  // can lead to a slowdown of up to 5% because the data is not in
  // adjacent memory.  The parameter data is not so much, so we will
  // copy between the KIM-published parameters and this internal data.
  struct ParamsZBL2 {
    // Two-body parameters.
    double Z_i; // XXX not needed    
    double Z_j; // XXX not needed     
    double ZBLcut;
    double ZBLexpscale;
    // Pre-computed.
    double a;       // 0.8854 * a0 / (Z_i^0.23 + Z_j^0.23)
    double premult; // Z_i * Z_j * e^2 / (4 * pi * epsilon0)
  };
  struct KIMParamsZBL {
    explicit KIMParamsZBL(int N) // Number of particle types
      : Z_i(N,N), Z_j(N,N),
        ZBLcut(N,N), ZBLexpscale(N,N) {}
    // Copy data from a Params array.
    void from_params(const Array2D<ParamsZBL2>& p2) {
      for (int i = 0; i < Z_i.extent(0); ++i)
        for (int j = 0; j < Z_i.extent(1); ++j) {
          
        }
    };
    // Copy data to a Params array.
    void to_params(Array2D<ParamsZBL2>& p2) const {
      for (int i = 0; i < Z_i.extent(0); ++i)
        for (int j = 0; j < Z_i.extent(1); ++j) {
          
        }
    }
    Array2D<double> Z_i, Z_j;
    Array2D<double> ZBLcut, ZBLexpscale;
    // The size of all parameter arrays. Needed for calls to
    // KIM::ModelDriverCreate.SetParameterPointer().
    //const int size2; -> can take this from PairTersoff::KIMParams!
  };

  PairTersoffZBL(const std::string& parameter_file,
                 int n_spec,
                 std::map<std::string,int> type_map,
                 // Conversion factors.
                 double energy_conv,
                 double inv_energy_conv,
                 double length_conv,
                 double inv_length_conv,
                 double charge_conv);
  virtual ~PairTersoffZBL();

  virtual void update_params(); // Copy from KIM-published parameters to internal.

  KIMParamsZBL kim_params_zbl; // ZBL parameters published to KIM, see
                               // above why we keep two copies.

 protected:
  Array2D<ParamsZBL2> params_zbl_2; // n_spec*n_spec array of ZBL parameters

  void read_params(std::istream&, std::map<std::string,int>,
                   double, double, double);
  void prepare_params();

  virtual double repulsive(double, double, double, int, int,
                           bool, double&) const;
  virtual double ters_fa(double, double, int, int) const;
  virtual double ters_fa_d(double, double, double, int, int) const;

 private:
  const double global_a_0;       // Bohr radius for Coulomb repulsion
  const double global_epsilon_0; // permittivity of vacuum for Coulomb repulsion
  const double global_e;         // proton charge (negative of electron charge)

  // Derived, cached values.
  const double global_e_sq;

  // Fermi-like smoothing from ZBL to Tersoff.
  double F_fermi(double, double, double) const;
  double F_fermi_d(double, double, double) const;
};

}

#endif
