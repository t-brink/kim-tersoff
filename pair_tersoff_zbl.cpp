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
   Contributing author: Aidan Thompson (SNL) - original Tersoff implementation
                        David Farrell (NWU) - ZBL addition

   Modified for use with KIM by Tobias Brink (2020).
------------------------------------------------------------------------- */

#include "pair_tersoff_zbl.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

using namespace model_driver_Tersoff;
using namespace std;

/* ---------------------------------------------------------------------- */

PairTersoffZBL::PairTersoffZBL(const string& parameter_file,
                               int n_spec,
                               map<string,int> type_map,
                               // Conversion factors.
                               double energy_conv,
                               double inv_energy_conv,
                               double length_conv,
                               double inv_length_conv,
                               double charge_conv)
  : PairTersoff(n_spec, type_map),
    kim_params_zbl(n_spec),
    params_zbl_2(n_spec, n_spec),
    global_a_0(0.529 * length_conv),
    global_epsilon_0(0.00552635
                     * charge_conv * charge_conv
                     * inv_energy_conv * inv_length_conv),
    global_e(1.0 * charge_conv),
    global_e_sq(global_e * global_e)
{
  // Read parameter file.
  std::fstream f(parameter_file.c_str(), std::ios_base::in);
  read_params(f, type_map, energy_conv, length_conv, inv_length_conv);
}

PairTersoffZBL::~PairTersoffZBL() {}

/* ---------------------------------------------------------------------- */

/* Read parameters.

   Reads the parameters from an input file and stores
   them. Automatically calls prepare_params() to check validity of
   parameters and pre-compute some values.

   This is a bit of a copy-paste from pair_tersoff.cpp, but this way
   is still easiest.
 */
void PairTersoffZBL::read_params(istream& infile, std::map<string,int> type_map,
                                 double energy_conv,
                                 double length_conv,
                                 double inv_length_conv) {
  // Strip comment lines.
  stringstream buffer;
  string line;
  while(getline(infile, line))
    buffer << line.substr(0, line.find('#'));
  // Read in parameters.
  Array3D<bool> got_interaction(n_spec,n_spec,n_spec);
  got_interaction = false;
  Params2 temp_params2;
  Params3 temp_params3;
  ParamsZBL2 temp_params_zbl_2;
  double m; // m is an integer but some input files use floating point
            // notation, so we need to read into a double variable,
            // otherwise the C++ standard library chokes on the input.
  string type_i, type_j, type_k;
  double c, d; // Those are only stored in the published KIM parameters
               // since they are not needed for computation, we use
               // precomputed c², d² and c²/d² instead.
  double Z_i, Z_j; // Equivalent for these.
  while (buffer >> type_i
                >> type_j
                >> type_k
                >> m
                >> temp_params3.gamma
                >> temp_params3.lam3
                >> c
                >> d
                >> temp_params3.h // costheta0
                >> temp_params2.n
                >> temp_params2.beta
                >> temp_params2.lam2
                >> temp_params2.B
                >> temp_params3.R
                >> temp_params3.D
                >> temp_params2.lam1
                >> temp_params2.A
                >> Z_i
                >> Z_j
                >> temp_params_zbl_2.ZBLcut
                >> temp_params_zbl_2.ZBLexpscale) {
    // Convert m to integer.
    temp_params3.m = m;
    if (abs(m - temp_params3.m) > 1e-8)
      throw runtime_error("m must be an integer");
    // Unit conversion.
    temp_params2.A *= energy_conv;
    temp_params2.B *= energy_conv;
    temp_params2.lam1 *= inv_length_conv;
    temp_params2.lam2 *= inv_length_conv;
    temp_params3.lam3 *= inv_length_conv;
    temp_params3.R *= length_conv;
    temp_params3.D *= length_conv;
    temp_params_zbl_2.ZBLcut *= length_conv;
    temp_params_zbl_2.ZBLexpscale *= inv_length_conv;
    // Get the atom type indices.
    std::map<std::string,int>::const_iterator it;
    int i,j,k;
    it = type_map.find(type_i);
    if (it != type_map.end())
      i = it->second;
    else
      throw runtime_error("Unknown species: " + type_i);
    it = type_map.find(type_j);
    if (it != type_map.end())
      j = it->second;
    else
      throw runtime_error("Unknown species: " + type_j);
    it = type_map.find(type_k);
    if (it != type_map.end())
      k = it->second;
    else
      throw runtime_error("Unknown species: " + type_k);
    // Check if parameters were given twice.
    if (got_interaction(i,j,k))
      throw runtime_error("Interaction " + type_i +
                          "-" + type_j +
                          "-" + type_k +
                          " is defined twice!");
    // All OK, store.
    got_interaction(i,j,k) = true;
    if (j == k) {
      temp_params2.R = temp_params3.R;
      temp_params2.D = temp_params3.D;
      params2(i,j) = temp_params2;
      params_zbl_2(i,j) = temp_params_zbl_2;
      kim_params_zbl.Z_i(i,j) = Z_i;
      kim_params_zbl.Z_j(i,j) = Z_j;
    }
    params3(i,j,k) = temp_params3;
    kim_params.c(i,j,k) = c;
    kim_params.d(i,j,k) = d;
  }
  if (!got_interaction.all())
    throw runtime_error("Not all interactions were set!");
  // Check parameters and pre-compute values.
  prepare_params();
  // Copy to KIM-published data.
  kim_params.from_params(params2, params3);
  kim_params_zbl.from_params(params_zbl_2);
}

void PairTersoffZBL::update_params() {
  kim_params.to_params(params2, params3);
  kim_params_zbl.to_params(params_zbl_2);
  prepare_params();
}

void PairTersoffZBL::prepare_params() {
  PairTersoff::prepare_params();

  for (int i = 0; i != n_spec; ++i) {
    string type_i = to_spec.at(i);
    for (int j = 0; j != n_spec; ++j) {
      string type_j = to_spec.at(j);
      ParamsZBL2& temp_params_zbl_2 = params_zbl_2(i,j);
      if (kim_params_zbl.Z_i(i,j) < 1)
        throw runtime_error("Parameter Z_i ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than one.");
      if (kim_params_zbl.Z_j(i,j) < 1)
        throw runtime_error("Parameter Z_j ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than one.");
      if (temp_params_zbl_2.ZBLcut < 0)
        throw runtime_error("Parameter ZBLcut ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than one.");
      if (temp_params_zbl_2.ZBLexpscale < 0)
        throw runtime_error("Parameter ZBLexpscale ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than one.");
      // Pre-compute values.
      temp_params_zbl_2.a =
        0.8854 * global_a_0 / (pow(kim_params_zbl.Z_i(i,j), 0.23)
                               + pow(kim_params_zbl.Z_j(i,j), 0.23));
      temp_params_zbl_2.premult =
        (kim_params_zbl.Z_i(i,j) * kim_params_zbl.Z_j(i,j) * global_e_sq)
        /
        (4.0 * pi * global_epsilon_0);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairTersoffZBL::repulsive(double r, double fc, double fc_d,
                                 int itype, int jtype,
                                 bool eflag, double &eng) const
{
  const double lam1 = params2(itype,jtype).lam1;
  const double A = params2(itype,jtype).A;

  // Tersoff repulsive portion

  const double tmp_exp = exp(-lam1 * r);
  const double eng_ters = fc * A * tmp_exp;
  const double fforce_ters = A * tmp_exp * (fc_d - fc*lam1) / r;

  // ZBL repulsive portion

  const double a_ij    = params_zbl_2(itype,jtype).a;
  const double premult = params_zbl_2(itype,jtype).premult;
  const double ZBLexpscale = params_zbl_2(itype,jtype).ZBLexpscale;
  const double ZBLcut = params_zbl_2(itype,jtype).ZBLcut;

  const double r_over_a = r / a_ij;
  const double exp1 = 0.1818  * exp(-3.2    * r_over_a);
  const double exp2 = 0.5099  * exp(-0.9423 * r_over_a);
  const double exp3 = 0.2802  * exp(-0.4029 * r_over_a);
  const double exp4 = 0.02817 * exp(-0.2016 * r_over_a);
  const double phi = exp1 + exp2 + exp3 + exp4;
  const double dphi = (1.0/a_ij) * (-3.2   * exp1 -
                                    0.9423 * exp2 -
                                    0.4029 * exp3 -
                                    0.2016 * exp4);
  const double eng_ZBL = premult * (1.0/r) * phi;
  const double fforce_ZBL = premult * -phi / pow2(r) + premult * dphi / r;

  // combine two parts with smoothing by Fermi-like function

  const double f_F = F_fermi(r, ZBLexpscale, ZBLcut);
  const double f_F_d = F_fermi_d(r, ZBLexpscale, ZBLcut);

  if (eflag) {
    eng = (1.0 - f_F) * eng_ZBL + f_F * eng_ters;
  }

  return -(-f_F_d * eng_ZBL +
           (1.0 - f_F) * fforce_ZBL +
           f_F_d * eng_ters + f_F * fforce_ters) / r;
}

/* ---------------------------------------------------------------------- */

double PairTersoffZBL::ters_fa(double r, double fc,
                               int itype, int jtype) const
{
  if (fc == 0.0) return 0.0;

  const double B = params2(itype,jtype).B;
  const double lam2 = params2(itype,jtype).lam2;
  const double ZBLexpscale = params_zbl_2(itype,jtype).ZBLexpscale;
  const double ZBLcut = params_zbl_2(itype,jtype).ZBLcut;

  return -B * exp(-lam2 * r) * fc * F_fermi(r, ZBLexpscale, ZBLcut);
}

/* ---------------------------------------------------------------------- */

double PairTersoffZBL::ters_fa_d(double r, double fc, double fc_d,
                                 int itype, int jtype) const
{
  if (fc == 0.0) return 0.0;

  const double B = params2(itype,jtype).B;
  const double lam2 = params2(itype,jtype).lam2;
  const double ZBLexpscale = params_zbl_2(itype,jtype).ZBLexpscale;
  const double ZBLcut = params_zbl_2(itype,jtype).ZBLcut;

  const double f_F = F_fermi(r, ZBLexpscale, ZBLcut);
  const double f_F_d = F_fermi_d(r, ZBLexpscale, ZBLcut);

  return B * exp(-lam2 * r) * (lam2 * fc * f_F - fc_d * f_F - fc * f_F_d);
}

/* ----------------------------------------------------------------------
   Fermi-like smoothing function
------------------------------------------------------------------------- */

double PairTersoffZBL::F_fermi(double r, double ZBLexpscale, double ZBLcut) const
{
  return 1.0 / (1.0 + exp(-ZBLexpscale * (r - ZBLcut)));
}

/* ----------------------------------------------------------------------
   Fermi-like smoothing function derivative with respect to r
------------------------------------------------------------------------- */

double PairTersoffZBL::F_fermi_d(double r, double ZBLexpscale, double ZBLcut) const
{
  return ZBLexpscale * exp(-ZBLexpscale * (r - ZBLcut)) /
         pow2(1.0 + exp(-ZBLexpscale * (r - ZBLcut)));
}
