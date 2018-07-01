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
   Contributing author: Aidan Thompson (SNL)

   Modified for use with KIM by Tobias Brink (2012,2013,2014,2017,2018).

   process_dEdr support added by Mingjian Wen (2018)
------------------------------------------------------------------------- */

#include "pair_tersoff.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

#include <KIM_API_status.h>

using namespace model_driver_Tersoff;
using namespace std;

/* ---------------------------------------------------------------------- */

PairTersoff::PairTersoff(string parameter_file,
                         int n_spec,
                         map<string,int> type_map,
                         // Conversion factors.
                         double energy_conv,
                         double length_conv,
                         double inv_length_conv,
                         // KIM indices.
                         const KimIndices& ki
                         )
  : kim_indices(ki), kim_params(n_spec), n_spec(n_spec),
    params2(n_spec, n_spec), params3(n_spec, n_spec, n_spec)
{
  // Prepare index -> element name mapping.
  for (map<string,int>::const_iterator i = type_map.begin();
       i != type_map.end();
       ++i) {
    to_spec[i->second] = i->first;
  }
  // Read parameter file.
  std::fstream f(parameter_file.c_str(), std::ios_base::in);
  read_params(f, type_map, energy_conv, length_conv, inv_length_conv);
}


PairTersoff::~PairTersoff()
{
}


/* Helper to check if an atom is a ghost atom --------------------------- */
bool is_ghost(KIM_API_model& kim_model, int j) {
  int num_neigh = 0;
  int jj;          // unused
  int* neighbors;  // unused
  double* distvec; // unused
  int error =
    kim_model.get_neigh(1, j, &jj, &num_neigh, &neighbors, &distvec);
  if (error < KIM_STATUS_OK) {
    kim_model.report_error(__LINE__, __FILE__,
                           "KIM_API_get_neigh",
                           error);
    throw runtime_error("compute: Error in KIM_API_get_neigh");
  }
  return num_neigh == 0;
}

/* ---------------------------------------------------------------------- */

void PairTersoff::compute(KIM_API_model& kim_model,
                          int n_atoms, // Actual number of atoms
                          const int* atom_types,
                          const Array2D<double>& atom_coords,
                          double* energy, double* atom_energy,
                          Array2D<double>* forces,
                          bool compute_process_dEdr) const
{
  int ii;          // Iteration over all atoms.
  int error;       // KIM error code.
  int n_neigh;     // Number of neighbors of i.
  int* neighbors;  // The indices of the neighbors.
  double* distvec; // Rij vectors. TODO: remove     
  const bool eflag = energy || atom_energy; // Calculate energy?

  // If requested, reset energy.
  if (energy)
    *energy = 0.0;
  if (atom_energy)
    for (int i = 0; i != n_atoms; ++i) {
      atom_energy[i] = 0.0;
    }

  // Reset forces.
  if (forces)
    for (int i = 0; i != n_atoms; ++i) {
      (*forces)(i, 0) = 0.0;
      (*forces)(i, 1) = 0.0;
      (*forces)(i, 2) = 0.0;
    }

  // loop over full neighbor list of my atoms

  ii = 0;
  while (ii != n_atoms) {
    int i;
    // Get neighbors.
    error =
      kim_model.get_neigh(1, // locator mode
                          // The central atom. Increment after.
                          ii++,
                          // Output.
                          &i, // The central atom.
                          &n_neigh, // Number of neighbors
                          &neighbors, // The neighbor indices
                          &distvec // Rij
                          );
    if (error < KIM_STATUS_OK) {
      kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_neigh",
                             error);
      throw runtime_error("compute: Error in KIM_API_get_neigh");
    }

    const int itype = atom_types[i];
    const double xtmp = atom_coords(i,0);
    const double ytmp = atom_coords(i,1);
    const double ztmp = atom_coords(i,2);

    for (int jj = 0; jj != n_neigh; ++jj) {
      int j = neighbors[jj];
      const int jtype = atom_types[j];

      const double delr_ij[3] = { atom_coords(j,0) - xtmp,
                                  atom_coords(j,1) - ytmp,
                                  atom_coords(j,2) - ztmp };
      const double rsq_ij = delr_ij[0]*delr_ij[0]
                          + delr_ij[1]*delr_ij[1]
                          + delr_ij[2]*delr_ij[2];

      const double cutsq = params2(itype,jtype).cutsq;
      if (rsq_ij > cutsq) continue;

      const double r_ij = sqrt(rsq_ij);

      const double R = params2(itype,jtype).R;
      const double D = params2(itype,jtype).D;
      const double fc_ij = ters_fc(r_ij, R, D); // Value of the cutoff function.
      const double dfc_ij = ters_fc_d(r_ij, R, D); // Derivative of fc_ij.

      // two-body interactions, skip half of them unless j is a ghost atom
      const bool j_ghost = is_ghost(kim_model, j);
      if (j_ghost || (i < j)) {
        const double lam1 = params2(itype,jtype).lam1;
        const double A = params2(itype,jtype).A;

        double evdwl; // Particle energy.
        const double fpair =
          repulsive(r_ij, fc_ij, dfc_ij, lam1, A, eflag, evdwl);

        const double half_prefactor = j_ghost ? 0.5 : 1.0;

        if (energy)
          *energy += half_prefactor * evdwl;

        if (atom_energy) {
          atom_energy[i] += 0.5 * evdwl;
          if (!j_ghost)
            atom_energy[j] += 0.5 * evdwl;
        }

        if (forces || compute_process_dEdr) {
          const double fx = delr_ij[0]*fpair;
          const double fy = delr_ij[1]*fpair;
          const double fz = delr_ij[2]*fpair;

          if (forces) {
            (*forces)(i,0) -= half_prefactor*fx;
            (*forces)(i,1) -= half_prefactor*fy;
            (*forces)(i,2) -= half_prefactor*fz;
            (*forces)(j,0) += half_prefactor*fx;
            (*forces)(j,1) += half_prefactor*fy;
            (*forces)(j,2) += half_prefactor*fz;
          }

          if (compute_process_dEdr) {
            run_process_dEdr(&kim_model,
                             -half_prefactor*fpair*r_ij, // dEdr
                             r_ij, delr_ij, i, j,
                             __LINE__, __FILE__);
          }
        }
      }

      // three-body interactions

      // accumulate bondorder zeta for each i-j interaction via loop over k

      double zeta_ij = 0.0;

      for (int kk = 0; kk != n_neigh; ++kk) {
        if (jj == kk) continue;
        int k = neighbors[kk];
        const int ktype = atom_types[k];

        const double delr_ik[3] = { atom_coords(k,0) - xtmp,
                                    atom_coords(k,1) - ytmp,
                                    atom_coords(k,2) - ztmp };
        const double rsq_ik = delr_ik[0]*delr_ik[0]
                            + delr_ik[1]*delr_ik[1]
                            + delr_ik[2]*delr_ik[2];
        const double cutsq = params3(itype,jtype,ktype).cutsq;
        if (rsq_ik > cutsq) continue;

        const double r_ik = sqrt(rsq_ik);

        const double R = params3(itype,jtype,ktype).R;
        const double D = params3(itype,jtype,ktype).D;
        const int m = params3(itype,jtype,ktype).m;
        const double lam3 = params3(itype,jtype,ktype).lam3;
        const double gamma = params3(itype,jtype,ktype).gamma;
        const double h = params3(itype,jtype,ktype).h;
        const double c2 = params3(itype,jtype,ktype).c2;
        const double d2 = params3(itype,jtype,ktype).d2;
        const double c2_d2 = params3(itype,jtype,ktype).c2_d2;

        zeta_ij += zeta(r_ij,r_ik,m,lam3,R,D,gamma,c2,d2,c2_d2,h,
                        delr_ij,delr_ik);
      }

      // pairwise force due to zeta

      const double B = params2(itype,jtype).B;
      const double lam2 = params2(itype,jtype).lam2;
      const double beta = params2(itype,jtype).beta;
      const double n = params2(itype,jtype).n;
      const double* n_precomp = params2(itype,jtype).n_precomp;

      double prefactor; // -0.5 * fa * ∇bij
      double evdwl; // Particle energy.
      const double fzeta =
        force_zeta(r_ij, fc_ij, dfc_ij, zeta_ij, B, lam2, beta, n, n_precomp,
                   prefactor, eflag, evdwl);

      if (energy)
        *energy += evdwl;

      if (atom_energy) {
        const double en = 0.5 * evdwl;
        atom_energy[i] += en;
        atom_energy[j] += en;
      }

      if (forces || compute_process_dEdr) {
        const double fx = delr_ij[0]*fzeta;
        const double fy = delr_ij[1]*fzeta;
        const double fz = delr_ij[2]*fzeta;

        if (forces) {
          (*forces)(i,0) += fx;
          (*forces)(i,1) += fy;
          (*forces)(i,2) += fz;
          (*forces)(j,0) -= fx;
          (*forces)(j,1) -= fy;
          (*forces)(j,2) -= fz;
        }

        if (compute_process_dEdr) {
          run_process_dEdr(&kim_model,
                           fzeta*r_ij, // dEdr
                           r_ij, delr_ij, i, j,
                           __LINE__, __FILE__);
        }
      }


      // attractive term via loop over k

      if (forces || compute_process_dEdr) {
        for (int kk = 0; kk != n_neigh; ++kk) {
          if (jj == kk) continue;
          int k = neighbors[kk];
          const int ktype = atom_types[k];

          const double delr_ik[3] = { atom_coords(k,0) - xtmp,
                                      atom_coords(k,1) - ytmp,
                                      atom_coords(k,2) - ztmp };
          const double rsq_ik = delr_ik[0]*delr_ik[0]
                              + delr_ik[1]*delr_ik[1]
                              + delr_ik[2]*delr_ik[2];
          const double cutsq = params3(itype,jtype,ktype).cutsq;
          if (rsq_ik > cutsq) continue;

          const double r_ik = sqrt(rsq_ik);

          // r_jk
          double delr_jk[3];
          delr_jk[0] = delr_ik[0] - delr_ij[0];
          delr_jk[1] = delr_ik[1] - delr_ij[1];
          delr_jk[2] = delr_ik[2] - delr_ij[2];
          const double rsq_jk = delr_jk[0]*delr_jk[0]
                              + delr_jk[1]*delr_jk[1]
                              + delr_jk[2]*delr_jk[2];
          const double r_jk = sqrt(rsq_jk);

          const double R = params3(itype,jtype,ktype).R;
          const double D = params3(itype,jtype,ktype).D;
          const int m = params3(itype,jtype,ktype).m;
          const double lam3 = params3(itype,jtype,ktype).lam3;
          const double gamma = params3(itype,jtype,ktype).gamma;
          const double h = params3(itype,jtype,ktype).h;
          const double c2 = params3(itype,jtype,ktype).c2;
          const double d2 = params3(itype,jtype,ktype).d2;
          const double c2_d2 = params3(itype,jtype,ktype).c2_d2;

          // dEdr_ij not normalized by r_ij, same for dEdr_ik and dEdr_jk
          // note prefactor = -0.5 * fa * ∇bij, as defined above
          double dEdr_ij, dEdr_ik, dEdr_jk;
          attractive(-prefactor,
                     R, D, m, lam3, gamma, c2, d2, c2_d2, h,
                     r_ij, r_ik, r_jk,
                     dEdr_ij, dEdr_ik, dEdr_jk);

          if (forces) {
            for (int dim = 0; dim < 3; ++dim) {
              const double pair_ij = dEdr_ij*delr_ij[dim]/r_ij;
              const double pair_ik = dEdr_ik*delr_ik[dim]/r_ik;
              const double pair_jk = dEdr_jk*delr_jk[dim]/r_jk;
              (*forces)(i,dim) +=  pair_ij + pair_ik;
              (*forces)(j,dim) += -pair_ij + pair_jk;
              (*forces)(k,dim) += -pair_ik - pair_jk;
            }
          }

          if (compute_process_dEdr) {
            run_process_dEdr(&kim_model, dEdr_ij, r_ij, delr_ij, i, j,
                             __LINE__, __FILE__);
            run_process_dEdr(&kim_model, dEdr_ik, r_ik, delr_ik, i, k,
                             __LINE__, __FILE__);
            run_process_dEdr(&kim_model, dEdr_jk, r_jk, delr_jk, j, k,
                             __LINE__, __FILE__);
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

/* Read parameters.

   Reads the parameters from an input file and stores
   them. Automatically calls prepare_params() to check validity of
   parameters and pre-compute some values.
 */
void PairTersoff::read_params(istream& infile, std::map<string,int> type_map,
                              double energy_conv,
                              double length_conv,
                              double inv_length_conv) {
  // Strip comments lines.
  stringstream buffer;
  string line;
  while(getline(infile, line))
    buffer << line.substr(0, line.find('#'));
  // Read in parameters.
  Array3D<bool> got_interaction(n_spec,n_spec,n_spec);
  got_interaction = false;
  Params2 temp_params2;
  Params3 temp_params3;
  double m; // m is an integer but some input files use floating point
            // notation, so we need to read into a double variable,
            // otherwise the C++ standard library chokes on the input.
  string type_i, type_j, type_k;
  double c, d; // Those are only stored in the published KIM parameters
               // since they are not needed for computation, we use
               // precomputed c², d² and c²/d² instead.
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
                >> temp_params2.A) {
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
}

void PairTersoff::update_params() {
  kim_params.to_params(params2, params3);
  prepare_params();
}

/* Check parameters and pre-compute values. */
void PairTersoff::prepare_params() {
  max_cutoff = 0.0;
  for (int i = 0; i != n_spec; ++i) {
    string type_i = to_spec.at(i);
    for (int j = 0; j != n_spec; ++j) {
      string type_j = to_spec.at(j);
      Params2& temp_params2 = params2(i,j);
      if (temp_params2.n < 0)
        throw runtime_error("Parameter n ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than zero.");
      if (temp_params2.beta < 0)
        throw runtime_error("Parameter beta ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than zero.");
      if (temp_params2.lam1 < 0)
        throw runtime_error("Parameter lambda1 ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than zero.");
      if (temp_params2.lam2 < 0)
        throw runtime_error("Parameter lambda2 ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than zero.");
      if (temp_params2.A < 0)
        throw runtime_error("Parameter A ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than zero.");
      if (temp_params2.B < 0)
        throw runtime_error("Parameter B ("
                            + type_i +
                            "-" + type_j +
                            ") may not be smaller than zero.");
      // Pre-compute values.
      const double n2 = 2.0 * temp_params2.n;
      const double n_r = -1.0/temp_params2.n;
      temp_params2.n_precomp[0] = pow(n2 * 1e-16, n_r);
      temp_params2.n_precomp[1] = pow(n2 * 1e-8, n_r);
      temp_params2.n_precomp[2] = 1.0 / temp_params2.n_precomp[1];
      temp_params2.n_precomp[3] = 1.0 / temp_params2.n_precomp[0];
      for (int k = 0; k != n_spec; ++k) {
        string type_k = to_spec.at(k);
        Params3& temp_params3 = params3(i,j,k);
        // Check values of parameters.
        if (kim_params.c(i,j,k) < 0)
          throw runtime_error("Parameter c ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (kim_params.d(i,j,k) < 0)
          throw runtime_error("Parameter d ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params3.m != 1 && temp_params3.m != 3)
          throw runtime_error("Parameter m ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") must be one or three.");
        if (temp_params3.R < 0)
          throw runtime_error("Parameter R ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params3.D < 0)
          throw runtime_error("Parameter D ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params3.D > temp_params3.R)
          throw runtime_error("Parameter D ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") must be smaller than R.");
        if (temp_params3.gamma < 0)
          throw runtime_error("Parameter gamma ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        // Cutoff.
        double cut =
          temp_params3.R + temp_params3.D; // max cutoff.
        temp_params3.cutsq = cut * cut; // for fast check if inside cutoff
        if (j == k)
          temp_params2.cutsq = temp_params3.cutsq;
        // Get the cutoff to pass to KIM, which is the biggest cutoff.
        if (cut > max_cutoff)
          max_cutoff = cut;
        // Pre-compute values.
        temp_params3.c2 = kim_params.c(i,j,k) * kim_params.c(i,j,k);
        temp_params3.d2 = kim_params.d(i,j,k) * kim_params.d(i,j,k);
        temp_params3.c2_d2 = temp_params3.c2 / temp_params3.d2;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairTersoff::repulsive(double r, double fc, double fc_d,
                              double lam1, double A,
                              bool eflag, double &eng) const
{
  const double tmp_exp = exp(-lam1 * r);
  if (eflag) eng = fc * A * tmp_exp;
  return -A * tmp_exp * (fc_d - fc*lam1) / r;
}

/* ---------------------------------------------------------------------- */

double PairTersoff::zeta(double rij, double rik,
                         int m, double lam3, double R, double D,
                         double gamma, double c2, double d2, double c2_d2,
                         double h,
                         const double* delrij, const double* delrik) const
{
  const double costheta = (delrij[0]*delrik[0]
                           + delrij[1]*delrik[1]
                           + delrij[2]*delrik[2]
                           ) / (rij * rik);

  double arg;
  if (m == 3) arg = pow(lam3 * (rij-rik), 3);
  else arg = lam3 * (rij-rik); // m == 1

  double ex_delr;
  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc(rik,R,D) * ters_gijk(costheta,gamma,c2,d2,c2_d2,h) * ex_delr;
}

/* ---------------------------------------------------------------------- */

double PairTersoff::force_zeta(double r, double fc, double fc_d, double zeta_ij,
                               double B, double lam2,
                               double beta, double n,
                               const double n_precomp[4],
                               double &prefactor,
                               bool eflag, double &eng) const
{
  const double fa = ters_fa(r, fc, B, lam2);
  const double fa_d = ters_fa_d(r, fc, fc_d, B, lam2);
  const double bij = ters_bij(zeta_ij, beta, n, n_precomp);
  prefactor = -0.5*fa * ters_bij_d(zeta_ij, beta, n, n_precomp);
  if (eflag) eng = 0.5*bij*fa;
  return 0.5*bij*fa_d / r;
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairTersoff::attractive(double prefactor,
                             double R, double D, int m, double lam3,
                             double gamma, double c2, double d2, double c2_d2,
                             double h,
                             double rij, double rik, double rjk,
                             double &drij, double &drik, double &drjk) const
{

  ters_zetaterm_d(prefactor,
                  R, D, m, lam3, gamma, c2, d2, c2_d2, h,
                  rij,rik,rjk,drij,drik,drjk);
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fc(double r, double R, double D) const
{
  if (r < R-D) return 1.0;
  if (r > R+D) return 0.0;
  return 0.5*(1.0 - sin(pi_2*(r - R)/D));
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fc_d(double r, double R, double D) const
{
  if (r < R-D) return 0.0;
  if (r > R+D) return 0.0;
  return -(pi_4/D) * cos(pi_2*(r - R)/D);
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fa(double r, double fc,
                            double B, double lam2) const
{
  if (fc == 0.0) return 0.0;
  return -B * exp(-lam2 * r) * fc;
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fa_d(double r, double fc, double fc_d,
                              double B, double lam2) const
{
  if (fc == 0.0) return 0.0;
  return B * exp(-lam2 * r) * (lam2 * fc - fc_d);
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_bij(double zeta, double beta, double n,
                             const double n_precomp[4]) const
{
  const double tmp = beta * zeta;
  if (tmp > n_precomp[0]) return 1.0/sqrt(tmp);
  if (tmp > n_precomp[1]) return (1.0 - pow(tmp,-n) / (2.0*n))/sqrt(tmp);
  if (tmp < n_precomp[3]) return 1.0;
  if (tmp < n_precomp[2]) return 1.0 - pow(tmp,n)/(2.0*n);
  return pow(1.0 + pow(tmp,n), -1.0/(2.0*n));
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_bij_d(double zeta, double beta, double n,
                               const double n_precomp[4]) const
{
  const double tmp = beta * zeta;
  if (tmp > n_precomp[0]) return beta * -0.5*pow(tmp,-1.5);
  if (tmp > n_precomp[1]) return beta * (-0.5*pow(tmp,-1.5) *
                                         (1.0 - (1.0 +  1.0/(2.0*n)) *
                                          pow(tmp,-n)));
  if (tmp < n_precomp[3]) return 0.0;
  if (tmp < n_precomp[2]) return -0.5*beta * pow(tmp, n-1.0);

  const double tmp_n = pow(tmp, n);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*n)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void PairTersoff::ters_zetaterm_d(double prefactor,
                                  double R, double D, int m, double lam3,
                                  double gamma, double c2, double d2,
                                  double c2_d2, double h,
                                  double rij, double rik, double rjk,
                                  double &drij, double &drik, double &drjk) const
{

  const double fc = ters_fc(rik,R,D);
  const double dfc = ters_fc_d(rik,R,D);

  double ex_delr,ex_delr_d,tmp;
  if (m == 3) tmp = pow(lam3 * (rij-rik), 3);
  else tmp = lam3 * (rij-rik); // m == 1

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (m == 3)
    ex_delr_d = 3.0*pow(lam3, 3) * pow(rij-rik, 2)*ex_delr;
  else ex_delr_d = lam3 * ex_delr; // m == 1

  double cos_theta, dcosdrij,dcosdrik,dcosdrjk;
  costheta_d(rij,rik,rjk, cos_theta, dcosdrij, dcosdrik, dcosdrjk);

  const double gijk = ters_gijk(cos_theta, gamma, c2, d2, c2_d2, h);
  const double gijk_d = ters_gijk_d(cos_theta, gamma, c2, d2, h);


  // compute the derivative wrt rij, rik, rjk
  drij = (fc*gijk_d*ex_delr*dcosdrij + fc*gijk*ex_delr_d) * prefactor;
  drik = (dfc*gijk*ex_delr + fc*gijk_d*ex_delr*dcosdrik - fc*gijk*ex_delr_d) * prefactor;
  drjk = (fc*gijk_d*ex_delr*dcosdrjk) * prefactor;

}

/* ---------------------------------------------------------------------- */

void PairTersoff::costheta_d(double rij, double rik, double rjk, double &cos_ijk,
                             double &drij, double &drik, double &drjk) const
{
  const double rijsq = rij*rij;
  const double riksq = rik*rik;
  const double rjksq = rjk*rjk;

  // cos_ijk, (i is the apex atom)
  cos_ijk = (rijsq + riksq - rjksq)/(2*rij*rik);
  // dcos_ijk/drij, dcos_ijk/drik, dcos_ijk/drjk
  drij = (rijsq - riksq + rjksq)/(2*rijsq*rik);
  drik = (riksq - rijsq + rjksq)/(2*rij*riksq);
  drjk = -rjk/(rij*rik);
}

