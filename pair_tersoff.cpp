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

   Modified for use with KIM by Tobias Brink (2012,2013,2014,2017,2018,2019,
   2020,2021).

   process_dEdr support added by Mingjian Wen (2018)
------------------------------------------------------------------------- */

#include "pair_tersoff.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace model_driver_Tersoff;
using namespace std;

/* ---------------------------------------------------------------------- */

PairTersoff::PairTersoff(const string& parameter_file,
                         int n_spec,
                         map<string,int> type_map,
                         // Conversion factors.
                         double energy_conv,
                         double, // unused inv_energy_conv
                         double length_conv,
                         double inv_length_conv,
                         double // unused charge_conv
                         )
  : kim_params(n_spec), n_spec(n_spec),
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

PairTersoff::PairTersoff(int n_spec,
                         map<string,int> type_map)
  : kim_params(n_spec), n_spec(n_spec),
    params2(n_spec, n_spec), params3(n_spec, n_spec, n_spec)
{
  // Prepare index -> element name mapping.
  for (map<string,int>::const_iterator i = type_map.begin();
       i != type_map.end();
       ++i) {
    to_spec[i->second] = i->first;
  }
}


PairTersoff::~PairTersoff()
{
}

/* ---------------------------------------------------------------------- */

void PairTersoff::compute(const KIM::ModelComputeArguments&
                                            model_compute_arguments,
                          int n_atoms, // Actual number of atoms
                          const int * const atom_types,
                          const int * const contributing,
                          const Array2D<const double>& atom_coords,
                          double* energy, double* atom_energy,
                          Array2D<double>* forces,
                          double* virial,
                          Array2D<double>* particle_virial,
                          bool compute_process_dEdr) const
{
  int error;       // KIM error code.
  int n_neigh;     // Number of neighbors of i.
  const int * neighbors;  // The indices of the neighbors.
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

  // Reset virial.
  if (virial)
    for (int i = 0; i != 6; ++i)
      virial[i] = 0.0;
  if (particle_virial)
    for (int i = 0; i != n_atoms; ++i)
      for (int j = 0; j != 6; ++j)
        (*particle_virial)(i, j) = 0.0;

  // loop over full neighbor list of my atoms

  for (int i = 0; i != n_atoms; ++i) {
    // Skip central ghost atoms.
    if (!contributing[i]) continue;

    // Get neighbors.
    error =
      model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
    if (error) {
      throw runtime_error("Error in "
                          "KIM::ModelComputeArguments.GetNeighborList");
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
      if (!contributing[j] || (i < j)) {
        double evdwl; // Particle energy.
        const double fpair =
          repulsive(r_ij, fc_ij, dfc_ij, itype, jtype, eflag, evdwl);

        const double half_prefactor = contributing[j] ? 1.0 : 0.5;

        if (energy)
          *energy += half_prefactor * evdwl;

        if (atom_energy) {
          atom_energy[i] += 0.5 * evdwl;
          if (contributing[j])
            atom_energy[j] += 0.5 * evdwl;
        }

        if (forces || virial || particle_virial || compute_process_dEdr) {
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

          if (virial || particle_virial) {
            const double vxx = delr_ij[0] * fx;
            const double vyy = delr_ij[1] * fy;
            const double vzz = delr_ij[2] * fz;
            const double vyz = delr_ij[1] * fz;
            const double vxz = delr_ij[0] * fz;
            const double vxy = delr_ij[0] * fy;

            if (virial) {
              virial[0] -= half_prefactor * vxx;
              virial[1] -= half_prefactor * vyy;
              virial[2] -= half_prefactor * vzz;
              virial[3] -= half_prefactor * vyz;
              virial[4] -= half_prefactor * vxz;
              virial[5] -= half_prefactor * vxy;
            }

            if (particle_virial) {
              (*particle_virial)(i,0) -= 0.5 * vxx;
              (*particle_virial)(i,1) -= 0.5 * vyy;
              (*particle_virial)(i,2) -= 0.5 * vzz;
              (*particle_virial)(i,3) -= 0.5 * vyz;
              (*particle_virial)(i,4) -= 0.5 * vxz;
              (*particle_virial)(i,5) -= 0.5 * vxy;
              if (contributing[j]) {
                (*particle_virial)(j,0) -= 0.5 * vxx;
                (*particle_virial)(j,1) -= 0.5 * vyy;
                (*particle_virial)(j,2) -= 0.5 * vzz;
                (*particle_virial)(j,3) -= 0.5 * vyz;
                (*particle_virial)(j,4) -= 0.5 * vxz;
                (*particle_virial)(j,5) -= 0.5 * vxy;
              }
            }
          }

          if (compute_process_dEdr) {
            run_process_dEdr(model_compute_arguments,
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

      double prefactor; // -0.5 * fa * ∇bij
      double evdwl; // Particle energy.
      const double fzeta =
        force_zeta(r_ij, fc_ij, dfc_ij, zeta_ij, itype, jtype,
                   prefactor, eflag, evdwl);

      if (energy)
        *energy += evdwl;

      if (atom_energy) {
        // Non-symmetric assignment of bond energy.
        //
        // There are basically two ways to distribute this three-body
        // term to particle energies. The first one, used in Tersoff's
        // papers, is to write the sum over atom pairs (i.e., sum over
        // i != j) and have a non-symmetric distribution. Brenner [in
        // PRB 42, 9458–9471 (1990)], among others, favored to use a
        // sum over bonds (i < j) and an averaged b_ij term so that
        // the energy is distributed symmetrically. This is useful for
        // Brenner's overbinding correction. Regarding the physics,
        // the choice is arbitrary and does not affect the forces.
        //
        // KIM requires that non-contributing particles have zero
        // energy. This is because these particles can be used to
        // implement boundaries where the non-contributing particle
        // does not correspond to a real particle somewhere else (as
        // in periodic boundary conditions) and the energy assigned to
        // it would be lost. Then, the sum of particle energies is not
        // equal to the total energy. We are therefore forced to use
        // the non-symmetric energy assignment.
	atom_energy[i] += evdwl;
      }

      if (forces || virial || particle_virial || compute_process_dEdr) {
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

        if (virial || particle_virial) {
          const double vxx = delr_ij[0] * fx;
          const double vyy = delr_ij[1] * fy;
          const double vzz = delr_ij[2] * fz;
          const double vyz = delr_ij[1] * fz; // yz
          const double vxz = delr_ij[0] * fz; // xz
          const double vxy = delr_ij[0] * fy; // xy

          if (virial) {
            virial[0] += vxx;
            virial[1] += vyy;
            virial[2] += vzz;
            virial[3] += vyz;
            virial[4] += vxz;
            virial[5] += vxy;
          }

          if (particle_virial) {
            const double vxx2 = 0.5 * vxx;
            const double vyy2 = 0.5 * vyy;
            const double vzz2 = 0.5 * vzz;
            const double vyz2 = 0.5 * vyz;
            const double vxz2 = 0.5 * vxz;
            const double vxy2 = 0.5 * vxy;

            (*particle_virial)(i,0) += vxx2;
            (*particle_virial)(i,1) += vyy2;
            (*particle_virial)(i,2) += vzz2;
            (*particle_virial)(i,3) += vyz2;
            (*particle_virial)(i,4) += vxz2;
            (*particle_virial)(i,5) += vxy2;

            (*particle_virial)(j,0) += vxx2;
            (*particle_virial)(j,1) += vyy2;
            (*particle_virial)(j,2) += vzz2;
            (*particle_virial)(j,3) += vyz2;
            (*particle_virial)(j,4) += vxz2;
            (*particle_virial)(j,5) += vxy2;
          }
        }

        if (compute_process_dEdr) {
          run_process_dEdr(model_compute_arguments,
                           fzeta*r_ij, // dEdr
                           r_ij, delr_ij, i, j,
                           __LINE__, __FILE__);
        }
      }


      // attractive term via loop over k

      if (forces || virial || particle_virial || compute_process_dEdr) {
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

          double delr_jk[3];
          double r_jk;
          double fi[3], fj[3], fk[3];
          // dEdr_ij not normalized by r_ij, same for dEdr_ik and dEdr_jk
          // note prefactor = -0.5 * fa * ∇bij, as defined above
          double dEdr_ij, dEdr_ik, dEdr_jk;
          attractive(prefactor,
                     R, D, m, lam3, gamma, c2, d2, c2_d2, h,
                     r_ij, r_ik, rsq_ij, rsq_ik, delr_ij, delr_ik,
                     fi, fj, fk,
                     r_jk, delr_jk,
                     dEdr_ij, dEdr_ik, dEdr_jk,
                     forces || virial || particle_virial,
                     compute_process_dEdr);

          if (forces) {
            for (int dim = 0; dim < 3; ++dim) {
              (*forces)(i,dim) += fi[dim];
              (*forces)(j,dim) += fj[dim];
              (*forces)(k,dim) += fk[dim];
            }
          }

          if (virial || particle_virial) {
            const double vxx = delr_ij[0]*fj[0] + delr_ik[0]*fk[0];
            const double vyy = delr_ij[1]*fj[1] + delr_ik[1]*fk[1];
            const double vzz = delr_ij[2]*fj[2] + delr_ik[2]*fk[2];
            const double vyz = delr_ij[1]*fj[2] + delr_ik[1]*fk[2];
            const double vxz = delr_ij[0]*fj[2] + delr_ik[0]*fk[2];
            const double vxy = delr_ij[0]*fj[1] + delr_ik[0]*fk[1];

            if (virial) {
              virial[0] -= vxx;
              virial[1] -= vyy;
              virial[2] -= vzz;
              virial[3] -= vyz;
              virial[4] -= vxz;
              virial[5] -= vxy;
            }

            if (particle_virial) {
              const double vxx3 = 1.0/3.0 * vxx;
              const double vyy3 = 1.0/3.0 * vyy;
              const double vzz3 = 1.0/3.0 * vzz;
              const double vyz3 = 1.0/3.0 * vyz;
              const double vxz3 = 1.0/3.0 * vxz;
              const double vxy3 = 1.0/3.0 * vxy;

              (*particle_virial)(i,0) -= vxx3;
              (*particle_virial)(i,1) -= vyy3;
              (*particle_virial)(i,2) -= vzz3;
              (*particle_virial)(i,3) -= vyz3;
              (*particle_virial)(i,4) -= vxz3;
              (*particle_virial)(i,5) -= vxy3;

              (*particle_virial)(j,0) -= vxx3;
              (*particle_virial)(j,1) -= vyy3;
              (*particle_virial)(j,2) -= vzz3;
              (*particle_virial)(j,3) -= vyz3;
              (*particle_virial)(j,4) -= vxz3;
              (*particle_virial)(j,5) -= vxy3;

              (*particle_virial)(k,0) -= vxx3;
              (*particle_virial)(k,1) -= vyy3;
              (*particle_virial)(k,2) -= vzz3;
              (*particle_virial)(k,3) -= vyz3;
              (*particle_virial)(k,4) -= vxz3;
              (*particle_virial)(k,5) -= vxy3;
            }
          }

          if (compute_process_dEdr) {
            run_process_dEdr(model_compute_arguments,
                             dEdr_ij, r_ij, delr_ij, i, j,
                             __LINE__, __FILE__);
            run_process_dEdr(model_compute_arguments,
                             dEdr_ik, r_ik, delr_ik, i, k,
                             __LINE__, __FILE__);
            run_process_dEdr(model_compute_arguments,
                             dEdr_jk, r_jk, delr_jk, j, k,
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
  // Strip comment lines.
  stringstream buffer;
  string line;
  while(getline(infile, line))
    buffer << line.substr(0, line.find('#')) << endl;
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
    const string type_i = to_spec.at(i);
    for (int j = 0; j != n_spec; ++j) {
      const string type_j = to_spec.at(j);
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
        const string type_k = to_spec.at(k);
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

void PairTersoff::write_params(ofstream& outfile) {
  // Set maximum precision.
  outfile << setprecision(16);
  // Write.
  for (int i = 0; i < n_spec; ++i) {
    const string type_i = to_spec.at(i);
    for (int j = 0; j < n_spec; ++j) {
      const string type_j = to_spec.at(j);
      for (int k = 0; k < n_spec; ++k) {
        const string type_k = to_spec.at(k);
        outfile << type_i << " " << type_j << " " << type_k << " ";
        outfile << kim_params.m(i,j,k) << " ";
        outfile << kim_params.gamma(i,j,k) << " ";
        outfile << kim_params.lam3(i,j,k) << " ";
        outfile << kim_params.c(i,j,k) << " ";
        outfile << kim_params.d(i,j,k) << " ";
        outfile << kim_params.h(i,j,k) << " ";
        if (j == k) {
          // Two-body parameters are only taken from j == k, so in
          // order to make that clear, we write zeros otherwise.
          outfile << kim_params.n(i,j) << " ";
          outfile << kim_params.beta(i,j) << " ";
          outfile << kim_params.lam2(i,j) << " ";
          outfile << kim_params.B(i,j) << " ";
        } else {
          outfile << "0 0 0 0 ";
        }
        outfile << kim_params.R(i,j,k) << " ";
        outfile << kim_params.D(i,j,k) << " ";
        if (j == k) {
          // Two-body parameters are only taken from j == k, so in
          // order to make that clear, we write zeros otherwise.
          outfile << kim_params.lam1(i,j) << " ";
          outfile << kim_params.A(i,j) << endl;
        } else {
          outfile << "0 0" << endl;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairTersoff::repulsive(double r, double fc, double fc_d,
                              int itype, int jtype,
                              bool eflag, double &eng) const
{
  // In this case, we look up the parameters inside the function,
  // since derived classes (e.g., for ZBL support) need more
  // parameters.
  const double lam1 = params2(itype,jtype).lam1;
  const double A = params2(itype,jtype).A;

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
  if (m == 3) arg = pow3(lam3 * (rij-rik));
  else arg = lam3 * (rij-rik); // m == 1

  double ex_delr;
  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc(rik,R,D) * ters_gijk(costheta,gamma,c2,d2,c2_d2,h) * ex_delr;
}

/* ---------------------------------------------------------------------- */

double PairTersoff::force_zeta(double r, double fc, double fc_d, double zeta_ij,
                               int itype, int jtype,
                               double &prefactor,
                               bool eflag, double &eng) const
{
  // In this case, we look up the parameters inside the function,
  // since derived classes (e.g., for ZBL support) need more
  // parameters.
  const double beta = params2(itype,jtype).beta;
  const double n = params2(itype,jtype).n;
  const double* n_precomp = params2(itype,jtype).n_precomp;

  const double fa = ters_fa(r, fc, itype, jtype);
  const double fa_d = ters_fa_d(r, fc, fc_d, itype, jtype);
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
                             double rij, double rik,
                             double rijsq, double riksq,
                             const double *delrij, const double *delrik,
                             // output for force/virial
                             double *fi, double *fj, double *fk,
                             // output for process_dEdr
                             double &rjk, double *delrjk,
                             double &drij, double &drik, double &drjk,
                             //
                             bool want_force_or_virial,
                             bool want_dEdr) const
{
  /* This routine is a little bit more complicated now, since we want
     to both support virial calculation as implemented/defined in
     LAMMPS and the process_dEdr() callback of KIM. For the former, we
     need derivatives with regard to the positions and for the latter
     with regard to the distances. So we have some booleans to toggle
     what we want to compute. Currently, there might be a slight
     performance penalty if process_dEdr is requested but not virial.
     This is to avoid two code paths for certain computations
     depending on what values were requested.

     Will calculate r_jk (all related quantities) only if want_dEdr is
     true.

     TODO: can we save some computations and combine both at least partially?   
  */

  // Values needed by either approach.
  const double fc = ters_fc(rik,R,D);
  const double dfc = ters_fc_d(rik,R,D);

  const double tmp =
    (m == 3)
    ? pow3(lam3 * (rij-rik))    // m == 3
    : lam3 * (rij-rik);         // m == 1

  const double ex_delr =
    (tmp > 69.0776) ? 1.0e30 : ( (tmp < -69.0776) ? 0.0 : exp(tmp) );

  const double ex_delr_d =
    (m == 3)
    ? 3.0*pow3(lam3) * pow2(rij-rik)*ex_delr      // m == 3
    : lam3 * ex_delr;                             // m == 1

  double rij_hat[3];
  vec3_scale(1.0/rij, delrij, rij_hat); // rij_hat = delrij / rij

  double rik_hat[3];
  vec3_scale(1.0/rik, delrik, rik_hat); // rik_hat = delrik / rik

  const double cos_theta = vec3_dot(rij_hat, rik_hat);

  const double gijk = ters_gijk(cos_theta, gamma, c2, d2, c2_d2, h);
  const double gijk_d = ters_gijk_d(cos_theta, gamma, c2, d2, h);

  if (want_force_or_virial) {
    ters_zetaterm_d_pos(prefactor, fc, dfc, ex_delr, ex_delr_d, cos_theta,
                        gijk, gijk_d,
                        rij_hat, rij,
                        rik_hat, rik,
                        fi, fj, fk);
  }

  if (want_dEdr) {
    // compute jk vector
    delrjk[0] = delrik[0] - delrij[0];
    delrjk[1] = delrik[1] - delrij[1];
    delrjk[2] = delrik[2] - delrij[2];
    const double rjksq = delrjk[0]*delrjk[0]
                         + delrjk[1]*delrjk[1]
                         + delrjk[2]*delrjk[2];
    rjk = sqrt(rjksq);

    ters_zetaterm_d_dist(-prefactor, fc, dfc, ex_delr, ex_delr_d,
                         gijk, gijk_d,
                         rij, rik, rjk,
                         rijsq, riksq, rjksq,
                         drij, drik, drjk);
  }
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
                            int itype, int jtype) const
{
  if (fc == 0.0) return 0.0;

  // In this case, we look up the parameters inside the function,
  // since derived classes (e.g., for ZBL support) need more
  // parameters.
  const double B = params2(itype,jtype).B;
  const double lam2 = params2(itype,jtype).lam2;

  return -B * exp(-lam2 * r) * fc;
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fa_d(double r, double fc, double fc_d,
                              int itype, int jtype) const
{
  if (fc == 0.0) return 0.0;

  // In this case, we look up the parameters inside the function,
  // since derived classes (e.g., for ZBL support) need more
  // parameters.
  const double B = params2(itype,jtype).B;
  const double lam2 = params2(itype,jtype).lam2;

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

void PairTersoff::ters_zetaterm_d_pos(double prefactor,
                                      double fc, double dfc,
                                      double ex_delr, double ex_delr_d,
                                      double cos_theta,
                                      double gijk, double gijk_d,
                                      double *rij_hat, double rij,
                                      double *rik_hat, double rik,
                                      double *dri, double *drj, double *drk) const
{
  double dcosdri[3], dcosdrj[3], dcosdrk[3];
  costheta_d(rij_hat,rij, rik_hat,rik, cos_theta, dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

void PairTersoff::ters_zetaterm_d_dist(double prefactor,
                                       double fc, double dfc,
                                       double ex_delr, double ex_delr_d,
                                       double gijk, double gijk_d,
                                       double rij, double rik, double rjk,
                                       double rijsq, double riksq, double rjksq,
                                       double &drij, double &drik, double &drjk) const
{
  // dcos_ijk/drij, dcos_ijk/drik, dcos_ijk/drjk
  const double dcosdrij = (rijsq - riksq + rjksq)/(2*rijsq*rik);
  const double dcosdrik = (riksq - rijsq + rjksq)/(2*rij*riksq);
  const double dcosdrjk = -rjk/(rij*rik);

  // compute the derivative wrt rij, rik, rjk
  drij = (fc*gijk_d*ex_delr*dcosdrij + fc*gijk*ex_delr_d) * prefactor;
  drik = (dfc*gijk*ex_delr + fc*gijk_d*ex_delr*dcosdrik - fc*gijk*ex_delr_d) * prefactor;
  drjk = (fc*gijk_d*ex_delr*dcosdrjk) * prefactor;
}

/* ---------------------------------------------------------------------- */

void PairTersoff::costheta_d(double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double cos_theta,
                             double *dri, double *drj, double *drk) const
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk
  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}

