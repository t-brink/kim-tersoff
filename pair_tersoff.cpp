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

   Modified for use with KIM by Tobias Brink (2012,2013).
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
  : kim_indices(ki), n_spec(n_spec), params(n_spec, n_spec, n_spec)
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

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairTersoff::~PairTersoff()
{
}

/* ---------------------------------------------------------------------- */

void PairTersoff::compute(KIM_API_model& kim_model,
                          bool use_neighbor_list, // false -> cluster/MI_OPBC
                          bool use_distvec, // use rij array from neighbor list?
                          KIM_IterLoca access_mode,
                          int n_atoms, // Actual number of atoms
                          const int* atom_types,
                          const Array2D<double>& atom_coords,
                          double* boxSideLengths, // If != NULL -> MI_OPBC.
                          double* energy, double* atom_energy,
                          Array2D<double>* forces,
                          double* virial,
                          Array2D<double>* particleVirial) const
{
  int ii;          // Iteration over all atoms.
  int error;       // KIM error code.
  int n_neigh;     // Number of neighbors of i.
  int* neighbors;  // The indices of the neighbors.
  double* distvec; // Rij vectors.
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
  if (particleVirial)
    for (int i = 0; i != n_atoms; ++i)
      for (int j = 0; j != 6; ++j)
        (*particleVirial)(i, j) = 0.0;

  // In iterator mode: reset the iterator.
  if (use_neighbor_list && access_mode == KIM_ITERATOR_MODE) {
    error =
      kim_model.get_neigh(access_mode,
                          0, // reset iterator.
                          &ii, // not set when resetting
                          &n_neigh, // not set when resetting
                          &neighbors, // not set when resetting
                          &distvec // not set when resetting
                          );
    if (error != KIM_STATUS_OK && error != KIM_STATUS_NEIGH_ITER_INIT_OK) {
      kim_model.report_error(__LINE__, __FILE__,
                             "KIM_API_get_neigh (reset iterator)",
                             error);
      throw runtime_error("compute: Error in KIM_API_get_neigh");
    }
  } else { // cluster
    n_neigh = n_atoms;
    distvec = NULL;
  }

  // loop over full neighbor list of my atoms

  // ii will only be changed in locator mode, otherwise it stays zero
  // and the loop will be broken by iterator access return code.
  ii = 0;
  while (ii != n_atoms) {
    int i;
    // Get neighbors.
    if (use_neighbor_list) {
      error =
        kim_model.get_neigh(access_mode, // locator or iterator
                            // The central atom or command to get next
                            // atom in iterator mode.
                            (access_mode == KIM_LOCATOR_MODE
                             ? ii++
                             : 1),
                            // Output.
                            &i, // The central atom.
                            &n_neigh, // Number of neighbors
                            &neighbors, // The neighbor indices
                            &distvec // Rij
                            );
      if (access_mode == KIM_ITERATOR_MODE &&
          error == KIM_STATUS_NEIGH_ITER_PAST_END)
        break; // Loop breaking condition is that we used up all central atoms.
      if (error != KIM_STATUS_OK) {
        kim_model.report_error(__LINE__, __FILE__, "KIM_API_get_neigh",
                               error);
        throw runtime_error("compute: Error in KIM_API_get_neigh");
      }
    } else {
      i = ii;
      ++ii;
    }

    const int itype = atom_types[i];
    double xtmp, ytmp, ztmp;
    if (!use_distvec) {
      xtmp = atom_coords(i,0);
      ytmp = atom_coords(i,1);
      ztmp = atom_coords(i,2);
    }

    for (int jj = 0; jj != n_neigh; ++jj) {
      int j = use_neighbor_list ? neighbors[jj] : jj;
      if (!use_neighbor_list && i == j) continue;
      const int jtype = atom_types[j];

      double delr_ij[3];
      if (use_distvec) {
        delr_ij[0] = distvec[jj*3 + 0];
        delr_ij[1] = distvec[jj*3 + 1];
        delr_ij[2] = distvec[jj*3 + 2];
      } else {
        delr_ij[0] = atom_coords(j,0) - xtmp;
        delr_ij[1] = atom_coords(j,1) - ytmp;
        delr_ij[2] = atom_coords(j,2) - ztmp;
        if (boxSideLengths)
          for (int d = 0; d != 3; ++d)
            if (abs(delr_ij[d]) > 0.5 * boxSideLengths[d])
              delr_ij[d] -= (delr_ij[d]/abs(delr_ij[d]))*boxSideLengths[d];
      }
      const double rsq_ij = delr_ij[0]*delr_ij[0]
                          + delr_ij[1]*delr_ij[1]
                          + delr_ij[2]*delr_ij[2];

      const double cutsq = params(itype,jtype,jtype).cutsq;
      if (rsq_ij > cutsq) continue;

      const double r_ij = sqrt(rsq_ij);

      const double R = params(itype,jtype,jtype).R;
      const double D = params(itype,jtype,jtype).D;
      const double fc_ij = ters_fc(r_ij, R, D); // Value of the cutoff function.
      const double dfc_ij = ters_fc_d(r_ij, R, D); // Derivative of fc_ij.

      // two-body interactions

      const double lam1 = params(itype,jtype,jtype).lam1;
      const double A = params(itype,jtype,jtype).A;

      double evdwl; // Particle energy.
      const double fpair =
        repulsive(r_ij, fc_ij, dfc_ij, lam1, A, eflag, evdwl);

      if (energy)
        *energy += 0.5 * evdwl;

      if (atom_energy) {
        atom_energy[i] += 0.5 * evdwl;
      }

      if (forces || virial || particleVirial) {
        const double fx = delr_ij[0]*fpair;
        const double fy = delr_ij[1]*fpair;
        const double fz = delr_ij[2]*fpair;

        if (forces) {
          (*forces)(i,0) -= fx;
          (*forces)(i,1) -= fy;
          (*forces)(i,2) -= fz;
        }

        if (virial || particleVirial) {
          const double vxx = 0.5 * delr_ij[0] * fx;
          const double vyy = 0.5 * delr_ij[1] * fy;
          const double vzz = 0.5 * delr_ij[2] * fz;
          const double vyz = 0.5 * delr_ij[1] * fz;
          const double vxz = 0.5 * delr_ij[0] * fz;
          const double vxy = 0.5 * delr_ij[0] * fy;

          if (virial) {
            virial[0] -= vxx;
            virial[1] -= vyy;
            virial[2] -= vzz;
            virial[3] -= vyz;
            virial[4] -= vxz;
            virial[5] -= vxy;
          }

          if (particleVirial) {
            (*particleVirial)(i,0) -= vxx;
            (*particleVirial)(i,1) -= vyy;
            (*particleVirial)(i,2) -= vzz;
            (*particleVirial)(i,3) -= vyz;
            (*particleVirial)(i,4) -= vxz;
            (*particleVirial)(i,5) -= vxy;
          }
        }
      }

      // three-body interactions

      // accumulate bondorder zeta for each i-j interaction via loop over k

      double zeta_ij = 0.0;

      for (int kk = 0; kk != n_neigh; ++kk) {
        if (jj == kk) continue;
        int k = use_neighbor_list ? neighbors[kk] : kk;
        if (!use_neighbor_list && i == k)
          // Needed in cluster mode, as there is no neighbor list
          // which will make sure i = k doesn't happen!
          continue;
        const int ktype = atom_types[k];

        double delr_ik[3];
        if (use_distvec) {
          delr_ik[0] = distvec[kk*3 + 0];
          delr_ik[1] = distvec[kk*3 + 1];
          delr_ik[2] = distvec[kk*3 + 2];
        } else {
          delr_ik[0] = atom_coords(k,0) - xtmp;
          delr_ik[1] = atom_coords(k,1) - ytmp;
          delr_ik[2] = atom_coords(k,2) - ztmp;
          if (boxSideLengths)
            for (int d = 0; d != 3; ++d)
              if (abs(delr_ik[d]) > 0.5 * boxSideLengths[d])
                delr_ik[d] -= (delr_ik[d]/abs(delr_ik[d]))*boxSideLengths[d];
        }
        const double rsq_ik = delr_ik[0]*delr_ik[0]
                            + delr_ik[1]*delr_ik[1]
                            + delr_ik[2]*delr_ik[2];
        const double cutsq = params(itype,jtype,ktype).cutsq;
        if (rsq_ik > cutsq) continue;

        const double r_ik = sqrt(rsq_ik);

        const int m = params(itype,jtype,ktype).m;
        const double lam3 = params(itype,jtype,ktype).lam3;
        const double R = params(itype,jtype,ktype).R;
        const double D = params(itype,jtype,ktype).D;
        const double gamma = params(itype,jtype,ktype).gamma;
        const double c2 = params(itype,jtype,ktype).c2;
        const double d2 = params(itype,jtype,ktype).d2;
        const double c2_d2 = params(itype,jtype,ktype).c2_d2;
        const double h = params(itype,jtype,ktype).h;

        zeta_ij += zeta(r_ij,r_ik,m,lam3,R,D,gamma,c2,d2,c2_d2,h,
                        delr_ij,delr_ik);
      }

      // pairwise force due to zeta

      const double B = params(itype,jtype,jtype).B;
      const double lam2 = params(itype,jtype,jtype).lam2;
      const double beta = params(itype,jtype,jtype).beta;
      const double n = params(itype,jtype,jtype).n;
      const double* n_precomp = params(itype,jtype,jtype).n_precomp;

      double prefactor; // -0.5 * fa * ∇bij
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

      if (forces || virial || particleVirial) {
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

        if (virial || particleVirial) {
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

          if (particleVirial) {
            const double vxx2 = 0.5 * vxx;
            const double vyy2 = 0.5 * vyy;
            const double vzz2 = 0.5 * vzz;
            const double vyz2 = 0.5 * vyz;
            const double vxz2 = 0.5 * vxz;
            const double vxy2 = 0.5 * vxy;

            (*particleVirial)(i,0) += vxx2;
            (*particleVirial)(i,1) += vyy2;
            (*particleVirial)(i,2) += vzz2;
            (*particleVirial)(i,3) += vyz2;
            (*particleVirial)(i,4) += vxz2;
            (*particleVirial)(i,5) += vxy2;

            (*particleVirial)(j,0) += vxx2;
            (*particleVirial)(j,1) += vyy2;
            (*particleVirial)(j,2) += vzz2;
            (*particleVirial)(j,3) += vyz2;
            (*particleVirial)(j,4) += vxz2;
            (*particleVirial)(j,5) += vxy2;
          }
        }
      }

      // attractive term via loop over k

      if (forces || virial || particleVirial) {
        for (int kk = 0; kk != n_neigh; ++kk) {
          if (jj == kk) continue;
          int k = use_neighbor_list ? neighbors[kk] : kk;
          if (!use_neighbor_list && i == k)
            // Needed in cluster mode, as there is no neighbor list
            // which will make sure i = k doesn't happen!
            continue;
          const int ktype = atom_types[k];

          double delr_ik[3];
          if (use_distvec) {
            delr_ik[0] = distvec[kk*3 + 0];
            delr_ik[1] = distvec[kk*3 + 1];
            delr_ik[2] = distvec[kk*3 + 2];
          } else {
            delr_ik[0] = atom_coords(k,0) - xtmp;
            delr_ik[1] = atom_coords(k,1) - ytmp;
            delr_ik[2] = atom_coords(k,2) - ztmp;
            if (boxSideLengths)
              for (int d = 0; d != 3; ++d)
                if (abs(delr_ik[d]) > 0.5 * boxSideLengths[d])
                  delr_ik[d] -= (delr_ik[d]/abs(delr_ik[d]))*boxSideLengths[d];
          }
          const double rsq_ik = delr_ik[0]*delr_ik[0]
            + delr_ik[1]*delr_ik[1]
            + delr_ik[2]*delr_ik[2];
          const double cutsq = params(itype,jtype,ktype).cutsq;
          if (rsq_ik > cutsq) continue;

          const double r_ik = sqrt(rsq_ik);

          const double R = params(itype,jtype,ktype).R;
          const double D = params(itype,jtype,ktype).D;
          const int m = params(itype,jtype,ktype).m;
          const double lam3 = params(itype,jtype,ktype).lam3;
          const double gamma = params(itype,jtype,ktype).gamma;
          const double c2 = params(itype,jtype,ktype).c2;
          const double d2 = params(itype,jtype,ktype).d2;
          const double c2_d2 = params(itype,jtype,ktype).c2_d2;
          const double h = params(itype,jtype,ktype).h;

          double fi[3], fj[3], fk[3];

          attractive(prefactor, r_ij, r_ik,
                     R, D, m, lam3, gamma, c2, d2, c2_d2, h,
                     delr_ij, delr_ik, fi, fj, fk);

          if (forces) {
            (*forces)(i,0) += fi[0];
            (*forces)(i,1) += fi[1];
            (*forces)(i,2) += fi[2];
            (*forces)(j,0) += fj[0];
            (*forces)(j,1) += fj[1];
            (*forces)(j,2) += fj[2];
            (*forces)(k,0) += fk[0];
            (*forces)(k,1) += fk[1];
            (*forces)(k,2) += fk[2];
          }

          if (virial || particleVirial) {
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

            if (particleVirial) {
              const double vxx3 = 1.0/3.0 * vxx;
              const double vyy3 = 1.0/3.0 * vyy;
              const double vzz3 = 1.0/3.0 * vzz;
              const double vyz3 = 1.0/3.0 * vyz;
              const double vxz3 = 1.0/3.0 * vxz;
              const double vxy3 = 1.0/3.0 * vxy;

              (*particleVirial)(i,0) -= vxx3;
              (*particleVirial)(i,1) -= vyy3;
              (*particleVirial)(i,2) -= vzz3;
              (*particleVirial)(i,3) -= vyz3;
              (*particleVirial)(i,4) -= vxz3;
              (*particleVirial)(i,5) -= vxy3;

              (*particleVirial)(j,0) -= vxx3;
              (*particleVirial)(j,1) -= vyy3;
              (*particleVirial)(j,2) -= vzz3;
              (*particleVirial)(j,3) -= vyz3;
              (*particleVirial)(j,4) -= vxz3;
              (*particleVirial)(j,5) -= vxy3;

              (*particleVirial)(k,0) -= vxx3;
              (*particleVirial)(k,1) -= vyy3;
              (*particleVirial)(k,2) -= vzz3;
              (*particleVirial)(k,3) -= vyz3;
              (*particleVirial)(k,4) -= vxz3;
              (*particleVirial)(k,5) -= vxy3;
            }
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
  Params temp_params;
  string type_i, type_j, type_k;
  while (buffer >> type_i
                >> type_j
                >> type_k
                >> temp_params.m
                >> temp_params.gamma
                >> temp_params.lam3
                >> temp_params.c
                >> temp_params.d
                >> temp_params.h // costheta0
                >> temp_params.n
                >> temp_params.beta
                >> temp_params.lam2
                >> temp_params.B
                >> temp_params.R
                >> temp_params.D
                >> temp_params.lam1
                >> temp_params.A) {
    // Unit conversion.
    temp_params.A *= energy_conv;
    temp_params.B *= energy_conv;
    temp_params.lam1 *= inv_length_conv;
    temp_params.lam2 *= inv_length_conv;
    temp_params.lam3 *= inv_length_conv;
    temp_params.R *= length_conv;
    temp_params.D *= length_conv;
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
    params(i,j,k) = temp_params;
  }
  if (!got_interaction.all())
    throw runtime_error("Not all interactions were set!");
  // Check parameters and pre-compute values.
  prepare_params();
}

/* Check parameters and pre-compute values. */
void PairTersoff::prepare_params() {
  max_cutoff = 0.0;
  for (int i = 0; i != n_spec; ++i) {
    string type_i = to_spec.at(i);
    for (int j = 0; j != n_spec; ++j) {
      string type_j = to_spec.at(j);
      for (int k = 0; k != n_spec; ++k) {
        string type_k = to_spec.at(k);
        Params& temp_params = params(i,j,k);
        // Check values of parameters.
        if (temp_params.c < 0)
          throw runtime_error("Parameter c ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.d < 0)
          throw runtime_error("Parameter d ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.n < 0)
          throw runtime_error("Parameter n ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.beta < 0)
          throw runtime_error("Parameter beta ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.lam1 < 0)
          throw runtime_error("Parameter lambda1 ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.lam2 < 0)
          throw runtime_error("Parameter lambda2 ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.A < 0)
          throw runtime_error("Parameter A ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.B < 0)
          throw runtime_error("Parameter B ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.m != 1 && temp_params.m != 3)
          throw runtime_error("Parameter m ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") must be one or three.");
        if (temp_params.R < 0)
          throw runtime_error("Parameter R ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.D < 0)
          throw runtime_error("Parameter D ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        if (temp_params.D > temp_params.R)
          throw runtime_error("Parameter D ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") must be smaller than R.");
        if (temp_params.gamma < 0)
          throw runtime_error("Parameter gamma ("
                              + type_i +
                              "-" + type_j +
                              "-" + type_k +
                              ") may not be smaller than zero.");
        // Cutoff.
        temp_params.cut =
          temp_params.R + temp_params.D; // max cutoff.
        temp_params.cutmin =
          temp_params.R - temp_params.D; // below this the cutoff value is 1.0
        temp_params.cutsq =
          temp_params.cut * temp_params.cut; // for fast check if inside cutoff
        // Get the cutoff to pass to KIM, which is the biggest cutoff.
        if (temp_params.cut > max_cutoff)
          max_cutoff = temp_params.cut;
        // Pre-compute values.
        const double n2 = 2.0 * temp_params.n;
        const double n_r = -1.0/temp_params.n;
        temp_params.n_precomp[0] = pow(n2 * 1e-16, n_r);
        temp_params.n_precomp[1] = pow(n2 * 1e-8, n_r);
        temp_params.n_precomp[2] = 1.0 / temp_params.n_precomp[1];
        temp_params.n_precomp[3] = 1.0 / temp_params.n_precomp[0];
        temp_params.c2 = temp_params.c * temp_params.c;
        temp_params.d2 = temp_params.d * temp_params.d;
        temp_params.c2_d2 = temp_params.c2 / temp_params.d2;
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
                         double* delrij, double* delrik) const
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
                             double rij, double rik,
                             double R, double D, int m, double lam3,
                             double gamma, double c2, double d2, double c2_d2,
                             double h,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk) const
{
  double rij_hat[3];
  vec3_scale(1.0/rij, delrij, rij_hat); // rij_hat = delrij / rij

  double rik_hat[3];
  vec3_scale(1.0/rik, delrik, rik_hat); // rik_hat = delrik / rik

  ters_zetaterm_d(prefactor,
                  R, D, m, lam3, gamma, c2, d2, c2_d2, h,
                  rij_hat,rij,rik_hat,rik,fi,fj,fk);
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
                               (1.0 - 0.5*(1.0 +  1.0/(2.0*n)) *
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
                                  double *rij_hat, double rij,
                                  double *rik_hat, double rik,
                                  double *dri, double *drj, double *drk) const
{
  double ex_delr,ex_delr_d,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  const double fc = ters_fc(rik,R,D);
  const double dfc = ters_fc_d(rik,R,D);
  if (m == 3) tmp = pow(lam3 * (rij-rik), 3);
  else tmp = lam3 * (rij-rik); // m == 1

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (m == 3)
    ex_delr_d = 3.0*pow(lam3, 3) * pow(rij-rik, 2)*ex_delr;
  else ex_delr_d = lam3 * ex_delr; // m == 1

  const double cos_theta = vec3_dot(rij_hat,rik_hat);
  const double gijk = ters_gijk(cos_theta, gamma, c2, d2, c2_d2, h);
  const double gijk_d = ters_gijk_d(cos_theta, gamma, c2, d2, h);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

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

void PairTersoff::costheta_d(double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double *dri, double *drj, double *drk) const
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  const double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}
