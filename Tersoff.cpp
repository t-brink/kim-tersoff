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

#include "Tersoff.hpp"

#include <cmath>
#include <stdexcept>
#include <sstream>

#include <KIM_API_status.h>

using namespace std;

// Implementation of the Tersoff class.
namespace model_driver_Tersoff {

  // HELPERS for computing forces and energy.

  static inline double VR(double r, double A, double lambda1) {
    return A * exp(-lambda1 * r);
  }
  static inline double VA(double r, double B, double lambda2) {
    return -B * exp(-lambda2 * r);
  }

  static inline Vec3D<double> del_VR(Vec3D<double> R_hat, double VR,
                                     double lambda1) {
    R_hat *= lambda1 * VR;
    return R_hat;
  }
  static inline Vec3D<double> del_VA(Vec3D<double> R_hat, double VA,
                                     double lambda2) {
    R_hat *= -lambda2 * VA;
    return R_hat;
  }

  static inline double fc(double r, double R, double D,
                          double cutmin, double cutoff) {
    if (r > cutoff)
      return 0.0;
    else if (r < cutmin)
      return 1.0;
    else
      return 0.5 - 0.5 * sin(M_PI_2 * (r - R) / D);
  }

  static inline void fc_del_fc(const Vec3D<double>& R_hat, double r,
                               double R, double D,
                               double cutmin, double cutoff,
                               double& fc, Vec3D<double>& del_fc) {
    if (r > cutoff) {
      fc = 0.0;
      del_fc[0] = 0.0; del_fc[1] = 0.0; del_fc[2] = 0.0;
    } else if (r < cutmin) {
      fc = 1.0;
      del_fc[0] = 0.0; del_fc[1] = 0.0; del_fc[2] = 0.0;
    } else {
      const double inner = M_PI_2 * (r - R) / D;
      fc = 0.5 - 0.5 * sin(inner);
      del_fc = R_hat * M_PI_4 / D * cos(inner);
    }
  }

  // Single addend of the zeta term (for each k).
  /*
  static inline double zeta_k(double rij, double rik, ) {
  }
  */

  static inline double costheta(double rij, const Vec3D<double>& Rij,
                                double rik, const Vec3D<double>& Rik) {
    return dot(Rij, Rik) / (rij * rik);
  }

  static inline
  Vec3D<double> del_costheta(double rij, Vec3D<double> Rij_hat,
                             Vec3D<double> Rij,
                             double rik, Vec3D<double> Rik_hat,
                             const Vec3D<double>& Rik,
                             double costheta) {
    // Rij_hat/rij * costheta + Rik_hat/rik * costheta - (Rij + Rik)/(rij*rik)
    Rij_hat *= costheta / rij;
    Rik_hat *= costheta / rik;
    Rij += Rik;  Rij /= rij*rik;
    Rij_hat += Rik_hat - Rij;
    return Rij_hat;
  }

  static inline double g(double costheta,
                         double gamma, double c, double d, double h) {
    return gamma * (1
                    + pow2(c) / pow2(d)
                    - pow2(c) / ( pow2(d) + pow2(costheta - h) )
                    );
  }

  static inline Vec3D<double> del_g(double rij,
                                    const Vec3D<double>& Rij_hat,
                                    const Vec3D<double>& Rij,
                                    double rik,
                                    const Vec3D<double>& Rik_hat,
                                    const Vec3D<double>& Rik,
                                    double costheta,
                                    double gamma,
                                    double c, double d, double h) {
    Vec3D<double> del_costheta_ijk = del_costheta(rij, Rij_hat, Rij,
                                                  rik, Rik_hat, Rik,
                                                  costheta);
    del_costheta_ijk *=
      2 * gamma * pow2(c) * (costheta - h)
      / pow2( pow2(d) + pow2(costheta - h) );
    return del_costheta_ijk;
  }

  static inline double b(double zeta, double beta, double n) {
    return pow(1 + pow(beta, n) * pow(zeta, n), -1.0 / (2.0 * n));
  }

  static inline Vec3D<double> del_b(Vec3D<double> del_zeta, double zeta,
                                    double beta, double n) {
    del_zeta *=
      -0.5 * pow(beta, n)
      * pow(1 + pow(beta, n) * pow(zeta, n), -1.0 / (2.0 * n) - 1)
      * pow(zeta, n - 1);
    return del_zeta;
  }


  static inline double pairE(double fc, double VR, double bij, double VA) {
    return fc * (VR + bij * VA);
  }

  // PUBLIC

  void Tersoff::compute_cluster(int n_atoms, const int* atom_types,
                                const Array2D<double>& atom_coords,
                                double* energy, double* atom_energy,
                                Array2D<double>* forces) const {
    if (energy)
      *energy = 0.0;
    for (int i = 0; i != n_atoms; ++i) { //TODO: i < j!!!!!!!!!!!!!!!
      if (atom_energy)
        atom_energy[i] = 0.0;
      if (forces) {
        (*forces)(i, 0) = 0.0;
        (*forces)(i, 1) = 0.0;
        (*forces)(i, 2) = 0.0;
      }
      for (int j = 0; j != n_atoms; ++j) {
        if (i == j)
          continue;

        // pos_j - pos_i; points from j to i.
        const Vec3D<double> Rij(atom_coords(j,0) - atom_coords(i,0),
                                atom_coords(j,1) - atom_coords(i,1),
                                atom_coords(j,2) - atom_coords(i,2));
        const double rij_sq = Rij.abssq();

        const int type_i = atom_types[i];
        const int type_j = atom_types[j];

        // Check if inside cutoff.
        const Params& params_ij = params(type_i, type_j, type_j);
        if (rij_sq > params_ij.cutsq)
          continue;

        const double rij = sqrt(rij_sq);
        Vec3D<double> Rij_hat;
        if (forces)
          Rij_hat = Rij / rij;
        // Value of the cutoff function.
        double fc_ij;
        Vec3D<double> del_fc_ij;
        if (forces)
          fc_del_fc(Rij_hat, rij,
                    params_ij.R, params_ij.D,
                    params_ij.cutmin, params_ij.cut,
                    fc_ij, del_fc_ij);
        else
          fc_ij = fc(rij, params_ij.R, params_ij.D,
                     params_ij.cutmin, params_ij.cut);
        // Repulsive term.
        const double VRij = VR(rij, params_ij.A, params_ij.lam1);
        Vec3D<double> del_VRij;
        if (forces)
          del_VRij = del_VR(Rij_hat, VRij, params_ij.lam1);
        // Attractive term.
        const double VAij = VA(rij, params_ij.B, params_ij.lam2);
        Vec3D<double> del_VAij;
        if (forces)
          del_VAij = del_VA(Rij_hat, VAij, params_ij.lam2);
        // Bond order term.
        double zeta = 0.0;
        Vec3D<double> del_zeta(0.0, 0.0, 0.0);
        for (int k = 0; k != n_atoms; ++k) {
          if (k == i || k == j)
            continue;
          // pos_k - pos_i; points from k to i.
          const Vec3D<double> Rik(atom_coords(k,0) - atom_coords(i,0),
                                  atom_coords(k,1) - atom_coords(i,1),
                                  atom_coords(k,2) - atom_coords(i,2));
          const double rik_sq = Rik.abssq();
          const int type_k = atom_types[k];
          // Check if inside cutoff.
          const Params& params_ijk = params(type_i, type_j, type_k);
          if (rik_sq > params_ijk.cutsq)
            continue;
          const double rik = sqrt(rik_sq);
          Vec3D<double> Rik_hat;
          if (forces)
            Rik_hat = Rik / rik;
          double fc_ik;
          Vec3D<double> del_fc_ik;
          if (forces)
            fc_del_fc(Rik_hat, rik,
                      params_ijk.R, params_ijk.D,
                      params_ijk.cutmin, params_ijk.cut,
                      fc_ik, del_fc_ik);
          else
            fc_ik = fc(rik, params_ijk.R, params_ijk.D,
                       params_ijk.cutmin, params_ijk.cut);
          const double exp_term = exp(pow(params_ijk.lam3, params_ijk.m) *
                                      pow(rij - rik, params_ijk.m));
          const double costheta = dot(Rij, Rik) / (rij * rik);
          const double g_ijk = g(costheta, params_ijk.gamma,
                                 params_ijk.c, params_ijk.d, params_ijk.h);
          Vec3D<double> del_g_ijk;
          if (forces)
            del_g_ijk = del_g(rij, Rij_hat, Rij, rik, Rik_hat, Rik, costheta,
                              params_ijk.gamma,
                              params_ijk.c, params_ijk.d, params_ijk.h);
          zeta += fc_ik * exp_term * g_ijk;
          if (forces)
            del_zeta +=
              del_fc_ik * exp_term * g_ijk
              + fc_ik * exp_term * del_g_ijk
              - (Rij_hat - Rik_hat) * double(params_ijk.m) * pow(rij - rik, params_ijk.m-1) * pow(params_ijk.lam3, params_ijk.m) * zeta;
        }
        const double bij = b(zeta, params_ij.beta, params_ij.n);
        Vec3D<double> del_bij;
        if (forces) {
          del_bij = del_b(del_zeta, zeta, params_ij.beta, params_ij.n);
          const Vec3D<double> F =
            -0.5 * (del_fc_ij * (VRij + bij * VAij)
                    + fc_ij * (del_VRij + del_bij * VAij + bij * del_VAij)
                    );
          (*forces)(i, 0) += F[0];
          (*forces)(i, 1) += F[1];
          (*forces)(i, 2) += F[2];
        }
        // Sum it up.
        if (energy || atom_energy) {
          const double Eij = 0.5 * pairE(fc_ij, VRij, bij, VAij);
          if (energy)
            *energy += Eij;
          if (atom_energy)
            atom_energy[i] += Eij;
        }
      }
    }
  }


  void Tersoff::compute_neigh_list(KIM_API_model& kim_model,
                                   KIM_IterLoca access_mode,
                                   int n_atoms, // Actual number of atoms
                                   const int* atom_types,
                                   const Array2D<double>& atom_coords,
                                   double* energy, double* atom_energy) const {
    int error;
    int ii; // central atom
    int n_neigh; // number of neighbors of i
    int* neighbors; // the indices of the neighbors
    double* distvec; // Rij vectors

    // If requested, reset energy.
    if (energy)
      *energy = 0.0;

    // In iterator mode: reset the iterator.
    if (access_mode == KIM_ITERATOR_MODE) {
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
        throw runtime_error("compute_neigh_list: Error in KIM_API_get_neigh");
      }
    }
    ii = 0; // Start with atom 0.

    // Actual loop.
    while (ii != n_atoms) { // ii will only be changed in locator mode,
      // otherwise it stays zero and the loop will
      // be broken by iterator access return code.
      int i;
      // Get neighbors.
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
        throw runtime_error("compute_neigh_list: Error in KIM_API_get_neigh");
      }

      // Actually compute.
      if (atom_energy)
        atom_energy[i] = 0.0;
      for (int jj = 0; jj != n_neigh; ++jj) {
        int j = neighbors[jj];
        if (i == j)
          continue;

        // pos_j - pos_i; points from i to j.
        Vec3D<double> Rij;
        if (distvec) {
          Rij[0] = distvec[jj*3 + 0];
          Rij[1] = distvec[jj*3 + 1];
          Rij[2] = distvec[jj*3 + 2];
        } else {
          Rij[0] = atom_coords(j,0) - atom_coords(i,0);
          Rij[1] = atom_coords(j,1) - atom_coords(i,1);
          Rij[2] = atom_coords(j,2) - atom_coords(i,2);
        }
        const double rij_sq = Rij.abssq();

        const int type_i = atom_types[i];
        const int type_j = atom_types[j];

        // Check if inside cutoff.
        const Params& params_ij = params(type_i, type_j, type_j);
        if (rij_sq > params_ij.cutsq)
          continue;

        const double rij = sqrt(rij_sq);
        // Value of the cutoff function.
        const double fc_ij = fc(rij, params_ij.R, params_ij.D,
                                params_ij.cutmin, params_ij.cut);
        // Repulsive term.
        const double VRij = VR(rij, params_ij.A, params_ij.lam1);
        // Attractive term.
        const double VAij = VA(rij, params_ij.B, params_ij.lam2);
        // Bond order term.
        double zeta = 0.0;
        for (int kk = 0; kk != n_neigh; ++kk) {
          int k = neighbors[kk];
          if (k == i || k == j)
            continue;
          // pos_k - pos_i; points from k to i.
          Vec3D<double> Rik;
          if (distvec) {
            Rik[0] = distvec[kk*3 + 0];
            Rik[1] = distvec[kk*3 + 1];
            Rik[2] = distvec[kk*3 + 2];
          } else {
            Rik[0] = atom_coords(k,0) - atom_coords(i,0);
            Rik[1] = atom_coords(k,1) - atom_coords(i,1);
            Rik[2] = atom_coords(k,2) - atom_coords(i,2);
          }
          const double rik_sq = Rik.abssq();
          const int type_k = atom_types[k];
          // Check if inside cutoff.
          const Params& params_ijk = params(type_i, type_j, type_k);
          if (rik_sq > params_ijk.cutsq)
            continue;
          const double rik = sqrt(rik_sq);
          const double fc_ik = fc(rik, params_ijk.R, params_ijk.D,
                                  params_ijk.cutmin, params_ijk.cut);
          const double exp_term = exp(pow(params_ijk.lam3, params_ijk.m) *
                                      pow(rij - rik, params_ijk.m));
          const double costheta = dot(Rij, Rik) / (rij * rik);
          const double g_ijk = g(costheta, params_ijk.gamma,
                                 params_ijk.c, params_ijk.d, params_ijk.h);
          zeta += fc_ik * exp_term * g_ijk;
        }
        const double bij = b(zeta, params_ij.beta, params_ij.n);
        // Sum it up.
        if (energy || atom_energy) {
          const double Eij = 0.5 * pairE(fc_ij, VRij, bij, VAij);
          if (energy)
            *energy += Eij;
          if (atom_energy)
            atom_energy[i] += Eij;
        }
      }
    }
  }


  // PRIVATE

  void Tersoff::read_params(istream& infile, map<string,int> type_map,
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
    max_cutoff = 0.0;
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
      // Precalculated parameters.
      temp_params.cut = temp_params.R + temp_params.D; // max cutoff.
      temp_params.cutmin =
        temp_params.R - temp_params.D; // below this the cutoff value is 1.0
      // for neighbor list.
      temp_params.cutsq = temp_params.cut * temp_params.cut;
      // Get the cutoff to pass to KIM, which is the biggest cutoff.
      if (temp_params.cut > max_cutoff)
        max_cutoff = temp_params.cut;
      // Get the atom type indices.
      map<string,int>::const_iterator it;
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
  }

}
