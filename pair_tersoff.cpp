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

   Modified for use with KIM by Tobias Brink (2012).
------------------------------------------------------------------------- */

#include "pair_tersoff.hpp"

//#include <cmath>
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
#include <fstream>
#include <sstream>
/* DEBUG * /
#include <iostream>
/ * END DEBUG */

#include <KIM_API_status.h>

using namespace model_driver_Tersoff;
using namespace std;

/* ---------------------------------------------------------------------- */

PairTersoff::PairTersoff(std::string parameter_file,
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
  //cout << "Reading from: " << parameter_file << endl;
  std::fstream f(parameter_file.c_str(), std::ios_base::in);
  read_params(f, type_map, energy_conv, length_conv, inv_length_conv);
  /*
  cout << "n_spec = " << n_spec << endl;
  for (int i = 0; i != n_spec; ++i)
  for (int j = 0; j != n_spec; ++j)
  for (int k = 0; k != n_spec; ++k) {
    cout << " A " << i << "," << j << "," << k << " = " << params(i,j,k).A
         << endl;
  }
  */
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairTersoff::~PairTersoff()
{
}

/* ---------------------------------------------------------------------- */

void PairTersoff::compute(KIM_API_model& kim_model,
                          bool use_neighbor_list, // false -> cluster
                          KIM_IterLoca access_mode,
                          int n_atoms, // Actual number of atoms
                          const int* atom_types,
                          const Array2D<double>& atom_coords,
                          double* energy, double* atom_energy,
                          Array2D<double>* forces)
{
  int ii;
  int itype,jtype,ktype; // TODO: move in-loop
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  //int *ilist,*jlist,*numneigh,**firstneigh;
  int error;
  int n_neigh; // number of neighbors of i
  int* neighbors; // the indices of the neighbors
  double* distvec; // Rij vectors

  //TODO: figure this thing out.
  const static int eflag = 1;

  // If requested, reset energy.
  if (energy)
    *energy = 0.0;
  if (atom_energy)
    for (int i = 0; i != n_atoms; ++i) {
      atom_energy[i] = 0.0;
    }

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
      throw runtime_error("compute_neigh_list: Error in KIM_API_get_neigh");
    }
  } else { // cluster
    n_neigh = n_atoms;
    distvec = NULL;
  }

  // Reset forces.
  if (forces)
    for (int i = 0; i != n_atoms; ++i) {
      (*forces)(i, 0) = 0.0;
      (*forces)(i, 1) = 0.0;
      (*forces)(i, 2) = 0.0;
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
        throw runtime_error("compute_neigh_list: Error in KIM_API_get_neigh");
      }
    } else {
      i = ii;
      ++ii;
    }

    itype = atom_types[i];
    if (!distvec) {
      xtmp = atom_coords(i,0);
      ytmp = atom_coords(i,1);
      ztmp = atom_coords(i,2);
    }

    // two-body interactions, skip half of them  (TODO: skip half???? --brink)

    // TODO: Skip half here!  The two-body things can work with half
    // neighbor lists, while the three-body terms can't, so they
    // optimize here!  Perhaps combine with the next loop and skip on
    // repulsive() and f[...] = ...; when i > j or something?  This
    // seems more efficient!

    for (int jj = 0; jj != n_neigh; ++jj) {
      int j = use_neighbor_list ? neighbors[jj] : jj;
      if (i >= j) // skip half
        continue;
      jtype = atom_types[j];

      if (distvec) {
        delx = -distvec[jj*3 + 0]; // KIM outputs x_j - x_i, this code wants
        dely = -distvec[jj*3 + 1]; // x_i - x_j.
        delz = -distvec[jj*3 + 2];
      } else {
        delx = xtmp - atom_coords(j,0);
        dely = ytmp - atom_coords(j,1);
        delz = ztmp - atom_coords(j,2);
      }
      rsq = delx*delx + dely*dely + delz*delz;

      const double cutsq = params(itype,jtype,jtype).cutsq;
      if (rsq > cutsq) continue;

      const double lam1 = params(itype,jtype,jtype).lam1;
      const double A = params(itype,jtype,jtype).A;
      const double R = params(itype,jtype,jtype).R;
      const double D = params(itype,jtype,jtype).D;
      repulsive(rsq,lam1,A,R,D,fpair,eflag,evdwl);

      if (energy)
        *energy += evdwl;

      if (atom_energy) {
        atom_energy[i] += 0.5 * evdwl;
        atom_energy[j] += 0.5 * evdwl;
      }

      if (forces) {
        (*forces)(i,0) += delx*fpair;
        (*forces)(i,1) += dely*fpair;
        (*forces)(i,2) += delz*fpair;
        (*forces)(j,0) -= delx*fpair;
        (*forces)(j,1) -= dely*fpair;
        (*forces)(j,2) -= delz*fpair;
      }

      /* TODO: what does this do???
      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
      */
    }

    // three-body interactions
    // skip immediately if I-J is not within cutoff

    for (int jj = 0; jj != n_neigh; ++jj) {
      int j = use_neighbor_list ? neighbors[jj] : jj;
      if (i == j) // TODO: do we need this???? YES!
        continue;
      jtype = atom_types[j];
      //iparam_ij = elem2param[itype][jtype][jtype];

      if (distvec) {
        delr1[0] = distvec[jj*3 + 0]; // Here delr = x_j - x_i as in KIM!?
        delr1[1] = distvec[jj*3 + 1]; // TODO: above this seems to not matter,
        delr1[2] = distvec[jj*3 + 2]; // so unify that.
      } else {
        delr1[0] = atom_coords(j,0) - xtmp;
        delr1[1] = atom_coords(j,1) - ytmp;
        delr1[2] = atom_coords(j,2) - ztmp;
      }
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      const double cutsq = params(itype,jtype,jtype).cutsq;
      if (rsq1 > cutsq) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      zeta_ij = 0.0;

      for (int kk = 0; kk != n_neigh; ++kk) {
        if (jj == kk) continue; // TODO: This should suffice, but does it?
        int k = use_neighbor_list ? neighbors[kk] : kk;
        ktype = atom_types[k];
        //iparam_ijk = elem2param[itype][jtype][ktype];

        if (distvec) {
          delr2[0] = distvec[kk*3 + 0];
          delr2[1] = distvec[kk*3 + 1];
          delr2[2] = distvec[kk*3 + 2];
        } else {
          delr2[0] = atom_coords(k,0) - xtmp;
          delr2[1] = atom_coords(k,1) - ytmp;
          delr2[2] = atom_coords(k,2) - ztmp;
        }
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        const double cutsq = params(itype,jtype,ktype).cutsq;
        if (rsq2 > cutsq) continue;

        const int m = params(itype,jtype,ktype).m;
        const double lam3 = params(itype,jtype,ktype).lam3;
        const double R = params(itype,jtype,ktype).R;
        const double D = params(itype,jtype,ktype).D;
        const double gamma = params(itype,jtype,ktype).gamma;
        const double c = params(itype,jtype,ktype).c;
        const double d = params(itype,jtype,ktype).d;
        const double h = params(itype,jtype,ktype).h;

        zeta_ij += zeta(rsq1,rsq2,m,lam3,R,D,gamma,c,d,h,delr1,delr2);
      }

      // pairwise force due to zeta

      const double B = params(itype,jtype,jtype).B;
      const double lam2 = params(itype,jtype,jtype).lam2;
      const double R = params(itype,jtype,jtype).R;
      const double D = params(itype,jtype,jtype).D;
      const double beta = params(itype,jtype,jtype).beta;
      const double n = params(itype,jtype,jtype).n;

      force_zeta(rsq1,zeta_ij,B,lam2,R,D,beta,n,fpair,prefactor,eflag,evdwl);

      if (energy)
        *energy += evdwl;

      if (atom_energy) {
        atom_energy[i] += evdwl;
      }

      if (forces) {
        (*forces)(i,0) += delr1[0]*fpair;
        (*forces)(i,1) += delr1[1]*fpair;
        (*forces)(i,2) += delr1[2]*fpair;
        (*forces)(j,0) -= delr1[0]*fpair;
        (*forces)(j,1) -= delr1[1]*fpair;
        (*forces)(j,2) -= delr1[2]*fpair;
      }

      /* TODO: What does that do????
      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);
      */

      // attractive term via loop over k

      // TODO: restart a loop again?? why?

      for (int kk = 0; kk != n_neigh; ++kk) {
        if (jj == kk) continue; // TODO: This should suffice, but does it?
        int k = use_neighbor_list ? neighbors[kk] : kk;
        ktype = atom_types[k];
        //iparam_ijk = elem2param[itype][jtype][ktype];

        if (distvec) {
          delr2[0] = distvec[kk*3 + 0];
          delr2[1] = distvec[kk*3 + 1];
          delr2[2] = distvec[kk*3 + 2];
        } else {
          delr2[0] = atom_coords(k,0) - xtmp;
          delr2[1] = atom_coords(k,1) - ytmp;
          delr2[2] = atom_coords(k,2) - ztmp;
        }
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        const double cutsq = params(itype,jtype,ktype).cutsq;
        if (rsq2 > cutsq) continue;

        const double R = params(itype,jtype,ktype).R;
        const double D = params(itype,jtype,ktype).D;
        const int m = params(itype,jtype,ktype).m;
        const double lam3 = params(itype,jtype,ktype).lam3;
        const double gamma = params(itype,jtype,ktype).gamma;
        const double c = params(itype,jtype,ktype).c;
        const double d = params(itype,jtype,ktype).d;
        const double h = params(itype,jtype,ktype).h;

        attractive(prefactor, rsq1, rsq2,
                   R, D, m, lam3, gamma, c, d, h,
                   delr1, delr2, fi, fj, fk);

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

        /* TODO: what is this??? do we need it?
        if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);
        */
      }
    }
  }

  /* TODO: Do we need that???  I guess no.
  if (vflag_fdotr) virial_fdotr_compute();
  */
}

/* ---------------------------------------------------------------------- */

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


/* ---------------------------------------------------------------------- */

void PairTersoff::setup()
{
  //int i,j,k,m,n;

  // set elem2param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

/*
  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }
*/

  // compute parameter values derived from inputs
/*
  for (int m = 0; m < nparams; m++) {
    params_[m].cut = params_[m].bigr + params_[m].bigd;
    params_[m].cutsq = params_[m].cut*params_[m].cut;

    params_[m].c1 = pow(2.0*params_[m].powern*1.0e-16,-1.0/params_[m].powern);
    params_[m].c2 = pow(2.0*params_[m].powern*1.0e-8,-1.0/params_[m].powern);
    params_[m].c3 = 1.0/params_[m].c2;
    params_[m].c4 = 1.0/params_[m].c1;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (int m = 0; m < nparams; m++)
    if (params_[m].cut > cutmax) cutmax = params_[m].cut;
*/
}

/* ---------------------------------------------------------------------- */

void PairTersoff::repulsive(double rsq,
                            double lam1, double A, double R, double D,
                            double &fforce,
                            int eflag, double &eng)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;

  r = sqrt(rsq);
  tmp_fc = ters_fc(r, R, D);
  tmp_fc_d = ters_fc_d(r, R, D);
  tmp_exp = exp(-lam1 * r);
  fforce = -A * tmp_exp * (tmp_fc_d - tmp_fc*lam1) / r;
  if (eflag) eng = tmp_fc * A * tmp_exp;
}

/* ---------------------------------------------------------------------- */

double PairTersoff::zeta(double rsqij, double rsqik,
                         int m, double lam3, double R, double D,
                         double gamma, double c, double d, double h,
                         double *delrij, double *delrik)
{
  double rij,rik,costheta,arg,ex_delr;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

  if (m == 3) arg = pow(lam3 * (rij-rik), 3.0);
  else arg = lam3 * (rij-rik); // m == 1

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc(rik,R,D) * ters_gijk(costheta,gamma,c,d,h) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairTersoff::force_zeta(double rsq, double zeta_ij,
                             double B, double lam2, double R, double D,
                             double beta, double n,
                             double &fforce, double &prefactor,
                             int eflag, double &eng)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  fa = ters_fa(r, B, lam2, R, D);
  fa_d = ters_fa_d(r, B, lam2, R, D);
  bij = ters_bij(zeta_ij, beta, n);
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * ters_bij_d(zeta_ij, beta, n);
  if (eflag) eng = 0.5*bij*fa;
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairTersoff::attractive(double prefactor,
                             double rsqij, double rsqik,
                             double R, double D, int m, double lam3,
                             double gamma, double c, double d, double h,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,
                  R, D, m, lam3, gamma, c, d, h,
                  rij_hat,rij,rik_hat,rik,fi,fj,fk);
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fc(double r, double R, double D)
{
  if (r < R-D) return 1.0;
  if (r > R+D) return 0.0;
  return 0.5*(1.0 - sin(pi_2*(r - R)/D));
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fc_d(double r, double R, double D)
{
  if (r < R-D) return 0.0;
  if (r > R+D) return 0.0;
  return -(pi_4/D) * cos(pi_2*(r - R)/D);
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fa(double r,
                            double B, double lam2, double R, double D)
{
  if (r > R + D) return 0.0;
  return -B * exp(-lam2 * r) * ters_fc(r,R,D);
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fa_d(double r,
                              double B, double lam2, double R, double D)
{
  if (r > R + D) return 0.0;
  return B * exp(-lam2 * r) *
    (lam2 * ters_fc(r, R, D) - ters_fc_d(r, R, D));
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_bij(double zeta, double beta, double n)
{
  const double c1 = pow(2.0 * n * 1e-16, -1.0/n);
  const double c2 = pow(2.0 * n * 1e-8, -1.0/n);
  const double c3 = 1.0 / c2;
  const double c4 = 1.0 / c1;
  double tmp = beta * zeta;
  if (tmp > c1) return 1.0/sqrt(tmp);
  if (tmp > c2) return (1.0 - pow(tmp,-n) / (2.0*n))/sqrt(tmp);
  if (tmp < c4) return 1.0;
  if (tmp < c3) return 1.0 - pow(tmp,n)/(2.0*n);
  return pow(1.0 + pow(tmp,n), -1.0/(2.0*n));
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_bij_d(double zeta, double beta, double n)
{
  const double c1 = pow(2.0 * n * 1e-16, -1.0/n);
  const double c2 = pow(2.0 * n * 1e-8, -1.0/n);
  const double c3 = 1.0 / c2;
  const double c4 = 1.0 / c1;
  double tmp = beta * zeta;
  if (tmp > c1) return beta * -0.5*pow(tmp,-1.5);
  if (tmp > c2) return beta * (-0.5*pow(tmp,-1.5) *
                               (1.0 - 0.5*(1.0 +  1.0/(2.0*n)) *
                                pow(tmp,-n)));
  if (tmp < c4) return 0.0;
  if (tmp < c3) return -0.5*beta * pow(tmp, n-1.0);

  double tmp_n = pow(tmp, n);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*n)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void PairTersoff::ters_zetaterm_d(double prefactor,
                                  double R, double D, int m, double lam3,
                                  double gamma, double c, double d, double h,
                                  double *rij_hat, double rij,
                                  double *rik_hat, double rik,
                                  double *dri, double *drj, double *drk)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,R,D);
  dfc = ters_fc_d(rik,R,D);
  if (m == 3) tmp = pow(lam3 * (rij-rik),3.0);
  else tmp = lam3 * (rij-rik); // m == 1

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (m == 3)
    ex_delr_d = 3.0*pow(lam3, 3.0) * pow(rij-rik,2.0)*ex_delr;
  else ex_delr_d = lam3 * ex_delr; // m == 1

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta, gamma, c, d, h);
  gijk_d = ters_gijk_d(cos_theta, gamma, c, d, h);
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
                             double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}
