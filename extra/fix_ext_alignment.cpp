/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include "fix_ext_alignment.h"
#include "atom.h"
#include "atom_masks.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "atom_vec_ellipsoid.h"
#include "neighbor.h"
#include "comm.h"
#include "math_const.h"
#include <iostream>
#include <typeinfo>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace std;

enum{NONE,CONSTANT,EQUAL,ATOM,ANGLE_F,ANGLE_V,ANGLE_U};
#define SMALL 0.001
/* ---------------------------------------------------------------------- */

FixExtAlignment::FixExtAlignment(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idregion(nullptr), region(nullptr), sforce(nullptr)

{
  if (narg < 8) error->all(FLERR,"Illegal fix extalignment command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  energy_global_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;

  xstr = ystr = zstr = nullptr;

  // int style_one;
  if (utils::strmatch(arg[3],"angle_f")) {
      style_one = utils::numeric(FLERR,"4",false,lmp);}
  else if (utils::strmatch(arg[3],"angle_v")) {
      style_one = utils::numeric(FLERR,"5",false,lmp);}
  else if (utils::strmatch(arg[3],"angle_u")) {
      style_one = utils::numeric(FLERR,"6",false,lmp);}
  else
    error->all(FLERR,"Illegal extalignment style");

  xvalue = utils::numeric(FLERR,arg[4],false,lmp);
  xstyle = CONSTANT;
  yvalue = utils::numeric(FLERR,arg[5],false,lmp);
  ystyle = CONSTANT;
  zvalue = utils::numeric(FLERR,arg[6],false,lmp);
  zstyle = CONSTANT;
  kalign = utils::numeric(FLERR,arg[7],false,lmp);

  // optional args

  int iarg = 6;
  force_flag = 0;
  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;

  maxatom = atom->nmax;
  memory->create(sforce,maxatom,4,"addforce:sforce");

  maxatom_energy = 0;
}

/* ---------------------------------------------------------------------- */

FixExtAlignment::~FixExtAlignment()
{
  delete[] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixExtAlignment::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::init()
{
  // check variables

  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR,"fix ext-alignment requires Ellipsoids section");

  if (xstyle == CONSTANT || ystyle == CONSTANT || zstyle == CONSTANT)
    varflag = CONSTANT;

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double v[6];
  double fx,fy,fz;
  int nlocal = atom->nlocal;

  double **torque = atom->torque;
  int *ellipsoid = atom->ellipsoid;
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int newton_bond = force->newton_bond;

  double u0[3],uI[3];
  double values[3], valuesnorm,directionfield[3],normdf,normuI,uIcrossdf[3],uIdotdf,f1[3], aa[4], ab[6];
  double fterm1[3],fterm2[3];
  double *quatI;
  double quatItmp[4];
  double U,theta1, theta2, phi,dtheta1, dtheta2, dphi;
  double inverse;
  double a0,a1,a2, dUda0, dUda1, dUda2;
  double dUdfI[3], dUdfJ[3];
  double torq[3], torqI[3], torqJ[3];
  double GI[3];




  if (update->ntimestep % nevery) return;

  // virial setup

  v_init(vflag);

  // update region if necessary

  if (region) region->prematch();

  // fsum[0] = "potential energy" for added force
  // fsum[123] = force on atoms before extra force added

  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
    double unwrap[3];

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);

        if (style_one == ANGLE_F){
          u0[0] = 1; u0[1] = 0; u0[2] = 0; //f
        }
        else if (style_one == ANGLE_V){
          u0[0] = 0; u0[1] = 1; u0[2] = 0; //v
        }
        else if (style_one == ANGLE_U){
          u0[0] = 0; u0[1] = 0; u0[2] = 1; //u
        }


        quatI = bonus[ellipsoid[i]].quat;
        quatItmp[0] = quatI[0]; quatItmp[1] = quatI[1]; quatItmp[2] = quatI[2]; quatItmp[3] = quatI[3];

        aa[0] = quatItmp[0]*quatItmp[0];
        aa[1] = quatItmp[1]*quatItmp[1];
        aa[2] = quatItmp[2]*quatItmp[2];
        aa[3] = quatItmp[3]*quatItmp[3];

        ab[0] = quatItmp[0]*quatItmp[1];
        ab[1] = quatItmp[0]*quatItmp[2];
        ab[2] = quatItmp[0]*quatItmp[3];
        ab[3] = quatItmp[1]*quatItmp[2];
        ab[4] = quatItmp[1]*quatItmp[3];
        ab[5] = quatItmp[2]*quatItmp[3];


        uI[0] = (aa[0]+aa[1]-aa[2]-aa[3])*u0[0]+
                  ((ab[3]-ab[2])*u0[1]+(ab[1]+ab[4])*u0[2])*2.0;
        uI[1] = (aa[0]-aa[1]+aa[2]-aa[3])*u0[1]+
                  ((ab[2]+ab[3])*u0[0]+(ab[5]-ab[0])*u0[2])*2.0;
        uI[2] = (aa[0]-aa[1]-aa[2]+aa[3])*u0[2]+
                  ((ab[4]-ab[1])*u0[0]+(ab[0]+ab[5])*u0[1])*2.0;

        normuI = sqrt(uI[0] * uI[0] + uI[1] * uI[1] + uI[2] * uI[2]);

        uI[0] /= normuI;
        uI[1] /= normuI;
        uI[2] /= normuI;

        uIdotdf = xvalue*uI[0] + yvalue*uI[1] + zvalue*uI[2];

        if (uIdotdf > 0){
          directionfield[0] = xvalue;
          directionfield[1] = yvalue;
          directionfield[2] = zvalue;

          normdf = sqrt(directionfield[0] * directionfield[0] + directionfield[1] * directionfield[1] + directionfield[2] * directionfield[2]);
        }
        else{
          directionfield[0] = -1*xvalue;
          directionfield[1] = -1*yvalue;
          directionfield[2] = -1*zvalue;

          normdf = sqrt(directionfield[0] * directionfield[0] + directionfield[1] * directionfield[1] + directionfield[2] * directionfield[2]);
        }

        directionfield[0] /= normdf;
        directionfield[1] /= normdf;
        directionfield[2] /= normdf;

        uIcrossdf[0] = uI[1]*directionfield[2] - uI[2]*directionfield[1];
        uIcrossdf[1] = uI[2]*directionfield[0] - uI[0]*directionfield[2];
        uIcrossdf[2] = uI[0]*directionfield[1] - uI[1]*directionfield[0];

        torque[i][0] += kalign*uIcrossdf[0];
        torque[i][1] += kalign*uIcrossdf[1];
        torque[i][2] += kalign*uIcrossdf[2];
      }

  }
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixExtAlignment::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(fsum,fsum_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return fsum_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixExtAlignment::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(fsum,fsum_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return fsum_all[n+1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */
double FixExtAlignment::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
