/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(extalignment,FixExtAlignment)
// clang-format on
#else

#ifndef LMP_FIX_EXTALIGNMENT_H
#define LMP_FIX_EXTALIGNMENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixExtAlignment : public Fix {
 public:
  FixExtAlignment(class LAMMPS *, int, char **);
  ~FixExtAlignment();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  double xvalue,yvalue,zvalue,style_one,kalign;
  int varflag,iregion;
  char *xstr,*ystr,*zstr,*estr;
  char *idregion;
  class Region *region;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle,estyle;
  double fsum[4],fsum_all[4];
  double cosphi;
  int force_flag;
  int ilevel_respa;

  int maxatom, maxatom_energy;
  double **sforce;
};

}

#endif
#endif
