//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file readin.cpp
//  \brief fscanf a profile

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

//array to store log_rho (400 lines profile)
static AthenaArray<Real> logrhotable;
static int nlogrho;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {

  nlogrho = 212;
  logrhotable.NewAthenaArray(nlogrho);

  FILE *flogrho; 
  if ( (flogrho=fopen("./logrho.txt","r"))==NULL ){
    printf("Open input file error logRho");
  }

  for (int i=0; i<nlogrho; i++){
    fscanf(flogrho, "%lf", &(logrhotable(i))); 
    //%lf is for double, fscanf documentation: https://cplusplus.com/reference/cstdio/fscanf/
    //fscanf(the_file_pointer, format_you_want, file_to_write)
  }

  fclose(flogrho);

  int x1size = mesh_size.nx1+2*NGHOST;

  for (int i =0; i<nlogrho; i++){
    printf("value is :%g\n", logrhotable(i));
  }

  return;
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{

  // free memory
  logrhotable.DeleteAthenaArray();

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

// void MeshBlock::ProblemGenerator(ParameterInput *pin) {

//   // //  Initialize density and momenta
//   // for (int k=ks; k<=ke; ++k) {
//   //   for (int j=js; j<=je; ++j) {
//   //     for (int i=is; i<=ie; ++i) {
//   //       phydro->u(IDN,k,j,i) = 

//   //       phydro->u(IM1,k,j,i) =
//   //       phydro->u(IM2,k,j,i) = 
//   //       phydro->u(IM3,k,j,i) = 
//   //       if (NON_BAROTROPIC_EOS) {
//   //         phydro->u(IEN,k,j,i) = 
//   //       }
//   //     }
//   //   }
//   // }

//   return;
// }

