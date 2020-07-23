//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESSFOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file wind_mt_driven.cpp: line driven wind for mass transfer in a binary system
//======================================================================================

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
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../bvals/bvals.hpp"
#include "../utils/utils.hpp"
#include "../globals.hpp"


// wind  bousaries
//void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
	//	  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void TwoPointMass(MeshBlock *pmb, const Real time, const Real dt,
                  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void ParticleAccels(Real (&xi)[3],Real (&vi)[3],Real (&ai)[3]);

void SumGasOnParticleAccels(Mesh *pm, Real (&xi)[3],Real (&ag1i)[3], Real (&ag2i)[3]);

void particle_step(Real dt,Real (&xi)[3],Real (&vi)[3],Real (&ai)[3]);
void kick(Real dt,Real (&xi)[3],Real (&vi)[3],Real (&ai)[3]);
void drift(Real dt,Real (&xi)[3],Real (&vi)[3],Real (&ai)[3]);
int RefinementCondition(MeshBlock *pmb);

void cross(Real (&A)[3],Real (&B)[3],Real (&AxB)[3]);

void WritePMTrackfile(Mesh *pm, ParameterInput *pin);

Real fspline(Real r, Real eps);

void ParticleAccrete(Mesh *pm, Real(&xi)[3],Real(&vi)[3],Real(&mdot), Real(&pdot)[3]);

// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void StreamingOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
			int is, int ie, int js, int je, int ks, int ke, int ngh);

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
		      int is, int ie, int js, int je, int ks, int ke, int ngh);

void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
		      int is, int ie, int js, int je, int ks, int ke, int ngh);


Real massfluxix1(MeshBlock *pmb, int iout);
Real massfluxox1(MeshBlock *pmb, int iout);
Real momr_tot(MeshBlock *pmb, int iout);
//Real momr_ix1(MeshBlock *pmb, int iout);
//Real momr_ox1(MeshBlock *pmb, int iout);
Real momr_source(MeshBlock *pmb, int iout);
//Real pgas_ox1(MeshBlock *pmb, int iout);
//Real pgas_ix1(MeshBlock *pmb, int iout);
Real divrhovv(MeshBlock *pmb, int iout);
Real divpgas(MeshBlock *pmb, int iout);

void KeplerianVel(const Real rad, const Real theta, const Real phi, Real &v1, Real &v2, Real &v3);


// problem parameters which are useful to make global to this file
Real r0, gamma_gas;
Real dfloor;

// global (to this file) problem parameters
Real da,pa; // ambient density, pressure


// companion parameters and parameters for corotating frame
Real GM2, GM1; // point masses
Real rsoft2; // softening length of PM 2
int  include_gas_backreaction, corotating_frame; // flags for output, gas backreaction on EOM, frame choice
int n_particle_substeps; // substepping of particle integration
Real xi[3], vi[3], agas1i[3], agas2i[3]; // cartesian positions/vels of the secondary object, gas->particle acceleration
Real Omega[3];  // vector rotation of the frame
Real ecc, sma;

int particle_accrete;
Real mdot, pdot[3]; // accretion parameters


// for particle output file
Real trackfile_next_time, trackfile_dt;
int  trackfile_number;
Real Ggrav;


// restart simulations
int is_restart, change_setup;




//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // read in some global params (to this file)

  // first non-mode-dependent settings
  //pressure and density of background medium
  pa   = pin->GetOrAddReal("problem","pamb",1.0);
  da   = pin->GetOrAddReal("problem","damb",1.0);
  gamma_gas = pin->GetReal("hydro","gamma");


  // Get parameters for gravitatonal potential of central point mass and companion
  Ggrav = pin->GetOrAddReal("problem","Ggrav",6.67408e-8);
  GM1 = pin->GetOrAddReal("problem","GM1",0.0);
  GM2 = pin->GetOrAddReal("problem","GM2",1.0);
  r0 = pin->GetOrAddReal("problem","r0",1.0);
  corotating_frame = pin->GetInteger("problem","corotating_frame");
  

  // softening of companion gravity
  rsoft2 = pin->GetOrAddReal("problem","rsoft2",0.1);

  // for tracking particle orbits
  trackfile_dt = pin->GetOrAddReal("problem","trackfile_dt",0.01);

  include_gas_backreaction = pin->GetInteger("problem","gas_backreaction");//does the gas changes orbit of companion
  n_particle_substeps = pin->GetInteger("problem","n_particle_substeps");//for the Leap-frog scheme
  particle_accrete = pin->GetInteger("problem","particle_accrete");//companion accretion

  // local vars
  sma = pin->GetOrAddReal("problem","sma",2.0);//semi-major axis
  ecc = pin->GetOrAddReal("problem","ecc",0.0);
  Real Omega_orb, vcirc;


  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));

  change_setup = pin->GetOrAddReal("problem", "change_setup", 0);


  // // enroll the BCs
  // if(mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
  //   EnrollUserBoundaryFunction(OUTER_X1, DiodeOuterX1);
  // }


  //trying steaming bc
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, OutflowInnerX1);
  }
  
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    if (change_setup==0){
      EnrollUserBoundaryFunction(BoundaryFace::outer_x1, StreamingOuterX1);
    }else{
      EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OutflowOuterX1);
    }
  }


  // Enroll a Source Function
  EnrollUserExplicitSourceFunction(TwoPointMass);

  // Enroll AMR
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);

  //Enroll history dump
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0,massfluxix1,"massfluxix1");//mass flux inner
  EnrollUserHistoryOutput(1,massfluxox1,"massfluxox1");//mass flux outer
  //EnrollUserHistoryOutput(2,momr_tot,"momr_tot");//total momentum
  //EnrollUserHistoryOutput(3,divrhovv,"divrhovv");
  //EnrollUserHistoryOutput(4,divpgas,"divpgas");
  //EnrollUserHistoryOutput(3,momr_ix1,"momr_ix1");//rho*v*v inner
  //EnrollUserHistoryOutput(4,momr_ox1,"momr_ox1");//rho*v*v outer
  //EnrollUserHistoryOutput(5,momr_source,"momr_source");//from twopointmass
  //EnrollUserHistoryOutput(6,pgas_ox1,"pgas_ox1");//pagasx1 outer
  //EnrollUserHistoryOutput(7,pgas_ix1,"pgas_ix1");//pagasx1 inner

  // always write at startup
  trackfile_next_time = time;
  trackfile_number = 0;

  
  // allocate MESH data for the particle pos/vel, Omega frame
  // extra ones stores gas source terms in TwoPointMass function
  int blocksizex1 = pin->GetOrAddInteger("meshblock", "nx1", 1);
  int blocksizex2 = pin->GetOrAddInteger("meshblock", "nx2", 1);
  int blocksizex3 = pin->GetOrAddInteger("meshblock", "nx3", 1);
  blocksizex1 += 2*(NGHOST);
  if (blocksizex2 >1) blocksizex2 += 2*(NGHOST);
  if (blocksizex3 >1) blocksizex3 += 2*(NGHOST);

  AllocateRealUserMeshDataField(4);
  ruser_mesh_data[0].NewAthenaArray(3);
  ruser_mesh_data[1].NewAthenaArray(3);
  ruser_mesh_data[2].NewAthenaArray(3);
  ruser_mesh_data[3].NewAthenaArray(10,blocksizex3, blocksizex2, blocksizex1);

  

  //ONLY enter ICs loop if this isn't a restart
  if(time==0){
    // Print out some info
    if (Globals::my_rank==0){
      std::cout << "*** Setting initial conditions for t=0 ***\n";
    }

    // set the initial conditions for the pos/vel of the secondary
    xi[0] = sma*(1.0 + ecc);  // apocenter
    xi[1] = 0.0;
    xi[2] = 0.0;

    vcirc = sqrt((GM1+GM2)/sma);
    Omega_orb = vcirc/sma;

    vi[0] = 0.0;
    vi[1]= sqrt( vcirc*vcirc*(1.0 - ecc)/(1.0 + ecc) ); //v_apocenter
    vi[2] = 0.0;


    // now set the initial condition for Omega
    Omega[0] = 0.0;
    Omega[1] = 0.0;
    Omega[2] = 0.0;

    // In the case of a corotating frame,
    // subtract off the frame velocity and set Omega
    if(corotating_frame == 1){
      Omega[2] = Omega_orb;
      vi[1] -=  Omega[2]*xi[0];
    }

    // save the ruser_mesh_data variables
    for(int i=0; i<3; i++){
      ruser_mesh_data[0](i)  = xi[i];
      ruser_mesh_data[1](i)  = vi[i];
      ruser_mesh_data[2](i)  = Omega[i];
    }


  }else{
    change_setup = pin->GetOrAddReal("problem", "change_setup", 0);
    is_restart=1;
  }

  // Print out some info
  if (Globals::my_rank==0){
    std::cout << "==========================================================\n";
    std::cout << "==========   SIMULATION INFO =============================\n";
    std::cout << "==========================================================\n";
    std::cout << "time =" << time << "\n";
    std::cout << "Ggrav = "<< Ggrav <<"\n";
    std::cout << "gamma = "<< gamma_gas <<"\n";
    std::cout << "GM1 = "<< GM1 <<"\n";
    std::cout << "GM2 = "<< GM2 <<"\n";
    std::cout << "Omega_orb="<< Omega_orb << "\n";
    std::cout << "a = "<< sma <<"\n";
    std::cout << "e = "<< ecc <<"\n";
    std::cout << "P = "<< 6.2832*sqrt(sma*sma*sma/(GM1+GM2)) << "\n";
    std::cout << "rsoft2 ="<<rsoft2<<"\n";
    std::cout << "corotating frame? = "<< corotating_frame<<"\n";
    std::cout << "gas backreaction? = "<< include_gas_backreaction<<"\n";
    std::cout << "particle substeping n="<<n_particle_substeps<<"\n";
    if(time==0){
      std::cout << "==========================================================\n";
      std::cout << "==========   Particle        =============================\n";
      std::cout << "==========================================================\n";
      std::cout << "x ="<<xi[0]<<"\n";
      std::cout << "y ="<<xi[1]<<"\n";
      std::cout << "z ="<<xi[2]<<"\n";
      std::cout << "vx ="<<vi[0]<<"\n";
      std::cout << "vy ="<<vi[1]<<"\n";
      std::cout << "vz ="<<vi[2]<<"\n";
      std::cout << "==========================================================\n";
    }
  }


	return;   // SS: does this have to be here???

} // end





int RefinementCondition(MeshBlock *pmb)
{
  Real mindist=1.e10;
  for(int k=pmb->ks; k<=pmb->ke; k++){

    Real ph= pmb->pcoord->x2v(k);
    Real sin_ph = sin(ph);
    Real cos_ph = cos(ph);

    for(int j=pmb->js; j<=pmb->je; j++) {

      Real z_cyl= pmb->pcoord->x3v(j);

      for(int i=pmb->is; i<=pmb->ie; i++) {

	Real r = pmb->pcoord->x1v(i);

	Real x = r*cos_ph;
	Real y = r*sin_ph;
	Real z_cart = z_cyl; // z in cartesian

	Real dist = std::sqrt(SQR(x-xi[0]) +
			      SQR(y-xi[1]) +
			      SQR(z_cart-xi[2]) );

	mindist = std::min(mindist,dist);
      }
    }
  }
  if(mindist >  3.0*rsoft2) return -1;
  if(mindist <= 3.0*rsoft2) return 1;
}





// Source Function for two point masses
void TwoPointMass(MeshBlock *pmb, const Real time, const Real dt,
		  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

  if(is_restart>0){
    // else this is a restart, read the current particle state
    for(int i=0; i<3; i++){
      xi[i]    = pmb->pmy_mesh->ruser_mesh_data[0](i);
      vi[i]    = pmb->pmy_mesh->ruser_mesh_data[1](i);
      Omega[i] = pmb->pmy_mesh->ruser_mesh_data[2](i);
    }
    // print some info
    if (Globals::my_rank==0){
      std::cout << "*** Setting initial conditions for t>0 ***\n";
      std::cout <<"xi="<<xi[0]<<" "<<xi[1]<<" "<<xi[2]<<"\n";
      std::cout <<"vi="<<vi[0]<<" "<<vi[1]<<" "<<vi[2]<<"\n";
      std::cout <<"Omega="<<Omega[0]<<" "<<Omega[1]<<" "<<Omega[2]<<"\n";
    }
    if (change_setup==1){
      Real vcirc, Omega_orb;
      xi[0] = sma*(1.0 + ecc);  // apocenter
      xi[1] = 0.0;
      xi[2] = 0.0;  
    
      vcirc = sqrt((GM1+GM2)/sma);
      Omega_orb = vcirc/sma;

      vi[0] = 0.0;
      vi[1]= sqrt( vcirc*vcirc*(1.0 - ecc)/(1.0 + ecc) ); //v_apocenter
      vi[2] = 0.0;

      Omega[0] = pmb->pmy_mesh->ruser_mesh_data[2](0);
      Omega[1] = pmb->pmy_mesh->ruser_mesh_data[2](1);
      Omega[2] = pmb->pmy_mesh->ruser_mesh_data[2](2);

      if(corotating_frame == 1){
	Omega[2] = Omega_orb;
	vi[1] -=  Omega[2]*xi[0];
      }

      if (Globals::my_rank==0){
	std::cout << "*** Change initial conditions for restart ***\n";
	std::cout <<"xi="<<xi[0]<<" "<<xi[1]<<" "<<xi[2]<<"\n";
	std::cout <<"vi="<<vi[0]<<" "<<vi[1]<<" "<<vi[2]<<"\n";
	std::cout <<"Omega="<<Omega[0]<<" "<<Omega[1]<<" "<<Omega[2]<<"\n";
      }
    }

    is_restart=0;
    change_setup=0;
  }


  // Gravitational acceleration from orbital motion
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {

	Real r = pmb->pcoord->x1v(i);
	Real ph= pmb->pcoord->x2v(j);
	Real z_cyl = pmb->pcoord->x3v(k);

	Real vr  = prim(IVX,k,j,i);
	Real vph = prim(IVY,k,j,i);
	Real vz_cyl  = prim(IVZ,k,j,i);

	//get some angles
	Real sin_ph = sin(ph);
	Real cos_ph = cos(ph);
	//Real z_r_ang = atan2(z_cyl,r); //XH atan2 inclues the sign
	Real cos_zr = r/sqrt(r*r+z_cyl*z_cyl);
	Real sin_zr = z_cyl/sqrt(r*r+z_cyl*z_cyl);


	// current position of the secondary
	Real x_2 = xi[0];
	Real y_2 = xi[1];
	Real z_2 = xi[2];
	Real d12c = pow(xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2], 1.5);

	// cylindrical coordinates, get local cartesian
	Real x = r*cos_ph;
	Real y = r*sin_ph;
	Real z_cart = z_cyl;


	Real d2  = sqrt(pow(x-x_2, 2) +
			pow(y-y_2, 2) +
			pow(z_cart-z_2, 2) );

	//
	//  COMPUTE ACCELERATIONS
	//
	// PM1
	//Real a_r1 = -GM1/pow(r,2), r^2=r_y^2 + z^2 ;
	// cell volume avg'd version, see pointmass.cpp sourceterm code.
  // SS: is this coord_src1_i_(i) the right one??
	//XH: I think it's fine to use it here, just think coord_src1_i_(i) as a cell-volume averaged 1/r, say we note it as <1/r>
	//    so in spherical polar, it's like GM1*<1/r>/r. But in cylindrical, I guess we need to use either
	// Real a_r1 = -GM1/(r*r+z*z); //for not using cell-volume averaged quantities <1/r>, just use r*r+z*z
	//or 
	//Real a_r1 = -GM1/(pow(1./pmb->pcoord->coord_src1_i_(i),2)+z_cyl*z_cyl); //use cell-volume averaged r, (1./<1/r>)^2+z^2, 
        Real a_r1 = -GM1*pmb->pcoord->coord_src1_i_(i)/r;

	// PM2 gravitational accels in cartesian coordinates
	Real a_x = - GM2 * fspline(d2,rsoft2) * (x-x_2);
	Real a_y = - GM2 * fspline(d2,rsoft2) * (y-y_2);
	Real a_z_cart = - GM2 * fspline(d2,rsoft2) * (z_cart-z_2);

	// add the correction for the orbiting frame (relative to the COM)
	a_x += -  GM2 / d12c * x_2;
	a_y += -  GM2 / d12c * y_2;
	a_z_cart += -  GM2 / d12c * z_2;

	//store net external acceleration to user variable
	Real a_r_net = cos_ph*a_x + sin_ph*a_y;
	Real a_ph_net = -sin_ph*a_x + cos_ph*a_y;
	Real a_z_net = a_z_cart;
	//add m1
	a_r_net += a_r1*cos_zr;
	a_z_net += a_r1*sin_zr;

	Real den = prim(IDN,k,j,i);

        pmb->pmy_mesh->ruser_mesh_data[3](7,k,j,i) = den*a_r_net;
	pmb->pmy_mesh->ruser_mesh_data[3](8,k,j,i) = den*a_ph_net;
	pmb->pmy_mesh->ruser_mesh_data[3](9,k,j,i) = den*a_z_net;

	if(corotating_frame == 1){
	  // distance from the origin in cartesian (vector)
	  Real rxyz[3];
	  rxyz[0] = x;
	  rxyz[1] = y;
	  rxyz[2] = z_cart;

	  // get the cartesian velocities from the cylindrical (vector)
	  Real vgas[3];
	  vgas[0] = cos_ph*vr - sin_ph*vph;
	  vgas[1] = sin_ph*vr + cos_ph*vph;
	  vgas[2] = vz_cyl;

	  // add the centrifugal and coriolis terms

	  // centrifugal
	  Real Omega_x_r[3], Omega_x_Omega_x_r[3];

	  
	  cross(Omega,rxyz,Omega_x_r);
	  cross(Omega,Omega_x_r,Omega_x_Omega_x_r);
	  
	  a_x += - Omega_x_Omega_x_r[0];
	  a_y += - Omega_x_Omega_x_r[1];
	  a_z_cart += - Omega_x_Omega_x_r[2];

	  // coriolis
	  Real Omega_x_v[3];
	  cross(Omega,vgas,Omega_x_v);
	  
	  a_x += -2.0*Omega_x_v[0];
	  a_y += -2.0*Omega_x_v[1];
	  a_z_cart += -2.0*Omega_x_v[2];
	}

	// add the gas acceleration of the frame of ref
	if(include_gas_backreaction == 1){
	  a_x += -agas1i[0];
	  a_y += -agas1i[1];
	  a_z_cart += -agas1i[2];
	}
	
	// convert back to cylindrical 
	Real a_r  = cos_ph*a_x + sin_ph*a_y;
	Real a_ph = -sin_ph*a_x + cos_ph*a_y;
	Real a_z_cyl  = a_z_cart;
	
	// add the PM1 accel
	a_r += a_r1*cos_zr;
	a_z_cyl += a_r1*sin_zr;

	//
	// ADD SOURCE TERMS TO THE GAS MOMENTA/ENERGY
	//


	Real src_1 = dt*den*a_r;
	Real src_2 = dt*den*a_ph;
	Real src_3 = dt*den*a_z_cyl;

	// add the source term to the momenta  (source = - rho * a)
	
	//store source terms
	pmb->pmy_mesh->ruser_mesh_data[3](0,k,j,i) = src_1/dt; //changes to momentum 
	pmb->pmy_mesh->ruser_mesh_data[3](1,k,j,i) = src_2/dt;
	pmb->pmy_mesh->ruser_mesh_data[3](2,k,j,i) = src_3/dt;

	cons(IM1,k,j,i) += src_1;
	cons(IM2,k,j,i) += src_2;
	cons(IM3,k,j,i) += src_3;

	
	pmb->pmy_mesh->ruser_mesh_data[3](3,k,j,i) = 0.0;
	if (NON_BAROTROPIC_EOS) {
	  // update the energy (source = - rho v dot a
	  pmb->pmy_mesh->ruser_mesh_data[3](3,k,j,i) += (src_1/den * 0.5*(pmb->phydro->flux[X1DIR](IDN,k,j,i) + pmb->phydro->flux[X1DIR](IDN,k,j,i+1)))/dt;
	  cons(IEN,k,j,i) += src_1/den * 0.5*(pmb->phydro->flux[X1DIR](IDN,k,j,i) + pmb->phydro->flux[X1DIR](IDN,k,j,i+1));
	  pmb->pmy_mesh->ruser_mesh_data[3](3,k,j,i) += (src_2*prim(IVY,k,j,i) + src_3*prim(IVZ,k,j,i))/dt;
	  cons(IEN,k,j,i) += src_2*prim(IVY,k,j,i) + src_3*prim(IVZ,k,j,i);
	}

	pmb->pmy_mesh->ruser_mesh_data[3](4,k,j,i) = src_1/dt/prim(IDN,k,j,i);
	pmb->pmy_mesh->ruser_mesh_data[3](5,k,j,i) = src_2/dt/prim(IDN,k,j,i);
	pmb->pmy_mesh->ruser_mesh_data[3](6,k,j,i) = src_3/dt/prim(IDN,k,j,i);

      }
    }
  } // end loop over cells


}

//========================================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in MeshBlock class.  Can also be
//  used to initialize variables which are global to other functions in this file.
//  Called in MeshBlock constructor before ProblemGenerator.
//========================================================================================
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  AllocateUserOutputVariables(14); //store two point mass function
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){
  for(int k=ks; k<=ke; k++){
    for(int j=js; j<=je; j++){
      for(int i=is; i<=ie; i++){
	//Storing everything we added to RHS of hydro eqs in the user defined source TwoPointMass
	user_out_var(0,k,j,i) = pmy_mesh->ruser_mesh_data[3](0,k,j,i);//cons, fext_r
	user_out_var(1,k,j,i) = pmy_mesh->ruser_mesh_data[3](1,k,j,i);//cons, fext_theta
	user_out_var(2,k,j,i) = pmy_mesh->ruser_mesh_data[3](2,k,j,i);//cons, fext_z
	user_out_var(3,k,j,i) = pmy_mesh->ruser_mesh_data[3](4,k,j,i);//prim, aext_r
	user_out_var(4,k,j,i) = pmy_mesh->ruser_mesh_data[3](5,k,j,i);//prim, aext_theta
	user_out_var(5,k,j,i) = pmy_mesh->ruser_mesh_data[3](6,k,j,i);//prim, aext_z
	user_out_var(6,k,j,i) = pmy_mesh->ruser_mesh_data[3](7,k,j,i);//prim, fext_r, no corotate
	user_out_var(7,k,j,i) = pmy_mesh->ruser_mesh_data[3](8,k,j,i);//prim, fext_theta, no corotate
	user_out_var(8,k,j,i) = pmy_mesh->ruser_mesh_data[3](9,k,j,i);//prim, fext_z, no corotate


	//for angular momentum budget, only output zonal averaged quantities in <>, and volume
	Real vcirc = sqrt((GM1+GM2)/sma);
	Real Omega_orb = vcirc/sma;
	Real rad = pcoord->x2v(j);
	Real v_kep = sqrt(GM1/rad);
	Real deltav = phydro->w(IVY,k,j,i) + rad*Omega_orb - v_kep;
	Real rho = phydro->u(IDN,k,j,i);
	//dAMdt = d(rho*R*deltav)/dt
	user_out_var(9,k,j,i) = rho*rad*deltav*pcoord->GetCellVolume(k,j,i); 
	//AM_Mdot = <R*rho*vr>*d(R*vk)/dR, only output <R*rho*vr>
	user_out_var(10,k,j,i) = rad*rho*phydro->w(IVX,k,j,i)*pcoord->GetCellVolume(k,j,i); 
	//AM_FH = -d(R*R*<rho*vr*deltav>)/dr
	user_out_var(11,k,j,i) = rho*phydro->w(IVX,k,j,i)*deltav*pcoord->GetCellVolume(k,j,i); 
	//T(R) = <R cross Fext> R
	user_out_var(12,k,j,i) = pmy_mesh->ruser_mesh_data[3](8,k,j,i)*rad*pcoord->GetCellVolume(k,j,i);
	//cell volume
	user_out_var(13,k,j,i) = pcoord->GetCellVolume(k,j,i);
      }
    }
  }

}


void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);

	// local vars for background - not implemented
  //Real den, pres;


	// Prepare index bounds including ghost cells
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
	  jl -= (NGHOST);
	  ju += (NGHOST);
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
	  kl -= (NGHOST);
	  ku += (NGHOST);
  }

  Real rho_floor = pin->GetReal("hydro", "dfloor");
  if (NON_BAROTROPIC_EOS){
    rho_floor = 1.0e-5;
  }
 
  Real press_init = 1.0e-4;

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // compute initial conditions in cylindrical coordinates       
	Real r_local = pcoord->x1v(i);
        phydro->u(IDN,k,j,i) = rho_floor; 
        phydro->u(IM1,k,j,i) = 0.0;
	phydro->u(IM2,k,j,i) = rho_floor*pcoord->x1v(i)*(sqrt(GM1/pow(pcoord->x1v(i),3))-sqrt((GM1+GM2)/pow(sma,3)));
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = press_init/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}



//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief Function called once every time step for user-defined work.
//======================================================================================

void Mesh::UserWorkInLoop(){

  ParameterInput *pin;
  // first let the particle accrete
  if(particle_accrete>0){
    ParticleAccrete(pblock->pmy_mesh,xi,vi,mdot,pdot);
  }


  Real ai[3];
  // ONLY ON THE FIRST CALL TO THIS FUNCTION
  // (NOTE: DOESN'T WORK WITH RESTARTS)
  if(ncycle==0){
    // kick the initial conditions back a half step (v^n-1/2)

    // first sum the gas accel if needed
    if(include_gas_backreaction == 1){
      SumGasOnParticleAccels(pblock->pmy_mesh, xi,agas1i,agas2i);
    }

    ParticleAccels(xi,vi,ai);
    kick(-0.5*dt,xi,vi,ai);
  }

  // EVOLVE THE ORBITAL POSITION OF THE SECONDARY
  // do this on rank zero, then broadcast
  if (Globals::my_rank == 0){
    for (int ii=1; ii<=n_particle_substeps; ii++) {
      // add the particle acceleration to ai
      ParticleAccels(xi,vi,ai);
      // advance the particle
      particle_step(dt/n_particle_substeps,xi,vi,ai);
    }
  }

#ifdef MPI_PARALLEL
  // broadcast the position update from proc zero
  MPI_Bcast(xi,3,MPI_ATHENA_REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(vi,3,MPI_ATHENA_REAL,0,MPI_COMM_WORLD);
#endif

  // update the ruser_mesh_data variables
  for(int i=0; i<3; i++){
    ruser_mesh_data[0](i)  = xi[i];
    ruser_mesh_data[1](i)  = vi[i];
    ruser_mesh_data[2](i)  = Omega[i];
  }

  // sum the gas->part accel for the next step
  if(include_gas_backreaction == 1){
    SumGasOnParticleAccels(pblock->pmy_mesh, xi,agas1i,agas2i);

  }else if(time >= trackfile_next_time){
    SumGasOnParticleAccels(pblock->pmy_mesh, xi,agas1i,agas2i);
  }


  // write the output to the trackfile
  if(time >= trackfile_next_time){
    WritePMTrackfile(pblock->pmy_mesh,pin);
  }

}


void WritePMTrackfile(Mesh *pm, ParameterInput *pin){

  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign("pm_trackfile.dat");

    // open file for output
    FILE *pfile;
    std::stringstream msg;
    if((pfile = fopen(fname.c_str(),"a")) == NULL){
      msg << "### FATAL ERROR in function [WritePMTrackfile]" << std::endl
          << "Output file '" << fname << "' could not be opened";
      throw std::runtime_error(msg.str().c_str());
    }

    if(trackfile_number==0){
      fprintf(pfile,"#    ncycle     ");
      fprintf(pfile,"time           ");
      fprintf(pfile,"dt             ");
      fprintf(pfile,"x              ");
      fprintf(pfile,"y              ");
      fprintf(pfile,"z              ");
      fprintf(pfile,"vx             ");
      fprintf(pfile,"vy             ");
      fprintf(pfile,"vz             ");
      fprintf(pfile,"agas1x          ");
      fprintf(pfile,"agas1y          ");
      fprintf(pfile,"agas1z          ");
      fprintf(pfile,"agas2x          ");
      fprintf(pfile,"agas2y          ");
      fprintf(pfile,"agas2z          ");
      fprintf(pfile,"mdot            ");
      fprintf(pfile,"pdotx           ");
      fprintf(pfile,"pdoty           ");
      fprintf(pfile,"pdotz           ");
      fprintf(pfile,"\n");
    }


    // write the data line
    fprintf(pfile,"%20i",pm->ncycle);
    fprintf(pfile,"%20.6e",pm->time);
    fprintf(pfile,"%20.6e",pm->dt);
    fprintf(pfile,"%20.6e",xi[0]);
    fprintf(pfile,"%20.6e",xi[1]);
    fprintf(pfile,"%20.6e",xi[2]);
    fprintf(pfile,"%20.6e",vi[0]);
    fprintf(pfile,"%20.6e",vi[1]);
    fprintf(pfile,"%20.6e",vi[2]);
    fprintf(pfile,"%20.6e",agas1i[0]);
    fprintf(pfile,"%20.6e",agas1i[1]);
    fprintf(pfile,"%20.6e",agas1i[2]);
    fprintf(pfile,"%20.6e",agas2i[0]);
    fprintf(pfile,"%20.6e",agas2i[1]);
    fprintf(pfile,"%20.6e",agas2i[2]);
    fprintf(pfile,"%20.6e",mdot);
    fprintf(pfile,"%20.6e",pdot[0]);
    fprintf(pfile,"%20.6e",pdot[1]);
    fprintf(pfile,"%20.6e",pdot[2]);
    fprintf(pfile,"\n");

    // close the file
    fclose(pfile);

  } // end rank==0

  // increment counters
  trackfile_number++;
  trackfile_next_time += trackfile_dt;



  return;
}



void particle_step(Real dt,Real (&xi)[3],Real (&vi)[3],Real (&ai)[3]){
  // Leapfrog algorithm (KDK)

  // kick a full step
  kick(dt,xi,vi,ai);

  // drift a full step
  drift(dt,xi,vi,ai);

}

// kick the velocities dt using the accelerations given in ai
void kick(Real dt,Real (&xi)[3],Real (&vi)[3],Real (&ai)[3]){
  for (int i = 0; i < 3; i++){
    vi[i] += dt*ai[i];
  }
}

// drift the velocities dt using the velocities given in vi
void drift(Real dt,Real (&xi)[3],Real (&vi)[3],Real (&ai)[3]){
  for (int i = 0; i < 3; i++){
    xi[i] += dt*vi[i];
  }
}

void ParticleAccels(Real (&xi)[3],Real (&vi)[3],Real (&ai)[3]){

  Real d = sqrt(xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2]);

  // fill in the accelerations for the orbiting frame
  for (int i = 0; i < 3; i++){
    ai[i] = - GM1/pow(d,3) * xi[i] - GM2/pow(d,3) * xi[i];
  }

  // IF WE'RE IN A ROTATING FRAME
  if(corotating_frame == 1){
    Real Omega_x_r[3],Omega_x_Omega_x_r[3], Omega_x_v[3];

    // compute cross products
    cross(Omega,xi,Omega_x_r);
    cross(Omega,Omega_x_r,Omega_x_Omega_x_r);

    cross(Omega,vi,Omega_x_v);

    // fill in the accelerations for the rotating frame
    for (int i = 0; i < 3; i++){
      ai[i] += -Omega_x_Omega_x_r[i];
      ai[i] += -2.0*Omega_x_v[i];
    }
  }

  // add the gas acceleration to ai
  if(include_gas_backreaction == 1){
    for (int i = 0; i < 3; i++){
      ai[i] += -agas1i[i]+agas2i[i];
    }
  }

}



Real fspline(Real r, Real eps){
  // Hernquist & Katz 1989 spline kernel F=-GM r f(r,e) EQ A2
  Real u = r/eps;
  Real u2 = u*u;

  if (u<1.0){
    return pow(eps,-3) * (4./3. - 1.2*pow(u,2) + 0.5*pow(u,3) );
  } else if(u<2.0){
    return pow(r,-3) * (-1./15. + 8./3.*pow(u,3) - 3.*pow(u,4) + 1.2*pow(u,5) - 1./6.*pow(u,6));
  } else{
    return pow(r,-3);
  }
  
}




void ParticleAccrete(Mesh *pm, Real(&xi)[3],Real(&vi)[3], Real(&mdot), Real(&pdot)[3] ){

  // start by setting accelerations / positions to zero
  mdot = 0.0;
  for (int ii = 0; ii < 3; ii++){
    pdot[ii] = 0.0;
  }

  Real mshell = 0.0;
  Real Vshell = 0.0;

  MeshBlock *pmb=pm->pblock;
  Real dt = pm->dt;
  AthenaArray<Real> vol;

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  vol.NewAthenaArray(ncells1);

  while (pmb != NULL) {
    Hydro *phyd = pmb->phydro;

    // Sum history variables over cells.  Note ghost cells are never included in sums
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
	pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);
	for (int i=pmb->is; i<=pmb->ie; ++i) {
	  //coordinates
	  Real r = pmb->pcoord->x1v(i);
	  Real th= pmb->pcoord->x2v(j);
	  Real ph= pmb->pcoord->x3v(k);

	  //get some angles
	  Real sin_th = sin(th);
	  Real cos_th = cos(th);
	  Real sin_ph = sin(ph);
	  Real cos_ph = cos(ph);

	  // spherical polar coordinates, get local cartesian
	  Real x = r*sin_th*cos_ph;
	  Real y = r*sin_th*sin_ph;
	  Real z = r*cos_th;

	  // current position of the secondary
	  Real x_2 = xi[0];
	  Real y_2 = xi[1];
	  Real z_2 = xi[2];

	  Real d2 = sqrt(pow(x-x_2, 2) +
			 pow(y-y_2, 2) +
			 pow(z-z_2, 2) );

	  // conditions just outside sink
	  if((d2>rsoft2) && (d2<2.0*rsoft2)){
	    // add to the local sum
	    mshell += vol(i) * phyd->u(IDN,k,j,i);
	    Vshell += vol(i);
	  }

	}
      }
    }//end loop over cells
  pmb=pmb->next;
  } // end loop over blocks

#ifdef MPI_PARALLEL
  // sum over all ranks
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &mshell, 1, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &Vshell, 1, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
   } else {
    MPI_Reduce(&mshell,&mshell,1, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(&Vshell,&Vshell,1, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
   }

  // broadcast the result
  MPI_Bcast(&mshell,1,MPI_ATHENA_REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(&Vshell,1,MPI_ATHENA_REAL,0,MPI_COMM_WORLD);
 #endif



  // SECOND LOOP: apply sink
  Real rho_sink = 0.1 * mshell/Vshell;
  Real pres_sink = 0.5*rho_sink;

  if (Globals::my_rank == 0) {
    std::cout<<"rho_sink="<<rho_sink<<"\n";
  }

  pmb=pm->pblock;
  while (pmb != NULL) {
    Hydro *phyd = pmb->phydro;

    // Sum history variables over cells.  Note ghost cells are never included in sums
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
	pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);
	for (int i=pmb->is; i<=pmb->ie; ++i) {
	  //coordinates
	  Real r = pmb->pcoord->x1v(i);
	  Real th= pmb->pcoord->x2v(j);
	  Real ph= pmb->pcoord->x3v(k);

	  //get some angles
	  Real sin_th = sin(th);
	  Real cos_th = cos(th);
	  Real sin_ph = sin(ph);
	  Real cos_ph = cos(ph);

	  // spherical polar coordinates, get local cartesian
	  Real x = r*sin_th*cos_ph;
	  Real y = r*sin_th*sin_ph;
	  Real z = r*cos_th;

	  // current position of the secondary
	  Real x_2 = xi[0];
	  Real y_2 = xi[1];
	  Real z_2 = xi[2];

	  Real d2 = sqrt(pow(x-x_2, 2) +
			 pow(y-y_2, 2) +
			 pow(z-z_2, 2) );

	  // ADD a sink BC near the secondary mass
	  if(d2<rsoft2){
	    // get the cartesian velocities from the spherical (vector)
	    Real vr  = phyd->u(IM1,k,j,i)/phyd->u(IDN,k,j,i);
	    Real vth = phyd->u(IM2,k,j,i)/phyd->u(IDN,k,j,i);
	    Real vph = phyd->u(IM3,k,j,i)/phyd->u(IDN,k,j,i);

	    Real vgas[3];
	    vgas[0] = sin_th*cos_ph*vr + cos_th*cos_ph*vth - sin_ph*vph;
	    vgas[1] = sin_th*sin_ph*vr + cos_th*sin_ph*vth + cos_ph*vph;
	    vgas[2] = cos_th*vr - sin_th*vth;

	    // cell mass dm
	    Real dm = vol(i) * (phyd->u(IDN,k,j,i)-rho_sink);
	    // accreted momentum depends on velocity difference between particle/gas (cartesian)
	    Real dp1 = dm*(vgas[0]-vi[0]);
	    Real dp2 = dm*(vgas[1]-vi[1]);
	    Real dp3 = dm*(vgas[2]-vi[2]);

	    // reset values within the "sink"
	    phyd->u(IDN,k,j,i) = rho_sink;
	    phyd->u(IPR,k,j,i) = pres_sink;
	    phyd->u(IM1,k,j,i) = 0.0;
	    phyd->u(IM2,k,j,i) = 0.0;
	    phyd->u(IM3,k,j,i) = 0.0;
	    phyd->u(IEN,k,j,i) = pres_sink/(gamma_gas-1);

	    // add to the local sums
	    mdot += dm;
	    pdot[0] += dp1;
	    pdot[1] += dp2;
	    pdot[2] += dp3;
	  }

	}
      }
    }//end loop over cells
  pmb=pmb->next;
  } // end loop over blocks

#ifdef MPI_PARALLEL
  // sum over all ranks
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &mdot, 1, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, pdot, 3, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&mdot,&mdot,1, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(pdot,pdot,3, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
  }

  // divide by dt, then broadcast the result
  mdot=mdot/dt;
  pdot[0]=pdot[0]/dt;
  pdot[1]=pdot[1]/dt;
  pdot[2]=pdot[2]/dt;

  MPI_Bcast(&mdot,1,MPI_ATHENA_REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(pdot,3,MPI_ATHENA_REAL,0,MPI_COMM_WORLD);
#endif


}


void SumGasOnParticleAccels(Mesh *pm, Real (&xi)[3],Real (&ag1i)[3],Real (&ag2i)[3]){

  // start by setting accelerations / positions to zero
  for (int ii = 0; ii < 3; ii++){
    ag1i[ii] = 0.0;
    ag2i[ii] = 0.0;
  }

  MeshBlock *pmb=pm->pblock;
  AthenaArray<Real> vol;

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  vol.NewAthenaArray(ncells1);

  while (pmb != NULL) {
    Hydro *phyd = pmb->phydro;

    // Sum history variables over cells.  Note ghost cells are never included in sums
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
	pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);
	for (int i=pmb->is; i<=pmb->ie; ++i) {
	  // cell mass dm
	  Real dm = vol(i) * phyd->u(IDN,k,j,i);

	  //coordinates
	  Real r = pmb->pcoord->x1v(i);
	  Real th= pmb->pcoord->x2v(j);
	  Real ph= pmb->pcoord->x3v(k);

	    //get some angles
	  Real sin_th = sin(th);
	  Real cos_th = cos(th);
	  Real sin_ph = sin(ph);
	  Real cos_ph = cos(ph);

	  // spherical polar coordinates, get local cartesian
	  Real x = r*sin_th*cos_ph;
	  Real y = r*sin_th*sin_ph;
	  Real z = r*cos_th;

	  // current position of the secondary
	  Real x_2 = xi[0];
	  Real y_2 = xi[1];
	  Real z_2 = xi[2];

	  Real d2 = sqrt(pow(x-x_2, 2) +
			 pow(y-y_2, 2) +
			 pow(z-z_2, 2) );

	  Real d1c = pow(r,3);

	   // gravitational accels in cartesian coordinates

	  ag1i[0] += Ggrav*dm/d1c * x;
	  ag1i[1] += Ggrav*dm/d1c * y;
	  ag1i[2] += Ggrav*dm/d1c * z;

	  ag2i[0] += Ggrav*dm * fspline(d2,rsoft2) * (x-x_2);
	  ag2i[1] += Ggrav*dm * fspline(d2,rsoft2) * (y-y_2);
	  ag2i[2] += Ggrav*dm * fspline(d2,rsoft2) * (z-z_2);

	}
      }
    }//end loop over cells
    pmb=pmb->next;
  }//end loop over meshblocks

#ifdef MPI_PARALLEL
  // sum over all ranks
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, ag1i, 3, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, ag2i, 3, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
  } else {
    MPI_Reduce(ag1i,ag1i,3, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(ag2i,ag2i,3, MPI_ATHENA_REAL, MPI_SUM, 0,MPI_COMM_WORLD);
  }

  // and broadcast the result
  MPI_Bcast(ag1i,3,MPI_ATHENA_REAL,0,MPI_COMM_WORLD);
  MPI_Bcast(ag2i,3,MPI_ATHENA_REAL,0,MPI_COMM_WORLD);
#endif


}



void cross(Real (&A)[3],Real (&B)[3],Real (&AxB)[3]){
  // set the vector AxB = A x B
  AxB[0] = A[1]*B[2] - A[2]*B[1];
  AxB[1] = A[2]*B[0] - A[0]*B[2];
  AxB[2] = A[0]*B[1] - A[1]*B[0];
}





// //--------------------------------------------------------------------------------------
// //! \fn void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
// //                         FaceField &b, Real time, Real dt,
// //                         int is, int ie, int js, int je, int ks, int ke)
// //  \brief OUTFLOW boundary conditions, outer x1 boundary
//
// void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
// 		  FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
// {
//   // copy hydro variables into ghost zones, don't allow inflow
//   for (int n=0; n<(NHYDRO); ++n) {
//     if (n==(IVX)) {
//       for (int k=ks; k<=ke; ++k) {
// 	for (int j=js; j<=je; ++j) {
// #pragma simd
// 	  for (int i=1; i<=(NGHOST); ++i) {
// 	    prim(IVX,k,j,ie+i) =  std::max( 0.0, prim(IVX,k,j,(ie-i+1)) );  // positive velocities only
// 	  }
// 	}}
//     } else {
//       for (int k=ks; k<=ke; ++k) {
// 	for (int j=js; j<=je; ++j) {
// #pragma simd
// 	  for (int i=1; i<=(NGHOST); ++i) {
// 	    prim(n,k,j,ie+i) = prim(n,k,j,(ie-i+1));
// 	  }
// 	}}
//     }
//   }
//
//
//   // copy face-centered magnetic fields into ghost zones
//   if (MAGNETIC_FIELDS_ENABLED) {
//     for (int k=ks; k<=ke; ++k) {
//       for (int j=js; j<=je; ++j) {
// #pragma simd
// 	for (int i=1; i<=(NGHOST); ++i) {
// 	  b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
// 	}
//       }}
//
//     for (int k=ks; k<=ke; ++k) {
//       for (int j=js; j<=je+1; ++j) {
// #pragma simd
// 	for (int i=1; i<=(NGHOST); ++i) {
// 	  b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
// 	}
//       }}
//
//     for (int k=ks; k<=ke+1; ++k) {
//       for (int j=js; j<=je; ++j) {
// #pragma simd
// 	for (int i=1; i<=(NGHOST); ++i) {
// 	  b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
// 	}
//       }}
//   }
//
//   return;
// }



void StreamingOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
		      int is, int ie, int js, int je, int ks, int ke, int ngh){
  int L1flag = 0;
  Real local_dens = 1.0;
  Real local_vr = -0.01;
  
  Real local_press = 0.01;
  Real local_cs = 0.1;				
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);

  Real vcirc = sqrt((GM1+GM2)/sma);
  Real Omega_orb = vcirc/sma;

  for (int k=ks; k<=ke; ++k) {//z
    for (int j=js; j<=je; ++j) {//phi
      Real phi_coord = pco->x2v(j);
      for (int i=1; i<=(NGHOST); ++i) {//R

	if (fabs(phi_coord)<=0.1){// if within L1 point

	  prim(IDN,k,j,ie+i) = local_dens;
	  prim(IVX,k,j,ie+i) = local_vr;
	  prim(IVY,k,j,ie+i) = 0.0; //pco->coord(ie+i)*sqrt(GM1/pow(sma,3)); //since we are in the corotating frame
	  prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
	  
	  if (NON_BAROTROPIC_EOS) {
	    prim(IPR,k,j,ie+i) = local_press;
	  }
	  
	}else{//one-direction outflow
	  prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
	  prim(IVX,k,j,ie+i) = std::max(0.0, prim(IVX,k,j,ie));
	  prim(IVY,k,j,ie+i) = pco->x1v(ie+i)*(sqrt(GM1/pow(pco->x1v(ie+i),3))-sqrt((GM1+GM2)/pow(sma,3)));
	  prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
	  if (NON_BAROTROPIC_EOS){
	    prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie);
	  }
	}
	


      }//end R
    }//end theta
  }//end Phi


}


void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
		    int is, int ie, int js, int je, int ks, int ke, int ngh){
  for (int k=ks; k<=ke; ++k) {//Phi
    for (int j=js; j<=je; ++j) {//theta
      for (int i=1; i<=(NGHOST); ++i) {//R
	prim(IDN,k,j,is-i) = prim(IDN,k,j,is);
	prim(IVX,k,j,is-i) = std::min(0.0, prim(IVX,k,j,is));
	prim(IVY,k,j,is-i) = pco->x1v(is-i)*(sqrt(GM1/pow(pco->x1v(is-i),3))-sqrt((GM1+GM2)/pow(sma,3)));
	prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
	if (NON_BAROTROPIC_EOS){
	  prim(IPR,k,j,is-i) = prim(IPR,k,j,is);
	}

      }//end R
    }//end theta
  }//end Phi


}


void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
		    int is, int ie, int js, int je, int ks, int ke, int ngh){
  for (int k=ks; k<=ke; ++k) {//z
    for (int j=js; j<=je; ++j) {//phi
      Real phi_coord = pco->x2v(j);
      for (int i=1; i<=(NGHOST); ++i) {//R
	  prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
	  prim(IVX,k,j,ie+i) = std::max(0.0, prim(IVX,k,j,ie));
	  prim(IVY,k,j,ie+i) = pco->x1v(ie+i)*(sqrt(GM1/pow(pco->x1v(ie+i),3))-sqrt((GM1+GM2)/pow(sma,3)));
	  prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
	  if (NON_BAROTROPIC_EOS){
	    prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie);
	  }
      }
    }
  }
}

Real massfluxix1(MeshBlock *pmb, int iout){
  Real massflux = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> face1;
  face1.NewAthenaArray((ie-is)+2*NGHOST+2);

  AthenaArray<Real> x1flux = pmb->phydro->flux[X1DIR];

  if (pmb->pbval->apply_bndry_fn_[BoundaryFace::inner_x1]){
    for (int k=ks; k<=ke; k++){
      for (int j=js; j<=je; j++){
	pmb->pcoord->Face1Area(k , j, is, ie, face1);
	for (int i=is; i<=is; i++){
	  massflux += face1(is)*x1flux(0,k,j,is);//x1flux(0) is the density flux, multiply by volume to get mass
        }

      }
    }
  }

  face1.DeleteAthenaArray();
  x1flux.DeleteAthenaArray();

  return massflux;


}

Real massfluxox1(MeshBlock *pmb, int iout){
  Real massflux = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> face1;
  face1.NewAthenaArray((ie-is)+2*NGHOST+2);

  AthenaArray<Real> x1flux = pmb->phydro->flux[X1DIR];

  if (pmb->pbval->apply_bndry_fn_[BoundaryFace::outer_x1]){
    for (int k=ks; k<=ke; k++){
      for (int j=js; j<=je; j++){
	pmb->pcoord->Face1Area(k , j, is, ie+1, face1);
	for (int i=ie; i<=ie; i++){
	  massflux += face1(ie+1)*x1flux(0,k,j,ie+1);//x1flux(0) is the density flux, multiply by volume to get mass
        }

      }
    }
  }

  face1.DeleteAthenaArray();
  x1flux.DeleteAthenaArray();

  return massflux;


}

/*
//total momentum
Real momr_tot(MeshBlock *pmb, int iout){
  Real px = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  //AthenaArray<Real> crsource = pmb->pcr->cr_gas_source;

  for (int k=ks; k<=ke; k++){
    for (int j=js; j<=je; j++){
      for (int i=is; i<=ie; i++){
	px += pmb->phydro->u(IM1,k,j,i)*pmb->pcoord->GetCellVolume(k,j,i);
      }
    }
  }
  return px;

}


//rho*v*v
Real divrhovv(MeshBlock *pmb, int iout){
  Real MomFlux = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> face1,face2,face2_p1;
  face1.NewAthenaArray((ie-is)+2*NGHOST+2);
  face2.NewAthenaArray((ie-is)+2*NGHOST+2);	
  face2_p1.NewAthenaArray((ie-is)+2*NGHOST+2);

  for (int k=ks; k<=ke; k++){
    for (int j=js; j<=je; j++){
      pmb->pcoord->Face1Area(k , j, is, ie, face1);
      pmb->pcoord->Face2Area(k , j, is, ie, face2);
      pmb->pcoord->Face2Area(k , j+1, is, ie, face2_p1);
      for (int i=is; i<=ie; i++){
	//Real term1 = pmb->phydro->u(IM1,k,j,i)*(pmb->phydro->w(IVY,k,j,i+1)*face1(i+1)-pmb->phydro->w(IVY,k,j,i+1)*face1(i+1));
	//Real term2 = (pmb->phydro->u(IM2,k,j,i)/pmb->pcoord->x1v(i))*(face2_p1(i)*pmb->phydro->w(IVY,k,j+1,i)-face2(i)*pmb->phydro->w(IVY,k,j,i));
	// Real term3 = -pmb->phydro->u(IM2,k,j,i)*pmb->phydro->w(IVY,k,j,i)/pmb->pcoord->x1v(i);
	//MomFlux += (term1+term2+term3);
	 MomFlux += (face1(i+1)*pmb->phydro->u(IDN,k,j,i+1)*pow(pmb->phydro->w(IVX,k,j,i+1),2)-face1(i)*pmb->phydro->u(IDN,k,j,i)*pow(pmb->phydro->w(IVX,k,j,i),2));
      }
    }
  }

  
  return MomFlux;
}


//from twopointmass
Real momr_source(MeshBlock *pmb, int iout){
  Real source_x1 = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
    for (int k=ks; k<=ke; k++){
      for (int j=js; j<=je; j++){
	for (int i=is; i<=is; i++){
	  source_x1 += pmb->pcoord->GetCellVolume(k,j,i)*pmb->pmy_mesh->ruser_mesh_data[3](0,k,j,i);

	}
      }
    }


}

//grad pgas
Real divpgas(MeshBlock *pmb, int iout){
  Real Pgas = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> face1;
  face1.NewAthenaArray((ie-is)+2*NGHOST+2);

  for (int k=ks; k<=ke; k++){
    for (int j=js; j<=je; j++){
      pmb->pcoord->Face1Area(k , j, is, ie, face1);
      for (int i=is; i<=ie; i++){
	Pgas += (face1(is)*pmb->phydro->w(IEN,k,j,is)-face1(ie)*pmb->phydro->w(IEN,k,j,ie));//(area*pgas-area*pgas)/vol * vol
      }
    }
  }

  
  return Pgas;
}
*/

