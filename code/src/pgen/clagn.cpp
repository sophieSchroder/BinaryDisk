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


// wind  boundaries
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

// NEW 7/13
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
void VelProfileCyl(const Real rad, const Real phi, const Real z,
                   Real &v1, Real &v2, Real &v3);


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



// void AGNDiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
// 			int is, int ie, int js, int je, int ks, int ke, int ngh);
//
// void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
// 		      int is, int ie, int js, int je, int ks, int ke, int ngh);
//
// void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
// 		      int is, int ie, int js, int je, int ks, int ke, int ngh);
//
// void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
// 		      int is, int ie, int js, int je, int ks, int ke, int ngh);
//
// void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
// 		      int is, int ie, int js, int je, int ks, int ke, int ngh);


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
Real dfloor, pfloor;

// global (to this file) problem parameters
Real da,pa; // ambient density, pressure


// companion parameters and parameters for corotating frame
Real GM2, GM1; // point masses
Real rsoft2; // softening length of PM 2
int  include_gas_backreaction, corotating_frame; // flags for output, gas backreaction on EOM, frame choice
int  n_particle_substeps; // substepping of particle integration
Real xi[3], vi[3], agas1i[3], agas2i[3]; // cartesian positions/vels of the secondary object, gas->particle acceleration
Real Omega[3];  // vector rotation of the frame
Real ecc, sma, incl; // incl is the angle of inclination of the secondary's orbit
Real eccentric_anomaly; 

int particle_accrete;
Real mdot, pdot[3]; // accretion parameters

// NEW 7/13
Real dslope, p0_over_r0, pslope, rho0;


// for particle output file
Real trackfile_next_time, trackfile_dt;
int  trackfile_number;
Real Ggrav;


// restart simulations
int is_restart; //change_setup;




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

  // NEW 7/13 for disk.cpp add-ins:
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);
  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }


  // local vars
  sma = pin->GetOrAddReal("problem","sma",2.0); //semi-major axis
  ecc = pin->GetOrAddReal("problem","ecc",0.0); //eccentricity
  incl = pin->GetOrAddReal("problem","incl",0.0); //inclination angle
  eccentric_anomaly = pin->GetOrAddReal("problem","eccentric_anomaly",1.571); //determines where orbit starts
  Real Omega_orb, v_circ, r_max, r_min, b, r_sep; //r_max is the maximum distance from the focus of the ellipse
  Real sma_direction_coord, b_direction_coord; //for equation of ellipse, instead of x and y
  Real v_tot, true_anomaly, v_sma_direction_coord, v_b_direction_coord; //to calculate component starting velocities

  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));
  pfloor=pin->GetReal("hydro","pfloor");

  //change_setup = pin->GetOrAddReal("problem", "change_setup", 0);

  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DiskInnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, DiskOuterX2);
  }
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, DiskInnerX3);
  }
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, DiskOuterX3);
  }


  //trying steaming bc
  // if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
  //   EnrollUserBoundaryFunction(BoundaryFace::inner_x1, OutflowInnerX1);
  // }
  //
  // if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
  //   //if (change_setup==0){
  //   //  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, AGNDiskOuterX1);
  //   //}else{
  //     EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OutflowOuterX1);
  //   //  }
  // }
  //
  // if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
  //   EnrollUserBoundaryFunction(BoundaryFace::inner_x3, OutflowInnerX3);
  // }
  //
  // if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
  //   EnrollUserBoundaryFunction(BoundaryFace::outer_x3, OutflowOuterX3);
  // }

  // Enroll (do) a Source Function
  // This is the function that makes the gravitational of SMBH and companion
  EnrollUserExplicitSourceFunction(TwoPointMass);

  // Enroll AMR
  if(adaptive==true) {
    EnrollUserRefinementCondition(RefinementCondition);
  }

  //Enroll history dump
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0,massfluxix1,"massfluxix1");//mass flux inner
  EnrollUserHistoryOutput(1,massfluxox1,"massfluxox1");//mass flux outer

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

  AllocateRealUserMeshDataField(7);
  ruser_mesh_data[0].NewAthenaArray(3);
  ruser_mesh_data[1].NewAthenaArray(3);
  ruser_mesh_data[2].NewAthenaArray(3);
  ruser_mesh_data[3].NewAthenaArray(10,blocksizex3, blocksizex2, blocksizex1);
  //check the gravitational potential
  ruser_mesh_data[4].NewAthenaArray(3,blocksizex3, blocksizex2, blocksizex1);
  ruser_mesh_data[5].NewAthenaArray(6,blocksizex3, blocksizex2, blocksizex1);
  //store the time averaged quantities
  ruser_mesh_data[6].NewAthenaArray(11,blocksizex3, blocksizex2, blocksizex1);


  //ONLY enter Initial conditions loop if this isn't a restart
 if(time==0){
   // Print out some info
   if (Globals::my_rank==0){
     std::cout << "*** Setting initial conditions for t=0 ***\n";
   }

   r_max = sma*(1 + ecc);
   r_min = sma*(1 - ecc);
   b = sqrt(r_max*r_min); // semi-minor axis

   v_circ = sqrt((GM1+GM2)/sma);
   Omega_orb = v_circ/sma;
   
   // in the orbital plane of the semi-major axis
   sma_direction_coord = sma * cos(eccentric_anomaly) + (sma-r_min);
   // in the orbital plane of the semi-minor axis
   b_direction_coord = b * sin(eccentric_anomaly);

   // transforms above ellipse into x,y,z in code with effects of inclination
   xi[0] = sma_direction_coord;
   xi[1] = b_direction_coord * cos(incl);
   xi[2] = b_direction_coord * sin(incl);

   // equations needed to calculate starting velocity:
   // calc distance to starting point from central object
   r_sep = sqrt(pow(xi[0],2) + pow(xi[1],2) + pow(xi[2],2));
   // calc overall velocity at starting point
   v_tot = sqrt((GM1+GM2)*(2/r_sep - 1/sma));
   // calc true anomaly from eccentric anomaly 
   true_anomaly = 2 * atan(sqrt((1 + ecc)/(1 - ecc)) * tan(eccentric_anomaly / 2.0));
   // split components of velocity to components along semi-major and minor axes
   v_sma_direction_coord = v_tot * (-sin(true_anomaly) / 
                                    sqrt(1 + SQR(ecc) + (2 * ecc * cos(true_anomaly))));
   v_b_direction_coord = v_tot * ((ecc + cos(true_anomaly)) / 
                                   sqrt(1 + SQR(ecc) + (2 * ecc * cos(true_anomaly))));
 
   // transforms ellipse velocity components to x,y,z in code with inclination
   vi[0] = v_sma_direction_coord;
   vi[1] = v_b_direction_coord * cos(incl);
   vi[2] = v_b_direction_coord * sin(incl);

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
   is_restart=1;
 }

  // Print out some info
  if (Globals::my_rank==0){
    std::cout << "==========================================================\n";
    std::cout << "==========   SIMULATION INFO =============================\n";
    std::cout << "==========================================================\n";
    std::cout << "time = " << time << "\n";
    std::cout << "Ggrav = "<< Ggrav <<"\n";
    std::cout << "gamma = "<< gamma_gas <<"\n";
    std::cout << "GM1 = "<< GM1 <<"\n";
    std::cout << "GM2 = "<< GM2 <<"\n";
    std::cout << "Omega_orb="<< Omega_orb << "\n";
    std::cout << "a = "<< sma <<"\n";
    std::cout << "e = "<< ecc <<"\n";
    std::cout << "P = "<< 6.2832*sqrt(sma*sma*sma/(GM1+GM2)) << "\n";
    std::cout << "rsoft2 = "<<rsoft2<<"\n";
    std::cout << "incl = "<<incl<<"\n";
    std::cout << "corotating frame? = "<< corotating_frame<<"\n";
    std::cout << "gas backreaction? = "<< include_gas_backreaction<<"\n";
    std::cout << "particle substepping n = "<<n_particle_substeps<<"\n";
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


	return;

} // end



// SD: If we want to use AMR for our setup need to modify this
// Softening condition for the companion
 int RefinementCondition(MeshBlock *pmb)
 {
   Real mindist=1.e10;
   for(int k=pmb->ks; k<=pmb->ke; k++){

     Real ph = pmb->pcoord->x2v(k);
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
   if(mindist <= 3.0*rsoft2) return 1;   // Do refinement
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

    is_restart=0;
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
  // Real a_r1 = -GM1*pmb->pcoord->coord_src1_i_(i)/r; // SD&SS: See about 3d version - need cyl + radial polar
  Real a_r1 = -GM1/(r*(1./pmb->pcoord->coord_src1_i_(i))+z_cyl*z_cyl); //use cell-volume averaged r (1/src_1) and avg r length, r
  Real a_x, a_y, a_z_cart;
  if(pmb->pmy_mesh->time < 6.28){ //let shock move out
    // PM2 gravitational accels in cartesian coordinates
    // before shock has moved out, don't include effect of companion. 2pi = 1 orbit
    // SD: confirm there isn't a better way to be doing this
    a_x = 0.0;
    a_y = 0.0;
    a_z_cart = 0.0;
  } else {
    a_x = - GM2 * fspline(d2,rsoft2) * (x-x_2);
    a_y = - GM2 * fspline(d2,rsoft2) * (y-y_2);
    a_z_cart = - GM2 * fspline(d2,rsoft2) * (z_cart-z_2);
  }
	// add the correction for the orbiting frame (relative to the COM)
	a_x += -  GM2 / d12c * x_2;
	a_y += -  GM2 / d12c * y_2;
	a_z_cart += -  GM2 / d12c * z_2;

	//store net external acceleration to user variable
	Real a_r_net = a_r1*cos_zr + cos_ph*a_x + sin_ph*a_y;
	Real a_ph_net = -sin_ph*a_x + cos_ph*a_y;
	Real a_z_net = a_z_cart + a_r1*sin_zr;

  // get density of cell
	Real den = prim(IDN,k,j,i);

  // return force density to ruser_mesh_data
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

  AllocateUserOutputVariables(38); //store two point mass function
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){

  AthenaArray<Real> &x1flux = phydro->flux[X1DIR];
  AthenaArray<Real> &x2flux = phydro->flux[X2DIR];
  //AthenaArray<Real> face1, face1_m, face1_p;

  //face1.NewAthenaArray((ie-is)+2*NGHOST+2);
  //face1_m.NewAthenaArray((ie-is)+2*NGHOST+2);
  //face1_p.NewAthenaArray((ie-is)+2*NGHOST+2);


  for(int k=ks; k<=ke; k++){
    for(int i=is; i<=ie; i++){

      //cell center variables
      Real rad_c = pcoord->x1v(i);
      Real vkep_c = sqrt(GM1/rad_c);

      //face i+1/2
      Real rad_p = pcoord->x1f(i+1);
      Real vkep_p = sqrt(GM1/rad_p);

      //face i-1/2
      Real rad_m = pcoord->x1f(i);
      Real vkep_m = sqrt(GM1/rad_m);

      for(int j=js; j<=je; j++){

	//Storing everything we added to RHS of hydro eqs in the user defined source TwoPointMass
	user_out_var(0,k,j,i) = pmy_mesh->ruser_mesh_data[3](0,k,j,i);//cons, fext_r
	user_out_var(1,k,j,i) = pmy_mesh->ruser_mesh_data[3](1,k,j,i);//cons, fext_theta
	user_out_var(2,k,j,i) = pmy_mesh->ruser_mesh_data[3](2,k,j,i);//cons, fext_z
	user_out_var(3,k,j,i) = pmy_mesh->ruser_mesh_data[3](7,k,j,i);//cons, fext_r, no corotate
	user_out_var(4,k,j,i) = pmy_mesh->ruser_mesh_data[3](8,k,j,i);//cons, fext_theta, no corotate
	user_out_var(5,k,j,i) = pmy_mesh->ruser_mesh_data[3](9,k,j,i);//cons, fext_z, no corotate

	//AM check, zonal average
        //cell center variables
	Real rho_c = phydro->u(IDN,k,j,i);
	Real vphi_c = phydro->w(IM2,k,j,i);
	Real vr_c = phydro->w(IM1,k,j,i);

	//changing to non-corotating frames
	Real Omega_orb = sqrt((GM1+GM2)/sma)/sma;
	Real v_rotate_c = Omega_orb*rad_c;
	Real v_rotate_p = Omega_orb*rad_p;
	Real v_rotate_m = Omega_orb*rad_m;


	//dAMdt, only show AM = rho*R*(v-vkep)
	user_out_var(6,k,j,i) = rho_c*rad_c*(vphi_c-vkep_c)*pcoord->GetCellVolume(k,j,i);
	//dAMdt, no vkep
	user_out_var(11,k,j,i) = rho_c*rad_c*(vphi_c)*pcoord->GetCellVolume(k,j,i);
	//dAMdt, non-corotating
	user_out_var(34,k,j,i) = rho_c*rad_c*(vphi_c-vkep_c+v_rotate_c)*pcoord->GetCellVolume(k,j,i);

	//AM Mdot
	Real AMMdot = -rad_c*x1flux(IDN,k,j,i)*(vkep_p*pcoord->GetFace1Area(k,j,i+1)- vkep_m*pcoord->GetFace1Area(k,j,i));
	user_out_var(7,k,j,i) = AMMdot;

	//AMTH
	//Riemann solver flux? rho*vr*vphi = x1flux(rho*vphi), rho*vr*vk=x1flux(rho*vr)*vk
	Real AMTH = rad_p*pcoord->GetFace1Area(k,j,i+1)*(x1flux(IM2,k,j,i+1)-x1flux(IDN,k,j,i+1)*vkep_p) - rad_m*pcoord->GetFace1Area(k,j,i)*(x1flux(IM2,k,j,i)-x1flux(IDN,k,j,i)*vkep_m);
	user_out_var(8,k,j,i) = -AMTH;
	//AMTH, no vkep
	Real AMTH_net = rad_p*pcoord->GetFace1Area(k,j,i+1)*x1flux(IM2,k,j,i+1) - rad_m*pcoord->GetFace1Area(k,j,i)*x1flux(IM2,k,j,i);
	user_out_var(12,k,j,i) = -AMTH_net;
	//AMTH, in the non-corotating frame, adding omega*r

	Real AMTH_noncorotate = rad_p*pcoord->GetFace1Area(k,j,i+1)*(x1flux(IM2,k,j,i+1)-x1flux(IDN,k,j,i+1)*vkep_p + x1flux(IDN,k,j,i+1)*v_rotate_p) - rad_m*pcoord->GetFace1Area(k,j,i)*(x1flux(IM2,k,j,i)-x1flux(IDN,k,j,i)*vkep_m +x1flux(IDN,k,j,i)*v_rotate_m);
	user_out_var(35,k,j,i) = AMTH_noncorotate;

	//Torque
	// fect includes Coriolis force
	Real Torque = rad_c*pmy_mesh->ruser_mesh_data[3](1,k,j,i)*pcoord->GetCellVolume(k,j,i);
	// fext has no Coriolis force component
	Real Torque_ =  rad_c*pmy_mesh->ruser_mesh_data[3](8,k,j,i)*pcoord->GetCellVolume(k,j,i);
	user_out_var(9,k,j,i) = Torque;
	user_out_var(10,k,j,i) = Torque_;


	//store grav_phi
	user_out_var(13,k,j,i) = pmy_mesh->ruser_mesh_data[4](0,k,j,i); //gravphi(R, phi)
	user_out_var(14,k,j,i) = pmy_mesh->ruser_mesh_data[4](1,k,j,i); //F_grav_r (R, phi)
	user_out_var(15,k,j,i) = pmy_mesh->ruser_mesh_data[4](2,k,j,i); //F_grav_ph (R, phi)

	//store the time integrated AMMdot, AMTH and Torque
	user_out_var(16,k,j,i) = pmy_mesh->ruser_mesh_data[5](0,k,j,i); //AMMdot
	user_out_var(17,k,j,i) = pmy_mesh->ruser_mesh_data[5](1,k,j,i); //AMTH
	user_out_var(18,k,j,i) = pmy_mesh->ruser_mesh_data[5](2,k,j,i); //AMTH_net
	user_out_var(37,k,j,i) = pmy_mesh->ruser_mesh_data[5](5,k,j,i); //AMTH non coratating
	user_out_var(19,k,j,i) = pmy_mesh->ruser_mesh_data[5](3,k,j,i); //torque
	user_out_var(20,k,j,i) = pmy_mesh->ruser_mesh_data[5](4,k,j,i); //net torque (no Coriolis)


	//store the time integrated quantities for radial profile
	user_out_var(21,k,j,i) = pmy_mesh->ruser_mesh_data[6](0,k,j,i); //Mdot
	user_out_var(22,k,j,i) = pmy_mesh->ruser_mesh_data[6](1,k,j,i); //alpha_eff
	user_out_var(23,k,j,i) = pmy_mesh->ruser_mesh_data[6](2,k,j,i); //3*PI*sigma*cs*H
	user_out_var(36,k,j,i) = pmy_mesh->ruser_mesh_data[6](10,k,j,i); //3*PI*sigma*cs*H, with omega=omega_kep
	user_out_var(24,k,j,i) = pmy_mesh->ruser_mesh_data[6](3,k,j,i); //surface density
	user_out_var(25,k,j,i) = pmy_mesh->ruser_mesh_data[6](4,k,j,i); //AM

	//grid related values
	user_out_var(26,k,j,i) = pcoord->dx2f(j); //dphi
	user_out_var(27,k,j,i) = pcoord->GetCellVolume(k,j,i); //vol

	//time integration of radial profile
	user_out_var(28,k,j,i) = pmy_mesh->ruser_mesh_data[6](5,k,j,i); //TR = <sigma*vr*dvphi*dphi>
	user_out_var(29,k,j,i) = pmy_mesh->ruser_mesh_data[6](6,k,j,i); // <RcrossFext*dphi>
	user_out_var(30,k,j,i) = pmy_mesh->ruser_mesh_data[6](7,k,j,i); // <RcrossFext*dphi>

	//mdot
	user_out_var(31,k,j,i) = pcoord->GetFace1Area(k,j,i+1)*x1flux(IDN,k,j,i+1); //check m dot at outer boundary

	//pressure
	user_out_var(32,k,j,i) = pmy_mesh->ruser_mesh_data[6](8,k,j,i);
	//mach number
	user_out_var(33,k,j,i) = pmy_mesh->ruser_mesh_data[6](9,k,j,i);

      }
    }
  }
}


void MeshBlock::ProblemGenerator(ParameterInput *pin) {

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

  Real press_init;

  // NEW 7/13
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // compute initial conditions in cylindrical coordinates

        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
        phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
        phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad,phi,z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);


        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//!\f transform to cylindrical coordinate
// NEW 7/13

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(i);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes density in cylindrical coordinates
// NEW 7/13

Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*std::pow(rad/r0,dslope);
  Real dentem = denmid*std::exp(GM1/p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));
  den = dentem;
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! \f  computes pressure/density in cylindrical coordinates
// NEW 7/13

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, pslope);
  return poverr;
}

//----------------------------------------------------------------------------------------
//! \f  computes rotational velocity in cylindrical coordinates
// NEW 7/13

void VelProfileCyl(const Real rad, const Real phi, const Real z,
                   Real &v1, Real &v2, Real &v3) {
  Real p_over_r = PoverR(rad, phi, z);
  Real vel = (dslope+pslope)*p_over_r/(GM1/rad) + (1.0+pslope)
             - pslope*rad/std::sqrt(rad*rad+z*z);
  vel = std::sqrt(GM1/rad)*std::sqrt(vel);
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    v1=0.0;
    v2=vel;
    v3=0.0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    v1=0.0;
    v2=0.0;
    v3=vel;
  }
  return;
}



//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief Function called once every time step for user-defined work.
//======================================================================================

// function that does particle step and accretion onto companion
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

  MeshBlock *pmb = pblock->pmy_mesh->pblock;
  Hydro *phydro = pmb->phydro;
  Coordinates *pcoord = pmb->pcoord;
  int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;


  AthenaArray<Real> &x1flux = phydro->flux[X1DIR];
  AthenaArray<Real> &x2flux = phydro->flux[X2DIR];

  for(int k=ks; k<=ke; k++){
    for(int i=is; i<=ie; i++){

      //cell center variables
      Real rad_c = pcoord->x1v(i);
      Real vkep_c = sqrt(GM1/rad_c);

      //face i+1/2
      Real rad_p = pcoord->x1f(i+1);
      Real vkep_p = sqrt(GM1/rad_p);

      //face i-1/2
      Real rad_m = pcoord->x1f(i);
      Real vkep_m = sqrt(GM1/rad_m);

      for(int j=js; j<=je; j++){

	for (int n=0; n<(NHYDRO);n++){//loop over 0-4 for IDN, IVX, IVY, IVZ, IPR
	  if (std::isnan(phydro->u(n,k,j,i))){//if anything is NAN
	      printf("block: %d, n: %d ,k: %d,j: %d,i: %d\n", pmb->gid, n, k,j,i);//print mb index, variable index (IDN...), cell index
	      printf("x1v: %g, x2v:%g, x3v:%g\n",pmb->pcoord->x1v(i), pmb->pcoord->x2v(j),pmb->pcoord->x3v(k));//coordinate
	      abort();
	    }	  
	  }//end NHYDRO

	//cell center variables
	Real rho_c = phydro->u(IDN,k,j,i);
	Real vphi_c = phydro->w(IVY,k,j,i);
	Real vr_c = phydro->w(IVX,k,j,i);

	//changing to non-corotating frames
	Real Omega_orb = sqrt((GM1+GM2)/sma)/sma;
	Real v_rotate_c = Omega_orb*rad_c;
	Real v_rotate_p = Omega_orb*rad_p;
	Real v_rotate_m = Omega_orb*rad_m;

	//int AMModt
	Real AMMdot = -rad_c*x1flux(IDN,k,j,i)*(vkep_p*pcoord->GetFace1Area(k,j,i+1)- vkep_m*pcoord->GetFace1Area(k,j,i));
	ruser_mesh_data[5](0,k,j,i) += AMMdot*dt;

	//int AMTH
	Real AMTH = rad_p*pcoord->GetFace1Area(k,j,i+1)*(x1flux(IM2,k,j,i+1)-x1flux(IDN,k,j,i+1)*vkep_p) - rad_m*pcoord->GetFace1Area(k,j,i)*(x1flux(IM2,k,j,i)-x1flux(IDN,k,j,i)*vkep_m);
	ruser_mesh_data[5](1,k,j,i) += -AMTH*dt;
	//AMTH, no vkep
	Real AMTH_net = rad_p*pcoord->GetFace1Area(k,j,i+1)*x1flux(IM2,k,j,i+1) - rad_m*pcoord->GetFace1Area(k,j,i)*x1flux(IM2,k,j,i);
        ruser_mesh_data[5](2,k,j,i) += -AMTH_net*dt;

	Real AMTH_noncorotate = rad_p*pcoord->GetFace1Area(k,j,i+1)*(x1flux(IM2,k,j,i+1)-x1flux(IDN,k,j,i+1)*vkep_p + x1flux(IDN,k,j,i+1)*v_rotate_p) - rad_m*pcoord->GetFace1Area(k,j,i)*(x1flux(IM2,k,j,i)-x1flux(IDN,k,j,i)*vkep_m +x1flux(IDN,k,j,i)*v_rotate_m);
	ruser_mesh_data[5](5,k,j,i) += -AMTH_noncorotate*dt;

	//Torque
	// fect includes Coriolis force
	Real Torque = rad_c*ruser_mesh_data[3](1,k,j,i)*pcoord->GetCellVolume(k,j,i);
	// fext has no Coriolis force component
	Real Torque_ =  rad_c*ruser_mesh_data[3](8,k,j,i)*pcoord->GetCellVolume(k,j,i);
        ruser_mesh_data[5](3,k,j,i) += Torque*dt;
        ruser_mesh_data[5](4,k,j,i) += Torque_*dt;

	//Time integrated radial profiles
	//Mdot = int(sigma*vr*rad)
	Real sigma_c= rho_c*pcoord->GetCellVolume(k,j,i)/pcoord->GetFace3Area(k,j,i);
	Real mdot  = sigma_c*rad_c*vr_c*pcoord->dx2f(j);
	ruser_mesh_data[6](0,k,j,i) += mdot*dt;

	//alpha_eff = Mdot/(3pi*sigma*H*cs)
	//Real Omega_orb = sqrt((GM1+GM2)/sma)/sma;
	Real vphi_local = vphi_c + Omega_orb*rad_c;
	Real vtot = sqrt(vphi_local*vphi_local+vr_c*vr_c);
	Real cs2_c = pmb->peos->GetGamma()*phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);
	Real omega_local = vphi_c/rad_c + Omega_orb;
	Real press_avg = pcoord->dx2f(j)*3*PI*sigma_c*cs2_c/omega_local;
	Real omega_kep_c = sqrt(GM1/pow(rad_c,3));
	Real press_avg_ = pcoord->dx2f(j)*3*PI*sigma_c*cs2_c/omega_kep_c;

	ruser_mesh_data[6](1,k,j,i) += dt*mdot/press_avg;
	ruser_mesh_data[6](2,k,j,i) += dt*press_avg;
	ruser_mesh_data[6](10,k,j,i) += dt*press_avg_;

	//surface density
	ruser_mesh_data[6](3,k,j,i) += dt*sigma_c*pcoord->dx2f(j);

	//am
	ruser_mesh_data[6](4,k,j,i) += dt*sigma_c*rad_c*vphi_local*pcoord->dx2f(j);

	//the Reynold stress term of mdot
	ruser_mesh_data[6](5,k,j,i) += dt*sigma_c*vr_c*(vphi_local-vkep_c)*pcoord->dx2f(j); //Tr = <sigma*deltavphi>
	//the torque term of mdot
	ruser_mesh_data[6](6,k,j,i) += dt*rad_c*ruser_mesh_data[3](1,k,j,i)*pcoord->dx2f(j); //includes Coriolis
	ruser_mesh_data[6](7,k,j,i) += dt*rad_c*ruser_mesh_data[3](8,k,j,i)*pcoord->dx2f(j); //net fext
	//the pressure term of mdot
	ruser_mesh_data[6](8,k,j,i) += dt*phydro->w(IPR,k,j,i)*pcoord->dx2f(j);
	//the mach number
	ruser_mesh_data[6](9,k,j,i) += dt*(vtot/sqrt(cs2_c))*pcoord->dx2f(j);

  // if density is low and we are above the disk, modify z-vel and density
  // get cylindrical coordinates
  Real rad(0.0), phi(0.0), z(0.0);
  GetCylCoord(pcoord,rad,phi,z,i,j,k);

  // wait for disk to settle. see if this should be updated.
  if (time < 14.0) {
    // the second if statement is the hardcoded equation for a line
    // that sits just above the disk (checked up to t=200)
    if ((rho_c <= (5.0*dfloor)) && (std::abs(z)>(0.6*rad+0.2))) {
      phydro->w(IVZ,k,j,i) = 0.0;
      phydro->u(IDN,k,j,i) = dfloor;
      if (NON_BAROTROPIC_EOS) {
          phydro->w(IPR,k,j,i) = pfloor;
      }
    }
  }

      }//end phi
    }
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
  // softening for the companion
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
    // SD: does not work for cylindrical coordinates
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
  // find acceleration on particles from gas
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

	  //cylindrical coordinates
	  Real r = pmb->pcoord->x1v(i);
	  Real ph = pmb->pcoord->x2v(j);
	  Real z_cyl = pmb->pcoord->x3v(k);

	  Real sin_ph = sin(ph);
	  Real cos_ph = cos(ph);

	  // spherical polar coordinates, get local cartesian
	  Real x = r*cos_ph;
	  Real y = r*sin_ph;
	  Real z_cart = z_cyl;

	  // current position of the secondary
	  Real x_2 = xi[0];
	  Real y_2 = xi[1];
	  Real z_2 = xi[2];

	  Real d2 = sqrt(pow(x-x_2, 2) +
			 pow(y-y_2, 2) +
			 pow(z_cart-z_2, 2) );

	  Real d1c = pow(r,3);

	   // gravitational accels in cartesian coordinates

	  ag1i[0] += Ggrav*dm/d1c * x;
	  ag1i[1] += Ggrav*dm/d1c * y;
	  ag1i[2] += Ggrav*dm/d1c * z_cart;

	  ag2i[0] += Ggrav*dm * fspline(d2,rsoft2) * (x-x_2);
	  ag2i[1] += Ggrav*dm * fspline(d2,rsoft2) * (y-y_2);
	  ag2i[2] += Ggrav*dm * fspline(d2,rsoft2) * (z_cart-z_2);

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

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: sets solution in ghost zones to initial values
//

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        GetCylCoord(pco,rad,phi,z,il-i,j,k);
        prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        // prim(IM1,k,j,il-i) = v1;
        prim(IM1,k,j,il-i) = std::min(0.0, prim(IM1,k,j,il));
        //prim(IM2,k,j,il-i) = v2;
        prim(IM2,k,j,il-i) = prim(IM2,k,j,il);
       //prim(IM3,k,j,il-i) = v3;
        prim(IM3,k,j,il-i) = prim(IM3,k,j,il);
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
      }
    }
  }
}

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        GetCylCoord(pco,rad,phi,z,iu+i,j,k);
        prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        prim(IM1,k,j,iu+i) = v1;
        prim(IM2,k,j,iu+i) = v2;
        prim(IM3,k,j,iu+i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
      }
    }
  }
}

void DiskInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        GetCylCoord(pco,rad,phi,z,i,jl-j,k);
        prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        prim(IM1,k,jl-j,i) = v1;
        prim(IM2,k,jl-j,i) = v2;
        prim(IM3,k,jl-j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
      }
    }
  }
}

void DiskOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        GetCylCoord(pco,rad,phi,z,i,ju+j,k);
        prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        prim(IM1,k,ju+j,i) = v1;
        prim(IM2,k,ju+j,i) = v2;
        prim(IM3,k,ju+j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
      }
    }
  }
}

void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        GetCylCoord(pco,rad,phi,z,i,j,kl-k);
        prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        //prim(IM1,kl-k,j,i) = v1;
        prim(IM1,kl-k,j,i) = prim(IM1,kl,j,i);
        //prim(IM2,kl-k,j,i) = v2;
        prim(IM2,kl-k,j,i) = prim(IM2,kl,j,i);
        //prim(IM3,kl-k,j,i) = v3;
        prim(IM3,kl-k,j,i) = std::min(0.0, prim(IM3,kl,j,i)); 
        if (NON_BAROTROPIC_EOS)
          prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
      }
    }
  }
}

void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        GetCylCoord(pco,rad,phi,z,i,j,ku+k);
        prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        //prim(IM1,ku+k,j,i) = v1;
        prim(IM1,ku+k,j,i) = prim(IM1,ku,j,i);
        //prim(IM2,ku+k,j,i) = v2;
        prim(IM2,ku+k,j,i) = prim(IM2,ku,j,i);
        // prim(IM3,ku+k,j,i) = v3; 
        prim(IM3,ku+k,j,i) = std::max(0.0, prim(IM3,ku,j,i)); 
        if (NON_BAROTROPIC_EOS)
          prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
      }
    }
  }
}



// void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,
// 		    int is, int ie, int js, int je, int ks, int ke, int ngh){
//   // testing temp ceiling - hardcoded for now
//   // 6/17: setting new max temp in ghost cells really low - 10^-3 ~ .001
//   Real temp_max = 0.001;
//   for (int k=ks; k<=ke; ++k) {//Phi
//     for (int j=js; j<=je; ++j) {//theta
//       for (int i=1; i<=(NGHOST); ++i) {//R
// 	prim(IDN,k,j,is-i) = prim(IDN,k,j,is);
// 	prim(IVX,k,j,is-i) = std::min(0.0, prim(IVX,k,j,is));
// 	prim(IVY,k,j,is-i) = prim(IVY,k,j,is);//pco->x1v(is-i)*(sqrt(GM1/pow(pco->x1v(is-i),3))-sqrt((GM1+GM2)/pow(sma,3)));
// 	prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
// 	if (NON_BAROTROPIC_EOS){
//     prim(IPR,k,j,is-i) = std::max(prim(IPR,k,j,is), prim(IDN,k,j,is-i)*temp_max);
// 	}
//
//       }//end R
//     }//end phi
//   }//end z
//
//
// }
//
//
// void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){
//   for (int k=ks; k<=ke; ++k) {//z
//     for (int j=js; j<=je; ++j) {//phi
//       Real phi_coord = pco->x2v(j);
//       for (int i=1; i<=(NGHOST); ++i) {//R
// 	  prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
// 	  prim(IVX,k,j,ie+i) = std::max(0.0, prim(IVX,k,j,ie));
// 	  prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
// 	  //if (corotating_frame==1){
// 	  //  // last term subtracts the angular velocity of the frame
// 	  //  prim(IVY,k,j,ie+i) = pco->x1v(ie+i)*(sqrt(GM1/pow(pco->x1v(ie+i),3))-sqrt((GM1+GM2)/pow(sma,3)));
// 	  //} else {
// 	  //  // SD: confirm this is correct below
// 	  //  prim(IVY,k,j,ie+i) = pco->x1v(ie+i)*(sqrt(GM1/pow(pco->x1v(ie+i),3)));
// 	  //}
// 	  prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
// 	  if (NON_BAROTROPIC_EOS){
// 	    prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie);
// 	  }
//       }
//     }
//   }
// }
//
// void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){
//   for (int k=1; k<=(NGHOST); ++k) {//z
//     for (int j=js; j<=je; ++j) {//phi
//       for (int i=is; i<=ie; ++i) {//R
// 	prim(IDN,ks-k,j,i) = prim(IDN,ks,j,i);
// 	prim(IVX,ks-k,j,i) = prim(IVX,ks,j,i);
// 	prim(IVY,ks-k,j,i) = prim(IVY,ks,j,i);
// 	prim(IVZ,ks-k,j,i) = std::min(0.0, prim(IVZ,ks,j,i));
// 	if (NON_BAROTROPIC_EOS){
// 	  prim(IPR,ks-k,j,i) = prim(IPR,ks,j,i);
// 	}
//
//       }//end R
//     }//end phi
//   }//end z
//
// }
//
// void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){
//   for (int k=1; k<=(NGHOST); ++k) {//z
//     for (int j=js; j<=je; ++j) {//phi
//       for (int i=is; i<=ie; ++i) {//R
// 	prim(IDN,ke+k,j,i) = prim(IDN,ke,j,i);
// 	prim(IVX,ke+k,j,i) = prim(IVX,ke,j,i);
// 	prim(IVY,ke+k,j,i) = prim(IVY,ke,j,i);
// 	prim(IVZ,ke+k,j,i) = std::max(0.0, prim(IVZ,ke,j,i));
// 	if (NON_BAROTROPIC_EOS){
// 	  prim(IPR,ke+k,j,i) = prim(IPR,ke,j,i);
// 	}
//
//       }//end R
//     }//end phi
//   }//end z
//
// }




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
