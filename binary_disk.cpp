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


Real massfluxix1(MeshBlock *pmb, int iout);


// disk parameters
namespace{
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
void VelProfileCyl(const Real rad, const Real phi, const Real z,Real &v1, Real &v2, Real &v3);
// problem parameters which are useful to make global to this file
Real r0, rho0, dslope, p0_over_r0, pslope, gamma_gas;
Real dfloor;
}

// global (to this file) problem parameters
Real da,pa; // ambient density, pressure


// companion parameters and parameters for corotating frame
Real GM2, GM1; // point masses
Real rsoft2; // softening length of PM 2
int  include_gas_backreaction, corotating_frame; // flags for output, gas backreaction on EOM, frame choice
int n_particle_substeps; // substepping of particle integration
Real xi[3], vi[3], agas1i[3], agas2i[3]; // cartesian positions/vels of the secondary object, gas->particle acceleration
Real Omega[3];  // vector rotation of the frame

int particle_accrete;
Real mdot, pdot[3]; // accretion parameters


// for particle output file
Real trackfile_next_time, trackfile_dt;
int  trackfile_number;
Real Ggrav;


// restart simulations
int is_restart;




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

  // Get parameters for initial density profile
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

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
  Real sma = pin->GetOrAddReal("problem","sma",2.0);//semi-major axis
  Real ecc = pin->GetOrAddReal("problem","ecc",0.0);
  Real Omega_orb, vcirc;


  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));



  // // enroll the BCs
  // if(mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
  //   EnrollUserBoundaryFunction(OUTER_X1, DiodeOuterX1);
  // }


	// disk bc
	// enroll user-defined boundary condition
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


  // Enroll a Source Function
  EnrollUserExplicitSourceFunction(TwoPointMass);

  // Enroll AMR
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);

  //Enroll history dump
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0,massfluxix1,"massfluxix1");

  // always write at startup
  trackfile_next_time = time;
  trackfile_number = 0;


  // allocate MESH data for the particle pos/vel, Omega frame
  AllocateRealUserMeshDataField(3);
  ruser_mesh_data[0].NewAthenaArray(3);
  ruser_mesh_data[1].NewAthenaArray(3);
  ruser_mesh_data[2].NewAthenaArray(3);

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

    Real ph= pmb->pcoord->x3v(k);
    Real sin_ph = sin(ph);
    Real cos_ph = cos(ph);

    for(int j=pmb->js; j<=pmb->je; j++) {

      Real th= pmb->pcoord->x2v(j);
      Real sin_th = sin(th);
      Real cos_th = cos(th);

      for(int i=pmb->is; i<=pmb->ie; i++) {

	Real r = pmb->pcoord->x1v(i);
	Real x = r*sin_th*cos_ph;
	Real y = r*sin_th*sin_ph;
	Real z = r*cos_th;

	Real dist = std::sqrt(SQR(x-xi[0]) +
			      SQR(y-xi[1]) +
			      SQR(z-xi[2]) );

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
    is_restart=0;
  }


  // Gravitational acceleration from orbital motion
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {

	Real r = pmb->pcoord->x1v(i);
	Real th= pmb->pcoord->x2v(j);
	Real ph= pmb->pcoord->x3v(k);

	Real vr  = prim(IVX,k,j,i);
	Real vth = prim(IVY,k,j,i);
	Real vph = prim(IVZ,k,j,i);

	//get some angles
	Real sin_th = sin(th);
	Real cos_th = cos(th);
	Real sin_ph = sin(ph);
	Real cos_ph = cos(ph);

	// current position of the secondary
	Real x_2 = xi[0];
	Real y_2 = xi[1];
	Real z_2 = xi[2];
	Real d12c = pow(xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2], 1.5);

	// spherical polar coordinates, get local cartesian
	Real x = r*sin_th*cos_ph;
	Real y = r*sin_th*sin_ph;
	Real z = r*cos_th;

	Real d2  = sqrt(pow(x-x_2, 2) +
			pow(y-y_2, 2) +
			pow(z-z_2, 2) );

	//
	//  COMPUTE ACCELERATIONS
	//
	// PM1
	//Real a_r1 = -GM1/pow(r,2);
	// cell volume avg'd version, see pointmass.cpp sourceterm code.
	Real a_r1 = -GM1*pmb->pcoord->coord_src1_i_(i)/r;


	// PM2 gravitational accels in cartesian coordinates
	Real a_x = - GM2 * fspline(d2,rsoft2) * (x-x_2);
	Real a_y = - GM2 * fspline(d2,rsoft2) * (y-y_2);
	Real a_z = - GM2 * fspline(d2,rsoft2) * (z-z_2);

	//if(corotating_frame == 1){//Xiaoshan: excluding the COM-centerstarframe transformation too
	  // add the correction for the orbiting frame (relative to the COM)
	  a_x += -  GM2 / d12c * x_2;
	  a_y += -  GM2 / d12c * y_2;
	  a_z += -  GM2 / d12c * z_2;

	if(corotating_frame == 1){
	  // distance from the origin in cartesian (vector)
	  Real rxyz[3];
	  rxyz[0] = x;
	  rxyz[1] = y;
	  rxyz[2] = z;

	  // get the cartesian velocities from the spherical (vector)
	  Real vgas[3];
	  vgas[0] = sin_th*cos_ph*vr + cos_th*cos_ph*vth - sin_ph*vph;
	  vgas[1] = sin_th*sin_ph*vr + cos_th*sin_ph*vth + cos_ph*vph;
	  vgas[2] = cos_th*vr - sin_th*vth;

	  // add the centrifugal and coriolis terms

	  // centrifugal
	  Real Omega_x_r[3], Omega_x_Omega_x_r[3];
	  cross(Omega,rxyz,Omega_x_r);
	  cross(Omega,Omega_x_r,Omega_x_Omega_x_r);

	  a_x += - Omega_x_Omega_x_r[0];
	  a_y += - Omega_x_Omega_x_r[1];
	  a_z += - Omega_x_Omega_x_r[2];

	  // coriolis
	  Real Omega_x_v[3];
	  cross(Omega,vgas,Omega_x_v);

	  a_x += -2.0*Omega_x_v[0];
	  a_y += -2.0*Omega_x_v[1];
	  a_z += -2.0*Omega_x_v[2];
	}

	// add the gas acceleration of the frame of ref
	if(include_gas_backreaction == 1){
	  a_x += -agas1i[0];
	  a_y += -agas1i[1];
	  a_z += -agas1i[2];
	}

	// convert back to spherical
	Real a_r  = sin_th*cos_ph*a_x + sin_th*sin_ph*a_y + cos_th*a_z;
	Real a_th = cos_th*cos_ph*a_x + cos_th*sin_ph*a_y - sin_th*a_z;
	Real a_ph = -sin_ph*a_x + cos_ph*a_y;

	// add the PM1 accel
	a_r += a_r1;

	//
	// ADD SOURCE TERMS TO THE GAS MOMENTA/ENERGY
	//
	Real den = prim(IDN,k,j,i);

	Real src_1 = dt*den*a_r;
	Real src_2 = dt*den*a_th;
	Real src_3 = dt*den*a_ph;

	// add the source term to the momenta  (source = - rho * a)
        
	cons(IM1,k,j,i) += src_1;
	cons(IM2,k,j,i) += src_2;
	cons(IM3,k,j,i) += src_3;

	if (NON_BAROTROPIC_EOS) {
	  // update the energy (source = - rho v dot a
	  cons(IEN,k,j,i) += src_1/den * 0.5*(pmb->phydro->flux[X1DIR](IDN,k,j,i) + pmb->phydro->flux[X1DIR](IDN,k,j,i+1));
	  cons(IEN,k,j,i) += src_2*prim(IVY,k,j,i) + src_3*prim(IVZ,k,j,i);
	}
	

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
  return;
}

// Disk initial setup
//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

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



  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
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





// Disk functions


namespace {
//----------------------------------------------------------------------------------------
//!\f transform to cylindrical coordinate

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

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, pslope);
  return poverr;
}

//----------------------------------------------------------------------------------------
//! \f  computes rotational velocity in cylindrical coordinates

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
} // namespace

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
        prim(IM1,k,j,il-i) = v1;
        prim(IM2,k,j,il-i) = v2;
        prim(IM3,k,j,il-i) = v3;
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
        prim(IM1,kl-k,j,i) = v1;
        prim(IM2,kl-k,j,i) = v2;
        prim(IM3,kl-k,j,i) = v3;
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
        prim(IM1,ku+k,j,i) = v1;
        prim(IM2,ku+k,j,i) = v2;
        prim(IM3,ku+k,j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
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
	  massflux += face1(is)*x1flux(0,k,j,is);//x1flux(0) is the density flux, do we want to time volume too?
        }

      }
    }
  }

  face1.DeleteAthenaArray();
  x1flux.DeleteAthenaArray();

  return massflux;


}
