#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

double du2dx(int i, int j, int k, double dx, double*** U) {
  return ((U[i+1][j][k] + U[i][j][k])*(U[i+1][j][k] + U[i][j][k]) - (U[i][j][k] + U[i-1][j][k])*(U[i-1][j][k] + U[i+1-1][j][k]))/(4*dx);
}
double duvdy(int i, int j, int k, double dy, double*** U, double*** V) {
  return ((U[i][j+1][k] + U[i][j][k])*(V[i+1][j][k] + V[i][j][k]) - (U[i][j][k] + U[i][j-1][k])*(V[i][j-1][k] + V[i+1][j-1][k]))/(4*dy);
}
double duwdz(int i, int j, int k, double dz, double*** U, double*** W) {
  return ((U[i][j][k+1] + U[i][j][k])*(W[i+1][j][k] + W[i][j][k]) - (U[i][j][k] + U[i][j][k-1])*(W[i][j][k-1] + W[i+1][j][k-1]))/(4*dz);
}
double dvudx(int i, int j, int k, double dx, double*** V, double*** U) {
  return ((V[i+1][j][k] + V[i][j][k])*(U[i][j+1][k] + U[i][j][k]) - (V[i][j][k] + V[i-1][j][k])*(U[i-1][j][k] + U[i-1][j+1][k]))/(4*dx);
}
double dv2dy(int i, int j, int k, double dy, double*** V) {
  return ((V[i][j+1][k] + V[i][j][k])*(V[i][j+1][k] + V[i][j][k]) - (V[i][j][k] + V[i][j-1][k])*(V[i][j-1][k] + V[i][j+1-1][k]))/(4*dy);
}
double dvwdz(int i, int j, int k, double dz, double*** V, double*** W) {
  return ((V[i][j][k+1] + V[i][j][k])*(W[i][j+1][k] + W[i][j][k]) - (V[i][j][k] + V[i][j][k-1])*(W[i][j][k-1] + W[i][j+1][k-1]))/(4*dz);
}
double dwudx(int i, int j, int k, double dx, double*** W, double*** U) {
  return ((W[i+1][j][k] + W[i][j][k])*(U[i][j][k+1] + U[i][j][k]) - (W[i][j][k] + W[i-1][j][k])*(U[i-1][j][k] + U[i-1][j][k+1]))/(4*dx);
}
double dwvdy(int i, int j, int k, double dy, double*** W, double*** V) {
  return ((W[i][j+1][k] + W[i][j][k])*(V[i][j][k+1] + V[i][j][k]) - (W[i][j][k] + W[i][j-1][k])*(V[i][j-1][k] + V[i][j-1][k+1]))/(4*dy);
}
double dw2dz(int i, int j, int k, double dz, double*** W) {
  return ((W[i][j][k+1] + W[i][j][k])*(W[i][j][k+1] + W[i][j][k]) - (W[i][j][k] + W[i][j][k-1])*(W[i][j][k-1] + W[i][j][k+1-1]))/(4*dz);
}
double d2udx2(int i, int j, int k, double dx, double*** U) {
  return (U[i+1][j][k]+U[i-1][j][k]-2*U[i][j][k])/(dx*dx);
}
double d2udy2(int i, int j, int k, double dy, double*** U) {
  return (U[i][j+1][k]+U[i][j-1][k]-2*U[i][j][k])/(dy*dy);
}
double d2udz2(int i, int j, int k, double dz, double*** U) {
  return (U[i][j][k+1]+U[i][j][k-1]-2*U[i][j][k])/(dz*dz);
}
double d2vdx2(int i, int j, int k, double dx, double*** V) {
  return (V[i+1][j][k]+V[i-1][j][k]-2*V[i][j][k])/(dx*dx);
}
double d2vdy2(int i, int j, int k, double dy, double*** V) {
  return (V[i][j+1][k]+V[i][j-1][k]-2*V[i][j][k])/(dy*dy);
}
double d2vdz2(int i, int j, int k, double dz, double*** V) {
  return (V[i][j][k+1]+V[i][j][k-1]-2*V[i][j][k])/(dz*dz);
}
double d2wdx2(int i, int j, int k, double dx, double*** W) {
  return (W[i+1][j][k]+W[i-1][j][k]-2*W[i][j][k])/(dx*dx);
}
double d2wdy2(int i, int j, int k, double dy, double*** W) {
  return (W[i][j+1][k]+W[i][j-1][k]-2*W[i][j][k])/(dy*dy);
}
double d2wdz2(int i, int j, int k, double dz, double*** W) {
  return (W[i][j][k+1]+W[i][j][k-1]-2*W[i][j][k])/(dz*dz);
}
//Initialization
void INIT(double ***U, double ***V, double ***W, double ***P, double ***PHI, double ***NEW_PHI, double ***F, double ***G, double ***H, double nx, double ny, double nz, double dx, double dy, double dz, int R, int Cx, int Cy){
  //Initialize U V
  int n = 0;
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      for(int k = 0; k < nz; k++){
	if(i == 0){U[i][j][k] = 0.5;}
	else{U[i][j][k] = 0;}
	V[i][j][k] = 0.0;
	W[i][j][k] = 0.0;
	PHI[i][j][k] = 0.0;
	NEW_PHI[i][j][k] = 0.0;
	P[i][j][k] = 0.0;
	F[i][j][k] = 0.0;
	G[i][j][k] = 0.0;
	H[i][j][k] = 0.0;
	n+= 1;
	//if(n < 100){std::cout << U[i][j][k] << V[i][j][k] << W[i][j][k] << PHI[i][j][k] << NEW_PHI[i][j][k] << P[i][j][k] << F[i][j][k] << G[i][j][k] << H[i][j][k] << std::endl;}
      }
    }
  }
  n = 0;
  std::cout << "comp" << std::endl;
  //Initialize W
  for(int i = R-Cx; i <= R+Cx; i++){
    for(int j = R-Cy; j <= R+Cy; j++){
      if((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy) < R*R){
	W[i][j][0] = 1.0;
	PHI[i][j][0] = 1.0;
	n += 1;
	//if(n < 100){std::cout << U[i][j][0] << V[i][j][0] << W[i][j][0] << PHI[i][j][0] << NEW_PHI[i][j][0] << P[i][j][0] << F[i][j][0] << G[i][j][0] << H[i][j][0] << std::endl;}
      }
    }
  }
}

void COMP_F(double*** F, double*** U, double*** V, double ***W, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double Re, double Gx){
  int imax = nx-2; int jmax = ny-2; int kmax = nz-2;
  int n = 0;
  for(int i = 1; i <= imax-1; i++){
    for(int j = 1; j <= jmax; j++){
      for(int k = 1; k <= kmax; k++){
	F[i][j][k] = U[i][j][k] + dt*((1/Re)*(d2udx2(i,j,k,dx,U)+d2udy2(i,j,k,dy,U)+d2udz2(i,j,k,dz,U))-(du2dx(i,j,k,dx,U)+duvdy(i,j,k,dy,U,V)+duwdz(i,j,k,dz,U,W))+Gx);
      }
    }
  }
  for(int j = 1; j <= jmax; j++){
    for(int k = 1; k <= kmax; k++){
      F[0][j][k] = U[0][j][k];
      F[imax][j][k] = U[imax][j][k];
    }
  }
}

void COMP_G(double*** G, double*** U, double*** V, double ***W, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double Re, double Gy){
  int imax = nx-2; int jmax = ny-2; int kmax = nz-2;
  int n = 0;
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax-1; j++){
      for(int k = 1; k <= kmax; k++){
	G[i][j][k] = V[i][j][k] + dt*((1/Re)*(d2vdx2(i,j,k,dx,V)+d2vdy2(i,j,k,dy,V)+d2vdz2(i,j,k,dz,V))-(dvudx(i,j,k,dx,V,U)+dv2dy(i,j,k,dy,V)+dvwdz(i,j,k,dz,V,W))+Gy);
      }
    }
  }
  for(int k = 1; k <= kmax; k++){
    for(int i = 1; i <= imax; i++){
      G[i][0][k] = V[i][0][k];
      G[i][jmax][k] = V[i][jmax][k];
    }
  }
}
void COMP_H(double*** H, double*** U, double*** V, double ***W, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double Re, double Gz){
  int imax = nx-2; int jmax = ny-2; int kmax = nz-2;
  int n = 0;
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax; j++){
      for(int k = 1; k <= kmax-1; k++){
	H[i][j][k] = W[i][j][k] + dt*((1/Re)*(d2wdx2(i,j,k,dx,W)+d2wdy2(i,j,k,dy,W)+d2wdz2(i,j,k,dz,W))-(dwudx(i,j,k,dx,W,U)+dwvdy(i,j,k,dy,W,V)+dw2dz(i,j,k,dz,W))+Gz);
      }
    }
  }
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax; j++){
      H[i][j][0] = W[i][j][0];
      H[i][j][kmax] = W[i][j][kmax];
    }
  }
}

void COMP_RHS(double*** RHS, double*** F, double*** G, double*** H, int nx, int ny, int nz, double dx, double dy, double dz, double dt){
  int n = 0;
  for(int i = 1; i <= nx-2; i++){
    for(int j = 1; j <= ny-2; j++){
      for(int k = 1; k <= nz-2; k++){
	RHS[i][j][k] = ((F[i][j][k] - F[i-1][j][k])/dx + \
			(G[i][j][k] - G[i][j-1][k])/dy + \
			(H[i][j][k] - H[i][j][k-1])/dz)/dt;
	//std::cout << RHS[i][j][k] << "\t";
      }
    }
  }
}

void POISSON(double*** P, double*** RHS, double***F, double***G, double***H, double SOR, int nx, int ny, int nz, double dx, double dy, double dz, double dt, int step){
  int limit = 30;
  if(step < 10){limit = 100;}
  for(int iter = 1; iter <= limit; ++iter){
    //SET BOUNDARY CONDITION
    for(int i = 0; i <= nx-1; i++){
      for(int j = 0; j <= ny-1; j++){
	P[i][j][0] = P[i][j][1];
	P[i][j][nz-1] = P[i][j][nz-2];
      }
    }
    for(int j = 0; j <= ny-1; j++){
      for(int k = 0; k <= nz-1; k++){
	P[0][j][k] = P[1][j][k];
	P[nx-1][j][k] = P[nx-2][j][k];
      }
    }
    for(int k = 0; k <= nz-1; k++){
      for(int i = 0; i <= nx-1; i++){
	P[i][0][k] = P[i][1][k];
	P[i][ny-1][k] = P[i][ny-2][k];
      }
    }
    //
    for(int i = 1; i <= nx-2; ++i){
      for(int j = 1; j <= ny-2; ++j){
	for(int k = 1; k <= nz-2; k++){
	  //P[i][j][k] = ITER_P(i,j,k,P,RHS,F,G,H,dx,dy,dz,dt, SOR);
	  P[i][j][k] = (1-SOR) * P[i][j][k] + SOR/((2.0)/(dx*dx) + (2.0)/(dy*dy) + (2.0)/(dz*dz)) *
	     ((P[i+1][j][k] + P[i-1][j][k])/(dx*dx) +  \
	      (P[i][j+1][k] + P[i][j-1][k])/(dy*dy) + \
	      (P[i][j][k+1] + P[i][j][k-1])/(dz*dz) - RHS[i][j][k]);
	  //std::cout << P[i][j][k] << "\t"
	}
      }
    }
    //
    if(iter % 5 == 0){
      double summ = 0;
      double max_norm = 0;
      double res_ijk = 0;
      for(int i = 1; i <= nx-2; i++){
	for(int j = 1; j <= ny-2; j++){
	  for(int k = 1; k <= nz-2; k++){
            res_ijk = (((P[i+1][j][k]-P[i][j][k])-(P[i][j][k]-P[i-1][j][k]))/(dx*dx) + \
	      	      ((P[i][j+1][k]-P[i][j][k])-(P[i][j][k]-P[i][j-1][k]))/(dy*dy) + \
		      ((P[i][j][k+1]-P[i][j][k])-(P[i][j][k]-P[i][j][k-1]))/(dz*dz) - \
		      RHS[i][j][k]);
	    summ += res_ijk*res_ijk;
	    if(res_ijk > max_norm){max_norm = res_ijk;}
	  }
	}
      }
      double norm = sqrt(summ/((nx-2)*(ny-2)*(nz-2)));
      std::cout << iter << "\t"  << norm << "\t" << max_norm << std::endl;
      if (iter >= 10 && norm <= 0.1){ break;}
    }
  }
  std::cout << "------------------\n";
}

void ADAP_UVW(double*** U, double*** V, double*** W, double*** F, double*** G, double***H, double***P, int nx, int ny, int nz, double dx, double dy, double dz, double dt, int R, int Cx, int Cy){
  int imax = nx-2; int jmax = ny-2; int kmax = nz-2;
  //==============================
  //======== U COMPONENT =========
  //==============================
  for(int i = 1; i <= imax-1; i++){
    for(int j = 1; j <= jmax; j++){
      for(int k = 1; k <= kmax; k++){
	U[i][j][k] = F[i][j][k] - (dt/dx)*(P[i+1][j][k]-P[i][j][k]);
      }
    }
  }
  //INFLOW CONDITION
  // 1 //
  // westward //
  for(int j = 1; j <= jmax; j++){
    for(int k = 1; k <= kmax; k++){
      U[0][j][k] = 0.5;
    }
  }
  //NO SLIP CONDITION
  // 1 //
  // bottomward noslip while topward are free slip//
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax; j++){
      U[i][j][0] = -1.0*U[i][j][1];
    }
  }
  //FREE SLIP CONDITION
  // eastward //
  for(int j = 1; j <= jmax; j++){
    for(int k = 1; k <= kmax; k++){
      U[imax][j][k] = U[imax-1][j][k];
      U[imax+1][j][k] = U[imax][j][k];
    }
  }
  // north and south //
  for(int i = 1; i <= imax; i++){
    for(int k = 1; k <= kmax; k++){
      U[i][0][k] = U[i][1][k];
      U[i][jmax+1][k] = U[i][jmax][k];
    }
  }
  // topward //
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax; j++){
      U[i][j][kmax+1] = U[i][j][kmax]; 
    }
  }
  //==================================
  //=================================
  //==============================
  //======== V COMPONENT =========
  //==============================
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax-1; j++){
      for(int k = 1; k <= kmax; k++){
	V[i][j][k] = G[i][j][k] - (dt/dy)*(P[i][j+1][k]-P[i][j][k]);
      }
    }
  }
  //INFLOW CONDITION
  // 1 //
  // westward and eastward//
  for(int j = 1; j <= jmax; j++){
    for(int k = 1; k <= kmax; k++){
      //No slip
      //Westward
      V[0][j][k] = V[1][j][k];
      //Eastward
      V[imax+1][j][k] = V[imax][j][k];
    }
  }
  // bottomward noslip while topward are free slip//
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax; j++){
      //Bottomward
      V[i][j][0] = -1.0*V[i][j][1];
      //Topward
      V[i][j][kmax+1] = V[i][j][kmax];
      
    }
  }
  // free slip in both north and south //
  for(int i = 1; i <= imax; i++){
    for(int k = 1; k <= kmax; k++){
      //Southward
      V[i][0][k] = V[i][1][k];
      //NOrthward
      V[i][jmax][k] = V[i][jmax-1][k];
      V[i][jmax+1][k] = V[i][jmax][k];
    }
  }
  
  //==================================
  //=================================
  //==============================
  //======== W COMPONENT =========
  //==============================
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax; j++){
      for(int k = 1; k <= kmax-1; k++){
	W[i][j][k] = H[i][j][k] - (dt/dz)*(P[i][j][k+1]-P[i][j][k]);
      }
    }
  }
  //INFLOW CONDITION
  // 1 //
  // westward and eastward//
  for(int j = 1; j <= jmax; j++){
    for(int k = 1; k <= kmax; k++){
      //No slip
      //Westward
      W[0][j][k] = -1.0*W[1][j][k];
      //Eastward
      W[imax+1][j][k] = W[imax][j][k];
    }
  }
  // bottomward noslip while topward are free slip//
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax; j++){
      //Bottomward
      if((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy) < R*R){
	W[i][j][0] = 1;
      }
      else{W[i][j][0] = 0;}//Topward
      W[i][j][kmax] = W[i][j][kmax-1];
      W[i][j][kmax+1] = W[i][j][kmax];
      
    }
  }
  // free slip in both north and south //
  for(int i = 1; i <= imax; i++){
    for(int k = 1; k <= kmax; k++){
      //Southward
      W[i][0][k] = W[i][1][k];
      //NOrthward
      W[i][jmax+1][k] = W[i][jmax][k];
    }
  }
}

void ADAP_PHI(double*** PHI, double*** NEW_PHI, double*** U, double*** V, double*** W, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double Re, int R, int Cx, int Cy){
  int imax = nx-2; int jmax = ny-2; int kmax = nz-2;
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax; j++){
      for(int k = 1; k <= kmax; k++){
	NEW_PHI[i][j][k] = PHI[i][j][k] + dt*(((PHI[i+1][j][k]+PHI[i-1][j][k]-2*PHI[i][j][k])/(dx*dx) + \
					       (PHI[i][j+1][k]+PHI[i][j-1][k]-2*PHI[i][j][k])/(dy*dy) + \
					       (PHI[i][j][k+1]+PHI[i][j][k-1]-2*PHI[i][j][k])/(dz*dz))/Re - \
					      ((U[i][j][k]+U[i-1][j][k])*(PHI[i+1][j][k]-PHI[i-1][j][k])/(4*dx) + \
					       (V[i][j][k]+V[i][j-1][k])*(PHI[i][j+1][k]-PHI[i][j-1][k])/(4*dy) + \
					       (W[i][j][k]+W[i][j][k-1])*(PHI[i][j][k+1]-PHI[i][j][k-1])/(4*dz)));
      }
    }
  }
  //INFLOW CONDITION
  // 1 //
  // westward and eastward//
  for(int j = 1; j <= jmax; j++){
    for(int k = 1; k <= kmax; k++){
      //No slip
      //Westward
      NEW_PHI[0][j][k] = NEW_PHI[1][j][k];
      //Eastward
      NEW_PHI[imax+1][j][k] = NEW_PHI[imax][j][k];
    }
  }
  // bottomward noslip while topward are free slip//
  for(int i = 1; i <= imax; i++){
    for(int j = 1; j <= jmax; j++){
      //Bottomward
      if((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy) < R*R){
	NEW_PHI[i][j][0] = 1;
      }
      else{NEW_PHI[i][j][0] = 0;}
      //Topward
      NEW_PHI[i][j][kmax+1] = NEW_PHI[i][j][kmax];
      
    }
  }
  // free slip in both north and south //
  for(int i = 1; i <= imax; i++){
    for(int k = 1; k <= kmax; k++){
      //Southward
      NEW_PHI[i][0][k] = NEW_PHI[i][1][k];
      //NOrthward
      NEW_PHI[i][jmax+1][k] = NEW_PHI[i][jmax][k];
    }
  }
}

void UPDATE_PHI(double*** PHI, double ***NEW_PHI, int nx, int ny, int nz){
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      for(int k = 0; k < nz; k++){
	PHI[i][j][k] = NEW_PHI[i][j][k];
      }
    }
  }
}
void PARAVIEW(std::string filename, double*** U, double ***V, double*** W, double*** P, double*** PHI, int nx, int ny, int nz, double dx, double dy, double dz, int precision)
{
  std::ofstream myfile;
  myfile.open(filename);
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
  myfile << "POINTS " << nx*ny*nz << " float\n";
  for(int k = 0; k <= nz-1; k++){
    for(int j = 0; j <= ny-1; j++){
      for(int i = 0; i <= nx-1; i++){
	myfile << dx*i << " " << dy*j << " " << dz*k << std::endl;
      }
    }
  }
  myfile << "\n";
  myfile << "POINT_DATA" << " " << nx*ny*nz << std::endl;
  /*myfile << std::endl;
  myfile << "VECTORS " << "Vel" << " float\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int k = 0; k <= nz-1; k++){
    for(int j = 0; j <= ny-1; j++){
      for(int i = 0; i <= nx-1; i++){
	myfile << std::setprecision(precision) << U[i][j][k] << " " << V[i][j][k] << " " << W[i][j][k] << std::endl;
      }
    }
    }*/
  myfile << std::endl;
  myfile << "SCALARS " << "P" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int k = 0; k <= nz-1; k++){
    for(int j = 0; j <= ny-1; j++){
      for(int i = 0; i <= nx-1; i++){
	myfile << std::setprecision(precision) << P[i][j][k] << std::endl;
      }
    }
  }
  myfile << std::endl;
  myfile << "SCALARS " << "PHI" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int k = 0; k <= nz-1; k++){
    for(int j = 0; j <= ny-1; j++){
      for(int i = 0; i <= nx-1; i++){
	myfile << std::setprecision(precision) << PHI[i][j][k] << std::endl;
      }
    }
  }
  myfile << std::endl;
  myfile << "SCALARS " << "U" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int k = 0; k <= nz-1; k++){
    for(int j = 0; j <= ny-1; j++){
      for(int i = 0; i <= nx-1; i++){
	myfile << std::setprecision(precision) << U[i][j][k] << std::endl;
      }
    }
  }
  myfile << std::endl;
  myfile << "SCALARS " << "V" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int k = 0; k <= nz-1; k++){
    for(int j = 0; j <= ny-1; j++){
      for(int i = 0; i <= nx-1; i++){
	myfile << std::setprecision(precision) << V[i][j][k] << std::endl;
      }
    }
  }
  myfile << std::endl;
  myfile << "SCALARS " << "W" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int k = 0; k <= nz-1; k++){
    for(int j = 0; j <= ny-1; j++){
      for(int i = 0; i <= nx-1; i++){
	myfile << std::setprecision(precision) << W[i][j][k] << std::endl;
      }
    }
  }
  myfile.close();
}
void SAVE_RESTARTFILE(double*** TENSOR, std::string name, int nx, int ny, int nz, int iteration, int precision){
  std::ofstream myfileO;
  myfileO.open(name+"_"+std::to_string(iteration)+".dat");
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      for(int k = 0; k <= nz-1; k++){
	myfileO << std::setprecision(precision) << TENSOR[i][j][k] << std::endl;
      }
    }
  }
  myfileO.close();
}

void READ_RESTARTFILE(double*** TENSOR, std::string name, int nx, int ny, int nz, int iteration){
  std::ifstream myfileI;
  myfileI.open(name+"_"+std::to_string(iteration)+".dat");
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      for(int k = 0; k <= nz-1; k++){
	myfileI >> TENSOR[i][j][k];
      }
    }
  }
  myfileI.close();
}

void SAVE_T(double t, std::string name){
  std::ofstream myfileO;
  myfileO.open(name+".dat");
  myfileO << t << std::endl;
  myfileO.close();
}
double READ_I(std::string name){
  std::ifstream myfileI;
  double t;
  myfileI.open(name+".dat");
  myfileI >> t;
  myfileI.close();
  return t;
}

double COMP_dt(double*** U, double*** V, double***W, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double tau)
{
  double umax = U[0][0][0]; double vmax = V[0][0][0]; double wmax = W[0][0][0];
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      for(int k = 0; k <= nz-1; k++){
	if(U[i][j][k] > umax){umax = U[i][j][k];}
	if(V[i][j][k] > vmax){vmax = V[i][j][k];}
	if(W[i][j][k] > wmax){wmax = W[i][j][k];}
      }
    }
  }
  double ARR[4];
  ARR[0] = (Re/2)*(1/(dx*dx)+1/(dy*dy)+1/(dz*dz));
  ARR[1] = dx/umax;
  ARR[2] = dy/vmax;
  ARR[3] = dz/wmax;
  double tmin = ARR[0];
  for(int i = 0; i < 4; i++){
    if(ARR[i] < tmin){tmin = ARR[i];}
  }
  return tau*tmin;
}
int main(){
  const int nx = 200;
  const int ny = 100;
  const int nz = 150;
  const double dx = 20.0/(double)nx;
  const double dy = 10.0/(double)ny;
  const double dz = 15.0/(double)nz;
  const int R = 5;
  const int Cx = 45;
  const int Cy = 45;
  const double Re = 1000.0;
  const double SOR = 1.7;
  bool checkpoint = false;
  const int timestep = 13000;
  const int precision = 4;
  const int save_precision = 8;
  double dt = 0.0001;
  double t = 0;
  double tau = 0.05;
  double Gx = 0.0;
  double Gy = 0.0;
  double Gz = 0.0;
  //DEclare Pointer
  double*** U = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    U[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      U[i][j] = (double *) malloc (nz * sizeof(double));
    }
  }
  double*** V = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    V[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      V[i][j]= (double *) malloc (nz * sizeof(double));
    }
  }
  
  double*** W = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    W[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      W[i][j]= (double *) malloc (nz * sizeof(double));
    }
  }
  double*** F = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    F[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      F[i][j]= (double *) malloc (nz * sizeof(double));
    }
  }
  double*** G = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    G[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      G[i][j]= (double *) malloc (nz * sizeof(double));
    }
  }
  double*** H = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    H[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      H[i][j]= (double *) malloc (nz * sizeof(double));
    }
  }
  double*** RHS = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    RHS[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      RHS[i][j]= (double *) malloc (nz * sizeof(double));
    }
  }
  double*** P = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    P[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      P[i][j]= (double *) malloc (nz * sizeof(double));
    }
  }
  double*** PHI = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    PHI[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      PHI[i][j]= (double *) malloc (nz * sizeof(double));
    }
  }
  double*** NEW_PHI = (double ***) malloc (nx * sizeof(double**));
  for(int i = 0; i < nx; i++){
    NEW_PHI[i] = (double **) malloc (ny * sizeof(double*));
    for(int j = 0; j < ny; j++){
      NEW_PHI[i][j]= (double *) malloc (nz * sizeof(double));
    }
  }
  std::cout << "init value all" << std::endl;
  std::string name = "JICF_";
  INIT(U, V, W, P, PHI, NEW_PHI, F, G, H, nx, ny, nz, dx, dy, dz, R, Cx, Cy);
  PARAVIEW(name+std::to_string(1)+".vtk", U, V, W, P, PHI, nx, ny, nz, dx, dy, dz, precision);
  int readAt = 0;
  std::cout << "start solving ......\n";
  for(int it = readAt; it <= 13000; ++it){
    dt = COMP_dt(U, V, W, nx, ny, nz, Re, dx, dy, dz, tau);
    std::cout << it << "\t" << dt << "\t" << t << "\t";
    std::printf("(tau = %1.2f, Re = %4.1f) \n", tau, Re);
    //std::cout << dt << std::endl;
    COMP_F(F, U, V, W, nx, ny,nz, dx, dy, dz, dt, Re, Gx);
    //std::cout << "comp F" << std::endl;
    COMP_G(G, U, V, W, nx, ny,nz, dx, dy, dz, dt, Re, Gy);
    //std::cout << "comp G" << std::endl;
    COMP_H(H, U, V, W, nx, ny,nz, dx, dy, dz, dt, Re, Gz);
    //std::cout << "comp H" << std::endl;
    COMP_RHS(RHS, F, G, H, nx, ny, nz, dx, dy, dz, dt);
    //std::cout << "comp RHS" << std::endl;
    POISSON(P, RHS, F, G, H, SOR, nx, ny, nz,  dx, dy,dz, dt, it);
    //std::cout << "comp POISSION" << std::endl;
    ADAP_UVW(U, V, W, F, G, H,P,nx,ny,nz, dx, dy, dz,dt, R,Cx,Cy);
    //std::cout << "comp ADAP_UVW" << std::endl;
    ADAP_PHI(PHI, NEW_PHI, U, V, W, nx, ny, nz, dx, dy, dz, dt, Re,R,Cx,Cy);
    //std::cout << "comp ADAP_PHI" << std::endl;
    UPDATE_PHI(PHI,NEW_PHI,nx, ny,nz);
    //std::cout << "comp UPDATE_PHI" << std::endl;
    t += dt;
    if(it % 200 == 0){
      PARAVIEW(name+std::to_string(it)+".vtk", U,V, W, P, PHI,nx,ny,nz,dx, dy,dz,precision);
    }
    if(it % 500 == 0){
      SAVE_RESTARTFILE(U, "U", nx, ny, nz, it, save_precision);
      SAVE_RESTARTFILE(V, "V", nx, ny, nz, it, save_precision);
      SAVE_RESTARTFILE(W, "W", nx, ny, nz, it, save_precision);
      SAVE_RESTARTFILE(P, "P", nx, ny, nz, it, save_precision);
      SAVE_RESTARTFILE(PHI, "PHI", nx, ny, nz, it, save_precision);
      SAVE_T(t, "savetime");
    }
  }
  return 0;
}