
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <string>

using namespace std;

void initialize(double ***u, double ***u_new, double ***v, double ***v_new, double ***F, double ***G, double ***P, double ***P_new, double ***Phi, double ***Phi_new, int nx, int ny, int nz){

  for(int i=0; i<=nx-1; i++){
    for(int k=0; k<=nz-1; k++){
      for(int j=0; j<=ny-1; j++){
	P[i][j][k] =	0;
	P_new[i][j][k] = 0;
      }
    }
  }
  for(int i=0; i<=nx-1; i++){
    for(int k=0; k<=nz-2; k++){
      for(int j=0; j<=ny-2; j++){
	G[i][j][k] = 0;
	v[i][j][k] = 0;
	v_new[i][j][k] = 0;
      }
    }
  }
  for(int i=0; i<=nx-2; i++){
    for(int k=0; k<=nz-1; k++){
      for(int j=0; j<=ny-1; j++){
	F[i][j][k] = 0;
	u[i][j][k] = 0;
	u_new[i][j][k] = 0;  
      }
    }
  }
  for(int i=0; i<=nx-2; i++){
    for(int k=0; k<=nz-2; k++){
      for(int j=0; j<=ny-2; j++){
	Phi[i][j][k] = 0;
	Phi_new[i][j][k] = 0;
      }
    }
  }

  for(int k=1; k<=nz-2; k++){
    for(int j=1; j<=ny-2; j++){
      u[0][j][k] = 1;
      //P[1][j][1] = 1;  
    }
  }
  // for(int j=1; j<=((ny-3)/2); j++){
  //   Phi[0][j] = 1;  
  // }
  
}


//step 1 : find F&G
void simulation_FG(double **u, double **u_new, double **v, double **v_new, double **F, double **G, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double dt){
  double RHS_F;
  double RHS_F1;
  double RHS_F2;
  double RHS_F3;

  double RHS_G;
  double RHS_G1;
  double RHS_G2;
  double RHS_G3;

  
  for(int i=1; i<=nx-3; i++){
    for(int j =1; j<=ny-2; j++){
        RHS_F1 = (((u[i+1][j]+2*u[i][j]+u[i-1][j])/pow(dx,2.))+((u[i][j+1]+2*u[i][j]+u[i][j-1])/pow(dy,2.)))/Re;
      RHS_F2 = (pow(((u[i][j]+u[i+1][j])/2),2.)-pow(((u[i-1][j]+u[i][j])/2),2.))/dx;
      RHS_F3 = (((v[i][j]+v[i+1][j])/2)*((u[i][j]+u[i][j+1])/2)-((v[i][j-1]+v[i+1][j-1])/2)*((u[i][j-1]+u[i][j])/2))/dy;
      RHS_F = RHS_F1  - RHS_F3 -RHS_F2;
      F[i][j] = u[i][j] + dt*RHS_F;
    }
  }
  for(int j =1; j<=ny-2; j++){
   F[0][j]=1;
  }


  for(int j=1; j<=ny-2; j++){
  F[nx-2][j]=F[nx-3][j];
  }
  //////////////////////
  for(int i=1; i<=nx-2; i++){
    for(int j =1; j<=ny-3; j++){
      RHS_G1 = (((v[i+1][j]+2*v[i][j]+v[i-1][j])/pow(dx,2.))+((v[i][j+1]+2*v[i][j]+v[i][j-1])/pow(dy,2.)))/Re;
      RHS_G3 = (pow(((v[i][j]+v[i][j+1])/2),2.)-pow(((v[i][j-1]+v[i][j])/2),2.))/dy;
      RHS_G2 = (((u[i][j]+u[i][j+1])/2)*((v[i][j]+v[i+1][j])/2)-((u[i-1][j]+u[i-1][j+1])/2)*((v[i-1][j]+v[i][j])/2))/dx;
      RHS_G = RHS_G1 - RHS_G2 - RHS_G3;
      G[i][j] = v[i][j] + dt*RHS_G;
      
    }
  }

  for(int j = 1; j<=ny-3; j++){
    G[nx-1][j]=G[nx-2][j];
  }


}



void visualize(double ***var, int nx, int ny, int nz){

  for(int i=0; i<=nx-1; i++){
    for(int k=0; k<=nz-1; k++){
      for(int j=0; j<=ny-1; j++){
	cout << var[i][j][k] << "\t";
      }
      cout << "\n";
    }
    cout << "\n";
  }
  cout << "\n";
  
}

void paraview(string fileName, double ***var, int nx, int nz, int ny, double dx, double dy, double dz){
  ofstream myfile;
  myfile.open(fileName);
  //------------------------------------------------------------//
    // Paraview header
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

    // Grid
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << nz << " " << ny << "\n";
  myfile << "POINTS " << nx*nz*ny << " float\n";
  for(int j = 0; j <= ny-1; j++){
    for(int k=0; k<=nz-1; k++){
      for(int i = 0; i <= nx-1; i++){
	myfile << dx*i << " " << dy*j << " " << dz*k <<"\n";
      }
    }
  }
  
  // Data
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*ny*nz << "\n";

  myfile << "\n";
  myfile << "SCALARS PHI float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int j = 0; j <= ny-1; j++){
    for(int k = 0; k <= nz-1; k++){
      for(int i = 0; i <= nx-1; i++){
	myfile << var[i][j][k] << "\n";
      }
    }
  }
  myfile.close();

}


int main(){

  int nx = 5;
  int ny = 5;
  int nz = 5;
  double Re = 170.;
  double dx = 1;
  double dy = 1;
  double dz =1;
  double dt = 0.005;
  string fileName;

  double ***u;
  u = (double ***) malloc ((nx-1) * sizeof(double));
  for(int i=0; i<nx-1; i++){
    u[i] = (double **) malloc((ny) * sizeof(double));
    for(int j=0; j<ny; j++){
      u[i][j] = (double *) malloc((nz) * sizeof(double));
    }
  }
  double ***v;
  v = (double ***) malloc ((nx) * sizeof(double));
  for(int i=0; i<nx; i++){
    v[i] = (double **) malloc((ny-1) * sizeof(double));
    for(int j=0; j<ny-1; j++){
      v[i][j] = (double *) malloc((nz-1) * sizeof(double));
    }
  }

  double ***u_new;
  u_new = (double ***) malloc ((nx-1) * sizeof(double));
  for(int i=0; i<nx-1; i++){
    u_new[i] = (double **) malloc((ny) * sizeof(double));
    for(int j=0; j<ny; j++){
      u_new[i][j] = (double *) malloc((nz) * sizeof(double));
    }
  }
  double ***v_new;
  v_new = (double ***) malloc ((nx) * sizeof(double));
  for(int i=0; i<nx; i++){
    v_new[i] = (double **) malloc((ny-1) * sizeof(double));
    for(int j=0; j<ny-1; j++){
      v_new[i][j] = (double *) malloc((nz-1) * sizeof(double));
    }
  }
  
  double ***F;
  F = (double ***) malloc ((nx-1) * sizeof(double));
  for(int i=0; i<nx-1; i++){
    F[i] = (double **) malloc((ny) * sizeof(double));
    for(int j=0; j<ny; j++){
      F[i][j] = (double *) malloc((nz) * sizeof(double));
    }
  }
  double ***G;
  G = (double ***) malloc ((nx) * sizeof(double));
  for(int i=0; i<nx; i++){
    G[i] = (double **) malloc((ny-1) * sizeof(double));
    for(int j=0; j<ny-1; j++){
      G[i][j] = (double *) malloc((nz-1) * sizeof(double));
    }
  }
  
  double ***P;
  P = (double ***) malloc ((nx) * sizeof(double));
  for(int i=0; i<nx; i++){
    P[i] = (double **) malloc((ny) * sizeof(double));
    for(int j=0; j<ny; j++){
      P[i][j] = (double *) malloc((nz) * sizeof(double));
    }
  }
  double ***P_old;
  P_old = (double ***) malloc ((nx) * sizeof(double));
  for(int i=0; i<nx; i++){
    P_old[i] = (double **) malloc((ny) * sizeof(double));
    for(int j=0; j<ny; j++){
      P_old[i][j] = (double *) malloc((nz) * sizeof(double));
    }
  }
  double ***P_new;
  P_new = (double ***) malloc ((nx) * sizeof(double));
  for(int i=0; i<nx; i++){
    P_new[i] = (double **) malloc((ny) * sizeof(double));
    for(int j=0; j<ny; j++){
      P_new[i][j] = (double *) malloc((nz) * sizeof(double));
    }
  }
  double ***Phi;
  Phi = (double ***) malloc ((nx-1) * sizeof(double));
  for(int i=0; i<nx-1; i++){
    Phi[i] = (double **) malloc((ny-1) * sizeof(double));
    for(int j=0; j<ny-1; j++){
      Phi[i][j] = (double *) malloc((nz-1) * sizeof(double));
    }
  }

  double ***Phi_new;
  Phi_new = (double ***) malloc ((nx-1) * sizeof(double));
  for(int i=0; i<nx-1; i++){
    Phi_new[i] = (double **) malloc((ny-1) * sizeof(double));
    for(int j=0; j<ny-1; j++){
      Phi_new[i][j] = (double *) malloc((nz-1) * sizeof(double));
    }
  }

  
  initialize(u,u_new,v,v_new,F,G,P,P_new,Phi,Phi_new,nx,ny,nz);
  //visualize(u,nx-1,ny,nz);
  visualize(v,nx,ny-1,nz-1);
  //visualize(P,nx,ny,nz);
  // paraview("fileName1.vtk",u,nx-1,ny,nz,dx,dy,dz);

}
