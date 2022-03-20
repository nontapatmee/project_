//solve momentum equation with no pressure in x-axis or F (page33-3.29)
//explicit euler in term t && central finite different interm dimention
//use gauss seidel iteration for solve Pressure

//findu,v 3.31 P33 &&


#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <string>

using namespace std;

void initialize(double **u, double **u_new, double **v, double **v_new, double **F, double **G, double **P, double **P_new, int nx, int ny){

  for(int i=0; i<=nx-1; i++){
    for(int j=0; j<=ny-1; j++){
      u[i][j] = 0;
      u_new[i][j] = 0;
      v[i][j] =	0;
      v_new[i][j] = 0;
      P[i][j] =	0;
      P_new[i][j] =	0;
      
    }
  }
  
  for(int i=0; i<=nx-1; i++){
    u[i][0] = 1;
    u[i][1] = 0.5;
    
  }
  
}

void visualize(double **var, int nx, int ny){

  for(int i=0; i<=nx-1; i++){
    for(int j=0; j<=ny-1; j++){
      cout << var[i][j] << "\t";
    }
    cout << "\n";
  }
  cout << "\n";
  
}


//step 1 : find F&G
void simulation_FG(double **u, double **u_new, double **v, double **v_new, double **F, double **G, int nx, int ny, double Re, double dx, double dy, double dt){
  double RHS_F;
  double RHS_F1;
  double RHS_F2;
  double RHS_F3;

  double RHS_G;
  double RHS_G1;
  double RHS_G2;
  double RHS_G3;

  
  for(int i=1; i<=nx-2; i++){
    for(int j =1; j<=nx-2; j++){
      RHS_F1 = (((u[i+1][j]+2*u[i][j]+u[i-1][j])/pow(dx,2.))+((u[i+1][j]+2*u[i][j]+u[i-1][j])/pow(dy,2.)))/Re;
      RHS_F2 = (pow(((u[i][j]+u[i+1][j])/2),2.)-pow(((u[i-1][j]+u[i][j])/2),2.))/dx;
      RHS_F3 = (((v[i][j]+v[i+1][j])/2)*((u[i][j]+u[i][j+1])/2)-((v[i][j-1]+v[i+1][j-1])/2)*((u[i][j-1]+u[i][j])/2))/dy;
      RHS_F = RHS_F1 - RHS_F2 - RHS_F3;
      cout << RHS_F1 << "\n";
      F[i][j] = u[i][j] + dt*RHS_F;
    }
  }

  for(int i=1; i<=nx-2; i++){
    F[i][nx-1]=F[i][nx-2];
  }
  //////////////////////
  for(int i=1; i<=nx-2; i++){
    for(int j =1; j<=nx-2; j++){
      RHS_G1 = (((v[i+1][j]+2*v[i][j]+v[i-1][j])/pow(dx,2.))+((v[i+1][j]+2*v[i][j]+v[i-1][j])/pow(dy,2.)))/Re;
      RHS_G3 = (pow(((v[i][j]+v[i][j+1])/2),2.)-pow(((v[i][j-1]+v[i][j])/2),2.))/dy;
      RHS_G2 = (((u[i][j]+u[i][j+1])/2)*((v[i][j]+v[i+1][j])/2)-((u[i-1][j]+u[i-1][j+1])/2)*((v[i-1][j]+v[i][j])/2))/dx;
      RHS_G = RHS_G1 - RHS_G2 - RHS_G3;
      G[i][j] = v[i][j] + dt*RHS_G;
    }
  }
  //boundary Neumann
  for(int i=1; i<=nx-2; i++){
    G[i][nx-1]=G[i][nx-2];
  }}


//step 2 find P[i,j] use explicit euler and central finite different !! gauss seidel iteration
void simulation_p(double **u, double **u_new, double **v, double **v_new, double **F, double **G, double **P, double **P_new, int nx, int ny, double Re, double dx, double dy, double dt){
  double aa, bb, cc;
  bool status = true;
  aa = 1/pow(dx,2.);
  bb = 1/pow(dy,2.);
  cc = (pow(dx,2.)*pow(dy,2.))/(2*(pow(dx,2.)+pow(dy,2.)));

  while(true){
    //for(int c=1; c<100; c++){
    if(status == false){break;}

    for(int i=1; i<=nx-2; i++){
      for(int j =1; j<=nx-2; j++){
	P[i][j]=P_new[i][j];
      }
    }


	
    for(int i=1; i<=nx-2; i++){
      for(int j =1; j<=nx-2; j++){
	P_new[i][j]=((((F[i][j]-F[i-1][j])/dx+(G[i][j]-G[i][j-1])/dy)*(-1/dt))+(aa*P[i+1][j]+aa*P_new[i-1][j]+bb*P[i][j+1]+bb*P_new[i][j-1]))*cc;
	status = false;
      }
    }

    for(int i=1; i<=nx-2; i++){
      for(int j =1; j<=nx-2; j++){
	if(abs(P[i][j] - P_new[i][j]) > 0.0005){
	  status = true;
	}
      }
    }
 
  }  
}


//step3 find uv last step ran-tata
void simulation_uv(double **u, double **u_new, double **v, double **v_new, double **F, double **G, double **P, double **P_new, int nx, int ny, double Re, double dx, double dy, double dt){

  for(int i=1; i<=nx-2; i++){
      for(int j =1; j<=nx-2; j++){
	u_new[i][j]=F[i][j]+(dt*P_new[i][j])/dx;
	v_new[i][j]=G[i][j]+(dt*P_new[i][j])/dy;
      }
    }

  for(int i=0; i<=nx-1; i++){
    u_new[i][0] = 1; 
  }
  
}

void update(double **var, double **var_new, int nx, int ny){
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){    
      var[i][j] = var_new[i][j];
    }
  }
}

int main(){

  int nx = 11;
  int ny = 11;
  double Re = 300.;
  double dx = 2;
  double dy = 1;
  double dt = 0.01;
  
  
  double **u;
  u = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    u[i] = (double *) malloc(ny * sizeof(double));
  }
  double **v;
  v = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    v[i] = (double *) malloc(ny * sizeof(double));
  }
  double **u_new;
  u_new = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    u_new[i] = (double *) malloc(ny * sizeof(double));
  }
  double **v_new;
  v_new = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    v_new[i] = (double *) malloc(ny * sizeof(double));
  }

  double **F;
  F = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    F[i] = (double *) malloc(ny * sizeof(double));
  }
  double **G;
  G = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    G[i] = (double *) malloc(ny * sizeof(double));
  }

  double **P;
  P = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    P[i] = (double *) malloc(ny * sizeof(double));
  }
  double **P_new;
  P_new = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    P_new[i] = (double *) malloc(ny * sizeof(double));
  }

  
  initialize(u,u_new,v,v_new,F,G,P,P_new,nx,ny);
  visualize(u,nx,ny);
  
for (int n = 1; n<=3; n++){  
  simulation_FG(u, u_new, v, v_new, F, G, nx, ny, Re, dx, dy, dt);
  simulation_p(u, u_new, v, v_new, F, G, P, P_new, nx, ny, Re, dx, dy, dt);
  simulation_uv(u, u_new, v, v_new, F, G, P, P_new, nx, ny, Re, dx, dy, dt);
  update(u,u_new,nx,ny);
  update(v,v_new,nx,ny);
  
  visualize(F,nx,ny);
  //visualize(G,nx,ny); 
  visualize(u,nx,ny);
  //visualize(v,nx,ny);
 }
  //visualize(P,nx,ny);
  //visualize(P_new,nx,ny);
  
  
  
  
}
