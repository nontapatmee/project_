//Solving 3D jet in cross flow

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <ctype.h>
#include <complex>
#include <math.h>


using namespace std;

void visualize(double ***var, int nx, int ny, int nz, char value){
  cout<<value<<"\n";
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

void initialize(double ***u, double ***u_new, double ***v, double ***v_new, double ***w, double ***w_new, double ***F, double ***G, double ***H, double ***P, double ***P_new, double ***Phi, double ***Phi_new, int nx, int ny, int nz){

  
  for(int i=0; i<=nx-1; i++){
    for(int k=0; k<=nz-1; k++){
      for(int j=0; j<=ny-1; j++){
	P[i][j][k] =	0;
	P_new[i][j][k] = 0;
      }
    }
  }
  for(int i=0; i<=nx-1; i++){
    for(int k=0; k<=nz-1; k++){
      for(int j=0; j<=ny-2; j++){
	G[i][j][k] = 0;
	v[i][j][k] = 0;
	v_new[i][j][k] = 0;
      }
    }
  }
  for(int i=0; i<=nx-1; i++){
    for(int k=0; k<=nz-2; k++){
      for(int j=0; j<=ny-1; j++){
	H[i][j][k] = 0;
	w[i][j][k] = 0;
	w_new[i][j][k] = 0;
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
  //set cross flow velocity = 0.5
  for(int k=1; k<=nz-2; k++){
    for(int j=1; j<=ny-2; j++){
      u[0][j][k] = 0.5;
      u_new[0][j][k] = 0.5;
      F[0][j][k] = 0.5;
    }
  }
  
  //Boundary condition u v w F G and H
  
  //side and top side (u v w F G and H) - neumann boundary condition 
  for(int i=1; i<=nx-3; i++){
    for(int j=1; j<=ny-2;j++){
      F[i][j][nz-1]=F[i][j][nz-2];
      u[i][j][nz-1]=u[i][j][nz-2];
      u_new[i][j][nz-1]=u_new[i][j][nz-2];
    }
    for(int k=1; k<=nz-2; k++){
      F[i][0][k]=F[i][1][k];
      u[i][0][k]=u[i][1][k];
      u_new[i][0][k]=u_new[i][1][k];

      F[i][ny-1][k]=F[i][ny-2][k];
      u[i][ny-1][k]=u[i][ny-2][k];
      u_new[i][ny-1][k]=u_new[i][ny-2][k];
    }
  }

  for(int i=1; i<=nx-2; i++){
    for(int j=1; j<=ny-3;j++){
      G[i][j][nz-1]=G[i][j][nz-2];
      v[i][j][nz-1]=v[i][j][nz-2];
      v_new[i][j][nz-1]=v_new[i][j][nz-2];
    }
    for(int k=1; k<=nz-2; k++){
      G[i][0][k]=G[i][1][k];
      v[i][0][k]=v[i][1][k];
      v_new[i][0][k]=v_new[i][1][k];

      G[i][ny-2][k]=G[i][ny-3][k];
      v[i][ny-2][k]=v[i][ny-3][k];
      v_new[i][ny-2][k]=v_new[i][ny-3][k];
    }
  }
  
  for(int i=1; i<=nx-2; i++){
    for(int j=1; j<=ny-2;j++){
      H[i][j][nz-2]=H[i][j][nz-3];
      w[i][j][nz-2]=w[i][j][nz-3];
      w_new[i][j][nz-2]=w_new[i][j][nz-3];
    }
    for(int k=1; k<=nz-3; k++){
      H[i][0][k]=H[i][1][k];
      w[i][0][k]=w[i][1][k];
      w_new[i][0][k]=w_new[i][1][k];

      H[i][ny-1][k]=H[i][ny-2][k];
      w[i][ny-1][k]=w[i][ny-2][k];
      w_new[i][ny-1][k]=w_new[i][ny-2][k];
    }
  }
  
  //outlet flow (u v w F G and H) - neumann boundary condition
  for(int j=1; j<=ny-2; j++){
    for(int k=1; k<=nz-2; k++){
      F[nx-2][j][k]=F[nx-3][j][k];
      u[nx-2][j][k]=u[nx-3][j][k];
      u_new[nx-2][j][k]=u_new[nx-3][j][k];
    }
  }
  for(int j=1; j<=ny-3; j++){
    for(int k=1; k<=nz-2; k++){
      G[nx-1][j][k]=G[nx-2][j][k];
      v[nx-1][j][k]=v[nx-2][j][k];
      v_new[nx-1][j][k]=v_new[nx-2][j][k];
    }
  }
  for(int j=1; j<=ny-2; j++){
    for(int k=1; k<=nz-3; k++){
      H[nx-1][j][k]=H[nx-2][j][k];
      w[nx-1][j][k]=w[nx-2][j][k];
      w_new[nx-1][j][k]=w_new[nx-2][j][k];
    }
  }
  
  //inlet flow (v w G and H) - dirichlet coundary condition
  for(int j=0; j<=ny-2; j++){
    for(int k=0; k<=nz-1; k++){
      G[0][j][k]=(-1)*G[1][j][k];
      v[0][j][k]=(-1)*v[1][j][k];
      v_new[0][j][k]=(-1)*v_new[1][j][k];
    }
  }
  for(int j=0; j<=ny-1; j++){
    for(int k=0; k<=nz-2; k++){
      H[0][j][k]=(-1)*H[1][j][k];
      w[0][j][k]=(-1)*w[1][j][k];
      w_new[0][j][k]=(-1)*w_new[1][j][k];
    }
  }
  
  //bottom side (u v w F G and H) - no slip or dirichlet boundary condition
  for(int i=0; i<=nx-1; i++){
    for(int j=0; j<=ny-2; j++){
	G[i][j][0]=(-1)*G[i][j][1];
	v[i][j][0]=(-1)*v[i][j][1];
	v_new[i][j][0]=(-1)*v_new[i][j][1];
    }
    for(int j=0; j<=ny-1; j++){
	H[i][j][0]=0;
	w[i][j][0]=0;
	w_new[i][j][0]=0;
    }
  }
  for(int i=0; i<=nx-2; i++){
    for(int j=0; j<=ny-1; j++){
	F[i][j][0]=(-1)*u[i][j][1];
	u[i][j][0]=(-1)*u[i][j][1];
	u_new[i][j][0]=(-1)*u_new[i][j][1];
    }
  }
  
  //set out flow jet at bottom of the model
  //velocity jet = 1
  //passsive scaler = 1
  double radius = (ny-1)/20.;
  double i_c = (ny-1)/2;
  double j_c = (ny-1)/2;
  for(int i = 0; i <= ny-1; i++){
    for(int j = 0; j <= ny-1; j++){
      if( sqrt( pow(i+1-i_c,2) + pow(j+1-j_c,2)  ) < radius ){
	
	Phi[i+1][j+1][0] = 1.;
	w[i+1][j+1][0] = 1.;
	H[i+1][j+1][0] = 1.;
	w_new[i+1][j+1][0] = 1.;

      }
    }
  }
  
}

void update(double ***var, double ***var_new, int nx, int ny, int nz){
  for(int i = 0; i <= nx-1; i++){
    for(int k = 0; k <= nz-1; k++){
      for(int j = 0; j <= ny-1; j++){
	var[i][j] = var_new[i][j];
      }
    }
  }
}


//step 1 : find compute F, G and H of half step momentum equation (runge-kutta method) 
void simulation_FG(double ***u, double ***u_new, double ***v, double ***v_new, double ***w, double ***w_new, double ***F, double ***G, double ***H, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double dt){
  double RHS_F;
  double RHS_G;
  double RHS_H;
  double du2_dx2;
  double du2_dy2;
  double du2_dz2;
  double dv2_dx2;
  double dv2_dy2;
  double dv2_dz2;
  double dw2_dx2;
  double dw2_dy2;
  double dw2_dz2;
  double du2_dx;
  double dv2_dy;
  double dw2_dz;
  double duv_dy;
  double duw_dz;
  double duv_dx;
  double dvw_dz;
  double duw_dx;
  double dwv_dy;

  //compute F
  for(int i=1; i<=nx-3; i++){
    for(int k=1; k<=nz-2; k++){
      for(int j =1; j<=ny-2; j++){

	du2_dx2=(u[i+1][j][k]-2*u[i][j][k]+u[i-1][j][k])/pow(dx,2);
	du2_dy2=(u[i][j+1][k]-2*u[i][j][k]+u[i][j-1][k])/pow(dy,2);
	du2_dz2=(u[i][j][k+1]-2*u[i][j][k]+u[i][j][k-1])/pow(dz,2);

	du2_dx=(pow((u[i][j][k]+u[i+1][j][k])/2,2)-pow((u[i-1][j][k]+u[i][j][k])/2,2))/dx ;

	duv_dy=(((v[i][j][k]+v[i+1][j][k])*(u[i][j][k]+u[i][j+1][k])/4)-((v[i][j-1][k]+v[i+1][j-1][k])*(u[i][j-1][k]+u[i][j][k])/4))/dy;

	duw_dz=(((w[i][j][k]+w[i+1][j][k])*(u[i][j][k]+u[i][j][k+1])/4)-((w[i][j][k-1]+w[i+1][j][k-1])*(u[i][j][k-1]+u[i][j][k])/4))/dz;

	RHS_F=(du2_dx2+du2_dy2+du2_dz2)/Re-du2_dx-duv_dy-duw_dz;
	F[i][j][k] = u[i][j][k]+(dt*RHS_F);
	
      }
    }
  }

  //compute G
  for(int i=1; i<=nx-2; i++){
    for(int k=1; k<=nz-2; k++){
      for(int j =1; j<=ny-3; j++){

	dv2_dx2=(v[i+1][j][k]-2*v[i][j][k]+v[i-1][j][k])/pow(dx,2);
	dv2_dy2=(v[i][j+1][k]-2*v[i][j][k]+v[i][j-1][k])/pow(dy,2);
	dv2_dz2=(v[i][j][k+1]-2*v[i][j][k]+v[i][j][k-1])/pow(dz,2);

	dv2_dy=(pow((v[i][j][k]+v[i][j+1][k])/2,2)-pow((v[i][j-1][k]+v[i][j][k])/2,2))/dy ;

	duv_dx=(((u[i][j][k]+u[i][j+1][k])*(v[i][j][k]+v[i+1][j][k])/4)-((u[i-1][j][k]+u[i-1][j+1][k])*(v[i-1][j][k]+v[i][j][k])/4))/dx;

	dvw_dz=(((w[i][j][k]+w[i][j+1][k])*(v[i][j][k]+v[i][j][k+1])/4)-((w[i][j][k-1]+w[i][j+1][k-1])*(v[i][j][k-1]+v[i][j][k])/4))/dz;

	RHS_G=(dv2_dx2+dv2_dy2+dv2_dz2)/Re-dv2_dy-duv_dx-dvw_dz;
	G[i][j][k] = v[i][j][k]+(dt*RHS_G);
	
      }
    }
  }

  //compute H
  for(int i=1; i<=nx-2; i++){
    for(int k=1; k<=nz-3; k++){
      for(int j =1; j<=ny-2; j++){

	dw2_dx2=(w[i+1][j][k]-2*w[i][j][k]+w[i-1][j][k])/pow(dx,2);
	dw2_dy2=(w[i][j+1][k]-2*w[i][j][k]+w[i][j-1][k])/pow(dy,2);
	dw2_dz2=(w[i][j][k+1]-2*w[i][j][k]+w[i][j][k-1])/pow(dz,2);

	dw2_dz=(pow((w[i][j][k]+w[i][j][k+1])/2,2)-pow((w[i][j][k-1]+w[i][j][k])/2,2))/dz ;

	duw_dx=(((u[i][j][k]+u[i][j][k+1])*(w[i][j][k]+w[i+1][j][k])/4)-((u[i-1][j][k]+u[i-1][j][k+1])*(w[i-1][j][k]+w[i][j][k])/4))/dx;
	  
	dwv_dy=(((v[i][j][k]+v[i][j][k+1])*(w[i][j][k]+w[i][j+1][k])/4)-((v[i][j-1][k]+v[i][j-1][k+1])*(w[i][j-1][k]+w[i][j][k])/4))/dy;
	  
	RHS_H=(dw2_dx2+dw2_dy2+dw2_dz2)/Re-dw2_dz-duw_dx-dwv_dy;
	H[i][j][k] = w[i][j][k]+(dt*RHS_H);
	
      }
    }
  }

  //Update G, F and H boundary conditions
  //outlet side, side and top side - neumann boundary condition
  for(int i=1; i<=nx-3; i++){
    for(int j=1; j<=ny-2;j++){
      F[i][j][nz-1]=F[i][j][nz-2];
    }
    for(int k=1; k<=nz-2; k++){
      F[i][0][k]=F[i][1][k];

      F[i][ny-1][k]=F[i][ny-2][k];
    }
  }
  for(int i=1; i<=nx-2; i++){
    for(int j=1; j<=ny-3;j++){
      G[i][j][nz-1]=G[i][j][nz-2];
    }
    for(int k=1; k<=nz-2; k++){
      G[i][0][k]=G[i][1][k];

      G[i][ny-2][k]=G[i][ny-3][k];
    }
  }
  
  for(int i=1; i<=nx-2; i++){
    for(int j=1; j<=ny-2;j++){
      H[i][j][nz-2]=H[i][j][nz-3];
    }
    for(int k=1; k<=nz-3; k++){
      H[i][0][k]=H[i][1][k];

      H[i][ny-1][k]=H[i][ny-2][k];
    }
  }

  for(int j=1; j<=ny-2; j++){
    for(int k=1; k<=nz-2; k++){
      F[nx-2][j][k]=F[nx-3][j][k];
    }
  }
  for(int j=1; j<=ny-3; j++){
    for(int k=1; k<=nz-2; k++){
      G[nx-1][j][k]=G[nx-2][j][k];
    }
  }
  for(int j=1; j<=ny-2; j++){
    for(int k=1; k<=nz-3; k++){
      H[nx-1][j][k]=H[nx-2][j][k];
    }
  }
  //inlet - inflow or dirichlet boundary condition
 for(int j=0; j<=ny-2; j++){
    for(int k=0; k<=nz-1; k++){
      G[0][j][k]=(-1)*G[1][j][k];
    }
  }
  for(int j=0; j<=ny-1; j++){
    for(int k=0; k<=nz-2; k++){
      H[0][j][k]=(-1)*H[1][j][k];
    }
  }
  //bottom - no slip or dirichlet boundary condition
  for(int i=0; i<=nx-1; i++){
    for(int j=0; j<=ny-2; j++){
      G[i][j][0]=(-1)*G[i][j][1];
    }      
  }
  for(int i=0; i<=nx-2; i++){
    for(int j=0; j<=ny-1; j++){
      F[i][j][0]=(-1)*u[i][j][1];
    }
  }
  
}

//step 2 find P : implicit in pressure => SOR iteration
void simulation_p(double ***u, double ***u_new, double ***v, double ***v_new, double ***w, double ***w_new, double ***F, double ***G, double ***H, double ***P, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double dt, int step){
  
  double RHS;
  int limit = 30;
  double max_norm = 0;
  double norm ;
  double coef_P;
  double dx_2, dy_2, dz_2;
  double rit=1.;
  double eps=0.1;
  double it_max;
  double SOR = 1.7;
  dx_2 = pow(dx,2.);
  dy_2 = pow(dy,2.);
  dz_2 = pow(dz,2.);
  coef_P = ((2)/dx_2 + (2)/dy_2 + (2)/dz_2) ;
  
  if (step < 10)
    {
      it_max = 100;
    }
  else{it_max = 50;}
      
  int it=0;
  while(true){
    rit = 0;
    // Set Pressure boundary condition - every side => neumann ,except outlet => dirichlet 
    for (int j = 0; j <= ny - 1; j++){
      for (int k = 0; k <= nz - 1; k++){
	P[0][j][k] = P[1][j][k];
	P[nx - 1][j][k] = 0;
      } 
      for (int i = 0; i <= nx - 1; i++){
	P[i][j][0] = P[i][j][1];
	P[i][j][nz - 1] = P[i][j][nz - 2];
      }
    }
    for (int k = 0; k <= nz - 1; k++){
      for (int i = 0; i <= nx - 1; i++){
	P[i][0][k] = P[i][1][k];
	P[i][ny - 1][k] = P[i][ny - 2][k];
      }
    }

    //compute Pressure
    for (int i = 1; i <= nx - 2; ++i){
      for (int j = 1; j <= ny - 2; ++j){
	for (int k = 1; k <= nz - 2; k++){
	  RHS = ((F[i][j][k]-F[i-1][j][k])/dx + (G[i][j][k]-G[i][j-1][k])/dy + (H[i][j][k]-H[i][j][k-1])/dz)/dt;
	  
	  P[i][j][k] = (1 - SOR) * P[i][j][k] + (SOR / coef_P) *((P[i + 1][j][k] + P[i - 1][j][k]) / dx_2 +(P[i][j + 1][k] + P[i][j - 1][k]) / dy_2 +(P[i][j][k + 1] + P[i][j][k - 1]) / dz_2 - RHS);
          
	  rit = rit +pow( (((P[i + 1][j][k] - P[i][j][k]) - (P[i][j][k] - P[i - 1][j][k])) / dx_2 + ((P[i][j + 1][k] - P[i][j][k]) - (P[i][j][k] - P[i][j - 1][k])) / dy_2 + ((P[i][j][k + 1] - P[i][j][k]) - (P[i][j][k] - P[i][j][k - 1])) / dz_2 - RHS), 2);
	}
      }
    }
      
    it++;
    norm = sqrt(rit/((nx - 2) * (ny - 2) * (nz - 2)));

    if ((it >= 1 && norm <= eps) || it>=it_max)
      {
	break;
      }
        
  }

  cout << it << "\t" << norm << "\n";
}

//step3 find u, v and w : explicit velocities
void simulation_uv(double ***u, double ***u_new, double ***v, double ***v_new, double ***w, double ***w_new, double ***F, double ***G, double ***H, double ***P, double ***P_new, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double dt){

  for(int i=1; i<=nx-3; i++){
    for(int k=1; k<=nz-2; k++){
      for(int j =1; j<=ny-2; j++){
	u_new[i][j][k]=F[i][j][k]-(dt*(P_new[i+1][j][k]-P_new[i][j][k]))/dx;
      }
    }
  }
  for(int i=1; i<=nx-2; i++){
    for(int k=1; k<=nz-2; k++){
      for(int j =1; j<=ny-3; j++){
	v_new[i][j][k]=G[i][j][k]-(dt*(P_new[i][j+1][k]-P_new[i][j][k]))/dy;
      }
    }
  }
  for(int i=1; i<=nx-2; i++){
    for(int k=1; k<=nz-3; k++){
      for(int j =1; j<=ny-2; j++){
	w_new[i][j][k]=H[i][j][k]-(dt*(P_new[i][j][k+1]-P_new[i][j][k]))/dz;
      }
    }
  }

  //update u,v,w boundary condition
  //side and top side => neumann boundary condition
  for(int i=1; i<=nx-3; i++){
    for(int j=1; j<=ny-2;j++){
      u[i][j][nz-1]=u[i][j][nz-2];
      u_new[i][j][nz-1]=u_new[i][j][nz-2];
    }
    for(int k=1; k<=nz-2; k++){
      u[i][0][k]=u[i][1][k];
      u_new[i][0][k]=u_new[i][1][k];

      u[i][ny-1][k]=u[i][ny-2][k];
      u_new[i][ny-1][k]=u_new[i][ny-2][k];
    }
  }

  for(int i=1; i<=nx-2; i++){
    for(int j=1; j<=ny-3;j++){
      v[i][j][nz-1]=v[i][j][nz-2];
      v_new[i][j][nz-1]=v_new[i][j][nz-2];
    }
    for(int k=1; k<=nz-2; k++){
      v[i][0][k]=v[i][1][k];
      v_new[i][0][k]=v_new[i][1][k];

      v[i][ny-2][k]=v[i][ny-3][k];
      v_new[i][ny-2][k]=v_new[i][ny-3][k];
    }
  }
  
  for(int i=1; i<=nx-2; i++){
    for(int j=1; j<=ny-2;j++){
      w[i][j][nz-2]=w[i][j][nz-3];
      w_new[i][j][nz-2]=w_new[i][j][nz-3];
    }
    for(int k=1; k<=nz-3; k++){
      w[i][0][k]=w[i][1][k];
      w_new[i][0][k]=w_new[i][1][k];

      w[i][ny-1][k]=w[i][ny-2][k];
      w_new[i][ny-1][k]=w_new[i][ny-2][k];
    }
  }
  
  //outlet => neumann boundary condition
  for(int j=1; j<=ny-2; j++){
    for(int k=1; k<=nz-2; k++){
      u[nx-2][j][k]=u[nx-3][j][k];
      u_new[nx-2][j][k]=u_new[nx-3][j][k];
    }
  }
  for(int j=1; j<=ny-3; j++){
    for(int k=1; k<=nz-2; k++){
      v[nx-1][j][k]=v[nx-2][j][k];
      v_new[nx-1][j][k]=v_new[nx-2][j][k];
    }
  }
  for(int j=1; j<=ny-2; j++){
    for(int k=1; k<=nz-3; k++){
      w[nx-1][j][k]=w[nx-2][j][k];
      w_new[nx-1][j][k]=w_new[nx-2][j][k];
    }
  }
  //inlet => inflow, dirichlet boundary condition 
  for(int j=0; j<=ny-2; j++){
    for(int k=0; k<=nz-1; k++){
      v[0][j][k]=(-1)*v[1][j][k];
      v_new[0][j][k]=(-1)*v_new[1][j][k];
    }
  }
  for(int j=0; j<=ny-1; j++){
    for(int k=0; k<=nz-2; k++){
      w[0][j][k]=(-1)*w[1][j][k];
      w_new[0][j][k]=(-1)*w_new[1][j][k];
    }
  }
  
  //bottom => no slip, dirichlet boundary condition
  for(int i=0; i<=nx-1; i++){
    for(int j=0; j<=ny-2; j++){
	v[i][j][0]=(-1)*v[i][j][1];
	v_new[i][j][0]=(-1)*v_new[i][j][1];
    }
  }
  for(int i=0; i<=nx-2; i++){
    for(int j=0; j<=ny-1; j++){
	u[i][j][0]=(-1)*u[i][j][1];
	u_new[i][j][0]=(-1)*u_new[i][j][1];
    }
  }
 
}


//Compute passive scalar : explicit
void simulation_passiveScalar(double ***u, double ***u_new, double ***v, double ***v_new, double ***w, double ***w_new, double ***F, double ***G, double ***H, double ***P, double ***P_new, double ***Phi,double ***Phi_new, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double dt){

  for(int i = 1; i <= nx-3; i++){
    for(int k = 1; k <= nz-3; k++){
      for(int j = 1; j <= ny-3; j++){    
	
	Phi_new[i][j][k] = Phi[i][j][k] + dt*(((Phi[i+1][j][k]+Phi[i-1][j][k]-2*Phi[i][j][k])/(pow(dx,2)) + (Phi[i][j+1][k]+Phi[i][j-1][k]-2*Phi[i][j][k])/(pow(dy,2)) + (Phi[i][j][k+1]+Phi[i][j][k-1]-2*Phi[i][j][k])/(pow(dz,2)))/Re - ((u[i][j][k]+u[i-1][j][k])*(Phi[i+1][j][k]-Phi[i-1][j][k])/(4*dx) + (v[i][j][k]+v[i][j-1][k])*(Phi[i][j+1][k]-Phi[i][j-1][k])/(4*dy) + (w[i][j][k]+w[i][j][k-1])*(Phi[i][j][k+1]-Phi[i][j][k-1])/(4*dz)));

      }
    }
  }
  
  //update Phi boundary
  //side & top => neumann
  for(int i=1; i<=nx-3; i++){
    for(int j=1; j<=ny-3;j++){
      Phi_new[i][j][nz-2]=Phi_new[i][j][nz-3];
    }
    for(int k=1; k<=nz-3; k++){
      Phi_new[i][0][k]=Phi_new[i][1][k];

      Phi_new[i][ny-2][k]=Phi_new[i][ny-3][k];
    }
  }

  //outlet => neumann
  for(int j=1; j<=ny-3; j++){
    for(int k=1; k<=nz-3; k++){
      Phi_new[nx-2][j][k]=Phi_new[nx-3][j][k];
    }
  }
 
  //inlet => dirichlet
  for(int j=0; j<=ny-2; j++){
    for(int k=0; k<=nz-2; k++){
      Phi_new[0][j][k]=0;
    }
  } 
}

//Compute kinetic energy from u v and w
void simulation_kineticE(double ***Ke, double ***u, double ***u_new, double ***v, double ***v_new, double ***w, double ***w_new, double ***F, double ***G, double ***H, double ***P, double ***P_new, double ***Phi,double ***Phi_new, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double dt){

  for(int i = 0; i <= nx-4; i++){
    for(int k = 0; k <= nz-4; k++){
      for(int j = 0; j <= ny-4; j++){
        Ke[i][j][k] = sqrt(pow(u[i+1][j+1][k+1],2)+pow(v[i+1][j+1][k+1],2)+pow(w[i+1][j+1][k+1],2));
      }
    }
  }
}


void paraview(string fileName, double ***var, int nx, int ny, int nz, double dx, double dy, double dz){
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
    for(int k = 0; k <= nz-1; k++){
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

  int nx = 300;
  int ny = 150;
  int nz = 225;
  double Re = 1000.;
  double dx = 20./real(nx-2);
  double dy = 10./real(ny-2);
  double dz = 15/real(nz-2);
  double dt = 0.008;
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
      v[i][j] = (double *) malloc((nz) * sizeof(double));
    }
  }
  double ***w;
  w = (double ***) malloc ((nx) * sizeof(double));
  for(int i=0; i<nx; i++){
    w[i] = (double **) malloc((ny) * sizeof(double));
    for(int j=0; j<ny; j++){
      w[i][j] = (double *) malloc((nz-1) * sizeof(double));
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
      v_new[i][j] = (double *) malloc((nz) * sizeof(double));
    }
  }
  double ***w_new;
  w_new = (double ***) malloc ((nx) * sizeof(double));
  for(int i=0; i<nx; i++){
    w_new[i] = (double **) malloc((ny) * sizeof(double));
    for(int j=0; j<ny; j++){
      w_new[i][j] = (double *) malloc((nz-1) * sizeof(double));
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
      G[i][j] = (double *) malloc((nz) * sizeof(double));
    }
  }
  double ***H;
  H = (double ***) malloc ((nx) * sizeof(double));
  for(int i=0; i<nx; i++){
    H[i] = (double **) malloc((ny) * sizeof(double));
    for(int j=0; j<ny; j++){
      H[i][j] = (double *) malloc((nz-1) * sizeof(double));
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

  double ***Ke;
  Ke = (double ***) malloc ((nx-3) * sizeof(double));
  for(int i=0; i<nx-3; i++){
    Ke[i] = (double **) malloc((ny-3) * sizeof(double));
    for(int j=0; j<ny-3; j++){
      Ke[i][j] = (double *) malloc((nz-3) * sizeof(double));
    }
  }

  
  initialize(u,u_new,v,v_new,w,w_new,F,G,H,P,P_new,Phi,Phi_new,nx,ny,nz);
  update(Phi_new,Phi,nx-1,ny-1,nz-1);
  
  for(int n=1;n<=80000;n++){
    simulation_FG(u, u_new, v, v_new, w, w_new, F, G, H, nx, ny, nz, Re, dx, dy, dz, dt);
    simulation_p(u, u_new, v, v_new, w, w_new, F, G, H, P_new, nx, ny, nz, Re, dx, dy, dz, dt, n);
    update(P,P_new,nx,ny,nz);
    simulation_uv(u, u_new, v, v_new, w, w_new, F, G, H, P, P_new, nx, ny, nz, Re, dx, dy, dz, dt);
    update(u,u_new,nx-1,ny,nz);
    update(v,v_new,nx,ny-1,nz);
    update(w,w_new,nx,ny,nz-1);
    simulation_passiveScalar(u, u_new, v, v_new, w, w_new, F, G, H, P, P_new, Phi, Phi_new, nx, ny, nz, Re, dx, dy, dz, dt);
    update(Phi,Phi_new,nx-1,ny-1,nz-1);
    simulation_kineticE(Ke,u, u_new, v, v_new, w, w_new, F, G, H, P, P_new, Phi, Phi_new, nx, ny, nz, Re, dx, dy, dz, dt);
    cout << "n = " << n << "\n";
    //visualize(u,nx-1,ny,nz,'u');
    //visualize(v,nx,ny-1,nz,'v');  
    //visualize(G,nx,ny-1,nz,'g');
    //visualize(w,nx,ny,nz-1,'w');
    //visualize(Phi,nx-1,ny-1,nz-1,'Q');
    //visualize(F,nx-1,ny,nz,'F');
    //visualize(G,nx,ny-1,nz,'G');
    //visualize(H,nx,ny,nz-1,'H');
    //visualize(P,nx,ny,nz,'P');
    //visualize(P_new,nx,ny,nz,'j');
    if(n%250==0){ 
      fileName = "u_3D" + to_string(n) + ".vtk";
      paraview(fileName,u,nx-1,ny,nz,dx,dy,dz);
      fileName = "Phi_3D" + to_string(n) + ".vtk";
      paraview(fileName,Phi,nx-1,ny-1,nz-1,dx,dy,dz);
      fileName = "w_3D" + to_string(n) + ".vtk";
      paraview(fileName,w,nx,ny,nz-1,dx,dy,dz);
      fileName = "v_3D" + to_string(n) + ".vtk";
      paraview(fileName,v,nx,ny-1,nz,dx,dy,dz);
      fileName = "P_3D" + to_string(n) + ".vtk";
      paraview(fileName,P,nx,ny,nz,dx,dy,dz);
      fileName = "Ke_3D" + to_string(n) + ".vtk";
      paraview(fileName,Ke,nx-3,ny-3,nz-3,dx,dy,dz);
    }
  }

}
