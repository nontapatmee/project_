
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <ctype.h>
#include <complex>
#include <math.h>
//#include <conio.h>


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

  for(int k=1; k<=nz-2; k++){
    for(int j=1; j<=ny-2; j++){
      u[0][j][k] = 0.5;
      u_new[0][j][k] = 0.5;
      F[0][j][k] = 0.5;
    }
  }
  for(int i=0; i<=nx-2; i++){
    for(int k=1; k<=nz-2;k++){
      F[i][0][k]=(-1)*F[i][1][k];
      F[i][ny-1][k]=(-1)*F[i][ny-2][k];
      u[i][0][k]=(-1)*u[i][1][k];
      u[i][ny-1][k]=(-1)*u[i][ny-2][k];
      u_new[i][0][k]=(-1)*u_new[i][1][k];
      u_new[i][ny-1][k]=(-1)*u_new[i][ny-2][k];
    }
    for(int j=1; j<=ny-2; j++){
      F[i][j][0]=(-1)*F[i][j][1];
      F[i][j][nz-1]=(-1)*F[i][j][nz-2];
      u[i][j][0]=(-1)*u[i][j][1];
      u[i][j][nz-1]=(-1)*u[i][j][nz-2];
      u_new[i][j][0]=(-1)*u_new[i][j][1];
      u_new[i][j][nz-1]=(-1)*u_new[i][j][nz-2];
    }
    u[i][0][0]=u[i][1][0];
    u[i][0][nz-1]=u[i][0][nz-2];
    u[i][ny-1][0]=u[i][ny-2][0];
    u[i][ny-1][nz-1]=u[i][ny-2][nz-1];

    u_new[i][0][0]=u_new[i][1][0];
    u_new[i][0][nz-1]=u_new[i][0][nz-2];
    u_new[i][ny-1][0]=u_new[i][ny-2][0];
    u_new[i][ny-1][nz-1]=u_new[i][ny-2][nz-1];
    
    F[i][0][0]=F[i][1][0];
    F[i][0][nz-1]=F[i][0][nz-2];
    F[i][ny-1][0]=F[i][ny-2][0];
    F[i][ny-1][nz-1]=F[i][ny-2][nz-1];
  }


  // for(int j=1; j<=((ny-3)/2); j++){
  //   Phi[0][j] = 1.;
  //   Phi_new[0][j] = 1.;
  // }

  //set phi
  double radius = (ny)/20.;
  double i_c = (ny)/2;
  double j_c = (ny)/2;
 for(int i = 0; i <= ny-1; i++){ 
   for(int j = 0; j <= ny-1; j++){
     //for(int k = 0; k<=ny-1;k++){
     if( sqrt( pow(i+1-i_c,2) + pow(j+1-j_c,2)  ) < radius ){
	// Phi[i][j][nz-2] = 1.;
	// w[i][j][nz-2] = 2.;
	// H[i][j][nz-2] = 2.;
	// w_new[i][j][nz-2] = 2.;

	Phi[i][j][0] = 1.;
	w[i][j][0] = 1.;
	H[i][j][0] = 1.;
	w_new[i][j][0] = 1.;

	////////////
	//u[0][j][k] = 1.;
	//Phi[0][j][k] = 1.;
	//F[0][j][k] = 1.;
	//u_new[0][j][k] = 1.;
      }
      //else { Phi[i][j][nz-2] = 0.;}
      
   }
  }
  //visualize(w,nx,ny,nz-1,'w');
  //visualize(Phi,nx-1,ny-1,nz-1,'Q');
  
  
  
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


//step 1 : find F&G
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
  ///
  // for(int k=1; k<=nz-3; k++){
  //   for(int j=1; j<=ny-3; j++){
  //     v[0][j][k] = 1.;
  //     w[0][j][k] = 1.;
      
  //   }}
  ///
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

  //update inlet, outlet G,F,H boundary
  //G,H inlet Dirichlet, outlet neumann
  for(int k = 1; k<=nz-2; k++){
    for(int j = 1; j<=ny-3; j++){
      G[0][j][k]=(-1)*G[1][j][k];
      G[nx-1][j][k]=G[nx-2][j][k];
    }
  }
  for(int k = 1; k<=nz-3; k++){
    for(int j = 1; j<=ny-2; j++){
      H[0][j][k]=(-1)*H[1][j][k];
      H[nx-1][j][k]=H[nx-2][j][k];
    }
  }
  //F outlet neumann, wall Dirichlet
  for(int k = 1; k<=nz-2; k++){
    for(int j =1; j<=ny-2; j++){
      F[nx-2][j][k]=F[nx-3][j][k];
    }
  }
  for(int i=0; i<=nx-2; i++){
    for(int k=1; k<=nz-2;k++){
      F[i][0][k]=(-1)*F[i][1][k];
      F[i][ny-1][k]=(-1)*F[i][ny-2][k];
    }
    for(int j=1; j<=ny-2; j++){
      F[i][j][0]=(-1)*F[i][j][1];
      F[i][j][nz-1]=(-1)*F[i][j][nz-2];
    }
    F[i][0][0]=F[i][1][0];
    F[i][0][nz-1]=F[i][0][nz-2];
    F[i][ny-1][0]=F[i][ny-2][0];
    F[i][ny-1][nz-1]=F[i][ny-2][nz-1];
  }


}


//step 2 find P[i,j] use explicit euler and central finite different !! gauss seidel iteration
void simulation_p(double ***u, double ***u_new, double ***v, double ***v_new, double ***w, double ***w_new, double ***F, double ***G, double ***H, double ***P, double ***P_new, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double dt){
  double dx_2, dy_2, dz_2;
  double ee, ew, en, es, ea, eb;
  double RHS, coef_P;
  bool status = true;
  dx_2 = pow(dx,2.);
  dy_2 = pow(dy,2.);
  dz_2 = pow(dz,2.);
  
  double value_w=1.;
  double rit=1.;
  double eps=0.01;
  int it, it_max;
  double SOR = 1.7;
  it_max=1000;
  
  //update bd condition F(0,j) <= u(0,j) ,v(i,0) => G(i,0)
  for(int k=0; k<=nz-1; k++){
    for(int j=0; j<=ny-1; j++){
      F[nx-2][j][k]=u[nx-2][j][k];
    }
  }
   it=0;

   while(true){
     //check convergence ,norm of the residual<eps and max iteration
     if( (rit<eps) || (it>it_max) ){break;}

     rit =0.;
    
    for(int i=1; i<=nx-2; i++){
      for(int k=1; k<=nz-2; k++){
	for(int j =1; j<=ny-2; j++){

	  if(i==1){ew=0;}
	  else{ew=1;}
	  if(i==nx-2){ee=0;}
	  else{ee=1;}
	  if(j==1){es=0;}
	  else{es=1;}
	  if(j==ny-2){en=0;}
	  else{en=1;}
	  if(k==1){eb=0;}
	  else{eb=1;}
	  if(k==nz-2){ea=0;}
	  else{ea=1;}
	  
	  RHS = ((F[i][j][k]-F[i-1][j][k])/dx + (G[i][j][k]-G[i][j-1][k])/dy + (H[i][j][k]-H[i][j][k-1])/dz)/dt;
	  coef_P = ((ee+ew)/dx_2 + (en+es)/dy_2 + (ea+eb)/dz_2) ;
	  P_new[i][j][k] = ( ( ((ee*P[i+1][j][k]+ew*P_new[i-1][j][k])/dx_2 + (en*P[i][j+1][k]+es*P_new[i][j-1][k])/dy_2 + (ea*P[i][j][k+1]+eb*P_new[i][j][k-1])/dz_2) - RHS )/coef_P );
	  //rit=rit+pow( ( (ee*(P[i+1][j][k]-P[i][j][k])-ew*(P[i][j][k]-P[i-1][j][k]))/dx_2 + (en*(P[i][j+1][k]-P[i][j][k])-es*(P[i][j][k]-P[i][j-1][k]))/dy_2 + (ea*(P[i][j][k+1]-P[i][j][k])-eb*(P[i][j][k]-P[i][j][k-1]))/dz_2 - RHS ) , 2);

	  //P[i][j][k] = (1 - SOR) * P[i][j][k] + SOR / ((2.0) / (dx * dx) + (2.0) / (dy * dy) + (2.0) / (dz * dz)) *((P[i + 1][j][k] + P[i - 1][j][k]) / (dx * dx) +(P[i][j + 1][k] + P[i][j - 1][k]) / (dy * dy) +(P[i][j][k + 1] + P[i][j][k - 1]) / (dz * dz) - RHS);
	  
	  if(k==int(nz-2)/2){rit=rit+pow( (ee*(P[i+1][j][k]-P[i][j][k])-ew*(P[i][j][k]-P[i-1][j][k]))/pow(dx,2) + (en*(P[i][j+1][k]-P[i][j][k])-es*(P[i][j][k]-P[i][j-1][k]))/pow(dy,2) - RHS,2);}
	  status = false;

	}
      }
    }
    //cal rit
    //rit=sqrt(rit/((nx-2)*(ny-2)*(nz-2)));
    rit=sqrt(rit/((nx-2)*(ny-2)));
    //cout<<rit<<"\n";
    //check convergence
    for(int i=1; i<=nx-2; i++){
      for(int k=1; k<=nz-2; k++){
	for(int j =1; j<=ny-2; j++){
	  if(abs(P[i][j] - P_new[i][j]) > 10e-50){
	    status = true;
	  }
	}
      }
    }

    //update P wall boundary
    for(int k = 1; k<=nz-1; k++){
      for(int j =1; j<=ny-1; j++){
	P_new[0][j][k] = P_new[1][j][k];
	P_new[nx-1][j][k]=P_new[nx-2][j][k];
      }
    }
    for(int i=0; i<=nx-1; i++){
      for(int k=1; k<=nz-2;k++){
	P_new[i][0][k]=P_new[i][1][k];
	P_new[i][ny-1][k]=P_new[i][ny-2][k];
      }
      for(int j=1; j<=ny-2; j++){
	P_new[i][j][0]=P_new[i][j][1];
	P_new[i][j][nz-1]=P_new[i][j][nz-2];
      }
      P_new[i][0][0]=P_new[i][1][0];
      P_new[i][0][nz-1]=P_new[i][0][nz-2];
      P_new[i][ny-1][0]=P_new[i][ny-2][0];
      P_new[i][ny-1][nz-1]=P_new[i][ny-2][nz-1];
    }
    
    //update P_new=>P
    update(P,P_new,nx,ny,nz);
    it++;
    //visualize(P,nx,ny,nz,'P');
    //visualize(P_new,nx,ny,nz,'b');

   }
   cout<<rit<<"\n";
 
}

void POISSON(double ***P, double ***F, double ***G, double ***H, double SOR, int nx, int ny, int nz, double dx, double dy, double dz, double dt, int step)
{   
    double RHS;
    int limit = 300;
    if (step < 10)
    {
        limit = 1000;
    }
    for (int iter = 1; iter <= limit; iter++)
    {   
        double summ = 0;
        double max_norm = 0;
        double res_ijk = 0;
        // SET BOUNDARY CONDITION
        for (int i = 0; i <= nx - 1; i++)
        {
            for (int j = 0; j <= ny - 1; j++)
            {
                P[i][j][0] = P[i][j][1];
                P[i][j][nz - 1] = P[i][j][nz - 2];
            }
        }
        for (int j = 0; j <= ny - 1; j++)
        {
            for (int k = 0; k <= nz - 1; k++)
            {
                P[0][j][k] = P[1][j][k];
                P[nx - 1][j][k] = P[nx - 2][j][k];
            }
        }
        for (int k = 0; k <= nz - 1; k++)
        {
            for (int i = 0; i <= nx - 1; i++)
            {
                P[i][0][k] = P[i][1][k];
                P[i][ny - 1][k] = P[i][ny - 2][k];
            }
        }
        //
        for (int i = 1; i <= nx - 2; ++i)
        {
            for (int j = 1; j <= ny - 2; ++j)
            {
                for (int k = 1; k <= nz - 2; k++)
                {
                    RHS = ((F[i][j][k]-F[i-1][j][k])/dx + (G[i][j][k]-G[i][j-1][k])/dy + (H[i][j][k]-H[i][j][k-1])/dz)/dt;
                    // P[i][j][k] = ITER_P(i,j,k,P,RHS,F,G,H,dx,dy,dz,dt, SOR);
                    P[i][j][k] = (1 - SOR) * P[i][j][k] + SOR / ((2.0) / (dx * dx) + (2.0) / (dy * dy) + (2.0) / (dz * dz)) *
                                                              ((P[i + 1][j][k] + P[i - 1][j][k]) / (dx * dx) +
                                                               (P[i][j + 1][k] + P[i][j - 1][k]) / (dy * dy) +
                                                               (P[i][j][k + 1] + P[i][j][k - 1]) / (dz * dz) - RHS);
                    // std::cout << P[i][j][k] << "\t"
                     res_ijk = (((P[i + 1][j][k] - P[i][j][k]) - (P[i][j][k] - P[i - 1][j][k])) / (dx * dx) +
                                   ((P[i][j + 1][k] - P[i][j][k]) - (P[i][j][k] - P[i][j - 1][k])) / (dy * dy) +
                                   ((P[i][j][k + 1] - P[i][j][k]) - (P[i][j][k] - P[i][j][k - 1])) / (dz * dz) -
                                   RHS);
                    summ += res_ijk * res_ijk;
                    if (res_ijk > max_norm)
                        {
                            max_norm = res_ijk;
                        }
                }
            }
        }
    
        double norm = sqrt(summ / ((nx - 2) * (ny - 2) * (nz - 2)));
        if (iter >= 10 && norm <= 0.1)
        {
            break;
        }
        
    }
    cout << iter << "\t" << norm << "\t" << max_norm << "\n";
}

// //step3 find uv
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

 //update u,v,w boundary
 for(int k = 1; k<=nz-2; k++){
    for(int j = 1; j<=ny-3; j++){
      v_new[0][j][k]=(-1)*v_new[1][j][k];
      v_new[nx-1][j][k]=v_new[nx-2][j][k];
    }
  }
 for(int k = 1; k<=nz-3; k++){
    for(int j = 1; j<=ny-2; j++){
      w_new[0][j][k]=(-1)*w_new[1][j][k];
      w_new[nx-1][j][k]=w_new[nx-2][j][k];
    }
  }
  //u,v,w outlet neumann, u wall Dirichlet
  for(int k = 1; k<=nz-2; k++){
    for(int j =1; j<=ny-2; j++){
      u_new[nx-2][j][k]=F[nx-3][j][k];
    }
  }
  for(int i=0; i<=nx-2; i++){
    for(int k=1; k<=nz-2;k++){
      u_new[i][0][k]=(-1)*u_new[i][1][k];
      u_new[i][ny-1][k]=(-1)*u_new[i][ny-2][k];
    }
    for(int j=1; j<=ny-2; j++){
      u_new[i][j][0]=(-1)*u_new[i][j][1];
      u_new[i][j][nz-1]=(-1)*u_new[i][j][nz-2];
    }
    u_new[i][0][0]=u_new[i][1][0];
    u_new[i][0][nz-1]=u_new[i][0][nz-2];
    u_new[i][ny-1][0]=u_new[i][ny-2][0];
    u_new[i][ny-1][nz-1]=u_new[i][ny-2][nz-1];
  }
 
}



void simulation_passiveScalar(double ***u, double ***u_new, double ***v, double ***v_new, double ***w, double ***w_new, double ***F, double ***G, double ***H, double ***P, double ***P_new, double ***Phi,double ***Phi_new, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double dt){

  for(int i = 1; i <= nx-3; i++){
    for(int k = 1; k <= nz-3; k++){
      for(int j = 1; j <= ny-3; j++){    
	Phi_new[i][j][k] = ((((Phi[i+1][j][k]-2*Phi[i][j][k]+Phi[i-1][j][k])/(pow(dx,2)*Re))-((u[i][j][k]/(2*dx))*(Phi[i+1][j][k]-Phi[i-1][j][k])))+(((Phi[i][j+1][k]-2*Phi[i][j][k]+Phi[i][j-1][k])/(pow(dy,2)*Re))-((v[i][j][k]/(2*dy))*(Phi[i][j+1][k]-Phi[i][j-1][k])))+(((Phi[i][j][k+1]-2*Phi[i][j][k]+Phi[i][j][k-1])/(pow(dz,2)*Re))-((v[i][j][k]/(2*dz))*(Phi[i][j][k+1]-Phi[i][j][k-1]))))*dt+Phi[i][j][k];
      }
    }
  }
  for(int k =1; k<=nz-3 ; k++){
    for(int j = 1; j <= (ny-3); j++){    
      Phi_new[nx-2][j][k] = Phi_new[nx-3][j][k];
    }
  }
 
  // for(int i=0; i<=nx-2; i++){
  //   Phi_new[i][0]=(-1)*Phi_new[i][1];
  //   Phi_new[i][ny-2]=(-1)*Phi_new[i][ny-3];
  // }
  for(int i=0; i<=nx-2; i++){
    for(int k=0; k<=nz-2;k++){
      Phi_new[i][0][k]=(-1)*Phi_new[i][1][k];
      Phi_new[i][ny-2][k]=(-1)*Phi_new[i][ny-3][k];
    }
    for(int j=0; j<=ny-2; j++){
      Phi_new[i][j][0]=(-1)*Phi_new[i][j][1];
      Phi_new[i][j][nz-2]=(-1)*Phi_new[i][j][nz-3];
    }
    // u_new[i][0][0]=u_new[i][1][0];
    // u_new[i][0][nz-1]=u_new[i][0][nz-2];
    // u_new[i][ny-1][0]=u_new[i][ny-2][0];
    // u_new[i][ny-1][ny-1]=u_new[i][ny-2][ny-1];
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

  int nx = 100;
  int ny = 60;
  int nz = 80;
  double Re = 300.;
  double dx = 20./real(nx-2);
  double dy = 10./real(ny-2);
  double dz = 15/real(nz-2);
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

  
  initialize(u,u_new,v,v_new,w,w_new,F,G,H,P,P_new,Phi,Phi_new,nx,ny,nz);
  update(Phi_new,Phi,nx-1,ny-1,nz-1);
  //visualize(Phi,nx-1,ny-1,nz-1,'Q');
  //visualize(w,nx,ny-1,nz-1,'w');
  
  for(int n=1;n<=3000;n++){

    simulation_FG(u, u_new, v, v_new, w, w_new, F, G, H, nx, ny, nz, Re, dx, dy, dz, dt);
    //simulation_p(u, u_new, v, v_new, w, w_new, F, G, H, P, P_new, nx, ny, nz, Re, dx, dy, dz, dt);
    POISSON(P_new, F, G, H, 1.7, nx, ny, nz, dx, dy, dz, dt, n);
    update(P,P_new,nx,ny,nz);
    simulation_uv(u, u_new, v, v_new, w, w_new, F, G, H, P, P_new, nx, ny, nz, Re, dx, dy, dz, dt);
    update(u,u_new,nx-1,ny,nz);
    update(v,v_new,nx,ny-1,nz);
    update(w,w_new,nx,ny,nz-1);
    simulation_passiveScalar(u, u_new, v, v_new, w, w_new, F, G, H, P, P_new, Phi, Phi_new, nx, ny, nz, Re, dx, dy, dz, dt);
    update(Phi,Phi_new,nx-1,ny-1,nz-1);
    cout << "n = " << n << "\n";
    //visualize(u,nx-1,ny,nz,'u');
    //visualize(v,nx,ny-1,nz,'v');
  
    //visualize(G,nx,ny-1,nz,'v');
    // visualize(w,nx,ny,nz-1,'w');
    // visualize(Phi,nx-1,ny-1,nz-1,'Q');
    //visualize(F,nx-1,ny,nz,'F');
    //visualize(G,nx,ny-1,nz,'G');
    //visualize(H,nx,ny,nz-1,'H');
    //visualize(P,nx,ny,nz,'P');
    fileName = "u_3D" + to_string(n) + ".vtk";
    paraview(fileName,u,nx-1,ny,nz,dx,dy,dz);
    fileName = "Phi_3D" + to_string(n) + ".vtk";
    paraview(fileName,Phi,nx-1,ny-1,nz-1,dx,dy,dz);
    fileName = "w_3D" + to_string(n) + ".vtk";
    paraview(fileName,w,nx,ny,nz-1,dx,dy,dz);
  }
  //visualize(u,nx-1,ny,nz,'u');
}
