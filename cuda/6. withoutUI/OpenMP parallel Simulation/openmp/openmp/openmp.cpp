//// openmp.cpp : Defines the entry point for the console application.
////
//
//#include "stdafx.h"
//
//
//int _tmain(int argc, _TCHAR* argv[])
//{
//	return 0;
//}

//ramneek: after u have gone through the tute for openmp, try to make cloth engine as a stand alone component and then parralelize it
#include "stdafx.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

void GetPosInXYZDirection(double** nodepos, double *xPos, double *yPos, double *zPos, int N);
void GetForceInXYZDirection(double** force, double *xForce, double *yForce, double *zForce, int N);
void GetVelInXYZDirection(double** velocity, double *xVel, double *yVel, double *zVel, int N);
void create_cloth(int N,double separation, double offset, double ballsize, double** velocity, double** force, double** oldforce, double** nodepos);
void c_compute_force(int N,int delta, double gravity, double separation, double fcon, double *xForce, double *yForce, double *zForce,
double *xPos, double *yPos, double *zPos);
void c_compute_force_Cudafied(int N,int delta, double gravity, double separation, double fcon, double *xForce, double *yForce, double *zForce,
double *xPos, double *yPos, double *zPos);
double mag(double x, double y, double z);

int main(int argc, char* argv[])
{
	//int N=3,np=1;
	int N=50,np=1;

	
	if (argc != 3) {
		printf(" %s Number_of_nodes Number_of_threads \n",argv[0]);
		//return -1;
	 }
	 else {
		N = atoi(argv[1]);
		np = atoi(argv[2]);
		if (N < 1){
		  printf("Error: Number_of_nodes (%i) < 1 \n",N);
		  return -1;
		}
	 }

	double myBallX = 0.0;
	double myBallY = 0.0;
	double myBallZ = 0.0;
	double myBallRadius = 3.0;

	double separation=1.0;
	double mass = 1.0;
	double fcon = 10.0;
	int interact = 2;
	double gravity = 1.0;
	double ballsize = 3.0;//this should match with myballRadius
	double offset = 0.0;
	double dt = 0.01;
	double update = 1;

	omp_set_num_threads(np);
	
	//just for declaring it i will use a more conventional access pattern. otherwise i will access using nx*N+ny
	double **velocity;
	velocity  = (double **)malloc(N * N * sizeof(double *));
	for(int i = 0; i < N * N ; i++)
	velocity[i] = (double *)malloc(3 * sizeof(double));

	//just for declaring it i will use a more conventional access pattern. otherwise i will access using nx*N+ny
	double **force;
	force  = (double **)malloc(N * N * sizeof(double *));
	for(int i = 0; i < N * N ; i++)
	force[i] = (double *)malloc(3 * sizeof(double));

	//just for declaring it i will use a more conventional access pattern. otherwise i will access using nx*N+ny
	double **oldforce;
	oldforce  = (double **)malloc(N * N * sizeof(double *));
	for(int i = 0; i < N * N ; i++)
	oldforce[i] = (double *)malloc(3 * sizeof(double));

	//just for declaring it i will use a more conventional access pattern. otherwise i will access using nx*N+ny
	double **nodepos;
	nodepos  = (double **)malloc(N * N * sizeof(double *));
	for(int i = 0; i < N * N ; i++)
	nodepos[i] = (double *)malloc(3 * sizeof(double));

	create_cloth(N,separation,offset,ballsize,velocity,force,oldforce,nodepos);

	double *xVel, *yVel,  *zVel;
	xVel = (double *)malloc(N*N*sizeof(double));
	yVel = (double *)malloc(N*N*sizeof(double));
	zVel = (double *)malloc(N*N*sizeof(double));
	GetVelInXYZDirection(velocity,xVel,yVel,zVel,N);
	//xPos, yPos, zPos = GetPosInXYZDirection(nodes)

	double *xPos, *yPos,  *zPos;
	xPos = (double *)malloc(N*N*sizeof(double));
	yPos = (double *)malloc(N*N*sizeof(double));
	zPos = (double *)malloc(N*N*sizeof(double));
	GetPosInXYZDirection(nodepos,xPos,yPos,zPos,N);
	//xPos, yPos, zPos = GetPosInXYZDirection(nodes)

	double *xForce, *yForce,  *zForce;
	xForce = (double *)malloc(N*N*sizeof(double));
	yForce = (double *)malloc(N*N*sizeof(double));
	zForce = (double *)malloc(N*N*sizeof(double));
	GetForceInXYZDirection(force,xForce,yForce,zForce,N);
	//xForce, yForce, zForce = GetForceInXYZDirection(nodes) 

	c_compute_force(N,interact,gravity,separation,fcon,xForce, yForce, zForce,xPos, yPos, zPos);
	//c_compute_force_Cudafied(N,interact,gravity,separation,fcon,xForce, yForce, zForce,xPos, yPos, zPos);

	clock_t clktime1 = clock();
	time_t RealTime1 = time(0);
	
	int iter=0;
	//while(iter<0)
	while(iter<3000)
	{
		iter=iter+1;

		for(int nx=0;nx<N;nx++)
		{
			for(int ny=0;ny<N;ny++)
			{
				nodepos[nx*N+ny][0] = nodepos[nx*N+ny][0] + dt*(velocity[nx*N+ny][0]+dt*xForce[nx*N+ny]*.5);
				nodepos[nx*N+ny][1] = nodepos[nx*N+ny][1] + dt*(velocity[nx*N+ny][1]+dt*yForce[nx*N+ny]*.5);
				nodepos[nx*N+ny][2] = nodepos[nx*N+ny][2] + dt*(velocity[nx*N+ny][2]+dt*zForce[nx*N+ny]*.5);

				oldforce[nx*N+ny][0] = xForce[nx*N+ny];
				oldforce[nx*N+ny][1] = yForce[nx*N+ny];
				oldforce[nx*N+ny][2] = zForce[nx*N+ny];

				// if(iter==300)
				// {
					// printf("%E %E %E\n",nodepos[nx*N+ny][0],nodepos[nx*N+ny][1],nodepos[nx*N+ny][2]);
					// // printf("%E %E %E\n",xForce[nx*N+ny],yForce[nx*N+ny],zForce[nx*N+ny]);
				// }
			}
		}
		
		for(int i=0;i<N*N;i++)
		{
			double distX = nodepos[i][0] - myBallX;
			double distY = nodepos[i][1] - myBallY;
			double distZ = nodepos[i][2] - myBallZ;

			double dist = mag(distX,distY,distZ);
			//dist = node.pos-vector(myball.x,myball.y,myball.z)

			if(dist<myBallRadius)
			{
				// printf("%E %E %E\n",nodepos[i][0] ,nodepos[i][1],nodepos[i][2]);
				double fvectorX = (distX/dist)*myBallRadius;
				double fvectorY = (distY/dist)*myBallRadius;
				double fvectorZ = (distZ/dist)*myBallRadius;
				//fvector=dist/dist.mag*myball.radius

				nodepos[i][0] = myBallX +fvectorX;
				nodepos[i][1] = myBallY +fvectorY;
				nodepos[i][2] = myBallZ +fvectorZ;
				//node.pos=vector(myball.x,myball.y,myball.z)+fvector	

				double fvectorMag = mag(fvectorX,fvectorY,fvectorZ);
				velocity[i][0] = velocity[i][0] - (velocity[i][0]*fvectorX/fvectorMag)*(fvectorX/fvectorMag);
				velocity[i][1] = velocity[i][1] - (velocity[i][1]*fvectorY/fvectorMag)*(fvectorY/fvectorMag);
				velocity[i][2] = velocity[i][2] - (velocity[i][2]*fvectorZ/fvectorMag)*(fvectorZ/fvectorMag);
				
				// printf("%E %E %E\n",fvectorX,fvectorY,fvectorZ);
				//node.velocity = node.velocity - (dot(node.velocity,fvector/fvector.mag))*(fvector/fvector.mag)
			}
		}
		
		GetPosInXYZDirection(nodepos,xPos,yPos,zPos,N);
		GetForceInXYZDirection(force,xForce,yForce,zForce,N);		
		c_compute_force(N,interact,gravity,separation,fcon,xForce, yForce, zForce,xPos, yPos, zPos);
		//c_compute_force_Cudafied(N,interact,gravity,separation,fcon,xForce, yForce, zForce,xPos, yPos, zPos);

		for(int nx=0;nx<N;nx++)
		{
			for(int ny=0;ny<N;ny++)
			{
				velocity[nx*N+ny][0]=velocity[nx*N+ny][0]+dt*(xForce[nx*N+ny] + oldforce[nx*N+ny][0])*0.5;
				velocity[nx*N+ny][1]=velocity[nx*N+ny][1]+dt*(yForce[nx*N+ny] + oldforce[nx*N+ny][1])*0.5;
				velocity[nx*N+ny][2]=velocity[nx*N+ny][2]+dt*(zForce[nx*N+ny] + oldforce[nx*N+ny][2])*0.5;
				//nodes[nx*N+ny].velocity+=dt*(vector(TotalForceEachDim[nx*N+ny][0],TotalForceEachDim[nx*N+ny][1],TotalForceEachDim[nx*N+ny][2])+nodes[nx*N+ny].oldforce)*0.5
				// if(iter==300)
				// {
				// printf("%E %E %E \n",velocity[nx*N+ny][0],velocity[nx*N+ny][1],velocity[nx*N+ny][2]);
				// }
			}
		}
		if(iter==2000)
		{
			clock_t clktime2 = clock();
			time_t RealTime2 = time(0);

			double diffClock = ((double)(clktime2-clktime1))/CLOCKS_PER_SEC;
			double diffSystem= difftime(RealTime2,RealTime1);					
			printf("CPU clock timetime:%f \n",diffClock);
			printf("Wall time:%f \n",diffSystem);
			printf("Force is %E %E %E\n",xForce[0],yForce[0],zForce[0]);
		}
	}
	getchar();
}

void GetPosInXYZDirection(double** nodepos, double *xPos, double *yPos, double *zPos, int N)
{
	//initialize the two dim array(matrix) to 
	for(int i = 0; i < N*N; i++)		
	xPos[i] = 0.0;
	for(int i = 0; i < N*N; i++)		
	yPos[i] = 0.0;
	for(int i = 0; i < N*N; i++)		
	zPos[i] = 0.0;

	for(int i = 0; i < N*N; i++)		
	{
		xPos[i] = nodepos[i][0];
		yPos[i] = nodepos[i][1];
		zPos[i] = nodepos[i][2];
	}
}

void GetForceInXYZDirection(double** force, double *xForce, double *yForce, double *zForce, int N)
{
	//initialize the two dim array(matrix) to 
	for(int i = 0; i < N*N; i++)		
	xForce[i] = 0.0;
	for(int i = 0; i < N*N; i++)		
	yForce[i] = 0.0;
	for(int i = 0; i < N*N; i++)		
	zForce[i] = 0.0;

	for(int i = 0; i < N*N; i++)		
	{
		xForce[i] = force[i][0];
		yForce[i] = force[i][1];
		zForce[i] = force[i][2];
	}
}

void GetVelInXYZDirection(double** velocity, double *xVel, double *yVel, double *zVel, int N)
{
	//initialize the two dim array(matrix) to 
	for(int i = 0; i < N*N; i++)		
	xVel[i] = 0.0;
	for(int i = 0; i < N*N; i++)		
	yVel[i] = 0.0;
	for(int i = 0; i < N*N; i++)		
	zVel[i] = 0.0;

	for(int i = 0; i < N*N; i++)		
	{
		xVel[i] = velocity[i][0];
		yVel[i] = velocity[i][1];
		zVel[i] = velocity[i][2];
	}
}

void create_cloth(int N,double separation, double offset, double ballsize, double** velocity, double** force, double** oldforce, double** nodepos)
{		
	//This is the conventional way to access two dim array. 
	//initialize the two dim array(matrix) to 
	// for(int i = 0; i < N*N; i++)
	// for(int j = 0; j < 3; j++)
	// velocity[i][j] = 0.0;

	//this is same but bit absurd way of doind same.
	// for(int nx=0;i<N;i++)
	// {
	// for(int ny=0;i<N;i++)
	// {
	// velocity[nx*N+ny][0] = 0.0;
	// }
	// }

	for(int nx=0;nx<N;nx++)
	{
		double x = nx*separation-(N-1)*separation*0.5+offset;
		for(int ny=0;ny<N;ny++)
		{
			double y = ny*separation-(N-1)*separation*0.5+offset;

			nodepos[nx*N+ny][0] = x;
			nodepos[nx*N+ny][1] = ballsize+1.0;
			nodepos[nx*N+ny][2] = y;

			velocity[nx*N+ny][0] = 0.0;
			velocity[nx*N+ny][1] = 0.0;
			velocity[nx*N+ny][2] = 0.0;

			force[nx*N+ny][0] = 0.0;
			force[nx*N+ny][1] = 0.0;
			force[nx*N+ny][2] = 0.0;

			oldforce[nx*N+ny][0] = 0.0;
			oldforce[nx*N+ny][1] = 0.0;
			oldforce[nx*N+ny][2] = 0.0;
		}
	}
}



int max(int a, int b) {
	if (a > b) {
		return a;
	} else {
		return b;
	}
}//end maxOnDevice

int min(int a, int b) {
	if (a > b) {
		return b;
	} else {
		return a;
	}
}//end minOnDevice

double mag(double x, double y, double z) {
	return sqrt(x*x + y*y +z*z);
}//end minOnDevice

double normx(double x, double y, double z) {
	return x/sqrt(x*x + y*y +z*z);
}//end minOnDevice

double normy(double x, double y, double z) {
	return y/sqrt(x*x + y*y +z*z);
}//end minOnDevice

double normz(double x, double y, double z) {
	return z/sqrt(x*x + y*y +z*z);
}//end minOnDevice



//todo: now remove all these parameters as i dont need them
void c_compute_force(int N,int delta, double gravity, double separation, double fcon, double *xForce, double *yForce, double *zForce,
double *xPos, double *yPos, double *zPos){

	// double r12X =0.0;
	// double r12Y =0.0;
	// double r12Z =0.0;
	//r12=vector(0.0,0.0,0.0)
	int nx,ny,dx,dy;
	//double PE=0.0;
	//double len=0.0;

	#pragma omp parallel for default(none) shared(N,delta,gravity,separation,fcon,xForce,yForce,zForce,xPos,yPos,zPos) private(nx,ny,dx,dy)
	for (nx=0; nx<N; nx++)
	{
		for (ny=0; ny<N; ny++)
		{
			xForce[nx*N+ny] = 0.0;
			yForce[nx*N+ny] = -gravity;
			zForce[nx*N+ny] = 0.0;

			int lowerValuedx = max(nx-delta,0);
			int upperValuedx=min(nx+delta+1,N);
			for(dx=lowerValuedx; dx<upperValuedx;dx++)
			{
				int lowerValuedy=max(ny-delta,0);
				int upperValuedy=min(ny+delta+1,N);
				for(dy=lowerValuedy; dy<upperValuedy;dy++)
				{
					double len=sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)) ) *separation;

					if (nx!=dx || ny!=dy)
					{
						double r12X = xPos[dx*N+dy] - xPos[nx*N+ny];
						double r12Y = yPos[dx*N+dy] - yPos[nx*N+ny];
						double r12Z = zPos[dx*N+dy] - zPos[nx*N+ny];
						//PE = PE + fcon*(mag(r12X,r12Y,r12Z)-len)*(mag(r12X,r12Y,r12Z)-len);
						xForce[nx*N+ny] = xForce[nx*N+ny] +fcon*normx(r12X,r12Y,r12Z)*(mag(r12X,r12Y,r12Z)-len);
						yForce[nx*N+ny]= yForce[nx*N+ny] +fcon*normy(r12X,r12Y,r12Z)*(mag(r12X,r12Y,r12Z)-len);
						zForce[nx*N+ny]= zForce[nx*N+ny] +fcon*normz(r12X,r12Y,r12Z)*(mag(r12X,r12Y,r12Z)-len);
					}
				}
			}

		}
	}
}