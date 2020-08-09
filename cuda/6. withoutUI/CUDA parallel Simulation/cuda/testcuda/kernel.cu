//#include "Python.h"
//
//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include "arrayobject.h"
//#include <math.h>
//#include <string.h>
//#include <stdio.h>
//
//
//#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
//#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
//
//static PyObject *py_cComputeForce(PyObject *self, PyObject *args);
//PyObject* c_compute_force(int N,int delta, double gravity, double separation, double fcon, double *xForce, double *yForce, double *zForce,
//	double *xPos, double *yPos, double *zPos);
//PyObject* c_compute_force_Cudafied(int N,int delta, double gravity, double separation, double fcon, double *xForce, double *yForce, double *zForce,
//	double *xPos, double *yPos, double *zPos);
//
//
//__device__ int maxOnDevice(int a, int b) {
//	if (a > b) {
//		return a;
//	} else {
//		return b;
//	}
//}//end max
//
//__device__ int minOnDevice(int a, int b) {
//	if (a > b) {
//		return b;
//	} else {
//		return a;
//	}
//}//end min
//
//__device__ double magOnDevice(double x, double y, double z) {
//	return sqrt(x*x + y*y +z*z);
//}//end min
//
//__device__ double normxOnDevice(double x, double y, double z) {
//	return x/sqrt(x*x + y*y +z*z);
//}//end min
//
//__device__ double normyOnDevice(double x, double y, double z) {
//	return y/sqrt(x*x + y*y +z*z);
//}//end min
//
//__device__ double normzOnDevice(double x, double y, double z) {
//	return z/sqrt(x*x + y*y +z*z);
//}//end min
//
//int maxOnHost(int a, int b) {
//	if (a > b) {
//		return a;
//	} else {
//		return b;
//	}
//}//end max
//
//int minOnHost(int a, int b) {
//	if (a > b) {
//		return b;
//	} else {
//		return a;
//	}
//}//end min
//
//double mag(double x, double y, double z) {
//	return sqrt(x*x + y*y +z*z);
//}//end min
//
//double normx(double x, double y, double z) {
//	return x/sqrt(x*x + y*y +z*z);
//}//end min
//
//double normy(double x, double y, double z) {
//	return y/sqrt(x*x + y*y +z*z);
//}//end min
//
//double normz(double x, double y, double z) {
//	return z/sqrt(x*x + y*y +z*z);
//}//end min
//
//static PyObject *py_cComputeForce(PyObject *self, PyObject *args){
//	int N,delta,ok;
//	double gravity , separation, fcon;
//
//	PyObject *xposarray, *yposarray, *zposarray;
//	PyObject *xforcearray, *yforcearray, *zforcearray;
//	double *xPos, *yPos, *zPos;
//	double *xForce, *yForce, *zForce;
//
//	//double pe;
//	PyObject *lst;
//
//	ok = PyArg_ParseTuple(args,"iidddOOOOOO",&N,&delta,&gravity,&separation,&fcon,&xforcearray,&yforcearray,&zforcearray,
//		&xposarray,&yposarray,&zposarray);
//
//	/*if (true){
//	fprintf(stdout,"N= %d, delta=%d gravity=%f, separation=%f fcon=%f \n",N,delta,gravity,separation,fcon);
//	exit(1);
//	}*/
//
//	if (!ok){
//		fprintf(stderr,"Error (cComputeForce) in parsing arguments\n");
//		exit(1);
//	}
//
//	xPos = DDATA(xposarray);
//	yPos = DDATA(yposarray);
//	zPos = DDATA(zposarray);
//	xForce = DDATA(xforcearray);
//	yForce = DDATA(yforcearray);
//	zForce = DDATA(zforcearray);
//
//	//pe = c_compute_force(N,delta,gravity,separation,fcon,xForce,yForce,zForce,xPos,yPos,zPos);
//	
//	lst = c_compute_force(N,delta,gravity,separation,fcon,xForce,yForce,zForce,xPos,yPos,zPos);
//	//lst = c_compute_force_Cudafied(N,delta,gravity,separation,fcon,xForce,yForce,zForce,xPos,yPos,zPos);
//
//	return lst;
//	//return Py_BuildValue("d",pe);
//}
//
////ramneek: trying to cuda'fy this code
//__global__ void MyKernel(int *Nptr,int *deltaptr, double *gravityptr, double *separationptr, double *fconptr, double *xForce, double *yForce, double *zForce,
//	double *xPos, double *yPos, double *zPos/*, PyObject *force_on_each_ball_list*/ )
//{
//	int N = *Nptr;
//	int delta= *deltaptr;
//	double gravity= *gravityptr;
//	double separation = *separationptr;
//	double fcon = *fconptr;
//
//	double len=0.0;
//	double r12X =0.0;
//	double r12Y =0.0;
//	double r12Z =0.0;
//	double PE=0.0;
//
//	
//	int nx = blockDim.x * blockIdx.x + threadIdx.x;//use this place of nx
//	int ny = blockDim.x * blockIdx.x + threadIdx.y;//use this place of ny
//	if(!(nx< N && ny <N))
//		return;
//
//	xForce[nx*N+ny] = 0.0;
//	yForce[nx*N+ny] = -gravity;
//	zForce[nx*N+ny] = 0.0;
//
//	int lowerValuedx = maxOnDevice(nx-delta,0);
//	int upperValuedx=minOnDevice(nx+delta+1,N);
//	for(int dx=lowerValuedx; dx<upperValuedx;dx++)
//	{
//		int lowerValuedy=maxOnDevice(ny-delta,0);
//		int upperValuedy=minOnDevice(ny+delta+1,N);
//		for(int dy=lowerValuedy; dy<upperValuedy;dy++)
//		{
//			len=sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)) ) *separation;
//
//			bool condition = ny!=dy;
//			bool condition1 = nx!=dx;
//			if (condition || condition1)
//			//if (nx!=dx || ny!=dy)
//			{
//				r12X = xPos[dx*N+dy] - xPos[nx*N+ny];
//				r12Y = yPos[dx*N+dy] - yPos[nx*N+ny];
//				r12Z = zPos[dx*N+dy] - zPos[nx*N+ny];
//				PE = PE + fcon*(magOnDevice(r12X,r12Y,r12Z)-len)*(magOnDevice(r12X,r12Y,r12Z)-len);
//				xForce[nx*N+ny] = xForce[nx*N+ny] +fcon*normxOnDevice(r12X,r12Y,r12Z)*(magOnDevice(r12X,r12Y,r12Z)-len);
//				yForce[nx*N+ny]= yForce[nx*N+ny] +fcon*normyOnDevice(r12X,r12Y,r12Z)*(magOnDevice(r12X,r12Y,r12Z)-len);
//				zForce[nx*N+ny]= zForce[nx*N+ny] +fcon*normzOnDevice(r12X,r12Y,r12Z)*(magOnDevice(r12X,r12Y,r12Z)-len);
//
//				//i tried to first get the item and modify it and set it back.
//				//but then i thought....why not directly set the new item that position if it is anyways going to overwrite it.
//				//PyObject *temp=PyList_GetItem(force_on_each_ball_list, nx*N+ny);
//				/*ok = PyArg_ParseTuple(temp,"ddd",&N,&delta,&gravity,&separation,&fcon,&xforcearray,&yforcearray,&zforcearray,
//				&xposarray,&yposarray,&zposarray);*/
//
//				//ramneek: get the items out of xForce, yForce and zForce in the host method and use the follwing statements there. 
//				/*PyObject *item = Py_BuildValue("(ddd)",xForce[nx*N+ny],yForce[nx*N+ny],zForce[nx*N+ny]);
//				PyList_SetItem(force_on_each_ball_list, nx*N+ny, item);*/
//			}
//		}
//	}
//	/*   int i = threadIdx.x;
//    c[i] = a[i] + b[i];*/
//}
//
//PyObject* c_compute_force_Cudafied(int N,int delta, double gravity, double separation, double fcon, double *xForce, double *yForce, double *zForce,
//	double *xPos, double *yPos, double *zPos){
//
//		double r12X =0.0;
//		double r12Y =0.0;
//		double r12Z =0.0;
//		//r12=vector(0.0,0.0,0.0)
//		int nx,ny,dx,dy;
//		double PE=0.0;
//		double len=0.0;
//
//
//		//now i m building the force afresh every time i call the method. So no need to pass the force arguments array to it.
//		//just return a force array(list) to python form here for it to consume it.
//		//now issue is force is list in python. i can make a list from here and return it.
//		//but i have to let some values remiain (0,0,0) in the list and some have to have force (fx,fy,fz).
//		//what i could do was make a list here.add (0,0,0) to it for every element by using PyList_Append.
//		//but when i run i main loop, i have to update the (fx,fy,fz) values for which interaction is taking place.
//		//so i have to get a element out of the list and update it. so find a method for updating which shoudl be similiar to PyList_Append.
//
//
//		PyObject *force_on_each_ball_list = Py_BuildValue("[]");
//		if (!force_on_each_ball_list)
//			return NULL;
//
//		for (nx=0; nx<N; nx++)
//		{
//			for (ny=0; ny<N; ny++)
//			{
//				PyObject *lc = Py_BuildValue("(ddd)",0.0,0.0,0.0);
//				PyList_Append(force_on_each_ball_list,lc);
//				Py_DECREF(lc);
//			}
//		}
//
//		//allocate memory on device here
//		int *dev_N =0;
//		int *dev_delta=0;
//		double *dev_gravity=0;
//		double *dev_separation=0;
//		double *dev_fcon=0;
//
//		double *dev_xForce = 0;
//		double *dev_yForce = 0;
//		double *dev_zForce = 0;
//		double *dev_xPos = 0;
//		double *dev_yPos = 0;
//		double *dev_zPos = 0;
//
//		cudaError_t cudaStatus;
//
//		// Choose which GPU to run on, change this on a multi-GPU system.
//		cudaStatus = cudaSetDevice(0);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
//			exit(0);
//			goto Error;
//		}
//
//		cudaStatus = cudaMalloc((void**)&dev_N, sizeof(int));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMalloc((void**)&dev_delta, sizeof(int));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMalloc((void**)&dev_gravity, sizeof(double));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMalloc((void**)&dev_separation, sizeof(double));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMalloc((void**)&dev_fcon, sizeof(double));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//		// Allocate GPU buffers for 6 vectors    .
//		cudaStatus = cudaMalloc((void**)&dev_xForce, N*N * sizeof(double));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMalloc((void**)&dev_yForce,  N*N * sizeof(double));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMalloc((void**)&dev_zForce,  N*N * sizeof(double));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//		cudaStatus = cudaMalloc((void**)&dev_xPos, N*N * sizeof(double));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMalloc((void**)&dev_yPos,  N*N * sizeof(double));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMalloc((void**)&dev_zPos,  N*N * sizeof(double));
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMalloc failed!");
//			goto Error;
//		}
//
//
//		cudaStatus = cudaMemcpy(dev_N, &N,  sizeof(int), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMemcpy(dev_delta, &delta,  sizeof(int), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMemcpy(dev_gravity, &gravity,  sizeof(double), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMemcpy(dev_separation, &separation,  sizeof(double), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMemcpy(dev_fcon, &fcon,  sizeof(double), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		// Copy input vectors from host memory to GPU buffers.
//		cudaStatus = cudaMemcpy(dev_xForce, xForce, N*N * sizeof(double), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMemcpy(dev_yForce, yForce, N*N * sizeof(double), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMemcpy(dev_zForce, zForce, N*N * sizeof(double), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMemcpy(dev_xPos, xPos, N*N * sizeof(double), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMemcpy(dev_yPos, yPos, N*N * sizeof(double), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		cudaStatus = cudaMemcpy(dev_zPos, zPos, N*N * sizeof(double), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		dim3 threadsPerBlock(16, 16); 
//		int blocksPerGrid = (N*N + (threadsPerBlock.x*threadsPerBlock.y) - 1) / (threadsPerBlock.x*threadsPerBlock.y);
//		MyKernel<<<blocksPerGrid, threadsPerBlock>>>( dev_N,dev_delta,dev_gravity,dev_separation,dev_fcon,dev_xForce,
//			dev_yForce,dev_zForce,dev_xPos,dev_yPos,dev_zPos); 
//
//		////launch the kernel here
//		//int numBlocks = 1; 
//		//dim3 threadsPerBlock(N, N); 
//		//MyKernel<<<numBlocks, threadsPerBlock>>>( dev_N,dev_delta,dev_gravity,dev_separation,dev_fcon,dev_xForce,
//		//	dev_yForce,dev_zForce,dev_xPos,dev_yPos,dev_zPos); 
//
//		// cudaDeviceSynchronize waits for the kernel to finish, and returns
//		// any errors encountered during the launch.
//		cudaStatus = cudaDeviceSynchronize();
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
//			goto Error;
//		}
//
//		// Copy output vector from GPU buffer to host memory.
//		cudaStatus = cudaMemcpy(xForce, dev_xForce, N*N * sizeof(double), cudaMemcpyDeviceToHost);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//		// Copy output vector from GPU buffer to host memory.
//		cudaStatus = cudaMemcpy(yForce, dev_yForce, N*N * sizeof(double), cudaMemcpyDeviceToHost);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//		// Copy output vector from GPU buffer to host memory.
//		cudaStatus = cudaMemcpy(zForce, dev_zForce, N*N * sizeof(double), cudaMemcpyDeviceToHost);
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaMemcpy failed!");
//			goto Error;
//		}
//
//		//populate the list here
//		for (nx=0; nx<N; nx++)
//		{
//			for (ny=0; ny<N; ny++)
//			{
//				int lowerValuedx = maxOnHost(nx-delta,0);
//				int upperValuedx=minOnHost(nx+delta+1,N);
//				for(dx=lowerValuedx; dx<upperValuedx;dx++)
//				{
//					int lowerValuedy=maxOnHost(ny-delta,0);
//					int upperValuedy=minOnHost(ny+delta+1,N);
//					for(dy=lowerValuedy; dy<upperValuedy;dy++)
//					{
//						if (nx!=dx || ny!=dy)
//						{
//							PyObject *item = Py_BuildValue("(ddd)",xForce[nx*N+ny],yForce[nx*N+ny],zForce[nx*N+ny]);
//							PyList_SetItem(force_on_each_ball_list, nx*N+ny, item);
//						}
//					}
//				}
//			}
//		}
//		
//Error:
//		cudaFree(dev_N);
//		cudaFree(dev_delta);
//		cudaFree(dev_gravity);
//		cudaFree(dev_separation);
//		cudaFree(dev_fcon);
//
//		cudaFree(dev_xForce);
//		cudaFree(dev_yForce);
//		cudaFree(dev_zForce);
//		cudaFree(dev_xPos);
//		cudaFree(dev_yPos);
//		cudaFree(dev_zPos);
//
//		//pe returned is correct. so i m not returning it now. instead i will return a list of with each element being (xForce,yForce,zForce)
//		//return PE;
//		return force_on_each_ball_list;
//}
//
////ramneek:end
//
////todo: now remove all these parameters as i dont need them
//PyObject* c_compute_force(int N,int delta, double gravity, double separation, double fcon, double *xForce, double *yForce, double *zForce,
//	double *xPos, double *yPos, double *zPos){
//
//		double r12X =0.0;
//		double r12Y =0.0;
//		double r12Z =0.0;
//		//r12=vector(0.0,0.0,0.0)
//		int nx,ny,dx,dy;
//		double PE=0.0;
//		double len=0.0;
//
//
//		//now i m building the force afresh every time i call the method. So no need to pass the force arguments array to it.
//		//just return a force array(list) to python form here for it to consume it.
//		//now issue is force is list in python. i can make a list from here and return it.
//		//but i have to let some values remiain (0,0,0) in the list and some have to have force (fx,fy,fz).
//		//what i could do was make a list here.add (0,0,0) to it for every element by using PyList_Append.
//		//but when i run i main loop, i have to update the (fx,fy,fz) values for which interaction is taking place.
//		//so i have to get a element out of the list and update it. so find a method for updating which shoudl be similiar to PyList_Append.
//
//
//		PyObject *force_on_each_ball_list = Py_BuildValue("[]");
//		if (!force_on_each_ball_list)
//			return NULL;
//
//		for (nx=0; nx<N; nx++)
//		{
//			for (ny=0; ny<N; ny++)
//			{
//				PyObject *lc = Py_BuildValue("(ddd)",0.0,0.0,0.0);
//				PyList_Append(force_on_each_ball_list,lc);
//				Py_DECREF(lc);
//			}
//		}
//
//		for (nx=0; nx<N; nx++)
//		{
//			for (ny=0; ny<N; ny++)
//			{
//				xForce[nx*N+ny] = 0.0;
//				yForce[nx*N+ny] = -gravity;
//				zForce[nx*N+ny] = 0.0;
//
//				int lowerValuedx = max(nx-delta,0);
//				int upperValuedx=min(nx+delta+1,N);
//				for(dx=lowerValuedx; dx<upperValuedx;dx++)
//				{
//					int lowerValuedy=max(ny-delta,0);
//					int upperValuedy=min(ny+delta+1,N);
//					for(dy=lowerValuedy; dy<upperValuedy;dy++)
//					{
//						len=sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)) ) *separation;
//
//						if (nx!=dx || ny!=dy)
//						{
//							r12X = xPos[dx*N+dy] - xPos[nx*N+ny];
//							r12Y = yPos[dx*N+dy] - yPos[nx*N+ny];
//							r12Z = zPos[dx*N+dy] - zPos[nx*N+ny];
//							PE = PE + fcon*(mag(r12X,r12Y,r12Z)-len)*(mag(r12X,r12Y,r12Z)-len);
//							xForce[nx*N+ny] = xForce[nx*N+ny] +fcon*normx(r12X,r12Y,r12Z)*(mag(r12X,r12Y,r12Z)-len);
//							yForce[nx*N+ny]= yForce[nx*N+ny] +fcon*normy(r12X,r12Y,r12Z)*(mag(r12X,r12Y,r12Z)-len);
//							zForce[nx*N+ny]= zForce[nx*N+ny] +fcon*normz(r12X,r12Y,r12Z)*(mag(r12X,r12Y,r12Z)-len);
//
//							//i tried to first get the item and modify it and set it back.
//							//but then i thought....why not directly set the new item that position if it is anyways going to overwrite it.
//							//PyObject *temp=PyList_GetItem(force_on_each_ball_list, nx*N+ny);
//							/*ok = PyArg_ParseTuple(temp,"ddd",&N,&delta,&gravity,&separation,&fcon,&xforcearray,&yforcearray,&zforcearray,
//							&xposarray,&yposarray,&zposarray);*/
//
//							PyObject *item = Py_BuildValue("(ddd)",xForce[nx*N+ny],yForce[nx*N+ny],zForce[nx*N+ny]);
//							PyList_SetItem(force_on_each_ball_list, nx*N+ny, item);
//						}
//					}
//				}
//
//			}
//		}
//		//pe returned is correct. so i m not returning it now. instead i will return a list of with each element being (xForce,yForce,zForce)
//		//return PE;
//		return force_on_each_ball_list;
//}
//
//static PyMethodDef ClothEngine_methods[] =
//{	
//	{"cComputeForce",py_cComputeForce,METH_VARARGS},
//	{NULL,NULL} /* Sentinel */
//};
//
//PyMODINIT_FUNC
//	initClothEngine(){
//		(void) Py_InitModule("ClothEngine",ClothEngine_methods);
//}
