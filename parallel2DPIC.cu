/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <fstream>
#include <curand_kernel.h>
#include <cufft.h>

float L,LL; int N, C,itera;

using namespace std;

// función Maxwelliana de la distribución de las partículas.
__device__ float distribution (float vb, float aleatorio, curandState *states)     //generador de distribución maxwelliana para la velocidad
{

  // Genera un valor random v
   float fmax = 0.5 * (1.0 + exp (-2.0 * vb * vb));
   float vmin = - 5.0 * vb;
   float vmax = + 5.0 * vb;
   float v;
   float f;
   float x;
   int Idx = blockIdx.x*blockDim.x + threadIdx.x;

   while(true){
	   v = vmin + ((vmax - vmin) * aleatorio);
	   f = 0.5 * (exp (-(v - vb) * (v - vb) / 2.0) +
			    exp (-(v + vb) * (v + vb) / 2.0));
	   x = fmax * aleatorio;
	   if(x > f) aleatorio = curand_uniform(states + Idx);
	   else return v;
   }

}
//Distribución aleatoria de las partículas.
__global__ void distribucionParticulas(float *rx,float *ry,float *vx,float *vy,int N,curandState *states,float vb,float L){
	int Idx = blockIdx.x*blockDim.x + threadIdx.x;

	unsigned int seed = (unsigned int) (clock() * Idx);
	curand_init(seed, 0, 0, states + Idx);

	if(Idx < N){
		 rx[Idx] = L*curand_uniform(states + Idx);    //inicializando la posicion aleatoria en x
		 ry[Idx] = L*curand_uniform(states + Idx);
		 vx[Idx] = distribution(vb,curand_uniform(states + Idx),states);//;L*curand_uniform_float(states + Idx);//distribution(vb,states);                          //inicializa la velocidad con una distribucion maxwelliana
		 vy[Idx] = distribution(vb,curand_uniform(states + Idx),states);//L*curand_uniform_float(states + Idx);//distribution(vb,states);                          //inicializa la velocidad con una distribucion maxwelliana

	}


}
// inicialización de la densidad.
__global__ void inicializacionDensidad(float *ne,int C){
	int Id=blockIdx.x*blockDim.x + threadIdx.x;
		if(Id<(C*C)){
			ne[Id]=0.0;
		}
 }


//Calculo de la densidad en cada celda.

__global__ void calculoDensidad(float *rx, float *ry, float *ne, int N, int C,float L){
	int Id=blockIdx.x*blockDim.x + threadIdx.x;
	 float dx = L / float (C);
	 float dxx=L/float(C*C);
	if(Id<N){

				int jx = int(rx[Id]/dx); //posicion en x de la particula
			    int jy = int(ry[Id]/dx); //posicion en y de la particula
			    float yx = (rx[Id]/dx) - (float)jx; //posicion exacta de la particula en x de la celda "j"
			    //float yy = (ry[Id]/dx) - (float)jy; //posicion exacta de la particula en y de la celda "j"
			    ne[(jy*C)+jx] += (1. - yx)/dxx;
			    if(jx+1==C) ne[(jy*C)] += yx/dxx;
			    else ne[(jy*C)+jx+1] += yx/dxx;

    }

}
//pasar de reales a complejos.

__global__ void real2complex (float *ne, cufftComplex *u, int C)
{
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	int idy = blockIdx.y*blockDim.y+threadIdx.y;
	int index =idy*C+idx;

	if ( idx < C && idy <C)
	 {
		u[index].x = ne[index];
		u[index].y = 0.0f;
	}

}
//__global__ void prueba (cufftComplex *vf, float *vr, int C){
//	int idx = blockIdx.x*blockDim.x+threadIdx.x;
//	int idy = blockIdx.y*blockDim.y+threadIdx.y;
//	int index =idy*C+idx;
//
//	if(idx<C && idy<C){
//
//		vr[index]= (vf[index].x)/((float)C*(float)C*(float)C*(float)C);
//		vr[index]= (vf[index].y)/((float)C*(float)C*(float)C*(float)C);
//
//	}
//}

__global__ void solve_Poisson(cufftComplex *vf, cufftComplex *v, int C,float L){

	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	int idy = blockIdx.y*blockDim.y+threadIdx.y;
	float dx = L / float (C);
	float i,W,Wm,Wn;
	i = (0.0,L);
	W = exp(2.0 * M_PI * i / float(C));
	Wm = L;
	Wn = L;
	if(idx<C && idy<C){
		int index = idy*C+idx;
		float denom;
		denom = 4.0;
		denom -= (Wm + (L / Wm) + Wn +( L / Wn));
		if (denom != 0.0){
				vf[index].x *= dx*dx/denom;
				vf[index].y *= dx*dx/denom;
				}
				Wn *= W;//se multiplica por la constante W
			}
			Wm *= W;
			if(idx<C && idy<C){
			int index = idx*C+idy;
			v[index].x=vf[index].x;
			v[index].y=vf[index].y;
			}

	}



__global__ void complex2real(cufftComplex *v, float *vr, int C){
	/* compute idx and idy, the location of the element in the original CxC array */
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	int idy = blockIdx.y*blockDim.y+threadIdx.y;
	if ( idx < C && idy <C)
		 {
		 int index = idy*C+idx;
		 vr[index] = v[index].x /((float)C*(float)C);
		 }
}



////////////////////////////////////////////////////////////////////////////////////////////////////
int main(){
	// Parametros
	L = 64.0;            // dominio de la solucion 0 <= x <= L (en longitudes de debye)
	//L=LL*LL;
	N = 10000;            // Numero de particulas
	C = 64;            // Número de celdas.
	float vb = 3.0;    // velocidad promedio de los electrones
	//float kappa = 2. * M_PI / (L);
	//float dt=0.1;    // delta tiempo (en frecuencias inversas del plasma)
	//float tmax=10000;  // cantidad de iteraciones. deben ser 100 mil segun el material
	//int skip = int (tmax / dt) / 10; //saltos del algoritmo para reportar datos
	//int itera=0;
	float salida=0.;
	 //float dx = L / float (C);

/////////////////////////////////////////////////////////////////////////////////////////////////////
//Inicializacion de la posición de las particulas en x, y y velocidad en vx,vy del host y dispositivo.
	float *rx_h,*ry_h,*vx_h,*vy_h;
	float *rx_d,*ry_d, *vx_d,*vy_d;
////////////////////////////////////////////////////////////////////////////////////////////////////
	// inicialización de las variables de densidad del host y dispositivo.
	float *ne_h;
	float *ne_d;
	float *vr_h;
	float *vr_d;
////////////////////////////////////////////////////////////////////////////////////////////////////
	//inicializacion tipo complex a real.

	cufftComplex *u_complex_d,*vf_complex_d,*v_complex_d ;
	cudaMalloc((void**)&u_complex_d,sizeof(cufftComplex)*C*C);
	cudaMalloc((void**)&vf_complex_d,sizeof(cufftComplex)*C*C);
	cudaMalloc((void**)&v_complex_d,sizeof(cufftComplex)*C*C);

////////////////////////////////////////////////////////////////////////////////////////////////////
	int size = N*sizeof(float);
	int size_ne=C*C*sizeof(float);

//////////////////////////////////////////////////////////////////////////////////////////////////////
	//reserva en memoria al host
	rx_h = (float *)malloc(size);
	ry_h = (float *)malloc(size);
	vx_h = (float *)malloc(size);
	vy_h = (float *)malloc(size);
	ne_h = (float *)malloc(size_ne);
	vr_h = (float *)malloc(size_ne);

//////////////////////////////////////////////////////////////////////////////////////////////////////
	//reserva de memoria del dispositivo.
	cudaMalloc((void **)&rx_d,size);
	cudaMalloc((void **)&ry_d,size);
	cudaMalloc((void **)&vx_d,size);
	cudaMalloc((void **)&vy_d,size);
	cudaMalloc((void **)&ne_d,size_ne);
	cudaMalloc((void **)&vr_d,size_ne);

////////////////////////////////////////////////////////////////////////////////////////////////////

	//valores aleatorios y tamaños de los vectores.
	curandState *devStates;
	cudaMalloc((void **) &devStates, N * sizeof(curandState));


	float blockSize = 1024;
	dim3 dimBlock (ceil(N/blockSize), 1, 1);
	dim3 dimBlock2 (ceil((C*C)/blockSize), 1, 1);
	dim3 dimGrid (blockSize, 1, 1);


	distribucionParticulas<<<blockSize,dimBlock>>>(rx_d,ry_d,vx_d,vy_d,N,devStates,vb,L);
	cudaDeviceSynchronize();

	inicializacionDensidad<<<blockSize,dimBlock2>>>(ne_d,C);
	cudaDeviceSynchronize();

	calculoDensidad<<<blockSize,dimBlock>>>(rx_d,ry_d,ne_d,N,C,L);
	cudaDeviceSynchronize();

	cufftHandle plan;
	cufftPlan2d(&plan, C, C, CUFFT_C2C);

	real2complex<<<blockSize,dimBlock2>>>(ne_d,u_complex_d,C);
	cudaDeviceSynchronize();

	cufftExecC2C (plan, u_complex_d, vf_complex_d, CUFFT_FORWARD);
	// dividir el resultado por C4
	//prueba<<<dimGrid, dimBlock2>>> (vf_complex_d,vr_d,C);

	v_complex_d[0].x=0.0;
	v_complex_d[0].y=0.0;


	solve_Poisson<<<dimGrid, dimBlock2>>> (vf_complex_d,v_complex_d,C,L);
	cudaDeviceSynchronize();

	cufftExecC2C (plan, v_complex_d, v_complex_d, CUFFT_INVERSE);

	complex2real<<<dimGrid, dimBlock2>>> (v_complex_d,vr_d,C);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//posición en x.
	cudaMemcpy(rx_h, rx_d, size, cudaMemcpyDeviceToHost);

	// posición en y.
	cudaMemcpy(ry_h, ry_d, size, cudaMemcpyDeviceToHost);

	// velocidad en x.
	cudaMemcpy(vx_h, vx_d, size, cudaMemcpyDeviceToHost);

	//velocidad en y.
	cudaMemcpy(vy_h, vy_d, size, cudaMemcpyDeviceToHost);
	//inicializacion densidades
	cudaMemcpy(ne_h, ne_d, size_ne, cudaMemcpyDeviceToHost);
	//calculo poisson
	cudaMemcpy (vr_h , vr_d, size_ne, cudaMemcpyDeviceToHost);

	///////////////////Imprimir los resultados en archivos//////////////////////
	ofstream init;
		init.open("distribucionInicial.txt");
		  		    for (int i = 0; i < N; i++){
		  		    	init<<rx_h[i]<<" "<<ry_h[i]<<" "<<vx_h[i]<<" "<<vy_h[i]<<endl;

		  		    }

		  		    init.close();


		init.open("salida_densidad3.txt");
					for (int i = 0; i < C*C; i++){
						init<<ne_h[i]<<endl;
						salida+=ne_h[i];
					}

					init.close();
					cout<<salida<<" "<<endl;

		init.open("entrada_poisson");
					for (int i = 0; i < C; i++){
						for (int j = 0; j < C; j++){
								init<<ne_h[(C*i)+j]<<" ";
						}
						init<<endl;
					}
					init.close();

		init.open("poisson");
				for (int i = 0; i < C; i++){
					for (int j = 0; j < C; j++){
						init<< vr_h[(C*j)+i]<<" ";
							}
							init<<endl;
						}
						init.close();

////////////////////Liberar memoria//////////////////////////
	free(rx_h);
	free(ry_h);
	free(vx_h);
	free(vy_h);
	free(ne_h);
	free(vr_h);
	cufftDestroy(plan);
	cudaFree(rx_d);
	cudaFree(ry_d);
	cudaFree(vx_d);
	cudaFree(vy_d);
	cudaFree(ne_d);
	cudaFree(vr_d);
	cudaFree(u_complex_d);
	cudaFree(vf_complex_d);
	cudaFree(v_complex_d);
	return (0);

}
