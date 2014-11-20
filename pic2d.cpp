//============================================================================
 // Name        : pic2d.cpp
// Author      : spiderzelda
// Version     :
// Copyright   : GNU
// Description : Particle in Cell 2D code
//============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <fstream>
#include <fftw.h>
#include <complex.h>



using namespace std;

//funciones
double distribution(double vb);
void Density3(vector<double>rx, vector<double>ry, vector<double>& n);
void Electric2 (vector<double> phi, vector<double>& Ex, vector<double>& Ey);
void Poisson_prueba1 (vector<double>& u, vector<double> v, double kappa);
void Load (vector<double> r, vector<double> v, vector<double>& y);
void UnLoad (vector<double> y, vector<double>& r, vector<double>& v);
void rk4_fixed2 (double& x, vector<double>& y1, vector<double>& y2, double h);
void eval (double t, vector<double> y1,vector<double> y2, vector<double>& dydt1,vector<double>& dydt2);
void Output2 (char* fn1, char* fn2, double t, vector<double> rx,vector<double> rv, vector<double> vx,vector<double> vy);



//parametros
double L,LL; int N, C,itera;


int main()
{
  // Parametros
  L =64.0;            // dominio de la solucion 0 <= x <= L (en longitudes de debye) LA CANTIDAD DE CELDAS EN 1D DEBEN SER IGUAL AL ESPACIO DE SIMULACION
  N =10000;            // Numero de particulas
  C = 64;            // Numero de celdas EN UNA DIMENSION, EL TOTAL DE CELDAS ES C*C
  double vb = 3.0;    // velocidad promedio de los electrones
  double dt=0.1;    // delta tiempo (en frecuencias inversas del plasma)
  double tmax=40;  // cantidad de iteraciones. deben ser 100 mil segun el material
  int skip = int (tmax / dt) / 10; //saltos del algoritmo para reportar datos


  vector<double> rx,ry,vx,vy, n(C*C); //r: posicion de las particulas, v: velocidad de particulas n: densidad de particulas por celda


  //PARA PARALELIZAR. se hace por cada part√¨cula.
  double t = 0.;
  int seed = time (NULL); srand (seed);
  for (int i = 0; i < N; i++)
  {

	  rx.push_back(L*double (rand ()) / double (RAND_MAX));    //inicializando la posicion aleatoria en x
	  ry.push_back(L*double (rand ()) / double (RAND_MAX));




      vx.push_back(distribution(vb));                          //inicializa la velocidad con una distribucion maxwelliana
      vy.push_back(distribution(vb));                          //inicializa la velocidad con una distribucion maxwelliana


  }




  char* phase[11]; char* data[11]; //archivos para almacenar los datos de salida
    phase[0] = "phase0.txt";phase[1] = "phase1.txt";phase[2] = "phase2.txt";
    phase[3] = "phase3.txt";phase[4] = "phase4.txt";phase[5] = "phase5.txt";
    phase[6] = "phase6.txt";phase[7] = "phase7.txt";phase[8] = "phase8.txt";
    phase[9] = "phase9.txt";phase[10] = "phase10.txt";data[0] = "data0.txt";
    data[1] = "data1.txt"; data[2] = "data2.txt"; data[3] = "data3.txt";
    data[4] = "data4.txt"; data[5] = "data5.txt"; data[6] = "data6.txt";
    data[7] = "data7.txt"; data[8] = "data8.txt"; data[9] = "data9.txt";
    data[10] = "data10.txt";

    clock_t t1 = clock();   // variable para medir el tiempo de ejecucion



    Output2 (phase[0], data[0], t, rx, ry, vx, vy); //inicializacion del algoritmo y reporte de salida de estados iniciales en phase0 y data0


  // iterando la solucion
      vector<double> y1(2*N), y2(2*N);  //y1 y y2 son vectores que recibe la funcion de runge-kutta, cada vector posee la posicion y la velocidad de las particulas en 1 dimension
      //y1 tiene la posicion y la velocidad en x, mientras que y2 tiene la posicion y la velocidad en "y"
      int d1, d2 = 0;
      Load (rx, vx, y1);    //carga los valores de la posicion y la velocidad en el vector "y"
      Load (ry, vy, y2);




      for (int k = 1; k <= 10; k++)
      {
          for (int kk = 0; kk < skip; kk++) //estas iteraciones garantizan que el ciclo se itera la cantidad de veces especificada, arrojando 10 salidas intermedias
            {
               // Take time-step

        	   rk4_fixed2(t, y1, y2, dt);//se ejecuta runge kutta4, el cual itera para hallar las nuevas velocidades y posiciones teniendo las antiguas velocidades y el campo electrico


               // asegurarse que todas las partiΩculas estan dentro del espacio de fase (reinyeccion)
               for (int i = 0; i < N; i++)
                 {
                   if (y1[i] < 0.) y1[i] += L;
                   if (y1[i] > L) y1[i] -= L;
                   if (y2[i] < 0.) y2[i] += L;
                   if (y2[i] > L) y2[i] -= L;
                 }


               //las siguienes instrucciones son para el control del tiempo de ejecucion del algoritmo.
               d1=((kk*10)/skip) + ((k-1)*10);
               if (d1>d2){
               d2=d1;
               cout<<((kk*10)/skip) + ((k-1)*10)<<"% completo"<<endl;

               }
            }

          cout<<k*10<<"% completo"<<endl;
          if(k*10==1) {
        	  clock_t t3 = clock();
        	  cout<<"faltan "<<((t3-t1)*99)/60<<"segundos "<<endl;
          }

          d1=0;
          d2=0;

          // Output data, se genera la salida 10 veces en cada ejecucion
          UnLoad (y1, rx, vx);
          UnLoad (y2, ry, vy);
          Output2(phase[k], data[k], t, rx,ry, vx,vy);

        }



// mas instrucciones para la medida del tiempo de ejecucion
      clock_t t2 = clock();
      double exectime= (t2-t1)/1000000.0;
      int horas = exectime/60/60;
      int minutos =  exectime/60 - horas;
      double segundos = exectime - int(minutos*60);


      cout<<"tiempo algoritmo con "<<N<<" particulas, "<<C<<" celdas, dt = "<<dt<<", iteraciones = "<<tmax<<endl;
      cout<<"tiempo total en segundos: "<<exectime<<endl;
      cout<<((t2-t1)/1000000)/86400<<" horas, "<<minutos<<" minutos, "<<segundos <<" segundos"<<endl;

      ofstream params;
      params.open("params.txt");
      params<<"N: "<<N<<endl<<"C: "<<C<<endl<<"L: "<<L<<endl<<"vb: "<<vb<<endl<<"dt: "<<dt<<endl<<"tmax: "<<tmax<<endl<<"iteraciones: "<<tmax/dt<<endl;
      params<<((t2-t1)/1000000)/86400<<" horas, "<<minutos<<" minutos, "<<segundos <<" segundos"<<endl;
      params.close();

      return 0;

}

double distribution (double vb)     //generador de distribucion maxwelliana para la velocidad
{
  // inicializa el generador aleatorio
  static int flag = 0;
  if (flag == 0)
    {
      int seed = time (NULL);
      srand (seed);
      flag = 1;
    }

  // Genera un valor random v
  double fmax = 0.5 * (1. + exp (-2. * vb * vb));
  double vmin = - 5. * vb;
  double vmax = + 5. * vb;
  double v = vmin + (vmax - vmin) * double (rand ()) / double (RAND_MAX);


  // Acceptar y reinyectar particulas
  double f = 0.5 * (exp (-(v - vb) * (v - vb) / 2.) +
		    exp (-(v + vb) * (v + vb) / 2.));
  double x = fmax * double (rand ()) / double (RAND_MAX);
  if (x > f) return distribution (vb);
  else
  {
	  return v;

  }

}

void Density3 (vector<double> rx, vector<double> ry, vector<double>& n) //calcula la densidad de carga de cada celda segun la cantidad de particulas presentes en ella.
{
	double salida=0.; //variable de control para assegurarse que la densidad de carga es consevativa
	double dx = L / double (C); //delta x, tamaño de cada celda
	double dxx = L / double (C*C);

	for (int i = 0; i < N; i++)
	{
		int jx = int(rx[i]/dx); //posicion en x de la particula
	    int jy = int(ry[i]/dx); //posicion en y de la particula
	    double yx = (rx[i]/dx) - (double)jx; //posicion exacta de la particula en x de la celda "j"
	    //double yy = (ry[i]/dx) - (double)jy; //posicion exacta de la particula en y de la celda "j"
	    n[(jy*C)+jx] += (1. - yx)/dxx; //se le carga a cada celda el valor correspondiente de densidad
	    if(jx+1==C) n[(jy*C)] += yx/dxx; // si es la ultima celda de la fila, el valor de la carga restante se le carga a la primera celda de esa fila
	    else n[(jy*C)+jx+1] += yx/dxx; // si no es la ultima celda de la fila, se le carga el valor de la densidad restante a la celda de la derecha
	}

	ofstream init;
	init.open("salida_densidad.txt"); //archivo para corroborar la salida de la densidad
	for (int i = 0; i < C*C; i++){
		init<<n[i]<<endl;
		salida+=n[i];
	}
	init.close();
	//cout<<salida<<" "<<dx<<endl;


}

// se usa la libreria fftw para calcular la transformada rapida de Fourier
// los vectores de entrada y salida son de tamaÔøΩ  o C
// se calcula la fft del vector f en los vectores Fr y Fi

void Poisson_prueba1 (vector<double>&u, vector<double> v, double kappa) // recibe el vector v de densidad, la constante kappa  y arroja el resultado de potencial electroestatico.
{
	double dx=L/double(C-1); //delta x, el tamaño de la malla

	vector<double> Vr(C*C), Vi(C*C), Ur(C*C), Ui(C*C); // crea el vector de imagiarios y reales hasta la cantida de celdas y es igual para los vectores del potencial electroestatico.

	fftw_complex U[C][C]; //matriz compleja que almacena los resultados
	fftw_complex ff[C], FF[C]; //vectores temporales para realizar las transformadas por filas y columnas
	//la variable ff va a servir como la entrada para la funcion que calcula la transformada, y FF sera la salida

/*
  //estas lineas corresponden a la carga puntual colocada para verificar el funcionamiento de la funcion
	  for(int x=0;x<C;x++){
		  for(int y=0;y<C;y++){
			  v[(x*C)+y]=0.0;
		        	}
		        }
	  v[(C*C/2)+C/2]=10.0;
*/


	ofstream init; //archivo que imprime la matriz de espacio de simulacion con la densidad de carga de cada celda
	init.open("entrada_poisson");
	for (int i = 0; i < C; i++){
		for (int j = 0; j < C; j++){
			init<<v[(C*i)+j]<<" ";
		}
		init<<endl;
	}
	init.close();


	for(int i=0;i<C;i++){
		for(int j=0;j<C;j++){
			c_re (ff[j]) = v[(C*i)+j];//asignando valores a la parte real de ff
			c_im (ff[j]) = 0.;
		}

		fftw_plan p = fftw_create_plan (C, FFTW_FORWARD, FFTW_ESTIMATE); //se crea el plan para realizar la transformada hacia adelante
		fftw_one (p, ff, FF); //se ejecuta el plan en 1D, primero por cada fila
		fftw_destroy_plan (p);//se libera memoria
		for(int x=0;x<C;x++){
			U[i][x].re=FF[x].re/double(C*C);//el resultado real de la transformada se pone en la matriz U
			U[i][x].im=FF[x].im/double(C*C);//al igual que la parte imaginaria
			//se descargan los valores de la salida FF en la matriz original U
			}
	}


	init.close();
	for(int i=0;i<C;i++){
		for(int j=0;j<C;j++){
			c_re (ff[j]) = U[j][i].re;//se toman los valores encontrados en la transformada en x y se disponen para realizar la transformada en la direccion "y"
			c_im (ff[j]) = U[j][i].im;
		}
		fftw_plan p = fftw_create_plan (C, FFTW_FORWARD, FFTW_ESTIMATE);//se crea el plan para realizar la transformada la direccion "y"
		fftw_one (p, ff, FF);//se ejecuta el plan en 1D, ahora por columnas, con los resultados adquiridos previamente en la transformada en la direccion "x"
		fftw_destroy_plan (p);//se libera memoria
		for(int j=0;j<C;j++){
			U[j][i].re=FF[j].re/double(C*C);//el resultado real de la transformada se pone en la matriz U
			U[j][i].im=FF[j].im/double(C*C);//al igual que la parte imaginaria
		}
	}



	//calculo de poisson:
	// ver documento para entender esquema de discretizacion

	U[0][0].re =0.;

	complex<double> i(0.0, L); //creamos una variable compleja para poder aplicar la discretizacion.
	complex<double> W = exp(2.0 * M_PI * i / double(C));
	complex<double> Wm = L, Wn = L;
	for (int m = 0; m < C; m++)
	{
		for (int n = 0; n < C; n++)
		{
			complex<double> denom = 4.0;
			denom -= Wm + L / Wm + Wn + L / Wn; //se calcula el denominador para cada celda, segun el equema de discretizacion
			if (denom != 0.0){
				U[m][n].re *= dx *dx / denom.real();
				U[m][n].im *= dx *dx / denom.imag();
			}
			Wn *= W;//se multiplica por la constante W
		}
		Wm *= W;
	}

	U[0][0].re=0.;
	U[0][0].im=0.;



	for(int x=0;x<C;x++){
		for(int y=0;y<C;y++) //se toman los resultados de la matriz U y se ponen en los vectores temporales Ur y Ui, los cuales se les aplicara la transformada inversa, para recuperar los valores de phi
		{
			Ur[(x*C)+y]= U[x][y].re;
			Ui[(x*C)+y]= U[x][y].im;
		 }
	}



	for(int i=1;i<C;i++){ //en este caso, el vector de entrada para la transformada es FF y la salida ff. HAY QUE HACERLO DESDE 1 PORQUE O SINO DA NAN.
	  	for(int j=0;j<C;j++){
	  		c_re (FF[j]) = Ur[(C*i)+j];
	  		c_im (FF[j]) = Ui[(C*i)+j];
	  	}

	  	fftw_plan q = fftw_create_plan (C, FFTW_BACKWARD, FFTW_ESTIMATE);
	  	fftw_one (q, FF, ff);//se calcula la transformada inversa en el eje x
	  	fftw_destroy_plan (q);
	  	for(int j=0;j<C;j++){
	  		U[i][j].re=ff[j].re; //se retornan los resultados a la matriz U
	  		U[i][j].im=ff[j].im;
	  	}
	}



	for(int i=0;i<C;i++){//el mismo prodecimiento anterior pero ahora en el eje y
		for(int j=0;j<C;j++){
			c_re (FF[j]) = U[j][i].re;
			c_im (FF[j]) = U[j][i].im;
		}
		fftw_plan q = fftw_create_plan (C, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_one (q, FF, ff);//se calcula la transformada
		fftw_destroy_plan (q);
		for(int j=0;j<C;j++){
			U[j][i].re=ff[j].re;
			U[j][i].im=ff[j].im;
		}
	}
	U[0][0].re =0.;


	for(int x=0;x<C;x++){
		for(int y=0;y<C;y++){
			u[(x*C)+y]=U[x][y].re/double(1.e+7);//en este caso, solo tomamos la parte real, que es la que nos interesa. obtenemos como resultado el potencial electrostatico.
		}
	}


	init.open("despues_inversa_y_salida");//se escribe un archivo de salida para analizar los datos. la salida corresponde al potencial electrostatico en cada celda conocido como phi.
	for (int i = 0; i < C; i++){
		for (int j = 0; j < C; j++){
			init<<u[(i*C)+j]<<" ";
		}
		init<<endl;
	}

	init.close();


}

// calculo del campo electrico el potencial de cada celda


void Electric2 (vector<double> phi, vector<double>& Ex, vector<double>& Ey) // recibe el potencial electroestatico calculado por la funcion poisson  y se calcula el campo electrico, tanto para X como para Y
{
  double dx = L / double (C); // el delta de x representa el tamano de la malla

  ofstream init; //archivo de control para analizar la entrada correcta de los datos. debe ser igual a "despues_inversa_y_salida"
  init.open("phi_en_campo.txt");
  for (int i = 0; i < C*C; i++)
	  init<<phi[i]<<endl;
  init.close();


  for (int x=0; x<C;x++){
	  Ex[x*C] = (phi[((x+1)*C)-1] - phi[(x*C)+1])/(2. * dx);// hallando el campo en x, en la primera columna
	  Ex[((x+1)*C)-1] = (phi[((x+1)*C)-2] - phi[(x*C)]) / (2. * dx);// hallando el campo en x, en la ultima columna
	  Ey[((C-1)*C)+x] = (phi[((C-2)*C)+x] - phi[x]) / (2. * dx); //hallando el campo en "y" para la ultima fila
	  Ey[x] = (phi[((C-1)*C)+x] - phi[x+C]) / (2. * dx);//hallando el campo para la primera fila y la ultima
  }

  for (int z=0;z<C;z++){//hallando el campo en x, el resto de las celdas en x
	  for (int w=1;w<C-1;w++){
		  Ex[w+(C*z)] = (phi[w-1] - phi[w+1]) / (2. * dx);
	  }
  }

  for(int i=C; i<(C*C)-C;i++){ //hallando el campo en "y" excepto para la primera y ultima fila
	  Ey[i] = (phi[i-C] - phi[i+C]) / (2. * dx);
  }

}


// se actualizan las nuevas posiciones y velocidades.
void Load (vector<double> r, vector<double> v, vector<double>& y)
{
  for (int i = 0; i < N; i++)
    {
      y[i] = r[i];
      y[N+i] = v[i];
    }
}

// Unload particle coordinates from solution vector

void UnLoad (vector<double> y, vector<double>& r, vector<double>& v)
{
  for (int i = 0; i < N; i++)
    {
      r[i] = y[i];
      v[i] = y[N+i];
    }
}



//arroja los resultados

void rk4_fixed2 (double &t, vector<double>& y1,vector<double>& y2, double h)
{
  // Y1 y Y2 son de tamano 2N, pues cada vector guarda posicion y velocidad en cada dimension
  // equations
  int n = y1.size();

  // arreglos para realizar el calculo de runge kutta
  vector<double> k1(n), k2(n), k3(n), k4(n), l1(n), l2(n), l3(n), l4(n),f1(n),f2(n), dydx1(n),dydx2(n),ny1(n),ny2(n);

  for(int i=0; i<n; i++){
	  dydx1[i]=.0;
	  dydx2[i]=.0;
  }

  itera+=1;//VARIABLE DE CONTROL PARA CONTAR LA CANTIDAD DE ITERACIONES Y DONDE PUEDEN EXISTIR POSIBLES CUELLOS DE BOTELLA EN MEMORIA


  for(int i=0; i<n; i++){
	  ny1[i]=y1[i];
	  ny2[i]=y2[i];
  }


  eval (t, ny1,ny2, dydx1, dydx2);


  for (int j = 0; j < n; j++)
    {
      k1[j] = h * dydx1[j];
      l1[j] = h * dydx2[j];
      f1[j] = ny1[j] + (k1[j] / 2.);
      f2[j] = ny2[j] + (l1[j] / 2.);
    }



  // primer paso intermedio
  itera+=1;
  eval(t + h / 2., f1,f2, dydx1,dydx2); //se evalua densidad, potencial y campo.



  //cout<<"iteraciones: "<<itera<<endl;
  for (int j = 0; j < n; j++)
    {
	  k2[j] = h * dydx1[j];
	  l2[j] = h * dydx2[j];
	  f1[j] = ny1[j] + k2[j] / 2.;
	  f2[j] = ny2[j] + l2[j] / 2.;
    }

  // segundo paso intermedio
  itera+=1;
  eval (t + h / 2., f1,f2, dydx1,dydx2);


  for (int j = 0; j < n; j++)
    {
	  k3[j] = h * dydx1[j];
	  l3[j] = h * dydx2[j];
	  f1[j] = ny1[j] + k3[j] / 2.;
	  f2[j] = ny2[j] + l3[j] / 2.;
    }


  itera+=1;
  // tercer paso intermedio
  eval (t + h, f1,f2, dydx1,dydx2);


  for (int j = 0; j < n; j++)
    {
      k4[j] = h * dydx1[j];
      l4[j] = h * dydx2[j];
    }

  // calculo del paso actual
  for (int j = 0; j < n; j++)
    {//calculo de runge-kutta.
      y1[j] += k1[j] / 6. + k2[j] / 3. + k3[j] / 3. + k4[j] / 6.;
      y2[j] += l1[j] / 6. + l2[j] / 3. + l3[j] / 3. + l4[j] / 6.;
    }
  t += h;

  return;
}



void Output2 (char* fn1, char* fn2, double t,
	     vector<double> rx, vector<double> ry ,vector<double> vx, vector<double> vy)
{

//itera el algoritmo para generar una salida en los archivos.

  ofstream phase;
  phase.open(fn1);
  //phase<<" pos x  "<<"  pos y "<<"  vel x  "<<"  vel y  "<<endl; FORMATO DE DATOS EN ARCHIVO
  for (int i = 0; i < N; i++)
	  phase<<rx[i]<<" "<<ry[i]<<" "<<vx[i]<<" "<<vy[i]<<endl;
  phase.close();



  //inicialización de parámetros
  vector<double> ne(C*C), n(C*C), phi(C*C), Ex(C*C), Ey(C*C);
  for (int i=0; i<n.size(); i++)
  {
	  ne[i]=0.;
	  n[i]=0.;
	  phi[i]=0.;
	  Ex[i]=0.;
	  Ey[i]=0.;

  }

  //cálculo de la densidad
  Density3 (rx,ry,ne);

  for (int j = 0; j < C*C; j++)
  {
	  n[j] = double (C*C) * ne[j] / double (N) - 1.; //hallando rho (ver esquema de normalizacion)
  }
  double kappa = 2. * M_PI / (L);


  Poisson_prueba1 (phi, n, kappa);//se resuelven las ecuaciones de poisson para el potencial electrostatico


  Electric2 (phi, Ex, Ey);//se calcula el campo electrico


  //se proceden a impromir los resultados de la iteración
  ofstream data;
  data.open(fn2);
  data<<" densidad: "<<" discretizada: "<<" campoX:   "<<" campoY:   "<<" phi:   "<<endl;//formato de la salida
  for (int j = 0; j < C*C; j++)
  {

	  data<<ne[j]<<" "<<n[j]<<" "<<Ex[j]<<" "<<Ey[j]<<" "<<phi[j]<<endl;
  }

  phase.close();

}

void eval (double t, vector<double> y1,vector<double> y2, vector<double>& dydt1,vector<double>& dydt2)

//esta función itera en cada uno de los 4 pasos de RK4.
//recibe la posicion y la velocidad en los vectores Y1 y Y2, extrae y los opera de la misma forma que la funciòn output
{

	vector<double> r1(N), v1(N), r2(N), v2(N), r1dot(N), v1dot(N),r2dot(N), v2dot(N), rx(N), ry(N); //inicializa los vectores.
	vector<double> ne(C*C), rho(C*C), phi(C*C), Ex(C*C),Ey(C*C); // ne el numero de electrones, rho promedio de densidades, phi otencial electroestatico , e campo electrico en cada celda.

	if(itera==1)
	          {
	        	  ofstream init;
	        	        init.open("eval.txt");
	        	        for (int i = 0; i < N; i++)
	        	        	init<<y1[i]<<" "<<y2[i]<<" "<<y1[N+i]<<" "<<y2[N+i]<<endl;
	        	        init.close();
	          }


	if(itera==2)
		          {
		        	  ofstream init;
		        	        init.open("eval2.txt");
		        	        for (int i = 0; i < N; i++)
		        	        	init<<y1[i]<<" "<<y2[i]<<" "<<y1[N+i]<<" "<<y2[N+i]<<endl;
		        	        init.close();
		          }
	UnLoad (y1, r1, v1);//se toman los valores de posicion y velocidad de la entrada para operarlos
	UnLoad (y2, r2, v2);




	  // se reinyenctan partículas que escapan del espacio de simulación
	  //limites de logitud por donde se mueven las particulas
	  rx = r1;
	  ry = r2;

	  for (int i = 0; i < N; i++)
	      {
	        if (rx[i] < 0.) rx[i] += L;
	        if (ry[i] < 0.) ry[i] += L;
	        if (rx[i] > L) rx[i] -= L;
	        if (ry[i] > L) ry[i] -= L;

	      }



	    //calculan la densidad con el numero de electrones que se encuentra dentro de la malla de estudio
	    //cout<<"antes de la densidad en eval"<<endl;

	    Density3 (rx,ry,ne); //se calcula la densidad de carga en la malla



	     double n0 = double (N) / L;
	     for (int j = 0; j < C*C; j++)
	       rho[j] = ne[j] / n0 - 1.;
	     double kappa = 2. * M_PI / L*L;
	     //calculo de poisson
	     Poisson_prueba1 (phi, rho, kappa);


	     //calculo del campo electrico
	     Electric2 (phi, Ex,Ey);


	     // Ecuaciones de movimiento
	     // al tener el potencial electroestatico se obtiene las ecuacioes de movimiento.
	     //se calcula el campo electrico en la posicion de la particula.
	     for (int i = 0; i < N; i++)
	       {
	         double dx = L / double (C);
	         int jx = int (rx[i] / dx);
	         int jy = int (ry[i] / dx);
	         double yx = rx[i] / dx - double (jx);
	         double yy = ry[i] / dx - double (jy);

	         double Efieldx;
	         double Efieldy;


	         if ((jx+1)%C == 0)
	            Efieldx = Ex[jx] * (1. - yx) + Ex[jx-(C-1)] * yx;
	         else
	            Efieldx = Ex[jx] * (1. - yx) + Ex[jx+1] * yx;
	         if ((jy+1)%C == 0)
	            Efieldy = Ey[jy] * (1. - yy) + Ey[jy-(C-1)] * yy;
	         else
	            Efieldy = Ey[jy] * (1. - yy) + Ey[jy+1] * yy;


	         // por el esquema de normalización:
	         //derivada de la posicion =velocidad
	         //derivada de la velocidad = campo electrico en la posicion de la particula.
	         r1dot[i] = v1[i];
	         v1dot[i] = - Efieldx;
	         r2dot[i] = v2[i];
	         v2dot[i] = - Efieldy;

	       }

	     // se vuelven a cargar los valores de la posicion y la velocidad para una nueva iteracion
	     Load (r1dot, v1dot, dydt1);
	     Load (r2dot, v2dot, dydt2);



	   }
