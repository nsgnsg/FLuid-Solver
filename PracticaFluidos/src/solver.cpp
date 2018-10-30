#include "solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void Solver::Init(unsigned N, float dt, float diff, float visc)
{
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;
	this->N = N;
	
}
/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
*/
void Solver::FreeData(void)
{
	if (u_prev!=NULL)
		free(u_prev);
	if (v_prev != NULL)
		free(v_prev);
	if (dens_prev != NULL)
		free(dens_prev);
	if (u != NULL)
		free(u);
	if (v != NULL)
		free(v);
	if (dens != NULL)
		free(dens);
	//free(buffer)
//TODO: Libera los buffers de memoria.
}

void Solver::ClearData(void)
{
	for (int i = 0; i < (N+2)*(N+2); ++i) {
		u_prev[i] = 0;
		v_prev[i] = 0;
		dens_prev[i] = 0;
		u[i] = 0;
		v[i] = 0;
		dens[i] = 0;
		gauss = true;
	}
	
//TODO: Borra todo el contenido de los buffers
}

bool Solver::AllocateData(void)
{
//TODO:
//Reservamos memoria, en caso de fallo devlvemos false.
//Antes de devolver true, hay que limpiar la memoria reservada con un ClearData().

	 
	 u_prev = (float*)malloc((N+2) * (N+2) * sizeof(float));
	 v_prev = (float*)malloc((N + 2) * (N + 2) * sizeof(float));
	 dens_prev = (float*)malloc((N + 2) * (N + 2) * sizeof(float));
	 u = (float*)malloc((N + 2) * (N + 2) * sizeof(float));
	 v = (float*)malloc((N + 2) * (N + 2) * sizeof(float));
	 dens = (float*)malloc((N + 2) * (N + 2) * sizeof(float));
	 if (u_prev == NULL || v_prev == NULL || dens_prev == NULL || u == NULL || v == NULL || dens == NULL) {
		 return false;

	 }

	 ClearData();
	 return true;
}

void Solver::ClearPrevData() 
{
	for (int i = 0; i < (N + 2)*(N + 2); ++i) {
		u_prev[i] = 0;
		v_prev[i] = 0;
		dens_prev[i] = 0;
		
	}
//TODO: Borra el contenido de los buffers _prev
}

void Solver::AddDensity(unsigned x, unsigned y, float source)
{
	dens_prev[XY_TO_ARRAY(x, y)] += source;
//TODO: Añade el valor de source al array de densidades. Sería interesante usar la macro: XY_TO_ARRAY

}

void Solver::AddVelocity(unsigned x, unsigned y, float forceX, float forceY)
{
	u_prev[XY_TO_ARRAY(x, y)] += forceX;
	v_prev[XY_TO_ARRAY(x, y)] += forceY;
//TODO: Añade el valor de fuerza a sus respectivos arrays. Sería interesante usar la macro: XY_TO_ARRAY
}

void Solver::Solve()
{
	VelStep();
	DensStep();
}

//Estos metodos se usan para cambair el metofo iterativo a usar para resolver el sistema
void Solver::changeMethodGauss()
{
	gauss = true;
}

void Solver::changeMethodJacobi()
{
	gauss = false;
}

void Solver::DensStep()
{
	AddSource(dens, dens_prev);			//Adding input density (dens_prev) to final density (dens).
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Diffuse(0, dens, dens_prev);		//Writing result in dens because we made the swap before. bi = dens_prev. The initial trash in dens matrix, doesnt matter, because it converges anyways.
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Advect(0, dens, dens_prev, u, v);	//Advect phase, result in dens.
}

void Solver::VelStep()
{
	AddSource(u, u_prev);
	AddSource(v, v_prev);
	SWAP (u_prev,u)			
	SWAP (v_prev, v)
	Diffuse(1, u, u_prev);  
	Diffuse(2, v, v_prev); 
	Project(u, v, u_prev, v_prev);		//Mass conserving.
	
	SWAP (u_prev,u)			
	SWAP (v_prev,v)
	Advect(1, u, u_prev, u_prev, v_prev);
	Advect(2, v, v_prev, u_prev, v_prev);
	Project(u, v, u_prev, v_prev);		//Mass conserving.
}

void Solver::AddSource(float * base, float * source)
{
	
	int i, j;
	FOR_EACH_CELL
		base[XY_TO_ARRAY(i, j)] += source[XY_TO_ARRAY(i, j)] * dt;
	END_FOR
	
	/*for (int i = 0; i < ((N+2) * (N +2)); ++i) {
		base[i] = source[i] * dt;
	}*/
//TODO: Teniendo en cuenta dt (Delta Time), incrementar el array base con nuestro source. Esto sirve tanto para añadir las nuevas densidades como las nuevas fuerzas.
}


void Solver::SetBounds(int b, float * x)
{
	
	if (b == 0) {
		for (int i = 1; i < N + 1; ++i) {
			x[XY_TO_ARRAY(0, i)] = x[1, i];

			x[XY_TO_ARRAY(i, 0)] = x[i, 1];

			x[XY_TO_ARRAY(i, N + 1)] = x[i, N];

			x[XY_TO_ARRAY(N + 1, i)] = x[N, i];
		}
	}

	if (b == 1) {
		for (int i = 1; i < N + 1; ++i) {
			x[XY_TO_ARRAY(0, i)] = -x[1, i];
		
			x[XY_TO_ARRAY(i, 0)] = x[i, 1];
		
			x[XY_TO_ARRAY(i, N + 1)] = x[i, N];
		
			x[XY_TO_ARRAY(N + 1, i)] = -x[N, i];
		}
	}
	if (b == 2) {
		for (int i = 1; i < N + 1; ++i) {
			x[XY_TO_ARRAY(0, i)] = x[1, i];

			x[XY_TO_ARRAY(i, 0)] = -x[i, 1];

			x[XY_TO_ARRAY(i, N + 1)] = -x[i, N];

			x[XY_TO_ARRAY(N + 1, i)] = x[N, i];
		}
	}

	x[XY_TO_ARRAY(0, 0)] = 0;

	x[XY_TO_ARRAY(0, N+1)] = 0;

	x[XY_TO_ARRAY(N+1, N + 1)] = 0;

	x[XY_TO_ARRAY(N + 1, 0)] = 0;
/*TODO:
Input b: 0, 1 or 2.
	0: borders = same value than the inner value.
	1: x axis borders inverted, y axis equal.
	2: y axis borders inverted, x axis equal.
	Corner values allways are mean value between associated edges.
*/
}

/*
https://www.youtube.com/watch?v=62_RUX_hrT4
https://es.wikipedia.org/wiki/M%C3%A9todo_de_Gauss-Seidel <- Solución de valores independientes.
Despreciando posibles valores de x no contiguos, se simplifica mucho. Mirar diapositivas y la solución de Gauss Seidel de términos independientes.
Gauss Seidel -> Matrix x and x0
*/
void Solver::LinSolve(int b, float * x, float * x0, float aij, float aii)
{
	int i, j;
	int maxIteraciones = 25;
	for (int k = 0; k < maxIteraciones; ++k) {
		FOR_EACH_CELL
			x[XY_TO_ARRAY(i, j)] = (x0[XY_TO_ARRAY(i, j)] - aij * (x[XY_TO_ARRAY(i + 1, j)] + x[XY_TO_ARRAY(i - 1, j)] +
				x[XY_TO_ARRAY(i, j + 1)] + x[XY_TO_ARRAY(i, j - 1)])) / (1 + aii);
				
		END_FOR
		SetBounds(b, x);
	}
//TODO: Se recomienda usar FOR_EACH_CELL, END_FOR y XY_TO_ARRAY.
}


//Metodo hecho por Jacobi
void Solver::LinSolve2(int b, float * x, float * x0, float aij, float aii)
{
	int i, j;
	int maxIteraciones = 25;
	for (int k = 0; k < maxIteraciones; ++k) {
		FOR_EACH_CELL
			//x[XY_TO_ARRAY(i, j)] = (x0[XY_TO_ARRAY(i, j)] - aij * (x[XY_TO_ARRAY(i + 1, j)] + x[XY_TO_ARRAY(i - 1, j)] +
				//x[XY_TO_ARRAY(i, j + 1)] + x[XY_TO_ARRAY(i, j - 1)])) / (1 + aii);
			x[XY_TO_ARRAY(i, j)] = (1/(1 + aii)) *( x0[XY_TO_ARRAY(i, j)] - aij * (x[XY_TO_ARRAY(i + 1, j)] + x[XY_TO_ARRAY(i - 1, j)] +
				x[XY_TO_ARRAY(i, j + 1)] + x[XY_TO_ARRAY(i, j - 1)]));
		END_FOR
			SetBounds(b, x);
	}
	//TODO: Se recomienda usar FOR_EACH_CELL, END_FOR y XY_TO_ARRAY.
}


/*
Nuestra función de difusión solo debe resolver el sistema de ecuaciones simplificado a las celdas contiguas de la casilla que queremos resolver,
por lo que solo con la entrada de dos valores, debemos poder obtener el resultado.
*/
void Solver::Diffuse(int b, float * x, float * x0)
{
	
	
	float a = dt*diff*N*N;
	float aii = 4 * a;
	//float aii2 = 1 + 4 * a;
	float aij = -a;
	
	
	
	//LinSolve(b,x,x0,aij,aii);
	if (gauss) {
		LinSolve(b, x, x0, aij, aii);
		
	}
	else {
		LinSolve2(b, x, x0, aij, aii);
		
		
	}

	
	

//TODO: Solo necesitaremos pasar dos parámetros a nuestro resolutor de sistemas de ecuaciones de Gauss Seidel. Calculamos dichos valores y llamamos a la resolución del sistema.
}

/*
d is overwrited with the initial d0 data and affected by the u & v vectorfield.
Hay que tener en cuenta que el centro de las casillas representa la posición entera dentro de la casilla, por lo que los bordes estan
en las posiciones x,5.
*/
void Solver::Advect(int b, float * d, float * d0, float * u, float * v)
{
	
	
	
	int i, j;
	FOR_EACH_CELL
		float x = i - dt * N * u[XY_TO_ARRAY(i, j)];
		float y = j - dt * N * v[XY_TO_ARRAY(i, j)];
		//Evitar problemas con el margen
		if (x < 0.5)
			x = 0.5;
		else if (x > N + 0.5)
			x = N + 0.5;
		if (y < 0.5)
			y = 0.5;
		else if (y > N + 0.5)
			y = N + 0.5;
		//Coger la parte entera de x e y
		int cx1 = (int)x;
		int cy1 = (int)y;
		//Conseguimos las cuatro casillas mas cercanas que ademas de xInt e yINt seran sus casillas siguientes
		int cx2 = cx1 + 1;
		int cy2 = cy1 + 1;

		//Coger la parte decimal de x e y
		float fx2 = x - cx1;
		float fx1 = 1 - fx2;
		float fy2 = y - cy1;
		float fy1 = 1 - fy2;

		//Calcular la formula
		d[XY_TO_ARRAY(i, j)] = fx1 * (fy1 * d0[XY_TO_ARRAY(cx1, cy1)] + fy2 * d0[XY_TO_ARRAY(cx1, cy2)]) +
			fx2 * (fy1 * d0[XY_TO_ARRAY(cx2, cy1)] + fy2 * d0[XY_TO_ARRAY(cx2, cy2)]);
			
		
	END_FOR
	SetBounds(b, d);
		
		

//TODO: Se aplica el campo vectorial realizando una interploación lineal entre las 4 casillas más cercanas donde caiga el nuevo valor.
	
}

/*
Se encarga de estabilizar el fluido y hacerlo conservativo de masa. Se usa solo en las matrices de velocidades.
No necesaria implementación por su complejidad.
*/
void Solver::Project(float * u, float * v, float * p, float * div)
{
	int i, j;

	FOR_EACH_CELL
		div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
		p[XY_TO_ARRAY(i, j)] = 0;
	END_FOR
	SetBounds(0, div);
	SetBounds(0, p);

	if (gauss)
		LinSolve(0, p, div, 1, 4);
	else {
		LinSolve2(0, p, div, 1, 4);
		
	}


	//LinSolve(0, p, div, 1, 4);

	//Aproximamos: Laplaciano de q a su gradiente.
	FOR_EACH_CELL
		u[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
		v[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
	END_FOR
	SetBounds(1, u);
	SetBounds(2, v);
}


