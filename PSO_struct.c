#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

//Ejemplo f(x,y,z) = 50 - (x-5)² - (y-5)² - (z-5)²

typedef struct PARTICULA {
	float *Xi;
	double *Vi;
	float *Pi;
	double Fxi;//Evaluación de Xi con fitness function
	double Fpi;//Evalucación de Pi con fitness function
}PARTICULA;

typedef struct ENJAMBRE {
	PARTICULA* Parts;
	int NumeroParticulas;//Configurable
	int NumeroIteraciones;
	unsigned char NumeroParametrosPorParticula;//Configurable
	float LimiteInferiorXi;//Puede ser un arreglo, uno para cada parametro
	float LimiteSuperiorXi;
	unsigned int Pgid;
	float LimVelInferior;
	float LimVelSuperior;
	char BuscarMaximo; //1->Maximo, 0->Mínimo
	float Weight, fW, minW;
	float c1, c2, chi;
}ENJAMBRE;

ENJAMBRE* CrearEnjambre(int NumeroParticulas, int NumeroIteraciones, unsigned char NumeroParametrosPorParticula, 
	float LimiteInferiorXi, float LimiteSuperiorXi, float LimVelInferior, float LimVelSuperior, float c1, float c2,
	float Weight, float fW, float minW, char BuscarMaximo, char UsarChi);
void InicializarEnjambre(ENJAMBRE *pEnj);
void EvaluarEnjambre(ENJAMBRE *pEnj);
void EvaluarEnjambreInicial(ENJAMBRE* pEnj);
void ActualizarMejorLocal(ENJAMBRE* pEnj, int index);
void CopiarXiaPi(ENJAMBRE* pEnj, int index);
void ActualizarMejorEnjambre(ENJAMBRE *pEnj);
void ActualizarVelocidad(ENJAMBRE *pEnj);
void ActualizarPosicion(ENJAMBRE *pEnj);
void ImprimeParticula(ENJAMBRE *pEnj, unsigned int index);
void ImprimeEnjambre(ENJAMBRE *pEnj);
double Funcion(ENJAMBRE *pEnj, int index);
void LiberarEnjambre(ENJAMBRE *pEnj);

//Funciones de prueba
float DosBinomiosCuadrados(ENJAMBRE* pEnj, int index);
float TresBinomiosCuadrados(ENJAMBRE* pEnj, int index);
float BinomiosCuadrados(ENJAMBRE* pEnj, int index);
double Esferica(ENJAMBRE* pEnj, int index);
float f2(ENJAMBRE* pEnj, int index);
float f3(ENJAMBRE* pEnj, int index);
float f4(ENJAMBRE* pEnj, int index);
float f6(ENJAMBRE* pEnj, int index);
float f8(ENJAMBRE* pEnj, int index);
float Rosenbrock(ENJAMBRE* pEnj, int index);
float Rastrigin(ENJAMBRE* pEnj, int index);
float Griewank(ENJAMBRE* pEnj, int index);

int main(int argc, char* argv[])
{
	unsigned int iteracion = 0;
	ENJAMBRE* Enjambre;

	srand(time(NULL));
	Enjambre = CrearEnjambre(/*NPart*/25, /*NIter*/50000, /*NParamPPart*/30, /*LimInf*/-100.0, /*LimSup*/100.0, 
							/*LimVelInf*/-2.5, /*LimVelSup*/2.5, /*c1*/2.05, /*c2*/2.05, /*Weight*/1.0,
							/*fW*/1.0, /*minWeight*/ 0.4, /*BuscarMaximo*/0, /*UsarChi*/1);
	InicializarEnjambre(Enjambre);
	EvaluarEnjambreInicial(Enjambre);
	//ActualizarMejorEnjambre(Enjambre);
	//ImprimeEnjambre(Enjambre);
	//getchar();

	while(iteracion < Enjambre->NumeroIteraciones && fabs(Enjambre->Parts[Enjambre->Pgid].Fpi) > 0.00000001)
	//while(fabs(50 - Enjambre[Pgid].Fpi) > 0.00001)
	{
		ActualizarVelocidad(Enjambre);
		ActualizarPosicion(Enjambre);
		EvaluarEnjambre(Enjambre);
		ActualizarMejorEnjambre(Enjambre);

		iteracion++;
		//_________________________-inercia
		if(Enjambre->Weight > Enjambre->minW)
			Enjambre->Weight *= Enjambre->fW;
		//ImprimeEnjambre(Enjambre);
		//getchar();
	}
	//Mostrar la mejor posición histórica apuntada por Pgid
	printf("Solucion (Iteracion): %d \n", iteracion-1);
	ImprimeParticula(Enjambre, Enjambre->Pgid);
	LiberarEnjambre(Enjambre);

	return 0;
}

double Funcion(ENJAMBRE *pEnj, int index)/////////////////Evaluación
{
	return Esferica(pEnj, index);
}

ENJAMBRE* CrearEnjambre(int NumeroParticulas, int NumeroIteraciones, unsigned char NumeroParametrosPorParticula, 
	float LimiteInferiorXi, float LimiteSuperiorXi, float LimVelInferior, float LimVelSuperior, float c1, float c2,
	float Weight, float fW, float minW, char BuscarMaximo, char UsarChi)
{
	float ro;

	ENJAMBRE* Enj;
	Enj = (ENJAMBRE*)malloc(sizeof(ENJAMBRE));
	Enj->Parts = (PARTICULA*)malloc(sizeof(PARTICULA)*NumeroParticulas);

	Enj->NumeroParticulas = NumeroParticulas;
	Enj->NumeroIteraciones = NumeroIteraciones;
	Enj->NumeroParametrosPorParticula = NumeroParametrosPorParticula;
	Enj->LimiteInferiorXi = LimiteInferiorXi;
	Enj->LimiteSuperiorXi = LimiteSuperiorXi;
	Enj->LimVelInferior = LimVelInferior;
	Enj->LimVelSuperior = LimVelSuperior;
	//Enj->LimVelInferior = (LimiteSuperiorXi - LimiteInferiorXi)*-0.15;//15% del rango
	//Enj->LimVelSuperior = (LimiteSuperiorXi - LimiteInferiorXi)*0.15;//15% del rango
	Enj->Pgid = 0;
	Enj->BuscarMaximo = BuscarMaximo;
	Enj->Weight = Weight;
	Enj->fW = fW;
	Enj->minW = minW;
	Enj->c1 = c1;
	Enj->c2 = c2;

	ro = c1 + c2;
	//______________________________-Calcular chi
	if(UsarChi)
		Enj->chi = 2.0 / (fabs(2.0 - ro - sqrtf(ro*ro - 4.0*ro)));
	else
		Enj->chi = 1.0;
	return Enj;
}


void InicializarEnjambre(ENJAMBRE *pEnj)
{
	unsigned int i, d;
	for (i = 0; i < pEnj->NumeroParticulas; i++)
	{
		pEnj->Parts[i].Xi = (float *)malloc(sizeof(float)*pEnj->NumeroParametrosPorParticula); 
		pEnj->Parts[i].Vi = (double *)malloc(sizeof(double)*pEnj->NumeroParametrosPorParticula);
		pEnj->Parts[i].Pi = (float *)malloc(sizeof(float)*pEnj->NumeroParametrosPorParticula);
	}

	for (i = 0; i < pEnj->NumeroParticulas; i++)
	{
		for(d = 0; d < pEnj->NumeroParametrosPorParticula; d++)
		{
			pEnj->Parts[i].Xi[d] = pEnj->LimiteInferiorXi + (double)rand()/RAND_MAX*(pEnj->LimiteSuperiorXi-pEnj->LimiteInferiorXi);//Generar un número aleatorio en el rango del espacio de busqueda
			pEnj->Parts[i].Vi[d] = 0.0;
			pEnj->Parts[i].Pi[d] = pEnj->Parts[i].Xi[d];
		}
	}
}

void EvaluarEnjambreInicial(ENJAMBRE* pEnj)
{
	unsigned int i, d;

	for(i = 0; i < pEnj->NumeroParticulas; i++)
	{
		pEnj->Parts[i].Fxi = Funcion(pEnj, i);

		pEnj->Parts[i].Fpi = pEnj->Parts[i].Fxi;
		for (d = 0; d < pEnj->NumeroParametrosPorParticula; d++)
		{
			pEnj->Parts[i].Pi[d] = pEnj->Parts[i].Xi[d];
		}

	}
}

void EvaluarEnjambre(ENJAMBRE *pEnj)
{
	unsigned int i, d;

	for(i = 0; i < pEnj->NumeroParticulas; i++)
	{
		pEnj->Parts[i].Fxi = Funcion(pEnj, i);
		ActualizarMejorLocal(pEnj, i);
	}
}

void ActualizarMejorLocal(ENJAMBRE* pEnj, int index)
{
	if(pEnj->BuscarMaximo == 1)
	{
		if(pEnj->Parts[index].Fxi > pEnj->Parts[index].Fpi)
		{
			CopiarXiaPi(pEnj, index);
		}
	}
	else
	{
		if(pEnj->Parts[index].Fxi < pEnj->Parts[index].Fpi)
		{
			CopiarXiaPi(pEnj, index);
		}
	}
}

void CopiarXiaPi(ENJAMBRE* pEnj, int index)
{
	unsigned int d;

	pEnj->Parts[index].Fpi = pEnj->Parts[index].Fxi;
	for (d = 0; d < pEnj->NumeroParametrosPorParticula; d++)
	{
		pEnj->Parts[index].Pi[d] = pEnj->Parts[index].Xi[d];
	}
}

void ActualizarMejorEnjambre(ENJAMBRE *pEnj)
{
	unsigned int i;

	for(i = 0; i < pEnj->NumeroParticulas; i++)
	{
		if(pEnj->BuscarMaximo == 1)
		{
			if(pEnj->Parts[i].Fpi > pEnj->Parts[pEnj->Pgid].Fpi)
			{
				pEnj->Pgid = i;
			}
		}
		else
		{
			if(pEnj->Parts[i].Fpi < pEnj->Parts[pEnj->Pgid].Fpi)
			{
				pEnj->Pgid = i;
			}
		}
	}
}

void ActualizarVelocidad(ENJAMBRE *pEnj)
{
	double psi1, psi2;
	unsigned int i, d;
	double aux;
	//________________________________________________Obtener mejor posición
	/*static float* Pg = NULL;
	if(Pg == NULL)
	{
		Pg = (float*)malloc(sizeof(float)*pEnj->NumeroParametrosPorParticula);
	}
	for (d = 0; d < pEnj->NumeroParametrosPorParticula; d++)
	{
		Pg[d] = pEnj->Parts[pEnj->Pgid].Xi[d];
	}*/

	//_____________________________________________actualizar velocidad
	for (i = 0; i < pEnj->NumeroParticulas; i++)
	{
		for (d = 0; d < pEnj->NumeroParametrosPorParticula; d++)
		{
			psi1 = (double)rand()/RAND_MAX;
			psi2 = (double)rand()/RAND_MAX;

			aux = pEnj->chi * (pEnj->Weight * pEnj->Parts[i].Vi[d] + pEnj->c1*psi1*(pEnj->Parts[i].Pi[d] - pEnj->Parts[i].Xi[d]) + pEnj->c2*psi2*(pEnj->Parts[pEnj->Pgid].Pi[d] - pEnj->Parts[i].Xi[d]));	
			//printf("\n\n%lf, %u, %u\n\n", aux, i, d);
			if(aux > pEnj->LimVelSuperior){
				pEnj->Parts[i].Vi[d] = pEnj->LimVelSuperior;
				continue;
			}else
			if(aux < pEnj->LimVelInferior){
				pEnj->Parts[i].Vi[d] = pEnj->LimVelInferior;
				continue;
			}
			pEnj->Parts[i].Vi[d] = aux;
			//printf("\n%f\n", pEnj->Parts[i].Vi[d]);
		}
	}
}

void ActualizarPosicion(ENJAMBRE *pEnj)
{
	unsigned int i, d;

	for (i = 0; i < pEnj->NumeroParticulas; i++)
	{
		for (d = 0; d < pEnj->NumeroParametrosPorParticula; d++)
		{
			pEnj->Parts[i].Xi[d] += pEnj->Parts[i].Vi[d];
		}
	}
}


void ImprimeParticula(ENJAMBRE *pEnj, unsigned int index)
{
	unsigned int d;

	printf("\nX%u: ", index);
	for (d = 0; d < pEnj->NumeroParametrosPorParticula; d++)
	{
		printf("%f, ", pEnj->Parts[index].Xi[d]);
	}
	printf("\nV%u: ", index);
	for (d = 0; d < pEnj->NumeroParametrosPorParticula; d++)
	{
		printf("%f, ", pEnj->Parts[index].Vi[d]);
	}
	printf("\nP%u: ", index);
	for (d = 0; d < pEnj->NumeroParametrosPorParticula; d++)
	{
		printf("%f, ", pEnj->Parts[index].Pi[d]);
	}
	printf("\nX%uFit: %f", index, pEnj->Parts[index].Fxi);
	printf("\nP%uFit: %f", index, pEnj->Parts[index].Fpi);

	printf("\n");

}


void ImprimeEnjambre(ENJAMBRE *pEnj)
{
	unsigned int i;

	printf("ENJAMBRE: \n");
	for(i = 0; i < pEnj->NumeroParticulas; i++)
	{
		ImprimeParticula(pEnj, i);
	}
	printf("Mejor: %d\n", pEnj->Pgid);
}

void LiberarEnjambre(ENJAMBRE *pEnj)
{
	unsigned int i;

	for(i=0; i<pEnj->NumeroParticulas; i++)
	{
		free(pEnj->Parts[i].Xi);
		free(pEnj->Parts[i].Vi);
		free(pEnj->Parts[i].Pi);
	}

	free(pEnj->Parts);
	free(pEnj);
}


//______________________________________________________Funciones de prueba
float DosBinomiosCuadrados(ENJAMBRE* pEnj, int index)
{
	return (50.0 - ((pEnj->Parts[index].Xi[0] - 5)*(pEnj->Parts[index].Xi[0] - 5) + (pEnj->Parts[index].Xi[1] - 5)*(pEnj->Parts[index].Xi[1] - 5)));
}

float TresBinomiosCuadrados(ENJAMBRE* pEnj, int index)
{
	return (50.0 - ((pEnj->Parts[index].Xi[0] - 5)*(pEnj->Parts[index].Xi[0] - 5) + (pEnj->Parts[index].Xi[1] - 5)*(pEnj->Parts[index].Xi[1] - 5) + (pEnj->Parts[index].Xi[2] - 5)*(pEnj->Parts[index].Xi[2] - 5)));
}

float BinomiosCuadrados(ENJAMBRE* pEnj, int index)
{
	unsigned int i;
	float suma = 0.0;

	for(i = 0; i < pEnj->NumeroParametrosPorParticula; i++)
	{
		suma += (pEnj->Parts[index].Xi[i] - 5.0)*(pEnj->Parts[index].Xi[i] - 5.0);
	}

	return suma;
}

double Esferica(ENJAMBRE* pEnj, int index)
{
	unsigned int i;
	double suma = 0.0;

	for(i = 0; i < pEnj->NumeroParametrosPorParticula; i++)
	{
		suma += (double)pEnj->Parts[index].Xi[i]*pEnj->Parts[index].Xi[i];
	}

	return suma;
}

float f2(ENJAMBRE* pEnj, int index)
{
	unsigned int i;
	float sum = 0.0, mult = 1.0;

	for (i = 0; i < pEnj->NumeroParametrosPorParticula; i++)
	{
		sum += fabs(pEnj->Parts[index].Xi[i]);
		mult *= pEnj->Parts[index].Xi[i];
	}

	return sum + mult;
}

float f3(ENJAMBRE* pEnj, int index)
{
	unsigned int i, j;
	float sum = 0.0, aux = 0.0;

	for (i = 0; i < pEnj->NumeroParametrosPorParticula; i++)
	{
		aux = 0.0;
		for (j = 0; j < i; j++)
		{
			aux += pEnj->Parts[index].Xi[j];
		}
		aux *= aux;
		sum += aux;
	}

	return sum;
}

float f4(ENJAMBRE* pEnj, int index)
{
	unsigned int i;
	float max=fabs(pEnj->Parts[index].Xi[0]);

	for (i = 1; i < pEnj->NumeroParametrosPorParticula; i++)
	{
		if(fabs(pEnj->Parts[index].Xi[i]) > max)
		{
			max = fabs(pEnj->Parts[index].Xi[i]);
		}
	}

	return max;
}


float Rosenbrock(ENJAMBRE* pEnj, int index)
{
	unsigned int i;
	float sum = 0.0;

	for (i = 0; i < pEnj->NumeroParametrosPorParticula - 1; i++)
	{
		sum += 100 * ((pEnj->Parts[index].Xi[i+1] - pEnj->Parts[index].Xi[i] * pEnj->Parts[index].Xi[i])
			* (pEnj->Parts[index].Xi[i+1] - pEnj->Parts[index].Xi[i] * pEnj->Parts[index].Xi[i]))
			+ (pEnj->Parts[index].Xi[i] - 1) * (pEnj->Parts[index].Xi[i] - 1);
	}

	return sum;
}

float Rastrigin(ENJAMBRE* pEnj, int index)
{
	unsigned int i;
	float sum = 0.0;

	for (i = 0; i < pEnj->NumeroParametrosPorParticula; i++)
	{
		sum += pEnj->Parts[index].Xi[i] * pEnj->Parts[index].Xi[i] - 10.0*cos(2.0*3.14159265*pEnj->Parts[index].Xi[i]) + 10;
	}

	return sum;
}

float Griewank(ENJAMBRE* pEnj, int index)
{
	unsigned int i;
	float sum1 = 0.0;
	float sum2 = 0.0;

	for (i = 0; i < pEnj->NumeroParametrosPorParticula; i++)
	{
		sum1 += pEnj->Parts[index].Xi[i] * pEnj->Parts[index].Xi[i];
		sum2 += cos(pEnj->Parts[index].Xi[i]/sqrtf(i+1));
	}

	return (float)1/4000 * sum1 + sum2 + 1;
}

float f6(ENJAMBRE* pEnj, int index)
{
	unsigned int i;
	float sum = 0.0;

	for (i = 0; i < pEnj->NumeroParametrosPorParticula; i++)
	{
		sum += fabs(pEnj->Parts[index].Xi[i] + 0.5) * fabs(pEnj->Parts[index].Xi[i] + 0.5);
	}

	return sum;
}

float f8(ENJAMBRE* pEnj, int index)
{
	unsigned int i;
	float sum = 0.0;

	for (i = 0; i < pEnj->NumeroParametrosPorParticula; i++)
	{
		sum -= pEnj->Parts[index].Xi[i] * sin(sqrtf(fabs(pEnj->Parts[index].Xi[i])));
	}

	return sum;
}