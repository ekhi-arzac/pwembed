/***************************************************
 AC - OpenMP -- SERIE
 fun_s.c
 rutinas que se utilizan en el modulo grupopal_s.c
****************************************************/
#include "defineg.h" // definiciones
#include <float.h>
#include <math.h>
#include <stdlib.h>

// double inline pow2(double n)
// {
//   return n * n;
// }

/*******************************************************************
 1 - Funcion para calcular la distancia euclidea entre dos vectores
 Entrada: 2 elementos con NDIM caracteristicas (por referencia)
 Salida:  distancia (double)
********************************************************************/
double gendist(float *vec1, float *vec2)
{
  double distancia = 0.0;
  int i;
  for (i = 0; i < NDIM; i++)
  {
    distancia += (vec2[i] - vec1[i])*(vec2[i] - vec1[i]);
  }
  return sqrt(distancia);
}

/***********************************************************************************
 2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas
cercano) Entrada: nvec  numero de vectores, int mvec  vectores, una matriz de
tamanno MAXV x NDIM, por referencia cent  centroides, una matriz de tamanno
ngrupos x NDIM, por referencia Salida:  popul grupo mas cercano a cada elemento,
vector de tamanno MAXV, por ref.
************************************************************************************/
void grupo_cercano(int nvec, float mvec[][NDIM], float cent[][NDIM],
                   int *popul)
{
  int i, j;
  #pragma omp parallel for private(i, j)
  for (i = 0; i < nvec; i++)
  {
    double min_dist = DBL_MAX; // Infinito
    int min_index = -1;
    for (j = 0; j < ngrupos; j++)
    {
      double distancia = gendist(mvec[i], cent[j]);
      if (distancia < min_dist)
      {
        min_dist = distancia;
        min_index = j;
      }
      
    }
    popul[i] = min_index;
  }
  
}

/***************************************************************************************
 3 - Funcion para calcular la calidad de la particion de clusteres.
     Ratio entre a y b.
     El termino a corresponde a la distancia intra-cluster.
     El termino b corresponde a la distancia inter-cluster.
 Entrada: mvec    vectores, una matriz de tamanno MAXV x NDIM, por referencia
          listag  vector de ngrupos structs (informacion de grupos generados),
por ref. cent    centroides, una matriz de tamanno ngrupos x NDIM, por
referencia Salida:  valor del CVI (double): calidad/bondad de la particion de
clusters
****************************************************************************************/
double silhouette_simple(float mvec[][NDIM], struct lista_grupos *listag,
                         float cent[][NDIM], float a[])
{

  int i, j, k;
  double distancia;
  float b[ngrupos];

  // PARA COMPLETAR

  // aproximar a[i] de cada cluster: calcular la densidad de los grupos;
  //		media de las distancia entre todos los elementos del grupo;
  //   	si el numero de elementos del grupo es 0 o 1, densidad = 0
  // ...
  int nvecg;
  #pragma omp parallel for private(i, j, k, distancia, nvecg) shared(a, mvec, listag)
  for (i = 0; i < ngrupos; i++)
  {
    nvecg = listag[i].nvecg;

    if (nvecg <= 1)
    {
      a[i] = 0;
    }
    else
    {
      distancia = 0;
      for (j = 0; j < nvecg; j++)
      {
        for (k = j + 1; k < nvecg; k++)
        {
          distancia +=
              gendist(mvec[listag[i].vecg[k]], mvec[listag[i].vecg[j]]);
        }
      }
      a[i] = distancia / (nvecg * (nvecg - 1) / 2);
    }
  }
  // aproximar b[i] de cada cluster
  // ...
  #pragma omp parallel for private(i, j, distancia) shared(b, cent)
  for (i = 0; i < ngrupos; i++)
  {
    distancia = 0;
    for (j = 0; j < ngrupos; j++)
    {
      if (i != j)
        distancia += gendist(cent[i], cent[j]);
    }
    b[i] = distancia / (ngrupos - 1);
  }
  // calcular el ratio s[i] de cada cluster
  // ...
  double s = 0;
  #pragma omp parallel for reduction(+ : s)
  for (i = 0; i < ngrupos; i++)
  {
    s += (b[i] - a[i]) / max(a[i], b[i]);
  }
  // promedio y devolver
  // ...
  s /= ngrupos;
  return s;
}

void swap(float *arr, int i, int j)
{
  float temp = arr[i];
  arr[i] = arr[j];
  arr[j] = temp;
}

// A function to implement bubble sort
void bubbleSort(float arr[], int n)
{
  int i, j;
  for (i = 0; i < n - 1; i++)
    for (j = 0; j < n - i - 1; j++)
      if (arr[j] > arr[j + 1])
        swap(arr, j, j + 1);
}

/********************************************************************************************
 4 - Funcion para relizar el analisis de campos UNESCO
 Entrada:  listag   vector de ngrupos structs (informacion de grupos generados),
por ref. mcam     campos, una matriz de tamaño MAXV x NCAM, por referencia
 Salida:   info_cam vector de NCAM structs (informacion del analisis realizado),
por ref.
*****************************************************************************************/
void analisis_campos(struct lista_grupos *listag, float mcam[][NCAM],
                     struct analisis *info_cam)
{
  // PARA COMPLETAR
  // Realizar el analisis de campos UNESCO en los grupos:
  //    mediana maxima y el grupo en el que se da este maximo (para cada campo)
  //    mediana minima y su grupo en el que se da este minimo (para cada campo)

  int i, j, k, grp_max, grp_min, dim;
  float mediana_max, mediana_min;
  #pragma omp parallel for private(i, j, k, grp_max, grp_min, dim, mediana_max, mediana_min)
  for (i = 0; i < NCAM; i++)
  {

    mediana_max = -1;
    mediana_min = FLT_MAX;
    grp_max = -1;
    grp_min = -1;

    for (j = 0; j < ngrupos; j++)
    {

      dim = listag[j].nvecg;
      if (dim != 0)
      {
        float campo_grp[dim];

        for (k = 0; k < dim; k++)
        {
          campo_grp[k] = mcam[listag[j].vecg[k]][i];
        }

        bubbleSort(campo_grp, dim);
        float mediana = campo_grp[dim / 2];

        if (mediana > mediana_max)
        {
          grp_max = j;
          mediana_max = mediana;
        }

        if (mediana < mediana_min)
        {
          grp_min = j;
          mediana_min = mediana;
        }
      }
    }
    info_cam[i].gmax = grp_max;
    info_cam[i].gmin = grp_min;
    info_cam[i].mmax = mediana_max;
    info_cam[i].mmin = mediana_min;
  }
}

/*************************************
   OTRAS FUNCIONES DE LA APLICACION
**************************************/
void inicializar_centroides(float cent[][NDIM])
{
  int i, j;
  float rand_val;
  srand(147);
  for (i = 0; i < ngrupos; i++)
    for (j = 0; j < NDIM / 2; j++)
    {
      rand_val = ((rand() % 10000) / 10000.0) * 2 - 1;
      cent[i][j] = rand_val;
      cent[i][j + (NDIM / 2)] = cent[i][j];
    }
}

int nuevos_centroides(float mvec[][NDIM], float cent[][NDIM], int popul[],
                      int nvec)
{
  int i, j, fin;
  double discent;
  double additions[ngrupos][NDIM + 1];
  float newcent[ngrupos][NDIM];


  for (i = 0; i < ngrupos; i++)
    for (j = 0; j < NDIM + 1; j++)
      additions[i][j] = 0.0;

  // acumular los valores de cada caracteristica; numero de elementos al final
  for (i = 0; i < nvec; i++)
  {
    for (j = 0; j < NDIM; j++)
      additions[popul[i]][j] += mvec[i][j];
    additions[popul[i]][NDIM]++;
  }

  // calcular los nuevos centroides y decidir si el proceso ha finalizado o no
  // (en funcion de DELTA)
  fin = 1;
  for (i = 0; i < ngrupos; i++)
  {
    if (additions[i][NDIM] > 0)
    { // ese grupo (cluster) no esta vacio
      // media de cada caracteristica
      for (j = 0; j < NDIM; j++)
        newcent[i][j] = (float)(additions[i][j] / additions[i][NDIM]);

      // decidir si el proceso ha finalizado
      discent = gendist(&newcent[i][0], &cent[i][0]);
      if (discent > DELTA1)
      {
        fin = 0; // en alguna centroide hay cambios; continuar
      }

      // copiar los nuevos centroides
      for (j = 0; j < NDIM; j++)
        cent[i][j] = newcent[i][j];
    }
  }
  return fin;
}
