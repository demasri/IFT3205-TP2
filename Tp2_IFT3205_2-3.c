/*-------------------------------------------------------------------*/
/* Prog    : Tp2_IFT3205-2-3.c                                       */
/* Auteur  : Daniel El-Masri (20096261) daniel.el-masri@umontreal.ca */
/* Date    : 19 /02/2010                                              */
/* version :                                                         */ 
/* langage : C                                                       */
/* labo    : DIRO                                                    */
/*-------------------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo2.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/   
/*------------------------------------------------*/
#define NAME_VISUALISER "display "
#define NAME_IMG_IN1  "UdM_1"
#define NAME_IMG_IN2 "UdM_2"

/*------------------------------------------------*/
/* PROTOTYPE DE FONCTIONS  -----------------------*/   
/*------------------------------------------------*/

float RotatedMatrixElement(float** m, int j, int i, int length, int width) {

    if(i < 0 || i >= length)
        i = (i + length) % length;
    
    if(j < 0 || j >= width)
        j = (j + width) % width;

    return m[i][j];
}


void RotationBilinearInterpolation(float** MatriceImg, int length, int width, float** MatriceImgPivotee, float theta)
{
  int x0,y0, x1, y1, x2, y2;

  float f_xy, f_x1y, f_xy1, f_x0y1, f_x0y, f_x1y1;

  int half_length = length / 2;
  int half_width = width / 2;
  
  for(y2=0; y2 < length; y2++)
  {
    for(x2=0;x2 < width; x2++)
    {
      // Calcul de la position du pixel dans l'image initiale
      x0 =  (cos(-theta) * (x2 - half_width)) + (sin(-theta) * (y2 - half_length)) + half_width;
      y0 =  -(sin(-theta) * (x2 - half_width)) + (cos(-theta) * (y2 - half_length)) + half_length;

      x1 = floor(x0);
      y1 = floor(y0);

      if( x1 < 0 || y1 < 0 || x1 >= width || y1 >= length)
      {
        MatriceImgPivotee[y2][x2] = 0.0;
      }
      else
      {        
        // Calcul des pixels encadrant f(x', y')
        f_xy1 = RotatedMatrixElement(MatriceImg, x1, y1 + 1, length, width); // f(x, y+1)
        f_x1y1 = RotatedMatrixElement(MatriceImg, x1 + 1, y1 + 1, length, width); // f(x+1, y+1)
        f_xy = RotatedMatrixElement(MatriceImg, x1, y1, length, width); // f(x,y)
        f_x1y = RotatedMatrixElement(MatriceImg, x1 + 1, y1, length, width); // f(x+1, y)

        f_x0y = f_xy + (x0 - x1) * (f_x1y - f_xy); // f(x', y)
        f_x0y1 = f_xy1 + (x0 - x1) * (f_x1y1 - f_xy1); // f(x', y+1)

        // Interpolation bilineaire suivant les formules de l'annexe
        MatriceImgPivotee[y2][x2] = f_x0y + (y0 - y1) * (f_x0y1 - f_x0y); // f(x', y')
      }
    }
  }
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/
int main(int argc,char **argv)
 {
  int i,j,k;
  int length,width;
  int x0,y0;
  char BufSystVisuImg[100];

  //Constante
  length=512;
  width=512;
  
  //Allocation Memoire 
  float** MatriceImgI1=fmatrix_allocate_2d(length,width);
  float** MatriceImgM1=fmatrix_allocate_2d(length,width);
  float** MatriceImgI2=fmatrix_allocate_2d(length,width);
  float** MatriceImgM2=fmatrix_allocate_2d(length,width);
  float** MatriceImgPivotee=fmatrix_allocate_2d(length, width);

  //Lecture Image 
  float** MatriceImg1=LoadImagePgm(NAME_IMG_IN1,&length,&width);
  float** MatriceImg2=LoadImagePgm(NAME_IMG_IN2, &length, &width);

  // Initialisation de matrice final de l'image pivotee
  for(i=0; i<length; i++)
  {
    for(j=0;j<width;j++)
    {
      MatriceImgPivotee[i][j] = 0.0;
    }
  }

  // Calcul du module du spectre de l'image 1
  CenterImg(MatriceImg1, length, width);

  FFTDD(MatriceImg1, MatriceImgI1, length, width);

  Mod(MatriceImgM1, MatriceImg1, MatriceImgI1, length, width);

  Recal(MatriceImgM1, length, width);

  // Calcul du module du spectre de l'image 2
  CenterImg(MatriceImg2, length, width);

  FFTDD(MatriceImg2, MatriceImgI2, length, width);

  Mod(MatriceImgM2, MatriceImg2, MatriceImgI2, length, width);

  Recal(MatriceImgM2, length, width);

  // Pivot de le |G(u,v)| sur l'intervalle d'angle demander
  float theta;
  float erreur = 0;
  float bestErreur = INFINITY;
  float bestTheta;
  float pi = M_PI;

  for(theta = -pi/16; theta < pi/16; theta += 0.005)
  {
    // Rotation de l'image de theta
    RotationBilinearInterpolation(MatriceImgM2, length, width, MatriceImgPivotee, theta);

    // Calcul de l'erreur
    for(i=0;  i<length; i++)
    {
      for(j=0; j<width; j++)
      {
        erreur += abs(abs(MatriceImgPivotee[i][j]) - abs(MatriceImgM1[i][j]));
      }
    }

    // Sauvegarde du meilleur angle avec la plus petite erreur
    if(erreur < bestErreur)
    {
      bestErreur = erreur;
      bestTheta = theta;
    }

    printf("[%f::%f] > ", theta, erreur);

    // Reset de l'erreur pour prochaine iteration
    erreur = 0;
  }

  printf("\n\n Angle >>> [%f]", bestTheta);
  //==End=========================================================

  //Liberation memoire 
  free_fmatrix_2d(MatriceImg1);
  free_fmatrix_2d(MatriceImgPivotee);

  //retour sans probleme
  printf("\n C'est fini ... \n\n");
  return 0;    
}


