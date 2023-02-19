/*-------------------------------------------------------------------*/
/* Prog    : Tp2_IFT3205-2-1.c                                       */
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
#define NAME_IMG_IN2  "UdM_2"
#define NAME_IMG_OUT1 "image-TpIFT3205-2-1a"
#define NAME_IMG_OUT2 "image-TpIFT3205-2-1b"

/*------------------------------------------------*/
/* PROTOTYPE DE FONCTIONS  -----------------------*/   
/*------------------------------------------------*/

// Re-Utilisation de la fonction creer dans le TP1 1.3 pour mettre en valeur les basses frequences.
void MultEnLogarithme(float** mat, int length,int width)
{
  int i,j;

  for(i=0;i<length;i++)
  {
    for(j=0;j<width;j++) 
    {
        mat[i][j] = log(1 + mat[i][j]);
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
  float Theta0;
  int x0,y0;
  char BufSystVisuImg[100];

  //Constante
  length=512;
  width=512;
  
  //Allocation Memoire 
  float** MatriceImgI1=fmatrix_allocate_2d(length,width);
  float** MatriceImgM1=fmatrix_allocate_2d(length,width);
  float** MatriceImgR1=fmatrix_allocate_2d(length,width);
  float** MatriceImgI2=fmatrix_allocate_2d(length,width);
  float** MatriceImgM2=fmatrix_allocate_2d(length,width);
  float** MatriceImgR2=fmatrix_allocate_2d(length,width);
  float** MatriceImgI3=fmatrix_allocate_2d(length,width);
  float** MatriceImgM3=fmatrix_allocate_2d(length,width);
  float** MatriceImgR3=fmatrix_allocate_2d(length,width);
  float** MatriceImg3=fmatrix_allocate_2d(length,width);

  //Lecture Image 
  float** MatriceImg1=LoadImagePgm(NAME_IMG_IN1,&length,&width);
  float** MatriceImg2=LoadImagePgm(NAME_IMG_IN2,&length,&width);

 
  // Calcul du module spectre de l'image 1 et mis au centre humain
  CenterImg(MatriceImg1, length, width);

  FFTDD(MatriceImg1, MatriceImgI1, length, width);

  Mod(MatriceImgM1, MatriceImg1, MatriceImgI1, length, width);

  CenterImg(MatriceImgM1, length, width);

  MultEnLogarithme(MatriceImgM1, length, width);

  Recal(MatriceImgM1, length, width);

  // Calcul du module du spectre de l'image 2 et mis au centre humain
  CenterImg(MatriceImg2, length, width);

  FFTDD(MatriceImg2, MatriceImgI2, length, width);

  Mod(MatriceImgM2, MatriceImg2, MatriceImgI2, length, width);

  CenterImg(MatriceImgM2, length, width);

  MultEnLogarithme(MatriceImgM2, length, width);

  Recal(MatriceImgM2, length, width);

  //Sauvegarde
  SaveImagePgm(NAME_IMG_OUT1,MatriceImgM1,length,width);
  SaveImagePgm(NAME_IMG_OUT2,MatriceImgM2,length,width);

  //Commande systeme: VISU
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_OUT1);
  strcat(BufSystVisuImg,".pgm&");
  printf(" %s",BufSystVisuImg);
  system(BufSystVisuImg);
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_OUT2);
  strcat(BufSystVisuImg,".pgm&");
  printf(" %s",BufSystVisuImg);
  system(BufSystVisuImg);


  //==End=========================================================

  //Liberation memoire 
  free_fmatrix_2d(MatriceImgR1);
  free_fmatrix_2d(MatriceImgI1);
  free_fmatrix_2d(MatriceImgM1);
  free_fmatrix_2d(MatriceImgR2);
  free_fmatrix_2d(MatriceImgI2);
  free_fmatrix_2d(MatriceImgM2);
  free_fmatrix_2d(MatriceImgR3);
  free_fmatrix_2d(MatriceImgI3);
  free_fmatrix_2d(MatriceImgM3);
  free_fmatrix_2d(MatriceImg1);
  free_fmatrix_2d(MatriceImg2);  
  free_fmatrix_2d(MatriceImg3);

  //retour sans probleme
  printf("\n C'est fini ... \n\n");
  return 0; 	 
}


