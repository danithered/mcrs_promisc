#include <mainheader.h>
#include <stdlib.h>
#include <randomgen.h>

int diffTM (struct Cella *matrix_f, double *enzim_f, int enzaktszam_f, int numcol_f, int target_f) {
  /* Troffuli-Marghulis szerinti elforgatas: random elforgat egy negyzetet 90 fokkal
   * 
   * matrix_f: alapmatrix
   * enzim_f: enzimaktivitasokat tartalmazo matrix
   * enzaktszam_f: hany fele enzimaktivitas van
   * numcol_f:az alapmatrix oszlopainak szama 
   * target_f: az elforgatott negyzet bal felso cellajanak szama
   */
	int *negyzet_f;
	int matrixmeret_f=numcol_f*numcol_f, enzakt_f, x_f=1, y_f=1, csereszam_f;
	double *saveD_f;
	
	saveD_f = (double *) calloc((3+enzaktszam_f), sizeof(double));
	negyzet_f= (int *) calloc (4, sizeof(int));
	
  /* negyzet_f: az elforgatanado cellak szamait tartalmazo matrix
   * matrixmeret_f: numcol_f*numcol_f
   * enzakt_f: szamlalo, vegigmegy az enzimaktivitasokon
   * x_f: ket erteke lehet: numcol_f vagy 1, ettol fugg, hogy merre megy majd a forgatas
   * y_f: ket erteke lehet: numcol_f vagy 1, ettol fugg, hogy merre megy majd a forgatas
   * csereszam_f: szamlalo 1tol 4ig, megmondja, hogy a negyzet_f altal meghatarotozott szamu cellat az eggyel nagyobb szamu cellaval csereljen
   * saveD_f: az elso cella ertekeit menti el double ertekkent
   */
	
	//Random irany
	if ((gsl_rng_uniform(r) < 0.5)) x_f = numcol_f;
	else y_f=numcol_f;
	
	//negyzet meghat
	*(negyzet_f+0)=target_f;
	*(negyzet_f+1)=torus(matrixmeret_f, target_f+x_f);
	*(negyzet_f+2)=torus(matrixmeret_f, target_f+x_f+y_f);
	*(negyzet_f+3)=torus(matrixmeret_f, target_f+y_f);
	
	//1. pozicio lementese
	*(saveD_f+0)=(*(matrix_f+target_f)).k;
	*(saveD_f+1)=(double) (*(matrix_f+target_f)).szerk;
	*(saveD_f+2)=(double) (*(matrix_f+target_f)).spec;
	for (enzakt_f=0; enzakt_f<enzaktszam_f;enzakt_f++) {
		*(saveD_f+3+enzakt_f)= *(enzim_f+target_f*enzaktszam_f+enzakt_f);
	}
	
	//poziciok csereje
	for(csereszam_f=0; csereszam_f<3; csereszam_f++) {
		copypaste(matrix_f, enzim_f, enzaktszam_f, *(negyzet_f+csereszam_f+1), *(negyzet_f+csereszam_f+0));
	}
	
	//mentes vissza
	(*(matrix_f+(*(negyzet_f+3)))).k = *(saveD_f+0);
	(*(matrix_f+(*(negyzet_f+3)))).szerk = (int) *(saveD_f+1);
	(*(matrix_f+(*(negyzet_f+3)))).spec = (int) *(saveD_f+2);
	for (enzakt_f=0; enzakt_f<enzaktszam_f;enzakt_f++) {
		*(enzim_f+(*(negyzet_f+3))*enzaktszam_f+enzakt_f) = *(saveD_f+3+enzakt_f);
	}
	
	free(negyzet_f);
	free(saveD_f);
	return(0);
}

void copypaste (struct Cella *matrix_f, double *enzim_f, int enzaktszam_f, int copy_f, int paste_f) {
  /* az egyik cella adatait masolja egy masikba
   * 
   * matrix_f:...
   * enzim_f:...
   * enzaktszam_f:...
   * copy_f: az ebben a cellaban levo replikator adatait masolja egy masikba
   * paste_f: ebbe a cellaba masolja
   */
		int enzakt_f;
  /*enzakt_f: szamlalo, az enzimaktivitasokon megy vegig
   */
		(*(matrix_f+paste_f)).k = (*(matrix_f+copy_f)).k;
		(*(matrix_f+paste_f)).szerk = (*(matrix_f+copy_f)).szerk;
		(*(matrix_f+paste_f)).spec = (*(matrix_f+copy_f)).spec;
		for (enzakt_f=0; enzakt_f<enzaktszam_f;enzakt_f++) {
			*(enzim_f+paste_f*enzaktszam_f+enzakt_f) = *(enzim_f+copy_f*enzaktszam_f+enzakt_f);
		}
}		

