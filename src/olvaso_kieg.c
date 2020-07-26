#include <olvaso.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double metabolizmus(struct Cella *matrix_f, double *enzim_f, int *met_szomszedsag_f, int szomsz_cellaszam_f, int enzimaktszam_f, double reciprocEnzimaktszam_f, int cella_f) {
//	printf("\nenzimaktszam=%d", enzimaktszam_f);
	int szomszed_f, nezett_f=0, enzakt_f, dominans_f;
	double *enzimsum_f, metab_f=1;
	
	/* szomszed_f: a target (cella_f) szomszedjanak sorszama/ hanyadik szomszedot nezzuk eppen
	 * nezett_f: hanyadik cellat nezzuk eppen
	 * enzakt_f: a nezett enzakt szamlaloja
	 * enzimsum_f: az adott enzimaktivitasok osszege
	 * metab_f: M, azaz a metabolizmus erteke az adott cellara
	 */
	
	enzimsum_f=(double*) calloc(enzimaktszam_f, sizeof(double));
	
//	printf("\ncella: %d (szerk, spec, eakt)\n", cella_f);
	for(szomszed_f=0; szomszed_f < szomsz_cellaszam_f; szomszed_f++) {
		nezett_f=*(met_szomszedsag_f+cella_f*szomsz_cellaszam_f+szomszed_f);
		if (!(*(matrix_f+nezett_f)).k) continue;
		//*(enzimsum_f+(*(matrix_f+nezett_f)).k) = *(enzim_f+nezett_f*enzimaktszam_f+(*(matrix_f+nezett_f)).k);
	//	TRANSITION (matrix_f, enzim_f, enzimaktszam_f, 0.5, nezett_f);
		dominans_f=(*(matrix_f+nezett_f)).szerk;
		*(enzimsum_f+dominans_f) += *(enzim_f+nezett_f*enzimaktszam_f+dominans_f);
//		printf("\nnezett=%d\tdominans=%d", nezett_f, dominans_f);
//		printf(" (%d, %d, %g)", (*(matrix_f+nezett_f)).szerk, (*(matrix_f+nezett_f)).spec, *(enzim_f+nezett_f*enzimaktszam_f+dominans_f));
	}
	for(enzakt_f=0; enzakt_f<enzimaktszam_f; enzakt_f++){
		metab_f *= *(enzimsum_f+enzakt_f);
//		printf("\nenzimsum= %g", *(enzimsum_f+enzakt_f));
	}
//	printf("\n\nmetab_f=%g\t1/enzimaktszam_f=%g", metab_f, (1/ (double) enzimaktszam_f));
//	printf("\nmetab_f=%g, reciprocEnzimaktszam_f=%g", metab_f, reciprocEnzimaktszam_f);
	
	if (metab_f>0) metab_f=pow(metab_f, reciprocEnzimaktszam_f);
	
	free (enzimsum_f);
//		printf("\nvegso metab_f=%g", metab_f);
	return (metab_f);
} 

int fajlbaMetab (FILE *output_f, struct Cella *matrix_f, double *enzim_f, double met_szomsz_f, int enzimaktszam_f, int ncol_f, int nrow_f) {
	int *met_szomsz_matrix_f, met_neigh_cellaszam_f = szomsz_meret(met_szomsz_f), meret_f = ncol_f * nrow_f, cella_f;
	double reciprocEnzimaktszam_f = 1 / (double)enzimaktszam_f;
//	int sorszam =0;
	
	met_szomsz_matrix_f = (int *)metNeighInic(meret_f, ncol_f, met_szomsz_f);
	if(!met_szomsz_matrix_f) printf("nem sikerult metab szomszedsag matrixot lefoglalni!\n");
//	else printf("matrix lefoglalva!\n");
	
	for(cella_f=0; cella_f < meret_f; cella_f++) {
//		printf("cella: %d / %d\n", cella_f, meret_f);
		fprintf(output_f, "%f", metabolizmus(matrix_f, enzim_f, met_szomsz_matrix_f, met_neigh_cellaszam_f, enzimaktszam_f, reciprocEnzimaktszam_f, cella_f));
		if( cella_f==0 || ((cella_f % ncol_f) != (ncol_f-1)) ) fprintf(output_f, ";");
		else {
			fprintf(output_f, "\n");
//			sorszam++;
		}
	}
	
//	printf("\n%d sor kiirva\n", sorszam);
	free(met_szomsz_matrix_f);
	return (0);
}

int fajlbaSzomszTipusok (FILE *output_f, struct Cella *matrix_f, double met_szomsz_f, int enzimaktszam_f, int ncol_f, int nrow_f) {
	int *met_szomsz_matrix_f, *szomszedok_f, meret_f = nrow_f * ncol_f, cella_f = 0, szomsz_f = 0, met_neigh_cellaszam_f = szomsz_meret(met_szomsz_f), noEnzimtipus_f = hanyfajta(enzimaktszam_f);
	
	/*meret_f: az alapmatrix merete
	 * cella_f: szamlalo, a matrix cellain megy vegig
	 * szomsz_f: szamlalo, a met. szomszedsagon megy vegig, kesobb a kiirasnal a kulonboz tipusu enzimtipusokon megy vegig
	 */
	
	met_szomsz_matrix_f = (int *)metNeighInic(meret_f, ncol_f, met_szomsz_f);
	for (cella_f = 0; cella_f < meret_f; cella_f++) {
		szomszedok_f = (int *) calloc(noEnzimtipus_f, sizeof(int));
		for(szomsz_f = 0; szomsz_f < met_neigh_cellaszam_f; szomsz_f++) {
			(*(szomszedok_f + (*(matrix_f + (*(met_szomsz_matrix_f + cella_f*met_neigh_cellaszam_f + szomsz_f )) )).spec ))++;
		}
		for(szomsz_f = 0; szomsz_f < (noEnzimtipus_f-1); szomsz_f++) { //utolsot kiveve kiir
			fprintf(output_f, "%d;", *(szomszedok_f + szomsz_f) );
		}
		fprintf(output_f, "%d\n", *(szomszedok_f + noEnzimtipus_f)); //utolso kiir, sorzaras
		free(szomszedok_f);
	}
	
	free(met_szomsz_matrix_f);
	
	return(0);
}

int fajlbaAsszocTabla (FILE *output_f, struct Cella *matrix_f, double met_szomsz_f, int enzimaktszam_f, int ncol_f, int nrow_f) {
	int *met_szomsz_matrix_f, *szomszedok_f, *asszocT_f, meret_f = nrow_f * ncol_f, cella_f = 0, szomsz_f = 0, szomsz2_f = 0, met_neigh_cellaszam_f = szomsz_meret(met_szomsz_f), noEnzimtipus_f = hanyfajta(enzimaktszam_f), asszocMeret_f = noEnzimtipus_f*noEnzimtipus_f;
	
	/*meret_f: az alapmatrix merete
	 * cella_f: szamlalo, a matrix cellain megy vegig
	 * szomsz_f: szamlalo, a met. szomszedsagon megy vegig, kesobb a kiirasnal a kulonboz tipusu enzimtipusokon megy vegig
	 */
	
	met_szomsz_matrix_f = (int *)metNeighInic(meret_f, ncol_f, met_szomsz_f);
	asszocT_f = (int *)calloc(asszocMeret_f, sizeof(int));
	if(!asszocT_f) printf("nem sikerult asszoc tablat lefoglalni!\n");
	
	for (cella_f = 0; cella_f < meret_f; cella_f++) {
		szomszedok_f = (int *) calloc(noEnzimtipus_f, sizeof(int));
		for(szomsz_f = 0; szomsz_f < met_neigh_cellaszam_f; szomsz_f++) {
			(*(szomszedok_f + (*(matrix_f + (*(met_szomsz_matrix_f + cella_f*met_neigh_cellaszam_f + szomsz_f )) )).spec ))++;
		}
		for(szomsz_f = 0; szomsz_f < noEnzimtipus_f; szomsz_f++) {
			if( *(szomszedok_f + szomsz_f) > 0) { //van a cellaban vki
				if( *(szomszedok_f + szomsz_f) > 1 ) { //onmagaval asszocialt
					(*(asszocT_f + szomsz_f*noEnzimtipus_f + szomsz_f))++;
				}
				if(szomsz_f < (noEnzimtipus_f-1) ) { //ha nem utolso
					for (szomsz2_f = szomsz_f+1; szomsz2_f < noEnzimtipus_f; szomsz2_f++) {
						if( *(szomszedok_f + szomsz2_f) > 0 ) { //asszocialt
							(*(asszocT_f + szomsz2_f*noEnzimtipus_f + szomsz_f))++;
							(*(asszocT_f + szomsz_f*noEnzimtipus_f + szomsz2_f))++;
						}
					}
				}
			}
		}
		free(szomszedok_f);
	}
	for(szomsz_f=0; szomsz_f < asszocMeret_f; szomsz_f++) { //vegigmegy asszov tablan
		if( ((szomsz_f + 1)%noEnzimtipus_f) == 0 ) fprintf(output_f, "%d\n", *(asszocT_f + szomsz_f));
		else fprintf(output_f, "%d;", *(asszocT_f + szomsz_f));
	}
	
	free(met_szomsz_matrix_f);
	free(asszocT_f);
	
	return(0);
}
