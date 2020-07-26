#include <adattipusok.h>
#include <stdlib.h>
#include <stdio.h>

void tipusadatok(double *eredmeny_f, struct Cella *matrix_f, double *enzim_f, int enzaktnum_f, int meret_f, int enzimtipus_num_f, int sorhossz_f) {
/* megnezi, hogy matrix adott enzimfajtanak milyen tulsagai vanna: db, k, e1, eN <- atlag, szoras
 * 
 * matrix_f:...
 * enzim_f:...
 * enzaktnum_f:...
 * meret_f: a sejtautomata cellainak szama
 * enzimtipus_num_f: hany fajta enzim van
 */
	int cella_f, enzakt_f, enzimtipus_f, db_f=0;
	double *adat_f;
	//sorhossz_f=1+2*(enzaktnum_f+1)
	
/**/	//printf("enzaktnum: %d, meret: %d, enzimtipus_num_f: %d, sorhossz_f: %d", enzaktnum_f, meret_f, enzimtipus_num_f, sorhossz_f);
	
	/* cella_f: szamlalo, hanyadik cellanal jarunk 
	 * sorhossz_f: milyen hosszu egy adatsor, azaz az egyfajta enzim adatainak szama
	 * enzakt_f: szamlalo, hanyadik enzimaktivitasnal jarunk
	 * enzimtipus_f: szamlalo, hanyadik enzimtipust nezzuk
	 * db_f: hany db minta volt (n)
	 * eredmeny_f: az eredmenyeket tartalmazo matrix:
	 * 	0.: db szam
	 * 	1.: k atlaga
	 * 	2 - (enzaktnum_f+2).: enzmaktivitasok atlaga
	 * 	1+(enzaktnum_f+2).: k szorasa
	 * 	2+(enzaktnum_f+2) - 2+(enzaktnum_f+2)+enzaktnum_f.: enzimaktivitasok szorasa
	 * adat_f: a kovetkezo modositando adat
	 */
	
	//nullazas
	for(cella_f=0; cella_f<sorhossz_f*enzimtipus_num_f; cella_f++) {*(eredmeny_f+cella_f)=0;}
	//eredmeny_f= (double *) calloc((enzimtipus_num_f*sorhossz_f), sizeof(double));
	//if(!eredmeny_f) printf("\neredmeny_f lefoglalasa nem sikerult (tipusadatok fuggveny)\n");
	
	//for(cella_f=0; cella_f<(enzimtipus_num_f*sorhossz_f); cella_f++) { *(eredmeny_f+cella_f) = cella_f;}
	
	for(cella_f=0; cella_f<meret_f; cella_f++) {
		if( (*(matrix_f+cella_f)).k ) {
			adat_f= eredmeny_f+sorhossz_f* (*(matrix_f+cella_f)).spec;
/**/			//printf("\nadat_f: %g\n", *adat_f);
			//db szam noveles
			(*adat_f)++;
			//k osszeg -> atlag
			*(++adat_f) += (*(matrix_f+cella_f)).k;
			//k negyzetosszeg -> szoras
			*(adat_f+enzaktnum_f+1) += (*(matrix_f+cella_f)).k * (*(matrix_f+cella_f)).k;
			//enzimaktivitasokon vegigmegy
			for (enzakt_f=0; enzakt_f < enzaktnum_f; enzakt_f++) {
				//enzimaktivitasok osszege -> atlagok
				*(++adat_f) += *(enzim_f+enzaktnum_f*cella_f+enzakt_f);
				//enzimaktivitasok negyzetosszege -> szorasok
				*(adat_f+enzaktnum_f+1)+= *(enzim_f+enzaktnum_f*cella_f+enzakt_f) * *(enzim_f+enzaktnum_f*cella_f+enzakt_f);
			}
		}
	}
	
	//vegleges kiszamolas
	for(enzimtipus_f=0; enzimtipus_f<enzimtipus_num_f;enzimtipus_f++) {
		db_f =(int) *(eredmeny_f+enzimtipus_f*sorhossz_f);
		adat_f= eredmeny_f+enzimtipus_f*sorhossz_f+1+enzaktnum_f+1; //a kiszamolando cellara mutato pointer
		//szorasok kiszamitasa
		for (cella_f=0; cella_f < (enzaktnum_f+1); cella_f++) {
			if(*adat_f) *adat_f= szorasD(db_f,*(adat_f-enzaktnum_f-1) , *adat_f);
			adat_f++;
		}
		//atlagok kiszamitasa
		adat_f= eredmeny_f+enzimtipus_f*sorhossz_f+1;
		for (cella_f=0; cella_f < (enzaktnum_f+1); cella_f++) {
			if(*adat_f) *adat_f /= db_f;
			adat_f++;
		}
		
	}
  
	//return(eredmeny_f);
}

void kimenetOszlopnevek (FILE *output_f,  int enzaktnum_f) {
	int enzakt_f;
	
	fprintf(output_f, "ciklus;enzimtipus;db;k_atlag");
	for (enzakt_f=0; enzakt_f<enzaktnum_f; enzakt_f++) {
		fprintf(output_f, ";E%d_atlag", enzakt_f);
	}
	fprintf(output_f, ";k_sd");
	for (enzakt_f=0; enzakt_f<enzaktnum_f; enzakt_f++) {
		fprintf(output_f, ";E%d_sd", enzakt_f);
	}
	fprintf(output_f, "\n");
}

void kimenetTipusadat (FILE *output_f, struct Cella *matrix_f, double *eredmeny_f, int enzimtipus_num_f, int sorhossz_f, int ciklusszam_f) {
	int adatszamlalo_f, sorszamlalo_f;
	
	for (sorszamlalo_f=0; sorszamlalo_f<enzimtipus_num_f; sorszamlalo_f++) {
		fprintf(output_f, "%d;\"%010d\"", ciklusszam_f, kiirbit( sorszamlalo_f , BITNUM));
		for (adatszamlalo_f=0; adatszamlalo_f<sorhossz_f; adatszamlalo_f++) {
			fprintf(output_f, ";%G", (double) *(eredmeny_f+sorszamlalo_f*sorhossz_f+adatszamlalo_f));
		}
		fprintf(output_f, "\n");
	}
}

void fajlbaCellaHeaders (FILE *output_f, int enzaktnum_f) {
	int enzakt_f;
	
	fprintf(output_f, "ciklus;cella;spec;fold;k");
	for(enzakt_f=0; enzakt_f<enzaktnum_f; enzakt_f++) {
		fprintf(output_f, ";E%d", enzakt_f);
	}
	fprintf(output_f, "\n");
	
}

void fajlbaCella (FILE *output_f, struct Cella *matrix_f, double *enzim_f, int enzaktnum_f, int meret_f, int ciklusszam_f) {
	int cella_f, enzakt_f;
	
//	printf("done\n");
	for (cella_f=0; cella_f<meret_f; cella_f++) {
		fprintf(output_f, "%d;%d;%d;%d;%6.5G", ciklusszam_f, cella_f, (*(matrix_f+cella_f)).spec, (*(matrix_f+cella_f)).szerk, (*(matrix_f+cella_f)).k);
		for(enzakt_f=0; enzakt_f<enzaktnum_f; enzakt_f++) {
			fprintf(output_f, ";%8.7G", *(enzim_f+cella_f*enzaktnum_f+enzakt_f));
		}
		fprintf(output_f, "\n");
	}
}