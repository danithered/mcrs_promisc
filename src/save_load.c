#include <olvaso.h>
#include <randomgen.h>
#include <stdio.h>
#include <string.h> 

int betoltesRng(FILE *rngLoad_s) {
	r = (gsl_rng *) gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_fread (rngLoad_s, r);
  
	return(0);
}

struct Save betoltes(char *fajl_s, int hanyadik_s) {
/*errors:
 * 	kimenet.version = -1 :fajl nem megnyithato
 * 	kimenet.version = -2 :nem ismert tipus (1 és 2 ismert)
 * 	kimenet.version = -3 :nem jo verzioju fajl
 * 	kimenet.version = -4 :nem lehet lefoglalni matrixot
 * 	kimenet.version = -5 :kert adat nem kiolvashato
 * 
 */
	struct Save header_s;
	FILE *ptrFajl_s;
	int matrixMeret_s;
	struct Cella *targetD_s; 
	double *targetE_s;
	struct SaveAll out_s;
	
	//megnyit fajlt
	ptrFajl_s = fopen(fajl_s, "rb");
	if(!ptrFajl_s) {
		printf("nem lehet megnyitni fajlt: %s\nBeolvasas megszakitva\n", fajl_s);
		header_s.version=-1;
		fclose(ptrFajl_s);
		return(header_s);
	}
	
	//beolvas enzimaktivitas matrix headert
	fread(&header_s, sizeof(struct Save), 1, ptrFajl_s);
	
	//megfelelo adathoz teker
	if(hanyadik_s != 1) {
		fseek(ptrFajl_s, 0L, SEEK_SET); //rewind: teker-nek ez kell
		if(hanyadik_s < 1) teker(ptrFajl_s, 0); //az utolsot olvassa be, oda teker
		if(header_s.num < hanyadik_s) {
			printf("Nincs %d db adatsor a fajlban, csak %d. Ha utolsot akarod beolvastatni, akkor irj <1 erteket!\nBeolvasas megszakitva (%s)!\n", hanyadik_s, header_s.num, *fajl_s);
			header_s.version=-5;
			fclose(ptrFajl_s);
			return(header_s);
		}
		else teker(ptrFajl_s, hanyadik_s-1); //ennyit teker elore
		fread(&header_s, sizeof(struct Save), 1, ptrFajl_s); //uj headert beolvas
	}
	
	//jo verzioju-e enzimmatrix?
	if(header_s.version==1) {
		//milyen adatok vannak tarolva?
		switch(header_s.tipus) {
		  case 1: //alapmatrix
			matrixMeret_s=header_s.ncol*header_s.nrow;
			//lefoglal matrix
			targetD_s=(struct Cella *) calloc(matrixMeret_s, sizeof(struct Cella));
			if(!targetD_s) {
				printf("nem lehet lefoglalni alapmatrixot");
				header_s.version=-4;
			}
			//ADATOK BEOLVASASA
			fread(targetD_s, sizeof(struct Cella), matrixMeret_s, ptrFajl_s);
			break;
		  case 2: //enzimmatrix
			matrixMeret_s=header_s.ncol*header_s.nrow*header_s.eakt;
			//lefoglal matrix
			targetE_s=(double *) calloc(matrixMeret_s, sizeof(double));
			if(!targetE_s) {
				printf("nem lehet lefoglalni enzimaktivitas matrixot");
				header_s.version=-4;
			}
			//ADATOK BEOLVASASA
			fread(targetE_s, sizeof(double), matrixMeret_s, ptrFajl_s);
			break;
		  case 3: //alap + enzimmatrix
			//alap beolvas
			matrixMeret_s=header_s.ncol*header_s.nrow;
			//lefoglal matrix
			targetD_s=(struct Cella *) calloc(matrixMeret_s, sizeof(struct Cella));
			if(!targetD_s) {
				printf("nem lehet lefoglalni alapmatrixot");
				header_s.version=-4;
			}
			//ADATOK BEOLVASASA
			fread(targetD_s, sizeof(struct Cella), matrixMeret_s, ptrFajl_s);
			//enzim beolvas
			matrixMeret_s*=header_s.eakt;
			//lefoglal matrix
			targetE_s=(double *) calloc(matrixMeret_s, sizeof(double));
			if(!targetE_s) {
				printf("nem lehet lefoglalni enzimaktivitas matrixot\n");
				header_s.version=-4;
			}
			//ADATOK BEOLVASASA
			fread(targetE_s, sizeof(double), matrixMeret_s, ptrFajl_s);
			break;
		  default: //ha nem ismert tipus
			printf("nem ismert adattipus:%s\nBeolvasas megszakitva\n", fajl_s);
			header_s.version=-2;
			break;
		}
	}
	else {//ha nem jo verzioju az enzimmatrix
		printf("matrix nem jo verzio:%s\nBeolvasas megszakitva\n", fajl_s);
		header_s.version=-3;
	}
		
	fclose(ptrFajl_s);
	
	return(header_s);
}

struct Save teker(FILE *file_s, int teker_s) {
/*mentesi fajlt tekeri 'teker_s'-sel elore vagy vissza
 * fontos: header elott alljon a FILE pointer!!!
 * teker_s:
 * 	0: utolsora teker (az elejerol porget a legvegere)
 * 	negativ: visszafele teker (ez esetben az elejerol porget vissza)
 * 	pozitiv: elore teker
 */
	int poz_s=0;
	struct Save header_s;
	
	if(teker_s <= 0) { //utolsora teker vagy visszafele teker
		if(teker_s < 0) { //visszaporgeteshez beolvas poziciot ( hova menjen: poz-(-target_s) )
			fread(&header_s, sizeof(struct Save), 1, file_s);
			poz_s = (header_s.num < 0) ? (1 - header_s.num) : 1;
		}
		fseek(file_s, 0L, SEEK_SET);
		fread(&header_s, sizeof(struct Save), 1, file_s);
		fseek(file_s, 0L, SEEK_SET);
		
		if(teker_s < 0) { //visszafele teker
			teker_s = poz_s + teker_s -1;
			if(teker_s < 0) {
				printf("Nincs ennyi adat a fajlban, az elso adat ele akart tekerni\n");
				header_s.tipus=0;
				return(header_s);
			}
		}
		else { //vegere megy
			teker_s = header_s.num;
		}
	}
	//elore teker
	for(;teker_s > 0; teker_s--) {
		fread(&header_s, sizeof(struct Save), 1, file_s);
		fseek(file_s, adatMeret(header_s.ncol, header_s.nrow, header_s.tipus), SEEK_CUR);
	}
	
	//beolvas headert es visszaporget
	fread(&header_s, sizeof(struct Save), 1, file_s);
	if (header_s.version != 1) {//ha nem jo verzioju az enzimmatrix
		printf("matrix nem jo verzio. Beolvasas megszakitva\n");
		header_s.version=-3;
	}
	fseek(file_s, -(sizeof(struct Save)), SEEK_CUR);
	
	return(header_s);
}

int adatMeret(int ncol_s, int nrow_s, int tipus_s) {
	switch (tipus_s) {
		case 1: //alapmatrix
			return(ncol_s * nrow_s * sizeof(struct Save));
		case 2: //enzimmatrix
			return(ncol_s * nrow_s * sizeof(double));
		case 3: //alap + enzim
			return(ncol_s * nrow_s * (sizeof(double) + sizeof(struct Save) ) );
		default:
			printf("Nem ismert tipus: %d\n", tipus_s);
			return(-1);
	}
}
