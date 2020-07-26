#include <adattipusok.h>
#include <randomgen.h>
#include <stdio.h>
#include <string.h>


int binHeader(FILE *file_s, struct Save header_s, int ciklus_s, int tipus_s) {
	if(header_s.version==1) {
		header_s.gen=ciklus_s;
		header_s.tipus=tipus_s;
		fwrite(&header_s, sizeof(struct Save), 1, file_s);
	}
	else {
		printf("nem jo verzioju save header!!\n");
		return(1);
	}
	return(0);
}

int mentes(struct Cella *matrix_s, double *enzim_s, double *adatok_s, int state_s, int meret_s, int enzimaktszam_s, FILE *rngSave_s, FILE *dataSave_s) {
//int mentes(FILE *rngSave_s, FILE *dataSave_s) {
	
	//A randomszam generator mentese
	gsl_rng_fwrite (rngSave_s, r);
	fprintf(dataSave_s, "%s %d %d\n", gsl_rng_name (r), meret_s, enzimaktszam_s);
	
	//Az adatok mentese
	fajlbaCella(dataSave_s, matrix_s, enzim_s, enzimaktszam_s, meret_s, state_s);
	
	return(0);
}

int mentesBin(struct Cella *matrix_s, double *enzim_s, struct Save header_s, int meret_s, int ciklus_s, FILE *rngSave_s, FILE *dataSave_s, FILE *enzSave_s) {
	
	//A randomszam generator mentese
	gsl_rng_fwrite (rngSave_s, r);
		
	//A alapmatrix mentese
	binHeader(dataSave_s, header_s, ciklus_s, 1);
	fwrite(matrix_s, sizeof(struct Cella), meret_s, dataSave_s);

	//Az enzimaktivitas matrix mentese
	binHeader(enzSave_s, header_s, ciklus_s, 2);
	fwrite(enzim_s, sizeof( double), meret_s*(header_s.eakt), enzSave_s);

	
	return(0);
}

int writetoBin(struct Cella *matrix_s, double *enzim_s, struct Save header_s, int meret_s, int ciklus_s, FILE *dataSave_s, FILE *enzSave_s) {
/*ERRORS:
 * 	1: vmi baj van a megnyitott binaris fajllal
 * 	2: nem jo tipusu a binaris fajl
 */
	struct Save firstheader_s;
	int compare_s=0;
	
//	printf("ciklus: %d, header.tipus: %d\n", ciklus_s, header_s.tipus);
  
	//leellenoriz fajl meret
	fseek(dataSave_s, 0L, SEEK_END);
	//ha nem ures
	if(ftell(dataSave_s)) {
		//fajl elejere megy -> header modosit
		fseek(dataSave_s, 0L, SEEK_SET);
		fread(&firstheader_s, sizeof(struct Save), 1, dataSave_s);
//		printf("firstheader_s.tipus: %d\n", firstheader_s.tipus);
		compare_s=compareHeaders(firstheader_s, header_s);
		if(!(compare_s==0 || compare_s==2 || compare_s==4)) {
			printf("baj van adat hozzadasara megnyitott binaris fajllal\n");
			return(1);
		}
		if( !( (header_s.tipus==3 && firstheader_s.tipus == 3) || ((header_s.tipus==1 || header_s.tipus==2) && firstheader_s.tipus==1 ) ) ) {
			printf("nem jo tipusu a szerkesztesre megnyitott bin fajl:%d\n", firstheader_s.tipus);
			return(2);
		}
		firstheader_s.num++;
		fseek(dataSave_s, 0L, SEEK_SET);
		binHeader(dataSave_s, firstheader_s, firstheader_s.gen,firstheader_s.tipus);
		//fajl vegere, header modosit
		fseek(dataSave_s, 0L, SEEK_END);
		header_s.num=1-firstheader_s.num;
	}
	//A alapmatrix mentese
	if(header_s.tipus==3) binHeader(dataSave_s, header_s, ciklus_s, 3);
	else binHeader(dataSave_s, header_s, ciklus_s, 1);
	fwrite(matrix_s, sizeof(struct Cella), meret_s, dataSave_s);

	//Az enzimaktivitas matrix mentese
	if(header_s.tipus!=3) {
		binHeader(enzSave_s, header_s, ciklus_s, 2);
		fwrite(enzim_s, sizeof( double), meret_s*(header_s.eakt), enzSave_s);
	}
	else fwrite(enzim_s, sizeof( double), meret_s*(header_s.eakt), dataSave_s);
	
	return(0);
}


int compareHeaders(struct Save header1_s, struct Save header2_s) {
/*header fajlokat hasonlit ossze, megnezi volt-e bennuk hiba
 * visszateresi ertekek:
 *	0: minden rendben, szarmazhatnak (nem biztos!!!!!) azonos futasbol es mas tipusuak
 *	-1: nem lett feltoltve az 1. matrix, mert nem lehetett megnyitni a fajlt
 *	-2: nem ismert tipusu az 1. matrix
 *	-3: nem lett feltoltve az 1. matrix, mert nem ismerte fel a betolto program a fajl verziojat
 *	-4: nem lett feltoltve az 1. matrix, mert nem lehetett lefoglalni
 *	-10: ismeretlen hiba az 1. matrixban
 *	-101: nem lett feltoltve a 2. matrix, mert nem lehetett megnyitni a fajlt
 *	-102: nem ismert tipusu a 2. matrix
 *	-103: nem lett feltoltve a 2. matrix, mert nem ismerte fel a betolto program a fajl verziojat
 *	-104: nem lett feltoltve a 2. matrix, mert nem lehetett lefoglalni
 *	-110: ismeretlen hiba a 2. matrixban
 * 	1: eltero verzioju a ket header
 * 	2: megegyezo tipusu a ket header
 * 	3: eltero randomszam generatort ir le a ket header
 * 	4: eltero generaciobol szarmazik a ket header
 * 	5
 * 	6
 * 	7
 * 	8
 * 	9
 * 	10
 * 	11
 * 	12
 * 	13
 * 	14
 * 	15
 * 	16
 * 	17
 * 	18
 * 	19
 */	
	//version
		//header1
		if(header1_s.version<1){
			switch(header1_s.version){
			  case -1:
				printf("nem lett feltoltve az 1. matrix, mert nem lehetett megnyitni a fajlt\n");
				return(-1);
			  case -2:
				printf("nem ismert tipusu az 1. matrix\n");
				return(-2);
			  case -3:
				printf("nem lett feltoltve az 1. matrix, mert nem ismerte fel a betolto program a fajl verziojat\n");
				return(-3);
			  case -4:
				printf("nem lett feltoltve az 1. matrix, mert nem lehetett lefoglalni\n");
				return(-4);
			  case -5:
				 printf("kert adat nem kiolvashato (1.matrix)\n");
				 return(-5);
			  default:
				printf("ismeretlen hiba az 1. matrixban\n");
				return(-10);
			}
		}
		//header2
		if(header2_s.version<1){
			switch(header2_s.version){
			  case -1:
				printf("nem lett feltoltve a 2. matrix, mert nem lehetett megnyitni a fajlt\n");
				return(-101);
			  case -2:
				printf("nem ismert tipusu a 2. matrix\n");
				return(-102);
			  case -3:
				printf("nem lett feltoltve a 2. matrix, mert nem ismerte fel a betolto program a fajl verziojat\n");
				return(-103);
			  case -4:
				printf("nem lett feltoltve a 2. matrix, mert nem lehetett lefoglalni\n");
				return(-104);
			  case -5:
				 printf("kert adat nem kiolvashato (2.matrix)\n");
				 return(-105);
			  default:
				printf("ismeretlen hiba a 2. matrixban\n");
				return(-110);
			}
		}
		//osszehasonlit
		if(header1_s.version!=header2_s.version){
			  printf("eltero verzioju a ket header! (%d, %d)\n", header1_s.version, header2_s.version);
			  return(1);
		}
	//tipus
	if(header1_s.tipus==header2_s.tipus){
		  //printf("megegyezo tipusu a ket header: %d\n", header1_s.tipus);
		  return(2);
	}
	
	//rng
	if(strcmp(header1_s.rng, header2_s.rng)){
		  printf("eltero randomszam generatort ir le a ket header:\n1.: %s\t2.: %s\n", header1_s.rng, header2_s.rng);
		  return(3);
	}
	
	//gen
	if(header1_s.gen!=header2_s.gen){
		  printf("eltero generaciobol szarmazik a ket header\n");
		  return(4);
	}

	//ncol, nrow
	if((header1_s.ncol*header1_s.nrow) != (header2_s.ncol*header2_s.nrow)){
		  printf("eltero cellaszamu a ket header\n");
		  return(5);
	}
	
	if((header1_s.ncol!=header2_s.ncol) || (header1_s.nrow!=header2_s.nrow)){
		  printf("eltero meretu a ket header\n");
		  return(6);
	}
	
	//eakt
	if(header1_s.eakt != header2_s.eakt){
		  printf("eltero enzimaktivitas szamu a ket header\n");
		  return(7);
	}
	
	//met
	if(header1_s.met != header2_s.met){
		  printf("eltero metab szomsz szam a ket headerben\n");
		  return(8);
	}
	
	//repl
	if(header1_s.repl != header2_s.repl){
		  printf("eltero repl szomsz szam a ket headerben\n");
		  return(9);
	}
	
	//emax
	if(header1_s.emax != header2_s.emax){
		  printf("eltero emax a ket headerben\n");
		  return(10);
	}
	
	//trEE
	if(header1_s.trEE != header2_s.trEE){
		  printf("eltero trEE a ket headerben\n");
		  return(11);
	}
	
	//trEK
	if(header1_s.trEK != header2_s.trEK){
		  printf("eltero trEK a ket headerben\n");
		  return(12);
	}
	
	//kmin
	if(header1_s.kmin != header2_s.kmin){
		  printf("eltero kmin a ket headerben\n");
		  return(13);
	}
	
	//kmax
	if(header1_s.kmax != header2_s.kmax){
		  printf("eltero kmax a ket headerben\n");
		  return(14);
	}
	
	//phalal
	if(header1_s.phalal != header2_s.phalal){
		  printf("eltero phalal a ket headerben\n");
		  return(15);
	}
	
	//claimE
	if(header1_s.claimE != header2_s.claimE){
		  printf("eltero claimE a ket headerben\n");
		  return(16);
	}
	
	//sd
	if(header1_s.sd != header2_s.sd){
		  printf("eltero sd a ket headerben\n");
		  return(17);
	}
	
	//d
	if(header1_s.d != header2_s.d){
		  printf("eltero d ertekek a ket headerben\n");
		  return(18);
	}
	
	//pmut
	if(header1_s.pmut != header2_s.pmut){
		  printf("eltero pmut a ket headerben\n");
		  return(19);
	}
		  
	return(0);
};