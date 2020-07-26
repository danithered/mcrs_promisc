#include <mainheader.h>
#include <math.h>

/*Alapmatrix kiirasa konzolra
 !!!: Kell neki a stuct Cella*/
void konzolraMatrixStruct(struct Cella *matrix_f, int meret_f) {
	int cellaszam_f, col_f, row_f, gyok_f;
	
	gyok_f=sqrt(meret_f);
	printf("\nk ertekek:\n");
	for(row_f=0; row_f<gyok_f; row_f++) {
		for(col_f=0; col_f<gyok_f; col_f++) {
			printf("%g\t", (*(matrix_f+(row_f*gyok_f+col_f))).k);
		}
		printf("\n");
	}	
	printf("\nstruktura ertekek:\n");
	for(row_f=0; row_f<gyok_f; row_f++) {
		for(col_f=0; col_f<gyok_f; col_f++) {
			printf("%d\t", (*(matrix_f+(row_f*gyok_f+col_f))).szerk);
		}
		printf("\n");
	}
	printf("\nenzim specifikacio tipusok:\n");
	for(row_f=0; row_f<gyok_f; row_f++) {
		for(col_f=0; col_f<gyok_f; col_f++) {
			printf("%d\t", (*(matrix_f+(row_f*gyok_f+col_f))).spec);
		}
		printf("\n");
	}
	
	
	return ;
}

void konzolraMatrix(int *matrix_f, int ncol_f, int nrow_f) {
	int cellaszam_f, col_f, row_f;
	
	for(row_f=0; row_f<nrow_f; row_f++) {
		for(col_f=0; col_f<ncol_f; col_f++) {
			printf("%d\t", *(matrix_f+(row_f*ncol_f+col_f)));
		}
		printf("\n");
	}	
		
	return ;
}

void konzolraMatrixD(double *matrix_f, int ncol_f, int nrow_f) {
	int cellaszam_f, col_f, row_f;
	
	for(row_f=0; row_f<nrow_f; row_f++) {
		for(col_f=0; col_f<ncol_f; col_f++) {
			printf("%g\t", *(matrix_f+(row_f*ncol_f+col_f)));
		}
		printf("\n");
	}	
		
	return ;
}

void konzolraCella (struct Cella *matrix_f, double *enzim_f, int enzaktnum_f, int meret_f) {
	int cella_f, enzakt_f;
	
	printf("\n");
	for (cella_f=0; cella_f<meret_f; cella_f++) {
		printf("%d.:\tk= %6.5g\tfold=%d\tspec=%d", cella_f, (*(matrix_f+cella_f)).k, (*(matrix_f+cella_f)).szerk, (*(matrix_f+cella_f)).spec);
		for(enzakt_f=0; enzakt_f<enzaktnum_f; enzakt_f++) {
			printf("\tE%d=%8.7G", enzakt_f, *(enzim_f+cella_f*enzaktnum_f+enzakt_f));
		}
		printf("\n");
	}
}

void konzolraSave (struct Save header_f) {
	printf("%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", header_f.version, header_f.tipus, header_f.num, header_f.rng, header_f.gen, header_f.ncol, header_f.nrow, header_f.eakt, header_f.met, header_f.repl, header_f.emax, header_f.trEE, header_f.trEK, header_f.kmin, header_f.kmax, header_f.phalal, header_f.claimE, header_f.sd, header_f.d, header_f.pmut);
}

void konzolraSaveHeader(void) {
	printf("version\ttipus\tnum\trng\tgen\tncol\tnrow\teakt\tmet\trepl\temax\ttrEE\ttrEK\tkmin\tkmax\tphalal\tclaimE\tsd\td\tpmut\n");
}