#include <mainheader.h>
#include <stdlib.h>
#include <math.h>
#include <randomgen.h>

int extinkcio(struct Cella *matrix_f, double *matrixE_f, int enzimaktszam_f, int target_f) {
	/* matrix_f: alapmatrix, van benne: k, szerk
	 * matrixE_f: enzimaktivitas matrix
	 * enzimaktszam_f: az enzimaktivitasok szama
	 * target_f: a torlendo cella szama a matrixban
	 */
	int enzakt_szamlalo_f=0;
	
	(*(matrix_f+target_f)).k=0;
	(*(matrix_f+target_f)).szerk=-1;
	(*(matrix_f+target_f)).spec=0;
	for(enzakt_szamlalo_f=0; enzakt_szamlalo_f<enzimaktszam_f; enzakt_szamlalo_f++) {
		*(matrixE_f+enzimaktszam_f*target_f+enzakt_szamlalo_f)=0;
	}
	
	return (0);
}

//claim

double metabolizmus(struct Cella *matrix_f, double *enzim_f, int *met_szomszedsag_f, int szomsz_cellaszam_f, int enzimaktszam_f, double reciprocEnzimaktszam_f, int cella_f) {
/**/		//printf("\nenzimaktszam=%d", enzimaktszam_f);
	int szomszed_f, nezett_f=0, enzakt_f, dominans_f;
	double *enzimsum_f, metab_f=1;
	
	/* szomszed_f: a target (cella_f) szomszedjanak sorszama/ hanyadik szomszedot nezzuk eppen
	 * nezett_f: hanyadik cellat nezzuk eppen
	 * enzakt_f: a nezett enzakt szamlaloja
	 * enzimsum_f: az adott enzimaktivitasok osszege
	 * metab_f: M, azaz a metabolizmus erteke az adott cellara
	 */
	
	enzimsum_f=(double*) calloc(enzimaktszam_f, sizeof(double));
	
/**/	//printf("\ncella: %d (szerk, spec, eakt)\n", cella_f);
	for(szomszed_f=0; szomszed_f < szomsz_cellaszam_f; szomszed_f++) {
		nezett_f=*(met_szomszedsag_f+cella_f*szomsz_cellaszam_f+szomszed_f);
		if (!(*(matrix_f+nezett_f)).k) continue;
		//*(enzimsum_f+(*(matrix_f+nezett_f)).k) = *(enzim_f+nezett_f*enzimaktszam_f+(*(matrix_f+nezett_f)).k);
		TRANSITION (matrix_f, enzim_f, enzimaktszam_f, 0.5, nezett_f);
		dominans_f=(*(matrix_f+nezett_f)).szerk;
		*(enzimsum_f+dominans_f) += *(enzim_f+nezett_f*enzimaktszam_f+dominans_f);
/**/		//printf("\nnezett=%d\tdominans=%d", nezett_f, dominans_f);
/**/		//printf(" (%d, %d, %g)", (*(matrix_f+nezett_f)).szerk, (*(matrix_f+nezett_f)).spec, *(enzim_f+nezett_f*enzimaktszam_f+dominans_f));
	}
	for(enzakt_f=0; enzakt_f<enzimaktszam_f; enzakt_f++){
		metab_f *= *(enzimsum_f+enzakt_f);
/**/		//printf("\nenzimsum= %g", *(enzimsum_f+enzakt_f));
	}
/**/	//printf("\n\nmetab_f=%g\t1/enzimaktszam_f=%g", metab_f, (1/ (double) enzimaktszam_f));
/**/	//printf("\nmetab_f=%g, reciprocEnzimaktszam_f=%g", metab_f, reciprocEnzimaktszam_f);
	
	if (metab_f>0) metab_f=pow(metab_f, reciprocEnzimaktszam_f);
	
	free (enzimsum_f);
/**/		//printf("\nvegso metab_f=%g", metab_f);
	return (metab_f);
}

int replikacio (struct Cella * matrix_f, double *enzimakt_f, int enzaktszam_f, double speclimit_f, int templat_f, int target_f) {
	if ((*(matrix_f+target_f)).k) return (-1); //ha esetleg megsem lenne ures
	
	int masolandoEnzimakt=0;
	
	//k masolasa
	(*(matrix_f+target_f)).k = (*(matrix_f+templat_f)).k;
	
	//szerk atmasolasa <- mutaciokor ugyis megvaltozik
	(*(matrix_f+target_f)).szerk = (*(matrix_f+templat_f)).szerk;
	
	//enzimakt atmasolasa
	for(masolandoEnzimakt=0; masolandoEnzimakt<enzaktszam_f; masolandoEnzimakt++) {
		*(enzimakt_f+target_f*enzaktszam_f+masolandoEnzimakt) = *(enzimakt_f+templat_f*enzaktszam_f+masolandoEnzimakt);
	}
	
	//spec szamitasa
	specszamit(matrix_f, enzimakt_f, enzaktszam_f, speclimit_f, target_f);
	
	return(0);
}

void transition1 (struct Cella * matrix_f, double *enzim_f, int enzaktszam_f, int cella_f) {
	int jelenlegi_f= (*(matrix_f+cella_f)).szerk, masik_f=0;
	double pTransition_f=0, e1_f, e2_f;
	
	//kivalaszt masik aktivitas
	do {
		masik_f=gsl_rng_uniform_int(r, enzaktszam_f);
	} while(masik_f==jelenlegi_f);
	e2_f= *(enzim_f + enzaktszam_f*cella_f +masik_f);
	e1_f= *(enzim_f + enzaktszam_f*cella_f +jelenlegi_f);
	
	//pTransition_f
	pTransition_f= 1 - fabs((e1_f-e2_f)/((e1_f>e2_f)?e1_f:e2_f));
	
	//dontes
	if (gsl_rng_uniform(r) < pTransition_f) (*(matrix_f+cella_f)).szerk = masik_f;
}

void transition2 (struct Cella * matrix_f, double *enzim_f, int enzaktszam_f, int cella_f) {
/*TOROTT PALCA
 */
	int szamlalo_f=0, enzakt_f=(*(matrix_f+cella_f)).szerk;
	double esum_f=0.0, random_f=0.0, valaszto_f=0.0;
	int *nezettseg_f;
	
	nezettseg_f= (int *) calloc(enzaktszam_f, sizeof(int));
	
	for (szamlalo_f=0; szamlalo_f<enzaktszam_f; szamlalo_f++) {
		esum_f += *(enzim_f+cella_f*enzaktszam_f+szamlalo_f);
	}
	
	random_f = gsl_rng_uniform(r); 
	szamlalo_f=0;
	while (szamlalo_f <enzaktszam_f) {
		enzakt_f = gsl_rng_uniform_int(r, enzaktszam_f);
		if(! *(nezettseg_f+enzakt_f)) {
			valaszto_f += *(enzim_f+cella_f*enzaktszam_f+enzakt_f)/esum_f;
			if(random_f < valaszto_f) break;
			*(nezettseg_f+enzakt_f)=1;
			szamlalo_f++;
		}
		
	}
	(*(matrix_f+cella_f)).szerk = enzakt_f;
	
	free(nezettseg_f);
}

void transition3 (struct Cella * matrix_f, double *enzim_f, int enzaktszam_f, int hatar_f, int cella_f) {
/*PROD(E)^(1/enzaktszam)/max(E)
 */
	int enzaktReal_f=0, szamlalo_f=0, enzakt_f=(*(matrix_f+cella_f)).szerk;
	double esum_f=0.0, eprod_f=1.0, maxE_f=0.0, random_f=0.0, valaszto_f=0.0, pTransition_f=0.0;
	int *nezettseg_f;
	
	nezettseg_f= (int *) calloc(enzaktszam_f, sizeof(int));
	
	//prod, szum, max, n kiszamitasa
	for (szamlalo_f=0; szamlalo_f<enzaktszam_f; szamlalo_f++) {
		if (*(enzim_f+cella_f*enzaktszam_f+szamlalo_f) > hatar_f) { /*ha van jelentos enzimaktivitas*/
			//a hatarerteknel nagyobb enzimaktivitasok szamanak meghatarozasa
			enzaktReal_f++;
			//maximum
			if(*(enzim_f+cella_f*enzaktszam_f+szamlalo_f) > maxE_f) maxE_f = *(enzim_f+cella_f*enzaktszam_f+szamlalo_f);
			//produktum
			eprod_f *= *(enzim_f+cella_f*enzaktszam_f+szamlalo_f);
			//nezettseg, sum <- ezeket csak akkor, ha nem aktualis enzimaktivitas 
			if (szamlalo_f != enzakt_f) {
				*(nezettseg_f+szamlalo_f)=1;
				esum_f += *(enzim_f+cella_f*enzaktszam_f+szamlalo_f);
			}
		}
	}
//	printf("enzaktReal_f: %d, maxE_f: %g, prod: %g, sum: %g\n", enzaktReal_f, maxE_f, eprod_f, esum_f);
	
	//pTransition_f meghatározasa
	if(maxE_f) pTransition_f = pow(eprod_f, 1/(double)enzaktReal_f)/maxE_f;
//	printf("pTransition_f: %g\n", pTransition_f);
	  
	//dontes -> kicsavarodik-e
	if(gsl_rng_uniform(r)<pTransition_f) {
		//kicsavarodik -> valaszt egyet a tobbi enzimaktivitas kozul
		random_f = gsl_rng_uniform(r); 
		if(*(enzim_f+cella_f*enzaktszam_f+enzakt_f) > hatar_f) enzaktReal_f--;
		szamlalo_f=0;
		while (szamlalo_f < enzaktReal_f) {
			enzakt_f = gsl_rng_uniform_int(r, enzaktszam_f);
			if(*(nezettseg_f+enzakt_f)) {
				valaszto_f += *(enzim_f+cella_f*enzaktszam_f+enzakt_f)/esum_f;
//				printf("valaszto_f: %g, random_f: %g\n", valaszto_f, random_f);
				if(random_f < valaszto_f) break;
				*(nezettseg_f+enzakt_f)=0;
				szamlalo_f++;
			}
		}
		(*(matrix_f+cella_f)).szerk = enzakt_f;
	}
	
	free(nezettseg_f);
}