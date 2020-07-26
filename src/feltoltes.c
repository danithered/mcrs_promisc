#include <mainheader.h>
#include <stdlib.h>
#include <randomgen.h>
#include <math.h>


/*Alapmatrixot feltolt random ertekekkel.
!!!: Kell neki a struct Cella*/
void feltoltM(struct Cella *matrixCella_f, double *matrixE_f, int meret_f, int enzimaktszam_f, double feltolt_f, double kmin_f, double kmax_f, double emax_f, double tradeoffEE_f, double tradeoffEK_f, double speclimit_f, double sd_f) {
	      /*
	  * matrixCella_f: a feltoltendo alapmatrix.
	  * matrixE_f: a feloltendo enzimakt matrix
	  * meret_f: a feltoltendo matrix cellamerete (alapmartix cellamerete*enzimaktivitasok szama)
	  * enzimaktszam_f: hany enzimakt van osszesen
	  * feltolt_f: hany szazalekban van feltoltve a matrix
	  * kmin_f: k minimum erteke
	  */
	
	int cellaszam_f;
	double reciprocTradeoffEE_f=1/tradeoffEE_f, reciprocTradeoffEK_f=1/tradeoffEK_f; //a feltotltE fuggveny kimeneti vektora. Elso erteke: kumulalalt enzimaktivitas, második: emax_f
	
//	printf("\n%g\n", feltolt_f);
	
	for(cellaszam_f=0; cellaszam_f<meret_f; cellaszam_f++) {
		if(gsl_rng_uniform(r)<feltolt_f) {
			//enzimaktivitasok feltoltese
			feltoltCella(matrixCella_f, matrixE_f, enzimaktszam_f, cellaszam_f, emax_f, kmin_f, kmax_f, tradeoffEE_f, reciprocTradeoffEE_f, tradeoffEK_f, reciprocTradeoffEK_f);
			//*(matrixE_f+cellaszam_f*enzimaktszam_f+(rand()%enzimaktszam_f))=((double)rand()/RAND_MAX)*emax_f;
			
			//k generalasa
			//(*(matrixCella_f+cellaszam_f)).k= kszamit(*(eki_f+0), *(eki_f+1), tradeoffEK_f, kmax_f, kmin_f) *gsl_ran_gaussian(r, sd_f);
			//(*(matrixCella_f+cellaszam_f)).k= (kmax_f>kmin_f)?(kszamit(*(eki_f+0), *(eki_f+1), tradeoffEK_f, reciprocTradeoffEK_f, kmax_f, kmin_f)*gsl_rng_uniform_pos(r)+kmin_f):kmin_f;
			//if ((*(matrixCella_f+cellaszam_f)).k > 0) (*(matrixCella_f+cellaszam_f)).k+=kmin_f;
			//else (*(matrixCella_f+cellaszam_f)).k = kmin_f;
//				printf("\nkszamit: %g\tvegleges: %g", kszamit(*(eki_f+0), *(eki_f+1), tradeoffEK_f, kmax_f, kmin_f), (*(matrixCella_f+cellaszam_f)).k);
//				printf("\n%g", kszamit(0.6, 0.9, tradeoffEK_f, kmax_f, kmin_f));
//				printf("\n%g\t%g", *(eki_f+0), *(eki_f+1));
			
			//random enzimakt fold kivalasztasa
			(*(matrixCella_f+cellaszam_f)).szerk= gsl_rng_uniform_int(r, enzimaktszam_f);
//			printf("\nrandom struktura: %d", (*(matrixCella_f+cellaszam_f)).szerk);
			
			//specificitas beallit
			specszamit(matrixCella_f, matrixE_f, enzimaktszam_f, speclimit_f, cellaszam_f);
			
		}
	}
    return ;
}


// double* feltoltE(double *matrix_f, int enzimaktszam_f, int hanyadik_f, double emax_f, double tradeoff_f, double reciprocTradeoff_f) {
// 	int *nezettseg_f, kivalasztott_f=0, enzimakt_f=1;
// 	double *e_f, arany_f=0, esum_f=0, e1max_f=0, e2max_f=0;
// 	nezettseg_f= (int *) calloc(enzimaktszam_f, sizeof(int));
// 	e_f= (double *) calloc(2, sizeof(double)); //elso erteke: kumulalalt enzimaktivitas, második: emax_f
// 	
// 	kivalasztott_f=gsl_rng_uniform_int(r, enzimaktszam_f);
// 	*(nezettseg_f+kivalasztott_f)=1;
// 	*e_f=*(matrix_f+hanyadik_f*enzimaktszam_f+kivalasztott_f)=gsl_rng_uniform(r) *emax_f;
// 	
// 	while(enzimakt_f<enzimaktszam_f){
// 		kivalasztott_f=gsl_rng_uniform_int(r, enzimaktszam_f);
// 		if(!*(nezettseg_f+kivalasztott_f)) {
// 			*(matrix_f+hanyadik_f*enzimaktszam_f+kivalasztott_f)=tradeoff(*e_f, emax_f, tradeoff_f, reciprocTradeoff_f)*gsl_rng_uniform(r);
// 			*(nezettseg_f+kivalasztott_f)=1;
// 			enzimakt_f++;
// 			
// 			//uj maximalis kivetitett enzimaktivitas kiszamitasa
// 			//CSAK itt: e_f=e1; *(matrix_f+hanyadik_f*enzimaktszam_f+kivalasztott_f)=e2
// 			if (e_f) arany_f=*(matrix_f+hanyadik_f*enzimaktszam_f+kivalasztott_f)/ *e_f; // arany, meredekseg: e2/e1
// 			else arany_f=0;
// 			esum_f=atfogoPit(*(matrix_f+hanyadik_f*enzimaktszam_f+kivalasztott_f), *e_f); //az enzimatikus aktivitasok metszespontja
// 			e1max_f= emax_f/pow((1+pow(arany_f, tradeoff_f)), reciprocTradeoff_f );
// 			e2max_f=tradeoff(e1max_f, emax_f, tradeoff_f, reciprocTradeoff_f);
// 			emax_f=atfogoPit(e1max_f, e2max_f);
// 			
// 			//elozo enzimakt megvaltoztat osszesitettre
// 			*e_f=esum_f;
// 		}	
// 	}
// 	
// 	free(nezettseg_f);
// 	*(e_f+1) = emax_f;
// /**/	//printf("\nemax_f:")
//     return (e_f);
// }


int feltoltCella(struct Cella *matrix_f, double *enzim_f, int enzimaktszam_f, int hanyadik_f, double emax_f, double kmin_f, double kmax_f, double tradeoffEE_f, double reciprocTradeoffEE_f, double tradeoffEK_f, double reciprocTradeoffEK_f) {
	int *nezettseg_f, kivalasztott_f=0, enzimakt_f=1, ksorszam= enzimaktszam_f + 1;
	double *e_f, arany_f=0, esum_f=0, e1max_f=0, e2max_f=0, Value_f=0;
	nezettseg_f= (int *) calloc(ksorszam, sizeof(int));
	e_f= (double *) calloc(2, sizeof(double)); //elso erteke: kumulalalt enzimaktivitas, második: emax_f
	
	kivalasztott_f=gsl_rng_uniform_int(r, ksorszam);
//	printf("%d", kivalasztott_f);
	*(nezettseg_f+kivalasztott_f)=1;
	
	*e_f=gsl_rng_uniform(r) *emax_f;
	if(kivalasztott_f != enzimaktszam_f) {
		*(enzim_f+hanyadik_f*enzimaktszam_f+kivalasztott_f) = *e_f;
	}
	else {
		(*(matrix_f + hanyadik_f)).k = kskalaz(*e_f, emax_f, kmax_f, kmin_f);
//		printf("%d.cella: e_f= %g, k= %g\n", hanyadik_f, *e_f, (*(matrix_f + hanyadik_f)).k);
	}
	
	while(enzimakt_f < ksorszam){
		kivalasztott_f=gsl_rng_uniform_int(r, ksorszam);
//		printf("%d", kivalasztott_f);
		if(!*(nezettseg_f+kivalasztott_f)) {
			if(kivalasztott_f != enzimaktszam_f) {
				Value_f = *(enzim_f+hanyadik_f*enzimaktszam_f+kivalasztott_f) = tradeoff(*e_f, emax_f, tradeoffEE_f, reciprocTradeoffEE_f) * gsl_rng_uniform(r);
			}
			else {
				Value_f = tradeoff(*e_f, emax_f, tradeoffEK_f, reciprocTradeoffEK_f) * gsl_rng_uniform(r);
				(*(matrix_f + hanyadik_f)).k = kskalaz(Value_f, emax_f, kmax_f, kmin_f);
//				printf("%d.cella: e_f= %g, k= %g\n", hanyadik_f, *e_f, (*(matrix_f + hanyadik_f)).k);
			}
			
			*(nezettseg_f+kivalasztott_f)=1;
			enzimakt_f++;
			
			//uj maximalis kivetitett enzimaktivitas kiszamitasa
			//CSAK itt: e_f=e1; *(enzim_f+hanyadik_f*enzimaktszam_f+kivalasztott_f)=e2
			if (e_f) arany_f= Value_f / *e_f; // arany, meredekseg: e2/e1
			else arany_f=0;
			esum_f=atfogoPit(Value_f, *e_f); //az enzimatikus aktivitasok metszespontja
			
			if(kivalasztott_f != enzimaktszam_f) {
				e1max_f= emax_f/pow((1+pow(arany_f, tradeoffEE_f)), reciprocTradeoffEE_f );
				e2max_f=tradeoff(e1max_f, emax_f, tradeoffEE_f, reciprocTradeoffEE_f);
			}
			else {
				e1max_f= emax_f/pow((1+pow(arany_f, tradeoffEK_f)), reciprocTradeoffEK_f );
				e2max_f=tradeoff(e1max_f, emax_f, tradeoffEK_f, reciprocTradeoffEK_f);
			}
			emax_f=atfogoPit(e1max_f, e2max_f);
			
			//elozo enzimakt megvaltoztat osszesitettre
			*e_f=esum_f;
		}	
	}
	
	free(nezettseg_f);
	*(e_f+1) = emax_f;
//	printf("\nemax_f:")
//	printf("\n");
    return (0);
}


double tradeoff(double ertek1_f, double maxertek_f, double tradeoff_f, double reciprocTradeoff_f) {
	if (ertek1_f < maxertek_f) return(pow(pow(maxertek_f, tradeoff_f)-pow(ertek1_f, tradeoff_f), reciprocTradeoff_f));
	else return(0.0);
}

double atfogoPit(double befogo1_f, double befogo2_f) {
	return(pow((pow(befogo1_f, 2)+pow(befogo2_f, 2)), 0.5));
}

double kszamit(double e_f, double emax_f, double tradeoff_f, double reciprocTradeoff_f, double kmax_f, double kmin_f) {
	if (kmax_f>kmin_f) return(tradeoff(e_f, emax_f, tradeoff_f, reciprocTradeoff_f) * ((kmax_f-kmin_f)/emax_f));
	else return(kmin_f);
}

double kskalaz(double k_f, double emax_f, double kmax_f, double kmin_f) {
//	printf("in= %g\temax= %g\tkmin= %g\tkmax= %g\n", k_f, emax_f, kmin_f, kmax_f);
	if (kmax_f>kmin_f && k_f) return( (k_f * ((kmax_f - kmin_f)/emax_f)) + kmin_f );
	else return(kmin_f);
}

void specszamit(struct Cella *matrix_f, double *enzim_f, int enzimaktszam_f, double speclimit_f, int target_f) {
	int enzimakt_f;
	
	(*(matrix_f+target_f)).spec=0;
	for(enzimakt_f=0; enzimakt_f<enzimaktszam_f; enzimakt_f++) {
		if(*(enzim_f+target_f*enzimaktszam_f+enzimakt_f) > speclimit_f) (*(matrix_f+target_f)).spec = ir((*(matrix_f+target_f)).spec, enzimakt_f, 1);
	}
}