#include <mainheader.h>
#include <stdlib.h>
#include <math.h>
#include <randomgen.h>

int mutacioE1(struct Cella *matrix_f, double *enzim_f, int enzimaktszam_f, double sd_f, int target_f, double emax_f, double tradeoffEE_f, double reciprocTradeoffEE_f, double tradeoffEK_f, double reciprocTradeoffEK_f, double kmax_f, double kmin_f, double speclimit_f) {
	int *nezettseg_f, kivalasztott_f=0, enzimakt_f=1;
	double e_f, arany_f=0, esum_f=0, e1max_f=0, e2max_f=0, e2uj_f=0, e2ujmax_f=0;
	
	nezettseg_f= (int *) calloc(enzimaktszam_f, sizeof(int));
	
	kivalasztott_f=gsl_rng_uniform_int(r, enzimaktszam_f);
	*(nezettseg_f+kivalasztott_f)=1;
	e_f=gsl_ran_gaussian(r, sd_f) + (*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f));
/**/	//printf("\ne_f= %g", e_f);
	if (e_f<0) e_f = 0;
	if (e_f>emax_f) e_f = emax_f;
	*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f) = e_f;
/**/	//printf("\ne_f=%g", e_f);
	
	while(enzimakt_f<enzimaktszam_f){
		kivalasztott_f=gsl_rng_uniform_int(r, enzimaktszam_f);
		if(!*(nezettseg_f+kivalasztott_f)) {
			/*do {
				e2uj_f= *(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)* gsl_ran_gaussian(r, sd_f);
			}while ((tradeoff(e_f, emax_f, tradeoffEE_f) < e2uj_f) || (e2uj_f < 0));
			*/
			e2uj_f= *(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)+ gsl_ran_gaussian(r, sd_f);
/**/			//printf("\te2uj_f= %g", e2uj_f);
			e2ujmax_f= tradeoff(e_f, emax_f, tradeoffEE_f, reciprocTradeoffEE_f);
			if(e2ujmax_f < e2uj_f) e2uj_f= e2ujmax_f;
			if(e2uj_f<0) e2uj_f=0.0;
			
			*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)=e2uj_f;
/**/			//printf(", %g", e2uj_f);
			*(nezettseg_f+kivalasztott_f)=1;
			enzimakt_f++;
			
			//uj maximalis kivetitett enzimaktivitas kiszamitasa
			//CSAK itt: e_f=e1; *(matrix_f+hanyadik_f*enzimaktszam_f+kivalasztott_f)=e2
			if(e_f) arany_f=*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)/ e_f; // arany, meredekseg: e2/e1
			else arany_f= 0.0; //0-val osztas gazos
			esum_f=atfogoPit(*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f), e_f); //az enzimatikus aktivitasok metszespontja
			e1max_f= emax_f/pow((1+pow(arany_f, tradeoffEE_f)), reciprocTradeoffEE_f );
			e2max_f=tradeoff(e1max_f, emax_f, tradeoffEE_f, reciprocTradeoffEE_f);
			emax_f=atfogoPit(e1max_f, e2max_f);
			
			//elozo enzimakt megvaltoztat osszesitettre
			e_f=esum_f;
		}	
	}
	
	(*(matrix_f+target_f)).k = (kmax_f>kmin_f)?(kszamit(e_f, emax_f, tradeoffEK_f, reciprocTradeoffEK_f, kmax_f, kmin_f)* gsl_rng_uniform_pos(r)+kmin_f):kmin_f;
/**/	//printf(", %lf, emax_f=%lf, e2uj_f=%g, kmax=%g, kmin=%g, k=%lf, %lf", e_f, emax_f, e2uj_f, kmax_f, kmin_f, kszamit(e_f, emax_f, tradeoffEK_f, reciprocTradeoffEK_f, kmax_f, kmin_f), (*(matrix_f+target_f)).k);
	
	//szerk valasztas
	(*(matrix_f+target_f)).szerk= gsl_rng_uniform_int(r, enzimaktszam_f);
	
	//specificitas kiszamitasa
	specszamit(matrix_f, enzim_f, enzimaktszam_f, speclimit_f, target_f);
	
	
	free(nezettseg_f);
	return (0);
}

int mutacioE2(struct Cella *matrix_f, double *enzim_f, int enzimaktszam_f, double sd_f, int target_f, double emax_f, double tradeoffEE_f, double reciprocTradeoffEE_f, double tradeoffEK_f, double reciprocTradeoffEK_f, double kmax_f, double kmin_f, double speclimit_f) {
	  if (enzimaktszam_f==2) {
		  int kivalasztott_f=0, fixalt_f=0, i_f=0;
		  double arany_f=0.0, esum_f=0.0, e1max_f=0.0, e2max_f=0.0, emax2_f=0.0, e1_f=0.0, e0_f=0.0, emax_elm_f=0.0;
		  
		  //enzimaktivitasok meghat
		  kivalasztott_f=gsl_rng_uniform_int(r, enzimaktszam_f);
		  for(i_f=0;i_f<2;i_f++) {
			fixalt_f=(kivalasztott_f==0)?1:0;
			emax_elm_f=tradeoff(*(enzim_f+target_f*enzimaktszam_f+fixalt_f), emax_f, tradeoffEE_f, reciprocTradeoffEE_f);
			*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)=emax_elm_f*gsl_rng_uniform(r);
			kivalasztott_f=(kivalasztott_f==0)?1:0;
		  }
		  
		  //k meghat
		  if(kmax_f>kmin_f) {
			e0_f=*(enzim_f+target_f*enzimaktszam_f+0);
			e1_f=*(enzim_f+target_f*enzimaktszam_f+1);
			if(e0_f) arany_f=e1_f/e0_f; // arany, meredekseg
			else arany_f= 0.0; //0-val osztas gazos
			esum_f=atfogoPit(e0_f, e1_f); //az enzimatikus aktivitasok metszespontja
			e1max_f= emax_f/pow((1+pow(arany_f, tradeoffEE_f)), reciprocTradeoffEE_f );
			e2max_f=tradeoff(e1max_f, emax_f, tradeoffEE_f, reciprocTradeoffEE_f);
			emax2_f=atfogoPit(e1max_f, e2max_f);
			(*(matrix_f+target_f)).k = kszamit(esum_f, emax2_f, tradeoffEK_f, reciprocTradeoffEK_f, kmax_f, kmin_f)* gsl_rng_uniform_pos(r)+kmin_f;
		  }
		  else (*(matrix_f+target_f)).k=kmin_f;
		  
		  //szerk valasztas
		  (*(matrix_f+target_f)).szerk= gsl_rng_uniform_int(r, enzimaktszam_f);
		  
		  //specificitas kiszamitasa
		  specszamit(matrix_f, enzim_f, enzimaktszam_f, speclimit_f, target_f);
		  
		  
		  return (0);
	  }
	  else {
		  printf("rossz mutacios mechanizmus! Ez csak ketfele enzimre hasznalhato!");
		  return (1);
	  }
}

int mutacioE3(struct Cella *matrix_f, double *enzim_f, int enzimaktszam_f, double sd_f, int target_f, double emax_f, double tradeoffEE_f, double reciprocTradeoffEE_f, double tradeoffEK_f, double reciprocTradeoffEK_f, double kmax_f, double kmin_f, double speclimit_f) {
/*Tobbre is mukszik, egyenletes eloszlas
 * kivalaszt egy eakt -> mutaltat -> kivalaszt masik ea -> tradeoff-mutacio -> exponalas -> kivalaszt meg egy -> tradeoff-mutacio -> stb.
 */
	int *nezettseg_f, kivalasztott_f=0, enzimakt_f=1;
	double e_f, arany_f=0, esum_f=0, e1max_f=0, e2max_f=0, e2uj_f=0, e2ujmax_f=0;
	
	nezettseg_f= (int *) calloc(enzimaktszam_f, sizeof(int));
	
	kivalasztott_f=gsl_rng_uniform_int(r, enzimaktszam_f);
	*(nezettseg_f+kivalasztott_f)=1;
	e_f=gsl_rng_uniform(r) * emax_f;
/**/	//printf("\ne_f= %g", e_f);
//	if (e_f<0) e_f = 0;
//	if (e_f>emax_f) e_f = emax_f;
	*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f) = e_f;
/**/	//printf("\ne_f=%g", e_f);
	
	while(enzimakt_f<enzimaktszam_f){
		kivalasztott_f=gsl_rng_uniform_int(r, enzimaktszam_f);
		if(!*(nezettseg_f+kivalasztott_f)) {
			/*do {
				e2uj_f= *(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)* gsl_ran_gaussian(r, sd_f);
			}while ((tradeoff(e_f, emax_f, tradeoffEE_f) < e2uj_f) || (e2uj_f < 0));
			*/
			e2ujmax_f= tradeoff(e_f, emax_f, tradeoffEE_f, reciprocTradeoffEE_f);
			e2uj_f= gsl_rng_uniform(r) * e2ujmax_f;
/**/			//printf("\te2uj_f= %g", e2uj_f);
//			if(e2ujmax_f < e2uj_f) e2uj_f= e2ujmax_f;
//			if(e2uj_f<0) e2uj_f=0.0;
			
			*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)=e2uj_f;
/**/			//printf(", %g", e2uj_f);
			*(nezettseg_f+kivalasztott_f)=1;
			enzimakt_f++;
			
			//uj maximalis kivetitett enzimaktivitas kiszamitasa
			//CSAK itt: e_f=e1; *(matrix_f+hanyadik_f*enzimaktszam_f+kivalasztott_f)=e2
			if(e_f) arany_f=*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)/ e_f; // arany, meredekseg: e2/e1
			else arany_f= 0.0; //0-val osztas gazos
			esum_f=atfogoPit(*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f), e_f); //az enzimatikus aktivitasok metszespontja
			e1max_f= emax_f/pow((1+pow(arany_f, tradeoffEE_f)), reciprocTradeoffEE_f );
			e2max_f=tradeoff(e1max_f, emax_f, tradeoffEE_f, reciprocTradeoffEE_f);
			emax_f=atfogoPit(e1max_f, e2max_f);
			
			//elozo enzimakt megvaltoztat osszesitettre
			e_f=esum_f;
		}	
	}
	
	(*(matrix_f+target_f)).k = (kmax_f>kmin_f)?(kszamit(e_f, emax_f, tradeoffEK_f, reciprocTradeoffEK_f, kmax_f, kmin_f)* gsl_rng_uniform_pos(r)+kmin_f):kmin_f;
/**/	//printf(", %lf, emax_f=%lf, e2uj_f=%g, kmax=%g, kmin=%g, k=%lf, %lf", e_f, emax_f, e2uj_f, kmax_f, kmin_f, kszamit(e_f, emax_f, tradeoffEK_f, reciprocTradeoffEK_f, kmax_f, kmin_f), (*(matrix_f+target_f)).k);
	
	//szerk valasztas
	(*(matrix_f+target_f)).szerk= gsl_rng_uniform_int(r, enzimaktszam_f);
	
	//specificitas kiszamitasa
	specszamit(matrix_f, enzim_f, enzimaktszam_f, speclimit_f, target_f);
	
	
	free(nezettseg_f);
	return (0);
}

int mutacioE4(struct Cella *matrix_f, double *enzim_f, int enzimaktszam_f, double sd_f, int target_f, double emax_f, double tradeoffEE_f, double reciprocTradeoffEE_f, double tradeoffEK_f, double reciprocTradeoffEK_f, double kmax_f, double kmin_f, double speclimit_f) {
/*Tobbre is mukszik, egyenletes eloszlas, k nem utolso
 * kivalaszt egy eakt -> mutaltat -> kivalaszt masik ea -> tradeoff-mutacio -> exponalas -> kivalaszt meg egy -> tradeoff-mutacio -> stb.
 */
	int *nezettseg_f, kivalasztott_f=0, enzimakt_f=1, ksorszam_f= enzimaktszam_f + 1;
	double e_f, arany_f=0, esum_f=0, e1max_f=0, e2max_f=0, e2uj_f=0, e2ujmax_f=0;
	
	nezettseg_f= (int *) calloc(ksorszam_f, sizeof(int));
	
	kivalasztott_f=gsl_rng_uniform_int(r, ksorszam_f);
	*(nezettseg_f+kivalasztott_f)=1;
	e_f=gsl_rng_uniform(r) * emax_f;
//	printf("\ne_f= %g", e_f);
//	if (e_f<0) e_f = 0;
//	if (e_f>emax_f) e_f = emax_f;

	if(kivalasztott_f != enzimaktszam_f) {
		*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f) = e_f;
	}
	else {
		(*(matrix_f + target_f)).k = kskalaz(e_f, emax_f, kmax_f, kmin_f);
	}
//	printf("\ne_f=%g", e_f);
	
	while(enzimakt_f < ksorszam_f){
		kivalasztott_f=gsl_rng_uniform_int(r, ksorszam_f);
		if(!*(nezettseg_f+kivalasztott_f)) {
			/*do {
				e2uj_f= *(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)* gsl_ran_gaussian(r, sd_f);
			}while ((tradeoff(e_f, emax_f, tradeoffEE_f) < e2uj_f) || (e2uj_f < 0));
			*/
			if(kivalasztott_f != enzimaktszam_f) {
				e2ujmax_f= tradeoff(e_f, emax_f, tradeoffEE_f, reciprocTradeoffEE_f);
				e2uj_f= gsl_rng_uniform(r) * e2ujmax_f;	*(enzim_f+target_f*enzimaktszam_f+kivalasztott_f)=e2uj_f;
			}
			else {
				e2ujmax_f= tradeoff(e_f, emax_f, tradeoffEK_f, reciprocTradeoffEK_f);
				e2uj_f= gsl_rng_uniform(r) * e2ujmax_f;
				(*(matrix_f + target_f)).k = kskalaz(e2uj_f, emax_f, kmax_f, kmin_f);
			}
//			printf("\te2uj_f= %g", e2uj_f);
//			if(e2ujmax_f < e2uj_f) e2uj_f= e2ujmax_f;
//			if(e2uj_f<0) e2uj_f=0.0;

//			printf(", %g", e2uj_f);
			*(nezettseg_f+kivalasztott_f)=1;
			enzimakt_f++;
			
			//uj maximalis kivetitett enzimaktivitas kiszamitasa
			//CSAK itt: e_f=e1; *(matrix_f+hanyadik_f*enzimaktszam_f+kivalasztott_f)=e2
			if(e_f) arany_f = e2uj_f / e_f; // arany, meredekseg: e2/e1
			else arany_f= 0.0; //0-val osztas gazos
			esum_f=atfogoPit(e2uj_f, e_f); //az enzimatikus aktivitasok metszespontja
			
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
			e_f=esum_f;
		}	
	}
	
//	printf(", %lf, emax_f=%lf, e2uj_f=%g, kmax=%g, kmin=%g, k=%lf, %lf", e_f, emax_f, e2uj_f, kmax_f, kmin_f, kszamit(e_f, emax_f, tradeoffEK_f, reciprocTradeoffEK_f, kmax_f, kmin_f), (*(matrix_f+target_f)).k);
	
	//szerk valasztas
	(*(matrix_f+target_f)).szerk= gsl_rng_uniform_int(r, enzimaktszam_f);
	
	//specificitas kiszamitasa
	specszamit(matrix_f, enzim_f, enzimaktszam_f, speclimit_f, target_f);
	
	
	free(nezettseg_f);
	return (0);
}
