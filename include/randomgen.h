#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

time_t timer;
gsl_rng * r;


/*A main-be a kovetkezo sorokat be kell tenni az elejere es a vegere:*/
/*
//randomszam generator inicializalasa
	r = (gsl_rng *) gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r, time(&timer));
//randomszam generator lezarasa
	gsl_rng_free(r);
*/