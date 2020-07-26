#include <stdio.h>
#include <adattipusok.h>
#define MUTACIO mutacioE4
#define TRANSITION transition3

//#define BITNUMFORMAT %(BITNUM)d
	
//feltoltes.c
void feltoltM(struct Cella *, double *, int, int, double, double, double, double, double, double, double, double);
int feltoltCella(struct Cella *, double *, int, int, double, double, double, double, double, double, double);

double tradeoff(double, double, double, double);
double atfogoPit(double, double);
double kszamit(double, double, double, double, double, double);
void specszamit(struct Cella *, double *, int, double, int);
double kskalaz(double, double, double, double);

/*//konzolra.c
void konzolraMatrixStruct(struct Cella *, int);
void konzolraMatrix(int *, int, int);
void konzolraMatrixD(double *, int, int);*/

//szomszed.c
int* metNeighInic(int, int, double);
int szomsz_meret(double);

//torus.c
int torus(int, int);

//metab.c
int extinkcio(struct Cella *, double *, int, int);
double metabolizmus(struct Cella *, double *, int *, int, int, double, int);
int replikacio (struct Cella *, double *, int, double, int, int);
void transition1 (struct Cella *, double *, int, int);
void transition2 (struct Cella *, double *, int, int);
void transition3 (struct Cella *, double *, int, int, int);

//mutacio.c
//int mutacioE1(struct Cella *, double *, int, double, int, double, double, double, double, double, double, double, double);
int mutacioE4(struct Cella *, double *, int, double, int, double, double, double, double, double, double, double, double);

//diffuzio.c
int diffTM (struct Cella *, double *, int, int, int);
void copypaste (struct Cella *, double *, int, int, int);

//bitmuveletek.c
int olvas (int, int);
int ir (int, int, int);
int kiirbit (int, int);

//statisztika.c
int alatt (int, int);
int faktorialis (int);
int hanyfajta (int);
double szorasD (int, double, double);

//kimenet.c
void tipusadatok(double *, struct Cella *, double *, int, int, int, int);
void kimenetOszlopnevek (FILE *,  int);
void kimenetTipusadat (FILE *, struct Cella *, double *, int, int, int);
void fajlbaCella (FILE *, struct Cella *, double *, int, int, int);
void fajlbaCellaHeaders (FILE *, int);

//save.c
int mentes(struct Cella *, double *, double *, int, int, int, FILE *, FILE *);
int mentesBin(struct Cella *, double *, struct Save, int, int, FILE *, FILE *, FILE *);
int binHeader(FILE *, struct Save, int, int);
int compareHeaders(struct Save, struct Save);
int writetoBin(struct Cella *, double *, struct Save, int, int, FILE *, FILE *);
