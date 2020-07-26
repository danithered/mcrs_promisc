#include <stdio.h>
#include <adattipusok.h>


//load.c
struct Save teker(FILE *, int);
struct SaveAll betoltes(char *, int);

//olvaso_kieg.c
double metabolizmus(struct Cella *, double *, int *, int, int, double, int);
int fajlbaMetab (FILE *, struct Cella *, double *, double, int, int, int);
int fajlbaSzomszTipusok (FILE *, struct Cella *, double, int, int, int);
int fajlbaAsszocTabla (FILE *, struct Cella *, double, int, int, int);

//szomszed.c
int szomsz_meret(double);
int* metNeighInic(int, int, double);

//statisztika.c
int hanyfajta(int);
