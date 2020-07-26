#include <math.h>

int alatt (int osszes_f, int kivalasztott_f) {
	return(faktorialis(osszes_f)/ (faktorialis(kivalasztott_f)*faktorialis(osszes_f-kivalasztott_f)) );
}

int faktorialis (int szam_f) {
/* egy szam faktorialisat szamolja ki
 * 
 * szam_f: a szam, aminek a faktorialisat szamitjuk
 * szorzat_f: az eredmeny, az eddigi szamok szorzata 
 * szorzo_f: az a szam, amivel eppen megszorozzuk a szorzatot
 */
  
	int szorzat_f=1, szorzo_f;
	for(szorzo_f=2; szorzo_f<=szam_f; szorzo_f++) {
		szorzat_f*=szorzo_f;
	}
	
	return(szorzat_f);
}

int hanyfajta(int enzaktszam_f) {
  /* azt szamolja ki, hogy enzimaktivitasok szerint osszesen hanyfajta enzim lehetseges: parazita, specialista es generalista osszesen
   * 
   * enzaktszam_f: hanyfajta enzimaktivitas van
   * enzakt_f: szamlalo, korbemegy az enzimaktivitasokon
   * eredmeny_f: a permutaciok osszege
   * szamitas_max_f: ameddig maximalisan megeri ciklusban kiszamitani a permutaciokat
   */
	int enzakt_f=0, eredmeny_f=0, szamitas_max_f;
	
	//ha igaz -> paros -> asszimetrikus, mert van, h 0 enzimakt van, ezert az enzaktszam fele minusz 1-ig eri meg
	szamitas_max_f= ((enzaktszam_f%2)==0)?((enzaktszam_f/2) -1):((int)((double)enzaktszam_f/2 -0.5));
		
	for (enzakt_f=0; enzakt_f <= szamitas_max_f; enzakt_f++) {
		eredmeny_f+=alatt(enzaktszam_f, enzakt_f);
	}
	eredmeny_f*=2;
	if ((enzaktszam_f%2)==0) eredmeny_f+=alatt(enzaktszam_f, (szamitas_max_f+1));
	
	return(eredmeny_f);
}

double szorasD (int db_f, double sum_f, double negyzetSum_f) {
	if(db_f>1) return(sqrt(fabs(db_f*negyzetSum_f-sum_f*sum_f)/(db_f*(db_f-1))));
	else return (0); 
}