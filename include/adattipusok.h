#define BITNUM 10

struct Cella {
		double k;
		int szerk;
		unsigned int spec:BITNUM;
	};

struct Save {
	/*version 2 save
	 * 
	 * version=1
	 * tipus: milyen tipusu adat lesz a header utan mentve
	 * 	1: alapmatrix (struct Cella)
	 * 	2: enzimaktivitas matrix (double)
	 * 	3: struct Cella + double (egymás után)
	 * num: hany ciklus adata van osszesen elmentve a fajlba. Ha negativ szam (-x), akkor az azt jelenti, hogy meg mennyi (x) adat van elotte -> ossz adatszam az elso headerben lesz 
	 * gen: hanyadik generacio (ciklusszam)
	 * ncol: az alapmatrix sorainak szama
	 * nrow: az alapmatrix oszlopainak szama
	 * eakt: enzimaktivitasok szama
	 * met: metabolikus szomszedsag merete (0: vonNeumann, 1: 1es Moore, stb.)
	 * repl: replikacios szomszedsag merete
	 * emax: maximalis enzimatikus aktivitas
	 * trEE: a tradeoff az enzimatikus aktivitasok kozott
	 * trEK: a tradeoff az enzimatikus aktivitas es a replikacios koeff (k) kozott
	 * kmin: minimalis replikacio
	 * kmax: maximalis replikacio
	 * phalal: a replikator eltunesenek valoszinusege
	 * claimE: az uresen maradas Claim-je
	 * sd: az eloszlas szorasa
	 * d: diffuzio gyakorisaga
	 * pmut: a mutacio valoszinusege
	 * rng: randomszam generator tipusa
	 */
	int version, tipus, num;
	char rng[50];
	int gen, ncol, nrow, eakt;
	double met, repl, emax, trEE, trEK, kmin, kmax, phalal, claimE, sd, d, pmut;
};
 
struct SaveAll {
	struct Save h;
	struct Cella *m;
	double *e;
};