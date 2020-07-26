/*
MCRS, egyed alapu, tobb enzimatikus aktivitas, mutacio
Ez egy sejtautomata, ami az MCRS-n alapul. Egyed alapu (nem szekvencia alapu) replikatorok, amiknek tobb meg nem nevezett enzimatikus aktivitasa van, ami mind kell.

Argomentumok: 21 (ncol nrow enzakt_num met_neigh_meret repl_neigh_meret ciklusszam mintavetel_gyak cellamintavetel_gyak alapfeltoltes emax tradeoffEE tradeoffEK kmin kmax phalal claimEmpty sd spec_limit diffuzioGyak pMutacio ment) Lsd.:Argomentumok beolvasasa resz
test: 10 10 2 0 0 10 2 2 0.8 1 0.2 0.4 2 4 0.1 0.001 0.1 0.01 0.01 0 4 nev

main visszateresi ertekek:
    0: lefutott
    1: matrix lefoglalasi hiba
    2: replikacio hiba
    3: tul keves argomentum
    4: ilyen nevu fajl mar van -> futtasd ujra
    5: meghalt
    6: mentett allapot betolese nem sikerult
*/

#include <mainheader.h>
#include <stdlib.h>
#include <randomgen.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <load.h>




int main(int argc, char *argv[]) {
	//randomszam generator inicializalasa
	r = (gsl_rng *) gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r, time(&timer));
	/*gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r, time(&timer));*/

	//Argomentumok beolvasasa
	int argnum=21, argkiir=0;
	if (argc < argnum) {
		printf("tul keves argomentum: %d szukseges\n", argnum);
		return(3);
	}
	int ncol=atoi(argv[1]), nrow=atoi(argv[2]), enzakt_num=atoi(argv[3]), ciklusszam=atoi(argv[6]),  mintavetel_gyak=atof (argv[7]), cellamintavetel_gyak=atof (argv[8]), ment=atoi(argv[21]);
	double met_neigh_meret=atof(argv[4]), repl_neigh_meret=atof(argv[5]), alapfeltoltes=atof (argv[9]), emax=atof (argv[10]), tradeoffEE=atof (argv[11]), tradeoffEK=atof (argv[12]), kmin=atof (argv[13]), kmax=atof (argv[14]), phalal=atof (argv[15]), claimEmpty=atof (argv[16]), sd=atof (argv[17]), spec_limit=atof (argv[18]), diffuzioGyak=atof (argv[19]), pMutacio=atof (argv[20]);

	/* ncol: alapmatrix oszlopainak szama
	 * nrow: alapmatrix sorainak szama
	 * enzakt_num: az enzimaktivitasok szama
	 * met_neigh_meret: a metabolikus szomsz. merete: ha 0 -> von Neumann; ha több -> Moore szomsz
	 * repl_neigh_meret: a replikacios szomsz. merete: ha 0 -> von Neumann; ha több -> Moore szomsz
	 * ciklusszam: hany iteracios ciklus fusson le
	 * mintavetel_gyak: milyen gyakran vegyen mintat. Ha 0, akkor nem vesz, ha 1 akkor mindig, ha pl 3, akkor minden harmadik ciklusban
	 * ment: van -e a vegen mentes, vagy nincs
	 * 	0:nincs mentes
	 * 	1:van mentes 2 csv-be + rng mentes
	 * 	2:van mentes matrix.bin -be es enzim.bin -be +rng mentes
	 *	3:van mentes save.bin -be +rng mentes
	 *	4:van mentes 1 bin-be, a cellainfok kiirasa is ugyanabba a binbe megy +rng mentes  
	 * alapfeltoltes: a matrix hany szazaleka legyen alapbol feltoltve (random)
	 * emax: a maximalis enzimaktivitas merteke
	 * tradeoffEE: enzimaktivitas es enzimaktivitas kozotti tradeoff parameter
	 * tradeoffEK: enzimaktivitas es replikacios koefficiens kozotti tradeoff parameter
	 * kmin: minimalis k ertek (replikacios koefficiens)
	 * kmax: maximalis k ertek (replikacios koefficiens)
	 * phalal: az extinkcio valoszinusege
	 * claimEmpty: az uresen maradas ~valsege
	 * sd: a Gaussi eloszlas szorasa
	 * spec_limit: a minimum hatar, amin felul egy enzimaktivitas szamit
	 * diffuzioGyak: 
	 * pMutacio: replikaciot koveto mutacio valsege
	 */
	int load=0;
//	printf("load= %d, azon: %s, ment= %d\n", load, argv[argnum+1], ment);
	struct SaveAll betolt;
	if(ment == -1) {
		if((argc-1)>(argnum+3)) load=2;
		else load = 1;
		ment = 4;
	}
	if (load){ //mentett allapot visszaallitasa
		//matrix visszaallitasa
		betolt = betoltes(argv[argnum+2], -1);
		ncol = betolt.h.ncol;
		nrow = betolt.h.nrow;
		enzakt_num = betolt.h.eakt;
		met_neigh_meret = betolt.h.met;
		repl_neigh_meret = betolt.h.repl;
		emax = betolt.h.emax;
		tradeoffEE = betolt.h.trEE;
		tradeoffEK = betolt.h.trEK;
		kmin = betolt.h.kmin;
		kmax = betolt.h.kmax;
		phalal = betolt.h.phalal;
		claimEmpty = betolt.h.claimE;
		sd = betolt.h.sd;
		diffuzioGyak = betolt.h.d;
		pMutacio = betolt.h.pmut;
	}
	
	
	//Valtozok deklaralasa
	int meret=ncol*nrow, ciklus=0, iter=0, cella=0, met_neigh_cellaszam=0, repl_neigh_cellaszam=0, num_repl_neigh=0, nezett=0, templat=0, enztipus_num=0, sorhossz=1+2*(enzakt_num+1), replikator_num=1, enzakt;
	double knezett=0, claimSum=0, valaszto=0, claimJelolt=0, diffuzio_szam=0, reciprocTradeoffEE=1/tradeoffEE, reciprocTradeoffEK=1/tradeoffEK, reciprocEnzakt_num=1/(double)enzakt_num;
	char mappa[]="OUT/\0", mentesmappa[]="OUT/save/\0", fajlnev[50]="\0", kezdet[30]="\0", csvname[50]={0}, cellafajlnev[50]="\0", savetoR[50]="\0", savetoRs[51]="\0", savetoData[50]="\0", savetoE[50]="\0", azon[10]="\0";
	struct Save save;
	
//	printf("reciprocenzakt_num= %g\n", reciprocEnzakt_num);


	/*
	 * meret= ncol*nrow, az alapmatrix cellaszama
	 * ciklus: hanyadik iteracios ciklusnal jarunk (egy cikluson belul meret szamu iteracios lepes van)
	 * iter: egy cikluson belul hanyadik iteracios lepesnel tartunk
	 * cella: az epp nezett cella szama
	 * met_neigh_cellaszam: hany cellabol all a metabolikus szomszedsag
	 * repl_neigh_cellaszam: hany cellabol all a replikacios szomszedsag
	 * num_repl_neigh: hanyadik replikacios szomszedot nezzuk
	 * nezett: a replikacios szomszed, aminek a claim-jet eppen szamoljuk
	 * templat:
	 * enztipus_num: hanyfajta enzim van
	 * sorhossz: az adatok matrix soranak hossza, megmondja, hogy hany adat tartozik egytipusu replikatorhoz (dbszam, k, e1...eN atlag, szoras)
	 * replikator_num: van e elo replikator a rendszerben
	 * enzakt: szamlalo, vegigmegy az enzimaktivitasokon
	 * knezett: az eppen aktualisan nezett replikator k-ja
	 * claimSum: a replikacios szomszedsagban levo replikatorok claim-jenek osszege + claimEmpty
	 * valaszto:
	 * claimJelolt:
	 * */
	
	//mentesi adatok megadasa
	if(ment>1) {
		save.version=1; save.gen=0, save.num=1; save.ncol=ncol; save.nrow=nrow; save.eakt=enzakt_num; save.met=met_neigh_meret; save.repl=repl_neigh_meret; save.emax=emax; save.trEE=tradeoffEE; save.trEK=tradeoffEK; save.kmin=kmin; save.kmax=kmax; save.phalal=phalal; save.claimE=claimEmpty; save.sd=sd; save.d=diffuzioGyak; save.pmut=pMutacio;
		memset(save.rng, '\0', sizeof(save.rng));
		strcpy(save.rng, gsl_rng_name (r));
		if(ment==2) save.tipus=1;
		else save.tipus=3;
	}
	 
	//Pointerek deklaralasa
	struct Cella *matrix;
	double *enzim, *claimek, *adatok;
	int *met_neigh, *repl_neigh;
	FILE *output, *fp, *cellak, *rngSave, *rngSaveStart, *dataSave, *enzSave, *mentettRng;
	struct tm *ido;
	struct stat st = {0}, stCsv = {0};


	/*
	* matrix: struct: a matrix cellaira vonatkozo adatok
	    * 		matrix.k: a cellaban levo enzim replikacios koefficiense. Ha 0 akkor nincs benne replikator
	    * 		matrix.szerk: a cellaban levo enzim szerkezetenek szama -> melyik enzimfunkciot vegzi el eppen
	* enzim: enzimatikus aktivitasokat tartalmazo vektor, meret*enzakt_num meretu. (Forma: enzim1(1.funkció), e1(2), ..., e1(enzakt_num), e2(1), ...)
	* claimek: a repl szomszedsagban levo replikatorok claim-je (C)
	* adatok: a kiszamolt adatok: tipusok dbszama, k, e0...eN ertekeinek atlaga, varinaciaja
	* met_neigh: metabolikus szomszédság mátrix
	* repl_neigh: replikacios szomszédság mátrix
	*/


	//valtozok kiszamitasa
	met_neigh_cellaszam=szomsz_meret(met_neigh_meret);
	repl_neigh_cellaszam=szomsz_meret(repl_neigh_meret);
	enztipus_num=hanyfajta(enzakt_num);
//	printf("\nenzimtipus: %d db\n", enztipus_num);
	

	//matrixok lefoglalasa
	matrix=(struct Cella*) calloc(meret, sizeof(struct Cella));
	enzim=(double *) calloc(meret*enzakt_num, sizeof(double));
	claimek=(double *) calloc((repl_neigh_cellaszam+1), sizeof(double));
	adatok= (double *) calloc((enztipus_num *sorhossz), sizeof(double));

	if(!matrix||!enzim ||!claimek ||!adatok) {
            printf("nem lehet lefoglalni egy matrix-ot(main)\n");
            return(1);
	}

	
	//szovegmuveletek
	ido = localtime( &timer );
	if (argc >= argnum+1) {
		strcpy(azon, argv[argnum+1]);
	}
	else {
		sprintf(azon, "%03d\0", gsl_rng_uniform_int(r, 1000));
	}
	sprintf(kezdet, "%02d.%02d_%02d:%02d:%02d(%s)\0", (*ido).tm_mon, (*ido).tm_mday, (*ido).tm_hour, (*ido).tm_min, (*ido).tm_sec, azon);
	sprintf(csvname, "%s.csv\0", kezdet);
	sprintf(fajlnev, "%s%s\0", mappa, csvname);
	if(load==2) sprintf(fajlnev, "%s\0", argv[argnum+4]);

	sprintf(savetoR, "%ssaveR%s.bin\0", mentesmappa, kezdet);
	sprintf(savetoRs, "%ssaveRs%s.bin\0", mentesmappa, kezdet);
	sprintf(cellafajlnev, "%scells%s\0", mappa, csvname);
	
	switch(ment)
	{
	  case 1:
		sprintf(savetoData, "%ssaveD%s\0", mentesmappa, csvname);
		break;
	  case 2:
		sprintf(savetoE, "%senzim%s.bin\0", mentesmappa, kezdet);
		sprintf(savetoData, "%smatrix%s.bin\0", mentesmappa, kezdet);	
		break;
	  case 3:
		sprintf(savetoData, "%ssaveD%s.bin", mentesmappa, kezdet);
		break;
	  case 4:
		sprintf(savetoData, "%sdata%s.bin", mappa, kezdet);
		break;
	}
	
	//log file kezd
	if (stat("OUT", &st) == -1) mkdir("OUT", 7777);
	if (stat("OUT/save", &st) == -1) mkdir("OUT/save", 7777);
//	printf("kimenet innentol atiranyitva log-ba\n");
	
	fp = freopen("OUT/log", "a", stdout);
	if(load) {
		printf("\n$%s\njob REstarted from generation %d with parameters:\n", fajlnev, betolt.h.gen);
		printf("%d %d %d %g %g %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %d ", ncol, nrow, enzakt_num, met_neigh_meret, repl_neigh_meret, ciklusszam, mintavetel_gyak, cellamintavetel_gyak, alapfeltoltes, emax, tradeoffEE, tradeoffEK, kmin, kmax, phalal, claimEmpty, sd, spec_limit, diffuzioGyak, pMutacio, -1);
		for (argkiir = argnum + 1; argkiir < argc; argkiir++) {
			printf("%s ", argv[argkiir]);
		}
	}
	else {
		printf("\n$%s\njob started with parameters:\n",csvname);
		for (argkiir=1; argkiir < argc; argkiir++) {
			printf("%s ", argv[argkiir]);
		}
	}
	printf("\n");
	
	//Kimenet fajlok nyitasa
	if ((stat(fajlnev, &stCsv) != -1) && load != 2) {
		printf("\n%s ilyen nevu fajl mar van!\n", fajlnev);
		return(4);
	}
	
	output = fopen(fajlnev, "a");
	
	switch(ment)
	{
	  case 0:
		cellak = fopen(cellafajlnev, "a");
	  case 1:
		cellak = fopen(cellafajlnev, "a");
		rngSave = fopen(savetoR, "a");
		rngSaveStart = fopen(savetoRs, "a");
		dataSave = fopen(savetoData, "a");
		break;
	  case 2:
		cellak = fopen(cellafajlnev, "a");
		rngSave = fopen(savetoR, "a");
		rngSaveStart = fopen(savetoRs, "a");
		dataSave = fopen(savetoData, "ab");
		enzSave =fopen(savetoE, "ab");
		break;
	  case 3:
		cellak = fopen(cellafajlnev, "a");
		rngSave = fopen(savetoR, "a");
		rngSaveStart = fopen(savetoRs, "a");
		dataSave = fopen(savetoData, "ab");
	  case 4:
		rngSave = fopen(savetoR, "a");
		rngSaveStart = fopen(savetoRs, "a");
		dataSave = fopen(savetoData, "wb+");
		break;
	}
//	printf("fajlok megnyitva\n");
//	printf("argnum=%d\n", argnum);

	
	/******************************
	****A LENYEG****
	*******************************/
	
	//Matrixok feltoltese
	if (load){ //mentett allapot visszaallitasa
		//randomszam generator betoltese
//		printf("rng betoltese fajlbol: %s\n", argv[argnum+3]);
		mentettRng = fopen(argv[argnum+3], "rb");
		if( betoltesRng(mentettRng) ) {
			printf("randomszam generator visszallittasa fajlbol nem sikerult!");
			return(6);
		}
//		printf("r fajlbol (%s) betoltve\n", argv[argnum+3]);
		
		//matrix visszaallitasa
		matrix = betolt.m;
		enzim = betolt.e;
//		printf("matrixok fajlbol (%s) betoltve\n", argv[argnum+2]);
	}
	else { //ha nem mentett allapot visszaallitasa
		feltoltM(matrix, enzim, meret, enzakt_num, alapfeltoltes, kmin, kmax, emax, tradeoffEE, tradeoffEK, spec_limit, sd);
	}
//	konzolraMatrixStruct(matrix, meret);
//	printf("\nenzimakt:\n");
//	konzolraMatrixD(enzim, enzakt_num, meret);
//	printf("\nmeret: %d\n", meret);
//	konzolraCella(matrix, enzim, enzakt_num, meret);
	
	//szomszedsag matrixok feltoltese
	met_neigh= (int *)metNeighInic(meret, ncol, met_neigh_meret);
	repl_neigh= (int *)metNeighInic(meret, ncol, repl_neigh_meret);
//	konzolraMatrix(met_neigh, met_neigh_cellaszam, meret);
	
	
	//Kimenet
	tipusadatok(adatok, matrix, enzim, enzakt_num, meret, enztipus_num, sorhossz); //kezdeti allapot leolvasasa
	if(load != 2) { //ha nem mar meglevo fajlba ment
		kimenetOszlopnevek (output, enzakt_num);
		kimenetTipusadat (output, matrix, adatok, enztipus_num, sorhossz, 0);
//		printf("nem letezo fajlba iras, headerrel: %s\n", fajlnev);
	}
	else {
		kimenetTipusadat (output, matrix, adatok, enztipus_num, sorhossz, betolt.h.gen); //ha mar meglevo fajlba ment
//		printf("letezo fajl vegere iras: %s\n", fajlnev);
	}
	
	if(ment==4) writetoBin(matrix, enzim, save, meret, 0, dataSave, enzSave);
	else {
		fajlbaCellaHeaders(cellak, enzakt_num);
		fajlbaCella(cellak, matrix, enzim, enzakt_num, meret, 0);
	}
	if(ment) gsl_rng_fwrite (rngSaveStart, r); //elmenteni kezdeti randomszam generator allapotot

	
	//Iteracio
	for(ciklus=0; ciklus<ciklusszam && replikator_num; ciklus++) {
		if(load && (ciklus==0)) {
			ciklus = betolt.h.gen;
		}
		diffuzio_szam=0;
		//REPLIKACIOS LEPES
		for(iter=0; iter<meret; iter++) {
			cella=gsl_rng_uniform_int(r, meret);
			if ((*(matrix+cella)).k != 0) { //ha van benne vmi
				if(gsl_rng_uniform(r) < phalal) extinkcio(matrix, enzim, enzakt_num, cella); //meghal
			}
			else { //ha nincs benne semmi
				*(claimek+0)=claimEmpty;
				claimSum=claimEmpty;
				for (num_repl_neigh=1; num_repl_neigh<repl_neigh_cellaszam; num_repl_neigh++) {
					nezett= *(repl_neigh+cella*repl_neigh_cellaszam+num_repl_neigh);
					knezett=((*(matrix+nezett)).k);
					if (knezett) {
						*(claimek+num_repl_neigh) = (metabolizmus(matrix, enzim, met_neigh, met_neigh_cellaszam, enzakt_num, reciprocEnzakt_num, nezett) * knezett);
						claimSum+= *(claimek+num_repl_neigh);
						//*(claimek+num_repl_neigh);
						//metabolizmus(matrix, enzim, met_neigh, met_neigh_cellaszam, enzakt_num, 13);
						//double xx= metabolizmus(matrix, enzim, met_neigh, met_neigh_cellaszam, enzakt_num, 13);
//						printf("\n%g", *(claimek+num_repl_neigh));
//						printf("\nreciprocenzakt_num=%g", reciprocEnzakt_num);
					}
					else *(claimek+num_repl_neigh) = 0;
				}
//				printf("\nClaimSum= %g", claimSum);
//				printf("\nnezett cella: %d", cella);
				valaszto= gsl_rng_uniform_pos(r);
				templat=-2; //ha ilyet latsz, akkor nagy baj van...
				claimJelolt=0;
				for (num_repl_neigh=0; num_repl_neigh<repl_neigh_cellaszam; num_repl_neigh++) {
					if(*(claimek+num_repl_neigh)) {claimJelolt += *(claimek+num_repl_neigh)/claimSum;}
//					printf("\nclaim jelolt (%d): %g", (*(matrix+*(repl_neigh+cella*repl_neigh_cellaszam+num_repl_neigh))).szerk, claimJelolt);
					if (claimJelolt > valaszto) {
					  if (num_repl_neigh==0) templat=-1; //ha ures marad
					  else templat=*(repl_neigh+cella*repl_neigh_cellaszam+num_repl_neigh);
					  break;
					}
				}
//				printf("\ntemplat: %d\tclaimJelolt=%g\tvalaszto=%g", templat, claimJelolt, valaszto);
				/*REPLIKACIO!!!!*/
				if(templat > -1) { //ha van masolando templat
					if(replikacio(matrix, enzim, enzakt_num, spec_limit, templat, cella)) {
						printf("gaz van a replikacional");
						return (2); 
					}
//					printf("\n\nreplikacio: %d -> %d", templat, cella);
					/*MUTACIO!!!*/
					if (pMutacio > gsl_rng_uniform_pos(r)) {
						MUTACIO (matrix, enzim, enzakt_num, sd, cella, emax, tradeoffEE, reciprocTradeoffEE, tradeoffEK, reciprocTradeoffEK, kmax, kmin, spec_limit);
					}
				}
			}
			//DIFFUZIO
			for(diffuzio_szam+=diffuzioGyak; diffuzio_szam>=1; diffuzio_szam--) {
				cella=gsl_rng_uniform_int(r, meret);
//				konzolraMatrixStruct(matrix, meret);
				diffTM (matrix, enzim, enzakt_num, ncol, cella);
//				printf("\n\n");
//				konzolraMatrixStruct(matrix, meret);
			}
		}
//		fclose(fp);
		//tranzicio
// 		//for(cella= 0; cella<meret; cella++) transition1(matrix, enzim, enzakt_num, cella);
		
		//mintavetel
		if (mintavetel_gyak) {
			if ((ciklus%mintavetel_gyak)==0) {
//				fprintf(output, "\nciklusszam:;%d\n", ciklus+1);
				tipusadatok(adatok, matrix, enzim, enzakt_num, meret, enztipus_num, sorhossz);
				replikator_num=0;
				for (enzakt = 0; enzakt < enztipus_num; enzakt++) {
					replikator_num += *(adatok+ enzakt*sorhossz);
				}
				if(!replikator_num) printf("\na rendszer meghalt (%d .ciklus)\n", ciklus);
				kimenetTipusadat (output, matrix, adatok, enztipus_num, sorhossz, ciklus+1);
			}
		}
		if (cellamintavetel_gyak && ((ciklus%cellamintavetel_gyak)==0)) {
			if(ment==4) writetoBin(matrix, enzim, save, meret, ciklus+1, dataSave, dataSave);
			else fajlbaCella(cellak, matrix, enzim, enzakt_num, meret, ciklus+1);
		}
//		konzolraCella(matrix, enzim, enzakt_num, meret);
//		konzolraMatrixD(enzim, enzakt_num, meret);
//		printf("\n");
	}
	
	
//	konzolraMatrixStruct(matrix, meret);
// 	printf("\nenzimakt:\n");
//	konzolraMatrixD(enzim, enzakt_num, meret);
		
//  	double xxx=metabolizmus(matrix, enzim, met_neigh, met_neigh_cellaszam, enzakt_num, 13);
// 	printf("\n\n\nM=%g\n", xxx);
	switch(ment)
	{
	  case 1:
		tipusadatok(adatok, matrix, enzim, enzakt_num, meret, enztipus_num, sorhossz);
		mentes(matrix, enzim, adatok, ciklus, meret, enzakt_num, rngSave, dataSave);
		break;
	  case 2:
		mentesBin(matrix, enzim, save, meret, ciklus, rngSave, dataSave, enzSave);
		break;
	  case 3:
		writetoBin(matrix, enzim, save, meret, ciklus+1, dataSave, dataSave);
		gsl_rng_fwrite (rngSave, r);
		break;
	  case 4:
		writetoBin(matrix, enzim, save, meret, ciklus+1, dataSave, dataSave);
		gsl_rng_fwrite (rngSave, r);
		break;
	}
	
	time(&timer);
	printf("\n%s stopped at %s\n", kezdet, ctime( &timer ));
	
	/******************************
	****A LENYEG VEGE****
	*******************************/
	//matrixok felszabaditasa
	free(matrix);
	free(enzim);
	free(met_neigh);
	free(repl_neigh);
	free(claimek);
	free(adatok);
//	printf("felszabaditva\n");

	fclose(output);
	fclose(fp);
	if (ment != 4) fclose(cellak);
	switch(ment) 
	{
	  case 0:
		break;
	  case 2:
		fclose(enzSave);
	  default:
		fclose(rngSave);
		fclose(rngSaveStart);
		fclose(dataSave);
		break;
	}
	
	//randomszam generator lezarasa
	gsl_rng_free(r);
	
	if(!replikator_num) return(5);
	else return(0);
}
