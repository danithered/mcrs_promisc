#include <olvaso.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
/*
 * argomentumok:
 * 	<FILE>	: beolvasando 1. file cime
 * 	<FILE>	: beolvasando 2. file cime
 * 	-h	: olvassa ki a headereket csak
 * 	-gen	: olvassa ki a generaciokat csak
 * 	-g xn	: melyik adatokat olvassa ki a fajlbol. (0 v minusz: az utolsot, pozitiv: adott szamut) Kihagyasa eseten az elsot olvassa be.
 * 	--help	: kiir helpet
 */
	int argoment=0, gKezdArg=0, noFiles=0, mode=0, noKiir=0, *kiir, nezettGen=0, gen=0, compare=0, file=0, metab=0;
//	char *a="abcdefg\0";
	char *filename1, *filename2;
	FILE *file1, *file2, *out;
	struct Save header1, header2;
	struct SaveAll beolvasas1, beolvasas2;
	double *enzim1, *enzim2;
	struct Cella *matrix1, *matrix2;
	char help[]="This is a text which explains this program and its functions. \nWhen I have time, I'll write something here.";
	
	out=stdout;
	/* argoment: szamlalo: hanyadik argomentumot nezzuk
	 * noFiles: hany fajlt adtak meg
	 * mode:
	 * 	0: nem csinal semmit
	 * 	1: kiir fajl1-et es fajl 2-t
	 * 	2: kiir es kombinal fajl 1-et(alap) es fajl 2-t(enzim)
	 * 	3: kiir es kombinal fajl 1-et(enzim) es fajl 2-t(alap)
	 * 	-1: kiir helpet
	 * 	-2: kiir generaciokat
	 * 	-3: kiir header-eket
	 * noKiir: hany adatsort ir ki. Ha 0: az elsot irja ki.
	 * kiir: melyik adatsorokat irja ki
	 * nezettGen: szamlalo: melyik generaciot nezzuk epp (N-ediket a kiir vektorban)
	 * gen: melyik generaciot nezzuk epp
	 * compare: compareHeaders kimenetelet tarolja
	 * filename1, filename2: a fajlok helyet tartalmazo stringek
	 * file1, file2: a fajlokra mutato pointerek
	 */
	
	//ADATOK BEOLVASASA
	for(argoment=1; argoment<argc; argoment++) {
		if(strlen(argv[argoment]) < 2) continue; //2 karakternel kisebb argomentumokat nem olvas be
		
		//kapcsolok
		if(!strncmp(argv[argoment], "-", 1)){ //innentol csak kapcsolokat nez
			for(; argoment<argc; argoment++) { //vegignez kapcsolokat
//				printf("\t%s\n", argv[argoment]+1);
				//mode kapcsolok
				if(!strcmp(argv[argoment]+1, "-help\0")) { //--help
					mode=-1;
					break;
				}
				if(!strcmp(argv[argoment]+1, "h\0")) { //-h
					mode=-3;
					break;
				}
				if(!strcmp(argv[argoment]+1, "gen\0")) { //-gen
					mode=-2;
					break;
				}
				if(!strcmp(argv[argoment]+1, "met\0")) { //-met
					metab=1;
					break;
				}
				if(!strcmp(argv[argoment]+1, "ntc\0")) { //-ntc
					metab=2;
					break;
				}
				if(!strcmp(argv[argoment]+1, "ast\0")) { //-ast
					metab=3;
					break;
				}
				//adatok sorszama
				if(!strcmp(argv[argoment]+1, "g\0") && noKiir == 0 && argoment < argc-1) { //-g
					argoment++; //megnez kovetkezo argomentumot
					gKezdArg = argoment;
//					printf("itt kezdodnek a gen-ek: %d\n", gKezdArg);
					for(; argoment < argc && strncmp(argv[argoment], "-", 1); argoment++) {
						noKiir++;
					}
					argoment--; //visszateker, hogy kovetkezo kapcsolot megnezhesse a kulso ciklus
//					printf("argoment= %d\n", argoment);
				}
			}
			break;
		}
		
		//file nev eltarol
		switch (noFiles) {
			case 0: 
				filename1 = argv[argoment];
				noFiles++;
				mode=1;
				break;
		
			case 1:
				filename2 = argv[argoment];
				noFiles++;
				break;
		}
	}
	//beolvas megnezendo generaciokat
	if(noKiir) {
//		printf("van gen: %d db:\n", noKiir);
		kiir = (int *) calloc(noKiir, sizeof(int));
		for(argoment=0; argoment < noKiir; argoment++) {
			*(kiir + argoment) = atoi(argv[gKezdArg + argoment]);
//			printf("\t%d\n", *(kiir + argoment));
		}
	}
	
	//ADATOK FELDOLGOZASA
	if (mode == -1) {
		printf("help text for %s program:\n%s\n", argv[0], help);
		return(0);
	}
	if (noFiles<1) printf("You didn't add a file to process. You must add the filename before the - tags!\n");
	for(nezettGen=0; (nezettGen< noKiir || (nezettGen == 0 && noKiir ==0) ) && (mode == 1 || mode == 2 || mode == 3); nezettGen++) {	
		if(noKiir) gen = *(kiir + nezettGen); //ha megvan adva '-g' argomentum
		else gen = 1; //ha nincs megadva '-g' argomentum => elsot nezi
//		printf("betoltott fajl: %s\n", filename1);
		beolvasas1 = betoltes(filename1, gen);
		header1 = beolvasas1.h;
		if(header1.version != 1) return (header1.version);
		if (header1.tipus != 1) enzim1 = beolvasas1.e;
		else enzim1 = (double*) calloc(header1.ncol * header1.nrow * header1.eakt, sizeof(double));
		if(header1.tipus != 2) matrix1 = beolvasas1.m;
		else matrix1 = (struct Cella*) calloc(header1.ncol * header1.nrow, sizeof(struct Cella));
//		konzolraSaveHeader();
//		konzolraSave(header1);
		if (noFiles == 2){
			beolvasas2 = betoltes(filename2, gen);
			header2 = beolvasas2.h;
			//enzim2 = beolvasas2.e;
			//matrix2 = beolvasas2.m;
			if(header2.version != 1) return (header2.version);
			if (header2.tipus != 1) enzim2 = beolvasas2.e;
			else enzim2 = (double*) calloc(header2.ncol * header2.nrow * header2.eakt, sizeof(double));
			if(header2.tipus != 2) matrix2 = beolvasas2.m;
			else matrix2 = (struct Cella*) calloc(header2.ncol * header2.nrow, sizeof(struct Cella));
			//osszehasonlit ketto
			compare = compareHeaders(header1, header2);
			if(!compare) { //ket eltero tipusu header
				if (header1.tipus == 1 && header2.tipus == 2) mode = 2;
				if (header1.tipus == 2 && header2.tipus == 1) mode = 3;
			}
			else {
				if(compare < 0 && compare > -100) return(compare);
				else {if(compare < -100) noFiles=1;}
			}
		}
		
		//kiiras
//		//konzolraMatrixStruct(matrix1, 10);
		if( (!nezettGen || (mode == 1 && noFiles == 2)) && !metab ) fajlbaCellaHeaders (out, header1.eakt);
		switch (mode) {
			case 1:
				switch(metab) {
					case 0:
						fajlbaCella (out, matrix1, enzim1, header1.eakt, header1.ncol * header1.nrow, header1.gen);
						if (noFiles == 2) {
							fajlbaCellaHeaders (out, header2.eakt);
							fajlbaCella (out, matrix2, enzim2, header2.eakt, header2.ncol * header2.nrow, header2.gen);
						}
						break;
					case 1:
						fprintf(out, "\n");
						fajlbaMetab (out, matrix1, enzim1, header1.met, header1.eakt, header1.ncol, header1.nrow);
						break;
					case 2:
						fprintf(out, "\n");
						fajlbaSzomszTipusok (out, matrix1, header1.met, header1.eakt, header1.ncol, header1.nrow);
						break;
					case 3:
						fprintf(out, "\n");
						fajlbaAsszocTabla (out, matrix1, header1.met, header1.eakt, header1.ncol, header1.nrow);
						break;
				}
				break;
			case 2:
				fajlbaCella (out, matrix1, enzim2, header1.eakt, header1.ncol * header1.nrow, header1.gen);
				break;
			case 3:
				fajlbaCella (out, matrix2, enzim1, header1.eakt, header1.ncol * header1.nrow, header1.gen);
				break;
		}
//		printf("%d", (*(matrix1)).spec );
	}

//	printf("fileok: %s %s\nmode: %d, metab: %d\n", filename1, filename2, mode, metab);
//	printf("%s\n", a+2);
	
	for(file=0; (mode == -2 || mode == -3) && file < noFiles; file ++){
		if (file == 0) file1 = fopen64(filename1, "rb");
		else file1 = fopen64(filename2, "rb");
		if(!file1) {
			printf("nem lehet megnyitni fajlt: %s\nBeolvasas megszakitva\n", filename1);
			return(-1);
		}
		if (file == 0) printf("%s\n", filename1);
		else printf("%s\n", filename2);
		fread(&header1, sizeof(struct Save), 1, file1);
		if(mode == -3) {
			konzolraSaveHeader();
			//konzolraSave(header1);
		}
		//else printf("%d\t", header1.gen);
		rewind(file1); //visszateker elejere
		for(nezettGen=1; nezettGen < header1.num; nezettGen++) {
			//header2 = teker(file1, 1);
			fread(&header2, sizeof(struct Save), 1, file1);
			if(mode == -3) konzolraSave(header2);
			else printf("%d\t", header2.gen);
			fseek(file1, adatMeret(header2.ncol, header2.nrow, header2.tipus, header2.eakt), SEEK_CUR);
			//rewind(file1);
		}
		printf("\n");
		fclose(file1);
	}
//	printf("\nmeret: %d\n", (int) 300*300* (2 * sizeof(double) + sizeof(struct Cella)) + sizeof(struct Save) );
	
	
	if(noKiir) free(kiir);
	if(mode == 1 || mode == 2 || mode == 3) {
		free(matrix1);
		free(enzim1);
		if(noFiles == 2) {
			free(matrix2);
			free(enzim2);
		}
	}
	return (0);
}

