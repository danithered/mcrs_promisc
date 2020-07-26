#include <mainheader.h>
#include <stdlib.h>
#include <math.h>

//feltolti a szomszedsagmatrixot: sajat pozicio, N, W, E, S
int* metNeighInic(int cellaszam_f, int ncol_f, double neigh_tipus_f) {
	/*
	 * cellaszam_f: az alapmatrix cellaszama
	 * neigh_tipus_f: a szomsz. merete: ha 0 -> von Neumann; ha több -> Moore szomsz
	 * */
	
	//Valtozok deklaralasa
	int metNeigh_matrixmeret_f=0, *matrix_f, alapMCellaszamlalo_f=0, colszamlalo_f=0, rowszamlalo_f=0, vonN_f = 0, szomsz_f=0, szomszedszam_f=0;
	/*
	 * metNeigh_matrixmeret_f: a szomsz. matrix cellaszama
	 * alapMCellaszamlalo_f: szamlalo: hanyadik alapmatrix cellanal jarunk a feltoltesben
	 * ncol_f: alapmatrix oszlopainak szama =sqrt(cellaszam_f)
	 * colszamlalo_f: hanyadik oszlopban tartunk
	 * rowszamlalo_f: hanyadik sorban tartunk
	 * szomsz_f: egy cella szomszedsagaba hany cella tartozik pl.: 2es szomsz eseten (2*2+1)^2=25
	 * *matrix_f: a szomsz. matrix
	 * */
	
	
	szomsz_f = szomsz_meret(neigh_tipus_f);
	if(neigh_tipus_f==0 || ( neigh_tipus_f - (int)neigh_tipus_f )){
		vonN_f = 1;
		neigh_tipus_f = (double)(int)neigh_tipus_f;
	}
		
	metNeigh_matrixmeret_f=cellaszam_f* szomsz_f; //szomsz. matrix cellameretenek meghatarozasa
	matrix_f=(int *) calloc(metNeigh_matrixmeret_f, sizeof(int)); //szomsz. matrix lefoglalasa
	for(alapMCellaszamlalo_f=0; alapMCellaszamlalo_f<cellaszam_f; alapMCellaszamlalo_f++) {
		szomszedszam_f=0;
		//*(matrix_f+ alapMCellaszamlalo_f*szomsz_f+szomszedszam_f)=alapMCellaszamlalo_f;
		
		for(rowszamlalo_f = 0 - neigh_tipus_f; rowszamlalo_f <= 0 + neigh_tipus_f; rowszamlalo_f++) {
			for(colszamlalo_f = 0 - neigh_tipus_f; colszamlalo_f <= 0 + neigh_tipus_f; colszamlalo_f++) {
				//szomszedszam_f++;
				if(!colszamlalo_f && !rowszamlalo_f) *(matrix_f+ alapMCellaszamlalo_f*szomsz_f)=torus(cellaszam_f, alapMCellaszamlalo_f+rowszamlalo_f*ncol_f+colszamlalo_f);
				else {
					*(matrix_f+ alapMCellaszamlalo_f*szomsz_f+szomszedszam_f+1)=torus(cellaszam_f, alapMCellaszamlalo_f+rowszamlalo_f*ncol_f+colszamlalo_f);
					szomszedszam_f++;
				}
			}
		}
		if(vonN_f) {
			*(matrix_f + alapMCellaszamlalo_f*szomsz_f + szomszedszam_f + 1) = torus(cellaszam_f, alapMCellaszamlalo_f - (1 + neigh_tipus_f)*ncol_f); //1: N
			*(matrix_f + alapMCellaszamlalo_f*szomsz_f + szomszedszam_f + 2) = torus(cellaszam_f, alapMCellaszamlalo_f - (1 + neigh_tipus_f)); //2: W
			*(matrix_f + alapMCellaszamlalo_f*szomsz_f + szomszedszam_f + 3) = torus(cellaszam_f, alapMCellaszamlalo_f + (1 + neigh_tipus_f)); //3: E
			*(matrix_f + alapMCellaszamlalo_f*szomsz_f + szomszedszam_f + 4) = torus(cellaszam_f, alapMCellaszamlalo_f + (1 + neigh_tipus_f)*ncol_f); //4: S
		}
	}

	return (matrix_f);
}

int szomsz_meret(double tipus_f) {
	if(tipus_f < 1) return(5); //vonNeumann
	else {
		if(tipus_f - (int)tipus_f ) return(szomsz_meret((double)(int)tipus_f) + 4);
		else return (pow(((int)tipus_f*2+1), 2)); //Moore szomszedsag: pl.: 2es szomsz eseten (2*2+1)^2=25 
	}
}