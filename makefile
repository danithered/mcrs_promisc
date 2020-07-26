IDIR =./include
CC=gcc
CFLAGS=-I$(IDIR) `pkg-config --cflags --libs gsl`

ODIR=src/obj
LDIR =./lib
SRCDIR=src

LIBS=-lm

_DEPS = mainheader.h randomgen.h adattipusok.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o feltoltes.o szomszed.o torus.o metab.o mutacio.o diffuzio.o bitmuveletek.o statisztika.o kimenet.o save.o load.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJo = olvaso.o kimenet.o save.o load.o statisztika.o bitmuveletek.o konzolra.o olvaso_kieg.o szomszed.o torus.o
OBJo = $(patsubst %,$(ODIR)/%,$(_OBJo))


$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

progi: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  

wall: $(OBJ)
	$(CC) -o progi $^ $(CFLAGS) $(LIBS) -Wall

olvaso: $(OBJo)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
