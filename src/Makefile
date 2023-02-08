CC      = gcc
ACRM    = $(HOME)
LOPT    = -L$(ACRM)/lib -I$(ACRM)/include
OFILES  = carcons.o 
LIBS    = -lbiop -lgen -lm
EXE     = carcons
COPT    = -g -Wall -ansi -pedantic $(LOPT) -DDEBUG
#COPT    = -O3 -std=c99 -Wall -ansi -pedantic $(LOPT)

$(EXE) : $(OFILES)
	$(CC) $(COPT) -o $@ $(OFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c $<

clean :
	\rm -f $(OFILES)

distclean : clean
	\rm $(EXE)
