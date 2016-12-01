CXX=g++
CFLAGS=-I.
DEPS =	utils.h \
	typedef.h \
	models.h \
	callers.h

OBJ  =	models.o \
	utils.o \
	callers.o

OBJM =	main.o 

OBJT =	test.o

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

fml_sol : $(OBJM) $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)

ptest : $(OBJT) $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)

install : fml_sol
	./install.perl

.PHONY: clean

clean:
	rm -f *.o *~
