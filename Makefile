CXX=g++
CFLAGS=-I.
DEPS =	utils.h \
	typedef.h \
	models.h

OBJ  =	main.o \
	models.o \
	utils.o

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

fml_sol : $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o *~
