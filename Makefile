#.PHONY all clean

PROG = test
CC = g++-11
PATHFLAGS = -I/opt/homebrew/Cellar/armadillo/10.8.2/include
LIBFLAGS = -L/opt/homebrew/Cellar/armadillo/10.8.2/lib -larmadillo
OBJS = main.o functions.o system.o
$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS) -g3 $(PATHFLAGS) $(LIBFLAGS)

main.o : main.cpp
	$(CC) $(PATHFLAGS) $(LIBFLAGS) -c main.cpp

functions.o : functions.cpp functions.h
	$(CC) $(PATHFLAGS) $(LIBFLAGS) -c functions.cpp
	
system.o : system.cpp system.h
	$(CC) $(PATHFLAGS) $(LIBFLAGS) -c system.cpp
clean:
	rm -f core $(PROG) $(OBJS)
