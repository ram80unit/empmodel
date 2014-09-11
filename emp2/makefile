CC=g++ -fopenmp -Wall -O4
OBJS=emp2d.o
SRC=emp2d-fork.cpp
LFLAGS= -lm -L /home/ram80/include -lgomp
CFLAGS= -I /home/ram80/include -lgomp
all: emp2d
emp2d: $(OBJS)
	$(CC) -o emp2d $(OBJS) $(CFLAGS)
emp2d.o: $(SRC)
	$(CC) -o emp2d.o -c $(SRC) $(LFLAGS)
clean:
	rm *.o emp2d
