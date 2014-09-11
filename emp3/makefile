CC=g++ -fopenmp -Wall -O4
SRC=emp3d-fork.cpp
OBJS=emp3d.o
LFLAGS= -lm -L /home/ram80/include -lgomp
CFLAGS= -I /home/ram80/include -lgomp
all: emp3d
emp3d: $(OBJS)
	$(CC) -o emp3d $(OBJS) $(CFLAGS)
emp3d.o: $(SRC)
	$(CC) -o emp3d.o -c $(SRC) $(LFLAGS)
clean:
	rm *.o emp3d
