CC=g++ -fopenmp -Wall -O4
OBJS=emp1d.o
LFLAGS= -lm -L /home/ram80/include -lgomp
CFLAGS= -I /home/ram80/include -lgomp
all: emp1d
emp1d: $(OBJS)
	$(CC) -o emp1d $(OBJS) $(CFLAGS)
emp1d.o: emp1d-1.0.cpp
	$(CC) -o emp1d.o -c emp1d-1.0.cpp $(LFLAGS)
clean:
	rm *.o emp1d
