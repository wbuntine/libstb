##  Linux versions
# -lrt needed for clock_gettime(), Linux only
#      but only used in check.c
#  assumes GSL stored in regular paths
CFLAGS = -I../lib/ -Wall -g # O5 -DNDEBUG 
LDFLAGS = -L../lib/ -lstb -lm 

## MAC OSX versions
#CFLAGS = -I/opt/local/include -I../lib/ -Wall -O5 -DNDEBUG 
#LDFLAGS = -L/opt/local/lib -L../lib/ -lstb -lm 

LIBRARY = libstb.a

all:    demo list


demo:   demo.o 
	cc $(CFLAGS)  -o demo demo.o $(LDFLAGS)

gen:   gen.o 
	cc $(CFLAGS)  -o gen gen.o $(LDFLAGS)

list:   list.o 
	cc $(CFLAGS)  -o list list.o $(LDFLAGS)

check:	check.o ../lib/$(LIBRARY)
	cc  $(CFLAGS)  -o check check.o $(LDFLAGS) -lrt

clean: 
	rm -f *.o demo list check precision_test $(LIBRARY)

precision_test:   precision_test.o 
	cc $(CFLAGS)  -o precision_test precision_test.o $(LDFLAGS)

depend:
	makedepend $(FILES) 

# DO NOT DELETE THIS LINE -- make depend depends on it. 
