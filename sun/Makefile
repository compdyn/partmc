
all: simple

SIMPLE_OBJS = simple.o simple_helper.o

simple: $(SIMPLE_OBJS)
	gfortran -arch i386 -o $@ $(SIMPLE_OBJS) -L/Users/mwest/opt/lib -lsundials_cvode -lsundials_nvecserial

simple.o: simple.f90
	gfortran -c -I/Users/mwest/opt/include -o $@ $<

simple_helper.o: simple_helper.c
	gcc -arch i386 -c -I/Users/mwest/opt/include -o $@ $<




simple_c: simple_c.o
	gcc -arch i386 -o $@ $< -L/Users/mwest/opt/lib -lsundials_cvode -lsundials_nvecserial

simple_c.o: simple_c.c
	gcc -arch i386 -c -I/Users/mwest/opt/include -o simple_c.o simple_c.c

simple_f: simple_f.o
	gfortran -o simple_f simple_f.o -L/Users/mwest/opt/lib -lsundials_cvode -lsundials_nvecserial

simple_f.o: simple_f.f90
	gfortran -c -I/Users/mwest/opt/include -o simple_f.o simple_f.f90

clean:
	rm -rf simple_c simple_c.o simple_f simple_f.o simple simple.o simple_helper.o
