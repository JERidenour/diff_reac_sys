################################
CXX = mpicc
################################
CXXFLAGS = -std=c99 -Wall -O2
################################
I = -ICSparse/Include
CS = CSparse/Lib/libcsparse.a
################################

all: lib gray_scott

lib: 
		( cd CSparse/Lib ; $(MAKE) )

gray_scott_dep = gray_scott.o
gray_scott: lib $(gray_scott_dep)
	$(CXX) $(CXXFLAGS) $(I) $(gray_scott_dep) $(CS) -o $@

gray_scott.o: gray_scott.c
	$(CXX) $(CXXFLAGS) -c -o $@ gray_scott.c

clean:
	$(RM) *.o *.out

realclean: clean
	$(RM) -r *.dSYM
	$(RM) *~
	$(RM) -r docs
