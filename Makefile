################################
CXX = mpicc
################################
CXXFLAGS = -std=c99 -Wall -ferror-limit=1 -O
################################
I = -ICSparse/Include
CS = CSparse/Lib/libcsparse.a
################################

all: lib lap

lib:
	( cd CSparse/Lib ; $(MAKE) )

lap_dep = lap.o
lap: lib $(lap_dep)
	$(CXX) $(CXXFLAGS) $(I) $(lap_dep) $(CS) -o $@

lap.o: lap.c
	$(CXX) $(CXXFLAGS) -c -o $@ lap.c

# mtsort_dep = mtsort.o
# mtsort: $(mtsort_dep)
# 	$(CXX) $(CXXFLAGS) $(mtsort_dep) -o $@

# mtsort.o: mtsort.c
# 	$(CXX) $(CXXFLAGS) -c -o $@ mtsort.c

# test_merge_dep = test_merge.o
# test_merge: $(test_merge_dep)
# 	$(CXX) $(CXXFLAGS) $(test_merge_dep) -o $@

# test_merge.o: test_merge.c
# 	$(CXX) $(CXXFLAGS) -c -o $@ test_merge.c

# test_merge2_dep = test_merge2.o
# test_merge2: $(test_merge2_dep)
# 	$(CXX) $(CXXFLAGS) $(test_merge2_dep) -o $@

# test_merge2.o: test_merge2.c
# 	$(CXX) $(CXXFLAGS) -c -o $@ test_merge2.c

clean:
	$(RM) *.o *.out

realclean: clean
	$(RM) -r *.dSYM
	$(RM) *~
	$(RM) -r docs
