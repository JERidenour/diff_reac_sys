################################
CXX = mpicc
################################
CXXFLAGS = -std=c99 -Wall -ferror-limit=1 -O2
################################
# I = -ICSparse/Include
# CS = CSparse/Lib/libcsparse.a
################################

all: gs_seq

gs_seq_dep = gs_seq.o
gs_seq: $(gs_seq_dep)
	$(CXX) $(CXXFLAGS) $(gs_seq_dep) -o $@

gs_seq.o: gs_seq.c
	$(CXX) $(CXXFLAGS) -c -o $@ gs_seq.c



# all: lib crank

# lib:
# 	( cd CSparse/Lib ; $(MAKE) )

# crank_dep = crank.o
# crank: lib $(crank_dep)
# 	$(CXX) $(CXXFLAGS) $(I) $(crank_dep) $(CS) -o $@

# crank.o: crank.c
# 	$(CXX) $(CXXFLAGS) -c -o $@ crank.c

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
