################################
CXX = mpicc
DEBUG = -g
################################

CXXFLAGS = -std=c99 -Wall -ferror-limit=1 $(DEBUG)
################################

################################

################################

all: crank

crank_dep = crank.o
crank: $(crank_dep)
	$(CXX) $(CXXFLAGS) $(crank_dep) -o $@

crank.o: crank.c
	$(CXX) $(CXXFLAGS) -c -o $@ crank.c

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
