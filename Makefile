CC  = mpicc
CXX = mpic++

CXXFLAGS = -O3 -Wall -fopenmp
LIBFLAGS = -fPIC

LFLAGS = -L.
LIBS = -lpsort

SRCS = check.cpp time.cpp
LIBSRCS = bsort.cpp rsort.cpp qsort.cpp msort.cpp sort.cpp common.cpp

OBJS = $(SRCS:.cpp=.o)
LIBOBJS = $(LIBSRCS:.cpp=.o)

LIBRARY = libpsort.so
TEST = check time
.PHONY: clean debug depend

$(LIBRARY): $(LIBOBJS)
	$(CXX) $(CXXFLAGS) -shared -o $(LIBRARY) $(LIBOBJS)
	for host in $$(tail -n +2 hostfile); do scp -o ConnectTimeout=4 $(LIBRARY) $$host:~/mpi/ ; done

all: $(LIBRARY) $(TEST)
debug: CXXFLAGS += -DDEBUG
debug: all

profiler: CXXFLAGS += -Dgprofiler
profiler: LFLAGS += -L/usr/local/lib
profiler: LIBS += -lprofiler
profiler: debug

send: $(LIBRARY) $(TEST)
	for host in $$(tail -n +2 hostfile); do scp -o ConnectTimeout=4 $(LIBRARY) $(TEST) $$host:~/mpi/ ; done

$(TEST): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LFLAGS) $(LIBS)
	for host in $$(tail -n +2 hostfile); do scp -o ConnectTimeout=4 $@ $$host:~/mpi/ ; done

$(LIBOBJS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LIBFLAGS) -c $< -o $@

clean:
	$(RM) *.o *~ $(TEST) $(LIBRARY)
kill:
	for host in $$(tail -n +2 hostfile); do ssh $$host "kill -9 -1"; done

depend: $(SRCS) $(LIBSRCS)
	makedepend -Y. $^

# DO NOT DELETE THIS LINE -- make depend needs it

check.o: common.h sort.h payloadsize.h
time.o: common.h sort.h payloadsize.h
bsort.o: sort.h payloadsize.h
rsort.o: sort.h payloadsize.h common.h
qsort.o: sort.h payloadsize.h
msort.o: sort.h payloadsize.h
sort.o: sort.h payloadsize.h common.h
common.o: common.h sort.h payloadsize.h
