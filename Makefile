CXX = mpicxx
# use -xmic-avx512 instead of -xHost for Intel Xeon Phi platforms
OPTFLAGS = -O3 -fopenmp -DCHECK_RESULTS -DUSE_MPI_P2P #-DPRINT_RESULTS #-DDEBUG_PRINTF

# Options for choosing the MPI variants
#-DUSE_MPI_NPP
#-DUSE_MPI_NRM
#-DUSE_MPI_P2P
#-DUSE_MPI_NCL
#-DUSE_MPI_RMA
#-DUSE_MPI_UPX, use -DREPLACE_UPX_WITH_RMA to replace UPCXX with MPI calls

# use export ASAN_OPTIONS=verbosity=1 to check ASAN output
SNTFLAGS = -std=c++17 -fsanitize=address -O1 -fno-omit-frame-pointer
CXXFLAGS = -std=c++17 -g $(OPTFLAGS)

# Make this 0 if you don't want to use UPCXX
ENABLE_UPCXX=0
ifeq ($(ENABLE_UPCXX),1)
    CXXFLAGS += -DUPCXX_ASSERT_ENABLED=1 -DUPCXX_BACKEND=1 -DUPCXX_BACKEND_GASNET_SEQ=1 -DUPCXX_MPSC_QUEUE_ATOMIC=1 -D_GNU_SOURCE=1 -DGASNET_SEQ -I/home/sg/builds/upcxx/gasnet.debug/include -I/home/sg/builds/upcxx/gasnet.debug/include/smp-conduit -I/home/sg/builds/upcxx/upcxx.debug.gasnet_seq.smp/include -std=c++14 -g -wd654 -wd1125 -wd279 -wd1572
    GASNET_CONDUIT=smp
    GASNET_INSTALL=/home/sg/builds/upcxx/gasnet.debug
    LDFLAGS=-g -wd177 -wd279 -wd1572
    LIBS=-L/home/sg/builds/upcxx/upcxx.debug.gasnet_seq.smp/lib -lupcxx -L/home/sg/builds/upcxx/gasnet.debug/lib -lgasnet-smp-seq -lrt -lm
endif

OBJ = main.o
TARGET = match

all: $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

$(TARGET):  $(OBJ)
	$(CXX) $^ $(OPTFLAGS) -o $@ $(LDFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -rf *~ $(OBJ) $(TARGET)
