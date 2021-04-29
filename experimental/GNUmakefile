
CHOMBO_HOME := /mnt/beegfs/paolini-notos/Chombo_4/example/..
DIM := 2
VARS_OLD := $(.VARIABLES)

SRC := $(wildcard *.cpp) 
SRC += $(wildcard ../src/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BaseTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BoxTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/AMRTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/CrunchFlow/*.cpp) 

NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))


VPATH=o. ../src  ../../common $(CHOMBO_HOME)/src/BaseTools $(CHOMBO_HOME)/src/BoxTools $(CHOMBO_HOME)/src/AMRTools $(CHOMBO_HOME)/src/CrunchFlow
LIBS = 

CXXALLFLAGS := -Wall -Wno-unused-but-set-variable -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-99 -Wno-un
CXXGPUFLAGS := --compiler-options='$(CXXALLFLAGS)'

SRC += $(wildcard ../../common/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/EBTools/*.cpp) 

VPATH +=  $(CHOMBO_HOME)/src/EBTools 

CXXALLFLAGS := -Wall -Wno-unused-but-set-variable -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-99 -Wno-un
CXXGPUFLAGS := --compiler-options='$(CXXALLFLAGS)' -DPROTO_CUDA

CPPFLAGS = -std=c++14 -DCXXSTD=14

HOME = /mnt/beegfs/paolini-notos
PROTO_HOME = $(HOME)/proto
HDF5_ROOT = /usr/lib/x86_64-linux-gnu/hdf5/serial
DEPXX=nvcc

export OMPI_CXX=nvcc

CXXALLFLAGS := #-std=c++11
CXXDBGFLAGS:=$(CXXALLFLAGS) -O0 -g -pg #-pedantic 
CXXOPTFLAGS:=$(CXXALLFLAGS) -O3 -use_fast_math
MPI_DIR=/usr/lib/openmpi
LAPACKLIBS = -llapacke -llapack
LIBS +=-L$(MPI_DIR)/lib -lmpi
LIBS +=$(LAPACKLIBS)
LIBS += $(CUDALIBS)

HDFINCFLAGS=-I$(HDF5_ROOT)/include/ 
HDFLIBFLAGS=-L$(HDF5_ROOT)/lib/ -lhdf5 

#end mach dependent

LIBS     += $(HDFLIBFLAGS)
 
CXXFLAGS += -DCH_USE_DOUBLE

CXXFLAGS+=$(CXXDBGFLAGS)
#CXXFLAGS+=$(CXXOPTFLAGS) 

CPPFLAGS += -I../../common
CPPFLAGS += -DCH_USE_HDF5 
CPPFLAGS += $(HDFINCFLAGS)
CPPFLAGS += -I$(CHOMBO_HOME)/src/EBTools
CPPFLAGS += -I$(CHOMBO_HOME)/src/EBProto
CPPFLAGS += -DDIM=$(DIM) -I$(PROTO_HOME)/include  -I$(PROTO_HOME)/EBProto/include -I. -I../src 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BaseTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BoxTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/AMRTools
CPPFLAGS += -I$(CHOMBO_HOME)/src/CrunchFlow
CPPFLAGS += -DCH_SPACEDIM=$(DIM) -DCH_USE_64 
CPPFLAGS += -DCH_LANG_CC
CPPFLAGS += -std=c++14 -DCXXSTD=14 
CPPFLAGS += -I$(MPI)/include/

#####
#
# Uncomment next three lines to run on GPU (device), comment to run on CPU (host)
#
CXX=nvcc -g
CPPFLAGS += -DPROTO_CUDA=1 -x cu --expt-extended-lambda -arch=compute_70
CUDALIBS = -lcusolver -lcudart -lnvToolsExt
#
#####

LINK=$(CXX) -pg
OBJDIR = o
DEPDIR = d
TARGET = main.exe

$(TARGET): $(OBJ)
	$(LINK)  $(OBJ) $(LIBS) -o $(TARGET)

$(OBJDIR): ; mkdir -p $@


$(DEPDIR): ; mkdir  -p $@

DEPFILES := $(NODIRSRC:%.cpp=$(DEPDIR)/%.d)

$(DEPFILES):

DEPFLAGS = -MT $@ -MF $(DEPDIR)/$*.d 

ifeq ($(CXX),nvcc)
LIBS+= -lnvToolsExt
DEPFLAGS += --generate-dependencies
else
DEPFLAGS += -MM
endif

$(OBJDIR)/%.o : %.cpp  $(DEPDIR)/%.d | $(DEPDIR) 
	@mkdir -p $(OBJDIR)
	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $< 
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(CPPFLAGS) 

all: $(TARGET)

-include $(DEPFILES)

.PHONY: clean realclean print dirs

print-%: ; @echo $* = $($*)


vars:
	$(foreach v,$(filter-out $(VARS_OLD) VARS_OLD,$(.VARIABLES)), $(info $(v) = $($(v))))

clean:
	rm -rf d o *.o *.exe *.d ../../common/*.o ../src/*.o ../../common/*.d ../src/*.d

realclean:
	rm -rf d o *.o *.exe *.d ../../common/*.o ../src/*.o ../../common/*.d ../src/*.d

pristine:
	rm -rf d o *.o *.exe *.table *.d ../../common/*.o ../src/*.o ../../common/*.d ../src/*.d pout.* *.time.table.* *.hdf5
run:
	./main.exe
	gprof main.exe gmon.out > analysis.txt
	python /usr/local/lib/python3.5/dist-packages/gprof2dot.py<analysis.txt|dot -Tpng -o output.png
