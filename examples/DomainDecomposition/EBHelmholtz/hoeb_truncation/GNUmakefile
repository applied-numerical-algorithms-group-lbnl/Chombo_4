
CHOMBO_HOME := /Users/danielwaters/chombo_4/example/..
DIM := 2
VARS_OLD := $(.VARIABLES)

SRC := $(wildcard *.cpp) 
SRC += $(wildcard ../src/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BaseTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BoxTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/AMRTools/*.cpp) 
SRC += $(wildcard ../../common/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/EBTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/EBProto/*.cpp) 

VPATH=o. ../src  ../../common $(CHOMBO_HOME)/src/BaseTools $(CHOMBO_HOME)/src/BoxTools $(CHOMBO_HOME)/src/AMRTools 
VPATH +=  $(CHOMBO_HOME)/src/EBTools 

CHOMBO_CXXALLFLAGS := 
CHOMBO_CXXGPUFLAGS := --compiler-options='$(CHOMBO_CXXALLFLAGS)'

NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))


LIBS = $(LAPACKLIBS)

CHOMBO_CXXALLFLAGS := 
CHOMBO_CXXGPUFLAGS := --compiler-options='$(CHOMBO_CXXALLFLAGS)'

SRC += $(wildcard ../../hoeb_common/*.cpp) 

VPATH +=  ../../hoeb_common/

CHOMBO_CPPFLAGS = -std=c++14

PROTO_HOME =$(HOME)/proto

export OMPI_CXX=nvcc

LIBS=-L$(MPI_ROOT) -lmpi_ibm -lmpiprofilesupport

HDFINCFLAGS=-I$(OLCF_HDF5_ROOT)/include
HDFLIBFLAGS=-L$(OLCF_HDF5_ROOT)/lib -lhdf5 -lz
#end mach dependent





CHOMBO_CXXFLAGS += -DCH_USE_DOUBLE

CHOMBO_CXXFLAGS+= -DPR_STACK_ALLOC=4294967296
CHOMBO_CXXFLAGS+=$(CXXDBGFLAGS) 

CHOMBO_CPPFLAGS +=-DDIM=$(DIM) -I$(PROTO_HOME)/include   -I. -I../src 
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/BaseTools 
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/BoxTools 
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/AMRTools
CHOMBO_CPPFLAGS += -DCH_SPACEDIM=$(DIM) -DCH_USE_64 
CHOMBO_CPPFLAGS += -I../../common
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/EBTools -I$(CHOMBO_HOME)/src/EBProto


#need to use depxx because petsc links with g++
LINK=$(DEPXX)
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
endif

ifeq ($(DEPXX),nvcc)
DEPFLAGS += --generate-dependencies
else
DEPFLAGS += -MM
endif

$(OBJDIR)/%.o : %.cpp  $(DEPDIR)/%.d | $(DEPDIR) 
	@mkdir -p $(OBJDIR)
	$(DEPXX) $(DEPFLAGS) $(CHOMBO_CPPFLAGS) $< 
	$(CXX) -c -o $@ $< $(CHOMBO_CXXFLAGS) $(CHOMBO_CPPFLAGS) 

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

CHOMBO_CPPFLAG S += -I../../hoeb_common 
