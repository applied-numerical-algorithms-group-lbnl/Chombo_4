
CHOMBO_HOME := /home/graves/_chombo4/Chombo_4/example/..
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

CXXALLFLAGS := 
CXXGPUFLAGS := --compiler-options='$(CXXALLFLAGS)'

NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))


LIBS = $(LAPACKLIBS)

CXXALLFLAGS := 
CXXGPUFLAGS := --compiler-options='$(CXXALLFLAGS)'

CPPFLAGS = -std=c++11
PROTO_HOME =$(HOME)/_proto/proto
EFENCEFLAGS = -L/usr/lib -lefence
PETSC_DIR   = $(HOME)
include $(PETSC_DIR)/lib/petsc/conf/variables
DEPXX=g++
CXX=g++

CXXALLFLAGS := -std=c++11 -Wno-sign-compare -Wno-misleading-indentation -Wno-unknown-pragmas -Wno-format
CXXDBGFLAGS:=$(CXXALLFLAGS) -g -pedantic -Wno-variadic-macros
CXXOPTFLAGS:=$(CXXALLFLAGS) -O3

#LAPACKLIBS =  -lz -llapack_atlas -lblas -llapack 
LAPACKLIBS  =  -lgfortran -lm  -llapack_atlas -lblas -llapack -lrt -lgfortran
HDFINCFLAGS =  -I/usr/include/hdf5/serial
HDFLIBFLAGS = -L/usr/lib/x86_64-linux-gnu  -lhdf5_serial
#HDFLIBFLAGS = -L/usr/lib/x86_64-linux-gnu -lhdf5_hl -lhdf5_serial -ldl

CPPFLAGS += -DCH_USE_HDF5 

CPPFLAGS += $(HDFINCFLAGS)

LIBS     += $(HDFLIBFLAGS)
 
CXXFLAGS += -DCH_USE_DOUBLE

CXXFLAGS+=$(CXXDBGFLAGS) 

CPPFLAGS +=-DDIM=$(DIM) -I$(PROTO_HOME)/include   -I. -I../src 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BaseTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BoxTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/AMRTools
CPPFLAGS += -DCH_SPACEDIM=$(DIM) -DCH_USE_64 
CPPFLAGS += -I../../common
CPPFLAGS += -I$(CHOMBO_HOME)/src/EBTools -I$(CHOMBO_HOME)/src/EBProto



LINK=$(CXX)
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
	$(DEPXX) $(DEPFLAGS) $(CPPFLAGS) $< 
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
