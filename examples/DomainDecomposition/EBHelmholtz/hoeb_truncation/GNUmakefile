
CHOMBO_HOME := /autofs/nccs-svm1_home1/dwaters/chombo_4/examples/..
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

SRC += $(wildcard ../../hoeb_common/*.cpp) 

VPATH +=  ../../hoeb_common/

CPPFLAGS = -std=c++11

PROTO_HOME =$(HOME)/proto

HDFINCFLAGS=-I$(OLCF_HDF5_ROOT)/include
HDFLIBFLAGS=-L$(OLCF_HDF5_ROOT)/lib -lhdf5 -lz
#end mach dependent





CPPFLAGS += -DCH_USE_HDF5 

CPPFLAGS += $(HDFINCFLAGS)

LIBS     += $(HDFLIBFLAGS)
 
CXXFLAGS += -DCH_USE_DOUBLE

CXXGPUFLAGS:= -G -DTHRUST_DEBUG -Xcompiler '$(CXXDBGFLAGS)'

CXXFLAGS+=$(CXXGPUFLAGS) 

CXX = nvcc 
 DEP=nvcc

LIBS+= -lnvToolsExt -lgomp

CXXFLAGS+=-DPROTO_CUDA -x cu --expt-relaxed-constexpr --expt-extended-lambda -arch=compute_70

CPPFLAGS+= -DPR_STACK_ALLOC=4294967296
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

CPPFLAGS += -I../../hoeb_common 
