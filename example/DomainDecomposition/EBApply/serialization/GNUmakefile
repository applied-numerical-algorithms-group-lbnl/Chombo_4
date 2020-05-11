
CHOMBO_HOME := /home/spencer/users/graves/_chombo4/Chombo_4/example/..
DIM := 3
VARS_OLD := $(.VARIABLES)

SRC := $(wildcard *.cpp) 
SRC += $(wildcard ../src/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BaseTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BoxTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/AMRTools/*.cpp) 

NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))


VPATH=o. ../src  ../../common $(CHOMBO_HOME)/src/BaseTools $(CHOMBO_HOME)/src/BoxTools $(CHOMBO_HOME)/src/AMRTools 
LIBS = 

CXXALLFLAGS := -Wall -Wno-unused-but-set-variable -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-99 -Wno-un
CXXGPUFLAGS := --compiler-options='$(CXXALLFLAGS)'

SRC += $(wildcard ../../common/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/EBTools/*.cpp) 

VPATH +=  $(CHOMBO_HOME)/src/EBTools 

CXXALLFLAGS := -Wall -Wno-unused-but-set-variable -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-99 -Wno-un
CXXGPUFLAGS := --compiler-options='$(CXXALLFLAGS)'

CPPFLAGS = -std=c++11

PROTO_HOME =$(HOME)/_proto/proto
DEPXX=mpiCC
CXX=mpiCC

CXXALLFLAGS := -Wall  -std=c++11 -march=native -Wno-unused-but-set-variable -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-99 -Wno-unused-local-typedefs -Wno-literal-suffix -Wno-invalid-offsetof -Wno-variadic-macros -Wno-unknown-pragmas      
CXXDBGFLAGS:= -g -pedantic -Wall $(CXXALLFLAGS)
CXXOPTFLAGS:= -O3 $(CXXALLFLAGS)

LAPACKLIBS =  -lgfortran -lm  -llapack_atlas -lblas -llapack 
HDFINCFLAGS=-I/usr/local/anag/pkg/hdf5-1.8.7.parallel/include
HDFLIBFLAGS=-L/usr/local/anag/pkg/hdf5-1.8.7.parallel/lib -lhdf5 -lz
#end mach dependent





CPPFLAGS += -DCH_MPI -DMPICH_SKIP_MPICXX -ULAM_WANT_MPI2CPP -DMPI_NO_CPPBIND
LIBS     += $(MPILIBS)
CPPFLAGS += -DCH_USE_HDF5 

CPPFLAGS += $(HDFINCFLAGS)

LIBS     += $(HDFLIBFLAGS)
 
CXXFLAGS += -DCH_USE_DOUBLE

CXXOPTFLAGS:= -O3 $(CXXALLFLAGS)

CXXFLAGS+=$(CXXOPTFLAGS) 

CPPFLAGS += -I../../common
CPPFLAGS += -I$(CHOMBO_HOME)/src/EBTools


CPPFLAGS +=-DDIM=$(DIM) -I$(PROTO_HOME)/include  -I$(PROTO_HOME)/EBProto/include -I. -I../src 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BaseTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BoxTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/AMRTools
CPPFLAGS += -DCH_SPACEDIM=$(DIM) -DCH_USE_64 



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
	rm -rf d o *.o *.exe *.d ../../common/*.o ../src/*.o ../../common/*.d ../src/*.d pout.* *.time.table.* *.hdf5
