
#configuration: dim=3,opt=TRUE,MPI=FALSE,HDF5=TRUE,cuda=FALSE,hip=FALSE,PETSc=FALSE,OpenMP=FALSE

CHOMBO_HOME := /home/graves/_chombo4/Chombo_4/example/..
DIM := 3
VARS_OLD := $(.VARIABLES)

CHOMBO_DEPFLAGS = 
CHOMBO_CXXALLFLAGS = 
CHOMBO_CXXOPTFLAGS = 
CHOMBO_CXXDEBFLAGS = 
CHOMBO_CXXFLAGS =

CHOMBO_SRC := $(wildcard *.cpp) 
CHOMBO_SRC += $(wildcard ../src/*.cpp) 
CHOMBO_SRC += $(wildcard $(CHOMBO_HOME)/src/BaseTools/*.cpp)
CHOMBO_SRC += $(wildcard $(CHOMBO_HOME)/src/LinearAlgebra/*.cpp)
CHOMBO_SRC += $(wildcard $(CHOMBO_HOME)/src/BoxTools/*.cpp)
CHOMBO_SRC += $(wildcard ../../common/*.cpp) 
CHOMBO_SRC += $(wildcard $(CHOMBO_HOME)/src/EBTools/*.cpp) 
CHOMBO_SRC += $(wildcard $(CHOMBO_HOME)/src/EBProto/*.cpp) 

VPATH=o. ../src  ../../common $(CHOMBO_HOME)/src/BaseTools $(CHOMBO_HOME)/src/BoxTools $(CHOMBO_HOME)/src/LinearAlgebra
VPATH += $(CHOMBO_HOME)/src/EBTools

CHOMBO_NODIRSRC = $(notdir $(CHOMBO_SRC))
CHOMBO_NODIROBJ = $(subst .cpp,.o,$(CHOMBO_NODIRSRC))
CHOMBO_OBJ=$(patsubst %.o,o/%.o,  $(CHOMBO_NODIROBJ))


CHOMBO_LIBS = $(CHOMBO_LAPACKLIBS)


CHOMBO_SRC += $(wildcard ../../hoeb_common/*.cpp) 

VPATH +=  ../../hoeb_common/

CHOMBO_CPPFLAGS = -std=c++14

PROTO_HOME = $(HOME)/_proto/proto
CHOMBO_DEPXX=g++
CHOMBO_CXX=g++

CHOMBO_CXXALLFLAGS := -std=c++11 -march=native -Wno-unused-but-set-variable -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-99 -Wno-unused-local-typedefs -Wno-literal-suffix -Wno-invalid-offsetof -Wno-variadic-macros -Wno-unknown-pragmas      
CHOMBO_CXXDBGFLAGS:=$(CHOMBO_CXXALLFLAGS) -g -pedantic 
CHOMBO_CXXOPTFLAGS:=$(CHOMBO_CXXALLFLAGS) -O3

#CHOMBO_LAPACKLIBS =  -lz -llapack_atlas -lblas -llapack 
#CHOMBO_LAPACKLIBS =  -lgfortran -lm  -llapack_atlas -lblas -llapack 
CHOMBO_HDFINCFLAGS=-I$(ANAG_HDF5_DIR)/include 
CHOMBO_HDFLIBFLAGS=-L$(ANAG_HDF5_DIR)/lib -lhdf5_hl -lhdf5 -lz -llapack -lblas
CHOMBO_EFENCEFLAGS = -L/usr/lib -lefence

CHOMBO_PETSC_DIR   = $(ANAG_PETSC_DIR)
CHOMBO_SLEPC_DIR   = /usr/local/anag/pkg/slepc-3.7.4
#end mach dependent





CHOMBO_CPPFLAGS += -DCH_USE_HDF5 

CHOMBO_CPPFLAGS += $(CHOMBO_HDFINCFLAGS)

CHOMBO_LIBS     += $(CHOMBO_HDFLIBFLAGS)
 
CHOMBO_CXXFLAGS += $(CHOMBO_CXXALLFLAGS) -DCH_USE_DOUBLE

CHOMBO_CXXFLAGS+=$(CHOMBO_CXXDBGFLAGS) 

CHOMBO_CXXOPTFLAGS:= -O2 $(CHOMBO_CXXALLFLAGS) -DNDEBUG

CHOMBO_CXXFLAGS+=$(CHOMBO_CXXOPTFLAGS) 

CHOMBO_CXXFLAGS+=$(CHOMBO_CXXOPTFLAGS) 

CHOMBO_CPPFLAGS +=-DDIM=$(DIM) -I$(PROTO_HOME)/include   -I. -I../src 
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/BaseTools
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/LinearAlgebra
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/BoxTools
CHOMBO_CPPFLAGS += -DCH_SPACEDIM=$(DIM) -DCH_USE_64 
CHOMBO_CPPFLAGS += -I../../common
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/EBTools 
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/EBProto 
CHOMBO_CPPFLAGS += -I$(CHOMBO_HOME)/src/EBCM


#need to use depxx because petsc links with g++
CHOMBO_LINK=$(CHOMBO_CXX)
CHOMBO_OBJDIR = o
CHOMBO_DEPDIR = d
CHOMBO_TARGET = main.exe

$(CHOMBO_TARGET): $(CHOMBO_OBJ)
	$(CHOMBO_LINK)  $(CHOMBO_OBJ) $(CHOMBO_LIBS) -o $(CHOMBO_TARGET)

$(CHOMBO_OBJDIR): ; mkdir -p $@


$(CHOMBO_DEPDIR): ; mkdir  -p $@

CHOMBO_DEPFILES := $(CHOMBO_NODIRSRC:%.cpp=$(CHOMBO_DEPDIR)/%.d)

$(CHOMBO_DEPFILES):

CHOMBO_DEPFLAGS += -MT $@ -MF $(CHOMBO_DEPDIR)/$*.d 

ifeq ($(CHOMBO_CXX),nvcc)
#broke the build
#LIBS+= -lnvToolsExt
endif

ifeq ($(CHOMBO_DEPXX),nvcc)
CHOMBO_DEPFLAGS += --generate-dependencies
else
CHOMBO_DEPFLAGS += -MM
endif

$(CHOMBO_OBJDIR)/%.o : %.cpp  $(CHOMBO_DEPDIR)/%.d | $(CHOMBO_DEPDIR) 
	@mkdir -p $(CHOMBO_OBJDIR)
	$(CHOMBO_DEPXX) $(CHOMBO_DEPFLAGS) $(CHOMBO_CPPFLAGS) $< 
	$(CHOMBO_CXX) -c -o $@ $< $(CHOMBO_CXXFLAGS) $(CHOMBO_CPPFLAGS) 

all: $(CHOMBO_TARGET)

-include $(CHOMBO_DEPFILES)

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

CHOMBO_CPPFLAGS += -I../../hoeb_common 
