
CHOMBO_HOME := /home/spencer/users/graves/_chombo4/Chombo_4/example/..
DIM := 2

SRC := $(wildcard *.cpp) 
SRC += $(wildcard ../src/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BaseTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BoxTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/AMRTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/EBTools/*.cpp) 

NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))

NODIRDEP = $(subst .cpp,.d,$(NODIRSRC))
DEPENDS=$(patsubst %.d,d/%.d, $(NODIRDEP))

VPATH=o. ../src  ../../common $(CHOMBO_HOME)/src/BaseTools $(CHOMBO_HOME)/src/BoxTools $(CHOMBO_HOME)/src/AMRTools $(CHOMBO_HOME)/src/EBTools 
LIBS = 





VARS_OLD := $(.VARIABLES)
PROTO_HOME = $(HOME)/_proto/proto
PROTO_HOME =$(HOME)/_proto/proto
DEPXX=g++
CXX=g++

CXXALLFLAGS := -Wall  -std=c++11 -march=native -Wno-unused-but-set-variable -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-99 -Wno-unused-local-typedefs -Wno-literal-suffix -Wno-invalid-offsetof -Wno-variadic-macros -Wno-unknown-pragmas      
CXXDBGFLAGS:=$(CXXALLFLAGS) -g -pedantic -Wall
CXXOPTFLAGS:=$(CXXALLFLAGS) -O3

#LAPACKLIBS =  -lz -llapack_atlas -lblas -llapack 
#LAPACKLIBS =  -lgfortran -lm  -llapack_atlas -lblas -llapack 
HDFINCFLAGS=-I/usr/local/anag/pkg/hdf5-1.8.7.serial/include 
HDFLIBFLAGS=-L/usr/local/anag/pkg/hdf5-1.8.7.serial/lib -lhdf5_hl -lhdf5 -lz 
EFENCEFLAGS = -L/usr/lib -lefence

PETSC_DIR   = /usr/local/anag/pkg/petsc-3.11.3
SLEPC_DIR   = /usr/local/anag/pkg/slepc-3.7.4

#end mach dependent





CPPFLAGS += -DCH_USE_HDF5
CPPFLAGS += $(HDFINCFLAGS)
LIBS     += $(HDFLIBFLAGS) 
CXXFLAGS += -DCH_USE_DOUBLE

CXXFLAGS+=$(CXXDBGFLAGS) 

CPPFLAGS +=-DDIM=$(DIM) -I$(PROTO_HOME)/include  -I$(PROTO_HOME)/EBProto/include -I. -I../src 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BaseTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BoxTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/AMRTools
CPPFLAGS += -I$(CHOMBO_HOME)/src/EBTools
CPPFLAGS += -DCH_SPACEDIM=$(DIM) -DCH_USE_64 

DEP=$(DEPXX) -MM $(CXXFLAGS)
ifeq ($(DEPXX),nvcc)
DEP=$(DEPXX) -M $(CXXFLAGS)
LIBS+= -lnvToolsExt
endif

LINK=$(CXX)
OBJDIR = o
DEPDIR = d
TARGET = main.exe

$(TARGET): $(OBJ)
	$(LINK)  $(OBJ) $(LIBS) -o $(TARGET)

$(OBJDIR):
	if [ ! -d "./$(OBJDIR)" ];then     \
		mkdir $(OBJDIR);           \
	fi

$(DEPDIR):
	if [ ! -d "./$(DEPDIR)" ];then     \
		mkdir $(DEPDIR);           \
	fi

#compile non-main sources
$(OBJDIR)/%.o : %.cpp   $(OBJDIR) $(DEPDIR)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(CPPFLAGS) 
#	$(DEP) $<   $(CPPFLAGS) > $(DEPDIR)/$*.d

all: $(TARGET)

.PHONY: clean realclean print dirs


print-%: ; @echo $* = $($*)

-include $(DEPENDS)

vars:
	$(foreach v,$(filter-out $(VARS_OLD) VARS_OLD,$(.VARIABLES)), $(info $(v) = $($(v))))

clean:
	rm -rf d o *.o *.exe *.d ../../common/*.o ../src/*.o ../../common/*.d ../src/*.d

realclean:
	rm -rf d o *.o *.exe *.d ../../common/*.o ../src/*.o ../../common/*.d ../src/*.d

pristine:
	rm -rf d o *.o *.exe *.d ../../common/*.o ../src/*.o ../../common/*.d ../src/*.d pout.* *.time.table.* *.hdf5
