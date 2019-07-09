
CHOMBO_HOME := /home/spencer/users/graves/_chombo4/Chombo_4/example/..
DIM := 2

SRC := $(wildcard *.cpp) 
SRC += $(wildcard ../src/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BaseTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/BoxTools/*.cpp) 
SRC += $(wildcard $(CHOMBO_HOME)/src/AMRTools/*.cpp) 

NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))

NODIRDEP = $(subst .cpp,.d,$(NODIRSRC))
DEPENDS=$(patsubst %.d,d/%.d, $(NODIRDEP))

VPATH=o. ../src  ../../common $(CHOMBO_HOME)/src/BaseTools $(CHOMBO_HOME)/src/BoxTools $(CHOMBO_HOME)/src/AMRTools 







PROTO_HOME =$(HOME)/_proto/proto
DEPXX=mpiCC
CXX=mpiCC
DIM=2
CXXALLFLAGS := -Wall  -std=c++11 -march=native -Wno-unused-but-set-variable -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-99 -Wno-unused-local-typedefs -Wno-literal-suffix -Wno-invalid-offsetof -Wno-variadic-macros -Wno-unknown-pragmas      
CXXDBGFLAGS:= -g -pedantic -Wall $(CXXALLFLAGS)
CXXOPTFLAGS:= -O3 $(CXXALLFLAGS)

LAPACKLIBS =  -lgfortran -lm  -llapack_atlas -lblas -llapack 
HDFINCFLAGS=-I/usr/local/anag/pkg/hdf5-1.8.7.parallel/include
HDFLIBFLAGS=-L/usr/local/anag/pkg/hdf5-1.8.7.parallel/lib -lhdf5 -lz
#end mach dependent





CPPFLAGS += -DCH_MPI -DMPICH_SKIP_MPICXX -ULAM_WANT_MPI2CPP -DMPI_NO_CPPBIND
CPPFLAGS += -DCH_USE_HDF5
CXXFLAGS += -DCH_USE_DOUBLE
CXXFLAGS+=$(CXXOPTFLAGS) 
CPPFLAGS +=-DDIM=$(DIM) -I$(PROTO_HOME)/include  -I$(PROTO_HOME)/EBProto/include -I. -I../src 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BaseTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/BoxTools 
CPPFLAGS += -I$(CHOMBO_HOME)/src/AMRTools
CPPFLAGS += -DCH_SPACEDIM=$(DIM) -DCH_USE_64 
CPPFLAGS += $(HDFINCFLAGS) 

LIBS     += $(HDFLIBFLAGS) 

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
	$(LINK)  $(OBJ) $(LIBS) $(LAPACKLIBS) -o $(TARGET)

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



