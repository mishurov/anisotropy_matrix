CC = g++-4.6
SOURCES = SOP_anisotropy_matrix.cc
DSONAME = SOP_anisotropy_matrix.so
HVER = 13.0.237
INCDIRS = -I/opt/hfs$(HVER)/toolkit/include/
LIBDIRS = -L/opt/hfs$(HVER)/toolkit/dsolib -L/usr/X11R6/lib64 -L/usr/X11R6/lib/
SHAREDFLAG = -shared
OMPFLAG = -fopenmp
CFLAGS = -DVERSION=\"$(HVER)\" -DDLLEXPORT="" \
-D_GNU_SOURCE -DLINUX -fPIC -DFBX_ENABLED=1 -DOPENCL_ENABLED=1 \
-DOPENVDB_ENABLED=1 -DSESI_LITTLE_ENDIAN -DENABLE_THREADS -DUSE_PTHREADS \
-D_REENTRANT -D_FILE_OFFSET_BITS=64 -c -DGCC4 -DGCC3 -fno-strict-aliasing
WFLAGS = -Wno-deprecated -Wall -W -Wno-parentheses -Wno-sign-compare \
-Wno-reorder -Wno-uninitialized -Wunused -Wno-unused-parameter
OPTIMIZER = -O2
ARCHFLAGS = -DAMD64 -m64 -DSIZEOF_VOID_P=8
SYSLIBS = -lGLU -lGL -lX11 -lXext -lXi -ldl -lpthread -lgomp
MATHLIBS= -llapack

LINK ?= $(CC)
OBJFLAGS = -c $(OMPFLAG) $(CFLAGS) $(ARCHFLAGS) \
$(INCDIRS) $(WFLAGS) $(OPTIMIZER)
DSOFLAGS = $(LIBDIRS) $(SYSLIBS) $(MATHLIBS)

OBJECTS = $(SOURCES:.cc=.o)

all: $(DSONAME)
	@echo Done!
$(DSONAME): $(OBJECTS)
	$(LINK) $(CXXFLAGS) $(SHAREDFLAG) $^ $(DSOFLAGS) -o $@
%.o: %.cc
	$(CC) $(OMPFLAG) $(CXXFLAGS) $(OBJFLAGS) -DMAKING_DSO $< -o $@

# Exuberant ctags generation for VIM autocompletion
tags:
	$(CC) -M $(INCDIRS) $(SOURCES) | \
	sed -e 's/[\\ ]/\n/g' | \
	sed -e '/^$$/d' -e '/\.o:[ \t]*$$/d' | \
	ctags -L - --c++-kinds=+p --fields=+iaS --extra=+q	
