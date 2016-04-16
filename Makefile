# ########################################################################
# Copyright 2013 Alexander Mishurov
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ########################################################################

CC = g++-4.8
SOURCES = SOP_anisotropy_matrix.cc
DSONAME = SOP_anisotropy_matrix.so
HVER = 15.0.416
INCDIRS = -I/opt/hfs$(HVER)/toolkit/include/
LIBDIRS = -L/opt/hfs$(HVER)/toolkit/dsolib -L/usr/X11R6/lib64 -L/usr/X11R6/lib/
SHAREDFLAG = -shared
CFLAGS = -DVERSION=\"$(HVER)\" -DDLLEXPORT="" \
-D_GNU_SOURCE -DLINUX -fPIC -DFBX_ENABLED=1 -DOPENCL_ENABLED=1 \
-DOPENVDB_ENABLED=1 -DSESI_LITTLE_ENDIAN -DENABLE_THREADS -DUSE_PTHREADS \
-D_REENTRANT -D_FILE_OFFSET_BITS=64 -c -DGCC4 -DGCC3 -fno-strict-aliasing
WFLAGS = -Wno-deprecated -Wall -W -Wno-parentheses -Wno-sign-compare \
-Wno-reorder -Wno-uninitialized -Wunused -Wno-unused-parameter
OPTIMIZER = -O2
ARCHFLAGS = -DAMD64 -m64 -DSIZEOF_VOID_P=8
SYSLIBS = -lGLU -lGL -lX11 -lXext -lXi -ldl -lpthread

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
	$(CC) $(CXXFLAGS) $(OBJFLAGS) -DMAKING_DSO $< -o $@

# Exuberant ctags generation for VIM autocompletion
tags:
	$(CC) -M $(INCDIRS) $(SOURCES) | \
	sed -e 's/[\\ ]/\n/g' | \
	sed -e '/^$$/d' -e '/\.o:[ \t]*$$/d' | \
	ctags -L - --c++-kinds=+p --fields=+iaS --extra=+q	
