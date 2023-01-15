NCCX    = g++ 
CFLAGS  = --std=c++11 -w
LDFLAGS = -lm -lGL -lglut -lGLU -g -fopenmp
#INCDIR  = ./inc
#SRCDIR  = ./src
#SOURCES = $(SRCDIR)/test.cpp $(SRCDIR)/frepshader.cpp $(SRCDIR)/helper.cpp
#HEADERS = $(INCDIR)/helper.hpp $(INCDIR)/frepshader.hpp 
SOURCES = test.cpp opengl_setup.cpp fluid_dynamics.cpp
HEADERS = opengl_setup.hpp fluid_dynamics.hpp 

all: test.out
	
test.out: $(HEADERS) $(SOURCES)
	$(NCCX) $(SOURCES) -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f test.out
