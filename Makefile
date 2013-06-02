#!/bin/make
SHELL=/bin/bash

#Define the virtual paths
vpath %.cpp src/
vpath %.hpp include/
vpath %.o obj/

#Set some of the compile options
CXX = g++
CXXFLAGS = -fPIC -g -std=c++11 -Wall $(CINCLUDEDIRS)
LDLIBS = -ldl -lpthread
CINCLUDEDIRS = -Iinclude
c++SrcSuf = cpp

#Set the name of the program to be compiled
PROGRAM = fitterSample
VERSION = $(shell git describe --abbrev=0 --tags)

#Define Objects
FITTINGO = fitterSample.o

#DefineSharedObject
SO_NAME = $(subst .o,.so,lib$(PIPESO))
SO_NAME_W_VERSION = $(addsuffix .$(VERSION),$(SO_NAME))
SO_NAME_FULL=$(addprefix $(OBJDIR)/,$(SO_NAME_W_VERSION))

#Define the name of the header for the SO
HEADER_NAME = $(subst .o,.hpp,$(PIPESO))

#Make the object list and prefix the object directory
OBJS = $(FITTINGO) $(FUNCO)
OBJDIR = obj
OBJS_W_DIR = $(addprefix $(OBJDIR)/,$(OBJS))

ROOTCONFIG   := root-config
CXXFLAGS     += $(shell $(ROOTCONFIG) --cflags)
LDFLAGS      += $(shell $(ROOTCONFIG) --ldflags)
LDLIBS       += $(shell $(ROOTCONFIG) --libs) -lRooFit -lRooFitCore

.SUFFIXES: .$(c++SrcSuf)

all: $(OBJS_W_DIR) $(PROGRAM)

$(OBJS_W_DIR): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

$(PROGRAM): $(OBJS_W_DIR)
	$(CXX) $(CXXFLAGS) $(LDLIBS) $^ -o $@

$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean so
clean: 
	@echo "Cleaning..."
	@rm -f $(OBJDIR)/* $(PROGRAM) *~ src/*~ include/*~

so: $(HEADER_NAME)
	g++ -shared $(OBJDIR)/$(PIPESO)  -o $(SO_NAME_FULL)
	cp $(SO_NAME_FULL) $(SO_INSTALL_PATH)
	ln -sf $(SO_INSTALL_PATH)/$(SO_NAME_W_VERSION) $(SO_INSTALL_PATH)/$(SO_NAME)
