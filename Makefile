PROG := tightbinding
SHELL := /bin/sh
CC := gcc
CXX := g++

COMPILER_OPTIONS := -m64 -Wall -Wextra -Wshadow -Werror -pedantic -Iinclude
CFLAGS := -std=c99 $(COMPILER_OPTIONS)
CXXFLAGS := -std=c++11 -Weffc++ $(COMPILER_OPTIONS)
LDFLAGS := -Wl,--no-as-needed -lm

DEBUGFLAGS := -g -O0 -D _DEBUG
RELEASEFLAGS := -O2 -D NDEBUG

ifeq ($(PREFIX),)
	PREFIX := /usr/local
endif

SRCDIR	:= src
INCDIR	:= include
OBJDIR 	:= obj
BINDIR	:= bin
TARDIR	:= lib

SOURCES := $(shell echo $(SRCDIR)/*.cpp)
HEADERS := $(shell echo $(INCDIR)/*.h)
COMMON  :=
OBJECTS := $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SOURCES:.cpp=.o))
HEADER 	:= $(PROG).h
TARGET  := lib$(PROG).a

TAR_INSTALLDIR = $(PREFIX)/lib
INC_INSTALLDIR = $(PREFIX)/include

all: $(TARGET)

$(TARGET): $(OBJECTS)
	ar ru $(TARDIR)/$@ $^
	ranlib $(TARDIR)/$@

test:
	cd test && $(MAKE)

release: $(SOURCES) $(HEADERS) $(COMMON)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(RELEASEFLAGS) -o $(TARGET) $(SOURCES)

install: $(TARGET)
	install -d $(TAR_INSTALLDIR)
	install -m 644 $(TARDIR)/$(TARGET) $(TAR_INSTALLDIR)
	install -d $(INC_INSTALLDIR)
	install -m 644 $(INCDIR)/$(HEADER) $(INC_INSTALLDIR)

zip:
	-zip $(PROG).zip $(HEADERS) $(SOURCES) Makefile

clean:
	-rm -f $(TARDIR)/$(TARGET) $(OBJECTS) $(PROG).zip
	-rm -f $(TAR_INSTALLDIR)/$(TARGET) $(INC_INSTALLDIR)/$(HEADER)

remake: clean all

reinstall: clean install

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) $(COMMON)
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) -c -o $@ $<

.PHONY : all release install zip clean test remake reinstall
