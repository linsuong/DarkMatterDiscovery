.PHONY: all

include ../../FlagsForMake

Idir = .. chep_crt plot service2 getmem dynamicME SLHAplus ntools num

override CFLAGS += -Iinclude $(patsubst %,-I../%/include ,$(Idir))

special = pdf_dummy.c lhapdf2pdt.c  lhapdf.c

OBJ = $(notdir $(subst .c,.o,$(filter-out $(special),$(wildcard *.c ))))

#VPATH = strfun

lib = $(CALCHEP)/lib
bin = ../../bin

all: $(lib)/num_c.a  $(bin)/lhapdf2pdt  $(lib)/lhapdf.so

$(lib)/num_c.a:$(lib)/num_c.a($(OBJ))
	$(RANLIB) $(lib)/num_c.a
	@rm -f $(OBJ)

$(lib)/lhapdf.so: ../../FlagsForMake lhapdf.c lhapdf.cc
ifeq ($(LHAPDF),)
	$(CC) $(CFLAGS)  -shared -o $(lib)/lhapdf.so lhapdf.c
else
	$(CXX) $(CXXFLAGS)   -shared  -o $(lib)/lhapdf.so   -I$(LHAPDF)/include lhapdf.cc $(LHAPDF)/lib/libLHAPDF.so 
endif


$(bin)/lhapdf2pdt: lhapdf2pdt.c
	$(CC) $(CFLAGS)  -o $(bin)/lhapdf2pdt lhapdf2pdt.c $(lib)/ntools.a -lm

