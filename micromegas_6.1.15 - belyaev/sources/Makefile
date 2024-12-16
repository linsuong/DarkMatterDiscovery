  
.PHONY: all clean

ifneq ($(MAKECMDGOALS),clean)
include ../CalcHEP_src/FlagsForMake
endif

CFLAGS+= -I../include
all:../lib/micromegas.a ../lib/dummy.a ../lib/maxGap.so

special=improveCS.c

OBJ = $(notdir $(subst .c,.o,$(filter-out $(special),$(wildcard *.c))))


#ifneq ($(FC),)
#   override OBJ += faux.o
#endif
../lib/micromegas.a(cs22.o):../Data/et_tab.inc

../lib/micromegas.a:../lib/micromegas.a($(OBJ)) 
	$(RANLIB) ../lib/micromegas.a 
	rm -f *.o

../lib/dummy.a:../lib/dummy.a(improveCS.o)
	rm -f improveCS.o
	
../lib/maxGap.so:
	$(MAKE) -C ../Packages/maxGap

clean:
	rm -f ../lib/*.a ../lib/*.so  *.o ._*  */._*
	rm -f maxGap/pathData.f
