#read compiler options ..

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif 

.PHONY: all clean

aLibOBJ = $(subst .c,.o, $(wildcard *.c)) 

   all:   aLib.a 

aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

clean: 
	rm -f *.o *.a  *.so *.dll *.exe *_path.h isajet_slha so_locations 
	$(MAKE) -C nngg12     clean