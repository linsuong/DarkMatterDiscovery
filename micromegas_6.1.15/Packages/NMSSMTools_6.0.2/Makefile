.PHONY: clean init

ifeq ($(MAKECMDGOALS),init)
init:  micromegas/NMSSM/lib/lpath.h
	$(MAKE) -C sources
micromegas/NMSSM/lib/lpath.h:
	echo \#define LPATH \"`pwd`\" >  ../../NMSSM/lib/lpath.h
else 

ifeq ($(MAKECMDGOALS),clean)
clean:
	$(MAKE) -C sources clean
	$(MAKE) -C main clean
	-find . -name '*.DS_Store' -print | xargs rm -rf
	-find . -name '*.dSYM' -print | xargs rm -rf
	-find . -name '._*' -print | xargs rm -rf
else

micrO =  ../../

ifeq ($(wildcard $(micrO)/CalcHEP_src/FlagsForMake),)
$(error Use '[g]make init' for initialization)
endif

all:
	$(MAKE) -C main

endif

endif
