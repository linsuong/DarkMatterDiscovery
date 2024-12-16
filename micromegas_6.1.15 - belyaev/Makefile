
.PHONY: all clean flags


all:include/microPath.h
	$(MAKE) -C CalcHEP_src  MICROMEGAS=MICROMEGAS
	$(MAKE) -C sources

include/microPath.h:
	echo \#define micrO \"$(CURDIR)\"  > include/microPath.h
   
clean:  
	rm -f  include/microPath.h
	rm -f ._*
	rm -f Packages/._*
	rm -f man/._*
	rm -f Data/._* Data/*/._* Data/*/*/._*
	rm -f include/._*
	rm -fr lib/*.a lib/*.so  lib/*.o lib/._*  lib/*.dSYM
	./clean
flags: 
	$(MAKE) -C CalcHEP_src flags