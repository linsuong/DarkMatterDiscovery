CALCHEP:=/Users/linusong/Documents/DMDiscovery/calchep_3.9.1
WORK:=/Users/linusong/Documents/DMDiscovery/Work
include $(CALCHEP)/FlagsForMake
include EXTLIBmake
ifneq ($(SONAME),)
   soname_VandP := $(SONAME)../../../VandP.so
else
   soname_VandP :=
endif

ifneq ($(LHAPDFPATH),)
  EXTLIB := -L$(LHAPDFPATH) -lLHAPDF $(EXTLIB)
endif

VandP.so:VandP.c $(DEPEND) models/extlib7.mdl
	$(CC) $(CFLAGS)  -shared -o VandP.so $(soname_VandP) VandP.c $(CALCHEP)/include/VandPgate.c  $(EXTLIB) $(CALCHEP)/lib/dynamic_me.a  $(CALCHEP)/lib/libSLHAplus.a $(CALCHEP)/lib/ntools.a $(CALCHEP)/lib/serv.a $(CALCHEP)/lib/dummy.a  -lm $(lQuad) -lpthread
