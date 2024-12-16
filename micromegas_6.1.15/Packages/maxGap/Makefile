include ../../CalcHEP_src/FlagsForMake

../../lib/maxGap.so:
ifneq ($(FC),)
	./makePath $(CURDIR) > pathData.f
	$(FC) $(FFLAGS)  -shared  -o $(CALCHEP)/../lib/maxGap.so -ffixed-line-length-none -fPIC *.f         
else
	$(CC) $(CFLAGS)  -shared  -o  $(CALCHEP)/../lib/maxGap.so -fPIC  maxGapDummy.c 
endif
