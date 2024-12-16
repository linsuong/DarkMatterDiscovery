ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif 

.PHONY: all clean

all:  suspect.exe 


ifeq ($(FC),)

suspect.exe:  suspect2.c ../lf2c/libf2c.a
	$(CC) -O -o suspect.exe $^ -lm

../lf2c/libf2c.a:
	make -C ../lf2c 

else 


suspect.exe:  suspect2_call.f bsg.f suspect2.f  twoloophiggs.f
	$(FC) -o suspect.exe $(FFLAGS) $^


endif

clean:
	rm -f suspect.exe
