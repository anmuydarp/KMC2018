CC=g++
CFLAGS=-Wall -pedantic
CFILE1=KMCmain.cpp Math.cpp Initialize.cpp CreateDefects.cpp Device.cpp TransitionTable.cpp TransitionRateCalc.cpp TransitionRateFunctions.cpp
CFILE2=KMCLoop.cpp KMCEvents.cpp Write_Calc.cpp Lifetime.cpp
GDB=-ggdb -gdwarf-3

quickmake:$(CFILE)
	$(CC) $(CFLAGS) $(GDB) -o KineticMC $(CFILE1) $(CFILE2)
