CC = g++
CXX = g++
DEBUG = -g
OPTFLAGS = -O2
LDFLAGS = -Wall $(DEBUG) $(OPTFLAGS)
LDLIBS = -lbbhutil
CFLAGS = -std=c++11 -c $(LDFLAGS) $(LDLIBS)
CXXFLAGS = -std=c++11 -c $(LDFLAGS) $(LDLIBS)

p2 : p2.o

p2.o : fda-fns.h fda-io.h

p2-ctest : p2-ctest.o

p2-ctest.o : fda-io.h

.PHONY : clean
clean :
	rm -f *.o *~
