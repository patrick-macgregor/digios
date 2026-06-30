CC=g++
CFLAG= -std=c++17 -O2 -pthread  -lgsl -lgslcblas -lm
# CFLAG= -std=c++17 -g -O0
# Note spaces here instead of tabs - VERY IMPORTANT!
ifeq ($(shell uname), Darwin)
    # Test if using homebrew or macports
    ifneq ($(shell command -v brew 2>/dev/null),)
        $(info Assuming GSL installed via homebrew...)
        CFLAG += -w -I/opt/homebrew/opt/gsl/include -L/opt/homebrew/opt/gsl/lib
    else ifneq ( $(shell command -v port 2>/dev/null),)
        $(info Assuming GSL installed via MacPorts...)
        CFLAG += -w -I/opt/local/include/ -L/opt/local/lib
    else
        $(error Could not find GSL library! Exiting...)
    endif
endif
ROOTFLAG=`root-config --cflags --glibs`

pLIBS = $(shell pkg-config --libs arrow parquet)
pINCLUDE = $(shell pkg-config --cflags arrow parquet)

all: EventBuilder EventBuilder_Q EventBuilder_S EventBuilder_A

EventBuilder: EventBuilder.cpp BinaryReader.h Hit.h
	$(CC) $(CFLAG) EventBuilder.cpp -o EventBuilder ${ROOTFLAG}

EventBuilder_Q: EventBuilder_Q.cpp BinaryReader.h Hit.h
	$(CC) $(CFLAG) EventBuilder_Q.cpp -o EventBuilder_Q ${ROOTFLAG}

EventBuilder_S: EventBuilder_S.cpp BinaryReader.h Hit.h ../working/GeneralSortMapping.h
	$(CC) EventBuilder_S.cpp -o EventBuilder_S $(CFLAG) ${ROOTFLAG}

EventBuilder_A: EventBuilder_A.cpp BinaryReader.h Hit.h ../working/GeneralSortMapping.h
	$(CC) EventBuilder_A.cpp -o EventBuilder_A $(CFLAG) ${ROOTFLAG}

clean:
	-rm EventBuilder EventBuilder_Q EventBuilder_S EventBuilder_A
