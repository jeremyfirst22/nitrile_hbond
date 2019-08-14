##Author: J. First -- 14 Aug 2019 
#
# This make file is intended to build the gmx analysis tool g_nitrile hbond
# The variable GMX_PATH *MUST* be manually edited to point to the location of the Gromacs (v4) installation
#
NAME=g_nitrile_hbond
GMX_PATH=/usr/local/gromacs

#what should be done by default
all: $(NAME)

CPPFLAGS=-I$(GMX_PATH)/include
LDFLAGS=-L$(GMX_PATH)/lib -lgmx

#generate a list of object (.o) files
OBJS=$(patsubst %.cpp,%.o,$(NAME).cpp $(EXTRA_SRC))
CC=g++

#main program depend on all objects, rest is done by implicit rules
$(NAME): $(OBJS)
	$(CC) $(LDFLAGS) $(CPPFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(LDFLAGS) $(CPPFLAGS) -c $<

test:
	bash test.sh 

clobber:
	rm -f $(wildcard *.o) $(wildcard *.xvg)

#clean up rule
clean:
	rm -f $(NAME) $(OBJS)

#all, clean are phony rules, e.g. they are always run
.PHONY: all 
