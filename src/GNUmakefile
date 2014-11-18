# ===========================================================================
#  Makefile iS                                     Chun Shen Apr. 9, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS= -O3

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	HBT.e
endif

SRC		=	main.cpp readindata.cpp HBT.cpp arsenal.cpp \
			gauss_quadrature.cpp ParameterReader.cpp Table.cpp

INC		= 	parameters.h arsenal.h HBT.h ParameterReader.h \
                  readindata.h Stopwatch.h Table.h gauss_quadrature.h

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	../

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
./main.cpp: readindata.h HBT.h arsenal.h parameters.h Stopwatch.h ParameterReader.h
./arsenal.cpp: arsenal.h gauss_quadrature.h
./gauss_quadrature.cpp: gauss_quadrature.h
./HBT.cpp: HBT.h parameters.h arsenal.h readindata.h
./readindata.cpp: parameters.h readindata.h ParameterReader.h
./ParameterReader.cpp: ParameterReader.h arsenal.h
