# This file is based on the examples/Makefile from qpOASES3.0-beta

# select your operating system here!
include ../src/make_linux.mk
#include ../src/make_windows.mk


##
##	flags
##

IDIR        =  ../include
SRCDIR      =  ../src
EXTRAS_DIR  =  extras

IFLAGS      =  -I. \
               -I${IDIR}

QPOASES_EXES = \
	main_aircraft${EXE} \


QPOASES_HEADERS = \
	${IDIR}/SQProblem.hpp \
	${IDIR}/QProblem.hpp \
	${IDIR}/QProblemB.hpp \
	${IDIR}/Bounds.hpp \
	${IDIR}/Constraints.hpp \
	${IDIR}/SubjectTo.hpp \
	${IDIR}/Indexlist.hpp \
	${IDIR}/CyclingManager.hpp \
	${IDIR}/Utils.hpp \
	${IDIR}/Constants.hpp \
	${IDIR}/Types.hpp \
	${IDIR}/MessageHandling.hpp \
	${IDIR}/${EXTRAS_DIR}/SolutionAnalysis.hpp \
	${IDIR}/${EXTRAS_DIR}/OQPinterface.hpp



##
##	targets
##

all: ${QPOASES_EXES}

lcmpc:
	make -C cmpc

main_aircraft${EXE}: main_aircraft.${OBJEXT} lcmpc ../src/libqpOASES.a ../src/libqpOASESextras.a ./cmpc/libcmpc.a
	${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LIB} ${QPOASES_EXTRAS_LIB}  ${LIB_LAPACK} ${LIB_BLAS} -Lcmpc -lcmpc

clean:
	${RM} -f *.${OBJEXT} ${QPOASES_EXES}

clobber: clean


%.${OBJEXT}: %.cpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} ${IFLAGS} ${CPPFLAGS} -c $<


##
##	end of file
##
