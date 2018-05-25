SOURCES:= ../../mdatom
LIB := ${SOURCES}/lib
BUILD:= ${SOURCES}/build
HELPERS:= ../../util


.PHONY: base
base: ${BUILD}/mdatom params.inp
	@${BUILD}/mdatom params.inp coords.inp > out.dat
	@jmol coords.traj

.PHONY: harmonic
harmonic: ${BUILD}/mdatomHarmonic params.inp
	${BUILD}/mdatomHarmonic params.inp coords.inp > out.dat
	@jmol coords.traj

.PHONY: shake
shake: ${BUILD}/mdatomShake params.inp
	${BUILD}/mdatomShake params.inp coords.inp > out.dat
	@jmol coords.traj


${BUILD}/mdatom: ${BUILD} ${wildcard ${LIB}/*} ${SOURCES}/main.cpp ${SOURCES}/CMakeLists.txt
	@cd ${BUILD}; cmake ..; make mdatom;

${BUILD}/mdatomHarmonic: ${BUILD} ${wildcard ${LIB}/*} ${SOURCES}/main.cpp ${SOURCES}/CMakeLists.txt
	@cd ${BUILD}; cmake ..; make mdatomHarmonic;

${BUILD}/mdatomShake: ${BUILD} ${wildcard ${LIB}/*} ${SOURCES}/main.cpp ${SOURCES}/CMakeLists.txt
	@cd ${BUILD}; cmake ..; make mdatomShake;


${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;


params.inp: mdpar.ini
	@python ${HELPERS}/pargen.py $< > $@


.PHONY: clean
clean:
	@rm *.final *.traj params.inp 2>/dev/null || true
	@rm -r ${BUILD} 2>/dev/null || true

