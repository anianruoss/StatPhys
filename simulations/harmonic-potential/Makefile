SOURCES:= ../../mdatom
BUILD:= ${SOURCES}/build
HELPERS:= ../../util

.PHONY: all
all: ${BUILD}/mdatom params.inp
	@${BUILD}/mdatom params.inp coords.inp
	@jmol coords.traj

${BUILD}/mdatom: ${BUILD}
	@cd ${BUILD}; cmake ..; make;

${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;

params.inp: mdpar.ini
	@python ${HELPERS}/pargen.py $< > $@

.PHONY: clean
clean:
	@rm *.final *.traj params.inp 2>/dev/null || true

