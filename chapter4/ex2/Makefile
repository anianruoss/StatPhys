SOURCES:= ../../mdatom
BUILD:= ${SOURCES}/build
HELPERS:= ../../util


.PHONY: all
all: ${BUILD}/mdatom 
	${BUILD}/mdatom params.inp coords.inp

${BUILD}/mdatom: ${BUILD}
	@cd ${BUILD}; cmake ..; make;

${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;

.PHONY: clean
clean:
	@rm *.final *.traj 2>/dev/null || true

