SOURCES:= ../../mdatom
BUILD:= ${SOURCES}/build
HELPERS:= ../../util

EX3:= ../ex3
EX4:= ../ex4


.PHONY: all
all: ${BUILD}/mdatom ${EX4}/coords.final
	$< params.inp ${EX4}/coords.final

${BUILD}/mdatom: ${BUILD}
	@cd ${BUILD}; cmake ..; make;

${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;

${EX4}/coords.final: ${EX3}/coords.final
	@cd ${EX4}; ${BUILD}/mdatom params.inp $<;

${EX3}/coords.final: 
	@cd ${EX3}; ${BUILD}/mdatom params.inp;

.PHONY: clean
clean:
	@rm coords.final outputfile radialdistrfile *.png 2>/dev/null || true
	@rm ${EX3}/coords.final ${EX4}/coords.final 2>/dev/null || true

