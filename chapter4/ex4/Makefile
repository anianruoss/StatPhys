SOURCES:= ../../mdatom
BUILD:= ${SOURCES}/build
HELPERS:= ../../util

EX3:= ../ex3


.PHONY: all
all: ${BUILD}/mdatom ${EX3}/coords.final plot.png 

${BUILD}/mdatom: ${BUILD}
	@cd ${BUILD}; cmake ..; make;

${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;

${EX3}/coords.final: 
	@cd ${EX3}; ${BUILD}/mdatom params.inp coords.inp > outputfile;

plot.png: radialdistrfile
	python3 ${HELPERS}/plot_gr.py $< $@

radialdistrfile: outputfile
	python3	${HELPERS}/gr.py $< > $@

outputfile: params.inp ${EX3}/coords.final
	${BUILD}/mdatom $^ > $@

.PHONY: clean
clean:
	@rm coords.final outputfile radialdistrfile *.png 2>/dev/null || true
	@rm ${EX3}/coords.final ${EX3}/outputfile 2>/dev/null || true

