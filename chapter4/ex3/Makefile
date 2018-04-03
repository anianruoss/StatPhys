SOURCES:= ../../mdatom
BUILD:= ${SOURCES}/build
HELPERS:= ../../util


.PHONY: all
all: ${BUILD}/mdatom plot.png

${BUILD}/mdatom: ${BUILD}
	@cd ${BUILD}; cmake ..; make;

${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;

plot.png: radialdistrfile
	python3 ${HELPERS}/plot_gr.py $< $@

radialdistrfile: outputfile
	python3	${HELPERS}/gr.py $< > $@

outputfile: params.inp
	${BUILD}/mdatom $< > $@

.PHONY: clean
clean:
	@rm coords.final outputfile radialdistrfile *.png 2>/dev/null || true

