SOURCES:= ../../mdatom
BUILD:= ${SOURCES}/build
HELPERS:= ../../util


.PHONY: all
all: ${BUILD}/mdatom mdatom
	./calculate.sh

${BUILD}/mdatom: ${BUILD}
	@cd ${BUILD}; cmake ..; make;

${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;

mdatom:
	@cp ${BUILD}/mdatom $@


.PHONY: clean
clean:
	@rm -r calculations_dir mdatom 2>/dev/null || true

