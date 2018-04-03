SOURCES:= ../../mdatom
BUILD:= ${SOURCES}/build
HELPERS:= ../../util


.PHONY: all
all: mdatom_float mdatom_long ${BUILD}/mdatom mdatom_double 
	./calculate.sh

${BUILD}/mdatom: ${BUILD}
	@cd ${BUILD}; cmake ..; make;

${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;

mdatom_float:
	./compile_float.sh ${SOURCES}

mdatom_long:
	./compile_long_double.sh ${SOURCES}

mdatom_double:
	@cp ${BUILD}/mdatom mdatom_double


.PHONY: clean
clean:
	@rm -r build_* src_* calculations_dir mdatom_* 2>/dev/null || true

