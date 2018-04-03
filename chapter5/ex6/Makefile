SOURCES:= ../../mdatom
BUILD:= ${SOURCES}/build
HELPERS:= ../../util


.PHONY: all
all: ${BUILD}/mdatom mdatom gr.py
	./calculate.sh

${BUILD}/mdatom: ${BUILD}
	@cd ${BUILD}; cmake ..; make;

${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;

mdatom:
	@cp ${BUILD}/mdatom $@

gr.py:
	@cp ${HELPERS}/gr.py $@
	@chmod +x $@


.PHONY: clean
clean:
	@rm -r calculations_dir mdatom gr.py 2>/dev/null || true

