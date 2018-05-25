SOURCES:= ../../mdatom
LIB := ${SOURCES}/lib
BUILD:= ${SOURCES}/build


.PHONY: all
all: ${BUILD}/mdatom* run_simulations.sh plot_numatoms.py plot_temperature.py
	./run_simulations.sh
	@python plot_numatoms.py
	@python plot_temperature.py

${BUILD}/mdatom*: ${BUILD} ${wildcard ${LIB}/*} ${SOURCES}/main.cpp ${SOURCES}/CMakeLists.txt
	@cd ${BUILD}; cmake ..; make;

${BUILD}:
	@cd ${SOURCES}; mkdir build 2>/dev/null || true;


.PHONY: clean
clean:
	@rm -r data *.pdf 2>/dev/null || true
	@rm -r ${BUILD} 2>/dev/null || true
