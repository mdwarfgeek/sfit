PYTHON?=python
LIB?=sfit*.so

SRCS=sfit.c sfit.h sincos.h sysinfo.c wrap.c

all: $(LIB)

$(LIB): $(SRCS)
	$(PYTHON) setup.py build_ext --inplace

test: $(LIB)
	$(PYTHON) test.py

clean:
	rm -f $(LIB)
	rm -rf build
