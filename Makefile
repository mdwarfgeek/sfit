PYTHON?=python
LIB?=sfit*.so

SRCS=sfit.c sfit.h sincos.h sysinfo.c vsincos_gen.c vsincos_sleef.c wrap.c

all: $(LIB)

$(LIB): $(SRCS)
	$(PYTHON) setup.py build_ext --inplace

test: $(LIB)
	$(PYTHON) test.py

clean:
	rm -f *.o
	rm -f _configtest.c _configtest
	rm -f test_vsincos_sleef_sse4 test_vsincos_sleef_avx
	rm -f $(LIB)
	rm -rf build
