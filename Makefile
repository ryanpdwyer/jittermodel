CFLAGS=-shared -O3 -Wall
PLIBS=-I/usr/local/Cellar/python/2.7.8_1/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/usr/local/Cellar/python/2.7.8_1/Frameworks/Python.framework/Versions/2.7/include/python2.7
NLIB=-I/usr/local/lib/python2.7/site-packages/numpy/core/include
PLDFLAGS=-L/usr/local/Cellar/python/2.7.8_1/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config -ldl -framework CoreFoundation -lpython2.7

.PHONY:
	cython clean test

cython:
	cython -a jittermodel/_sim.pyx
	$(CC) $(CFLAGS) $(PLIBS) $(NLIB) $(PLDFLAGS) \
		 -o jittermodel/_sim.so \
		jittermodel/_sim.c

test: cython
	python jittermodel/tests/test_verification/cython.py

clean:
	rm jittermodel/_sim.so
	rm jittermodel/_sim.c
	rm -rf *.pyc