all	: mpihydro

mpihydro:
	$(MAKE) -C src
	cp -f src/mpihydro ./

clean:
	$(MAKE) -C src clean
	rm -f mpihydro

distclean:
	$(MAKE) -C src distclean
	rm -f mpihydro
