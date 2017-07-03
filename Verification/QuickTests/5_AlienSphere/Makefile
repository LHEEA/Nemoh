SHELL=/bin/bash

bindir=../../../bin/
logfile=log

.PHONY: run silent_run test clean

# Build and run Sphere reference test case.
run:
	$(bindir)/preProc .
	$(bindir)/solver .
	$(bindir)/postProc .

silent_run:
	@test -f $(logfile) || echo '' > $(logfile)
	@$(bindir)/preProc  . >> $(logfile)
	@$(bindir)/solver   . >> $(logfile)
	@$(bindir)/postProc . >> $(logfile)

test: silent_run
	@echo "Forces.dat"
	@cat <(echo -n "    Last result") results/Forces.dat
	@cat <(echo -n "    Nemoh2     ") reference_results/Forces.dat
	@echo "Line 20 of freesurface.00001.dat"
	@cat <(echo -n "    Last result") <(head -n 20 results/freesurface.00001.dat | tail -1)
	@cat <(echo -n "    Nemoh2     ") <(head -n 20 "reference_results/freesurface.    1.dat" | tail -1)

# Delete results of the computations
clean:
	rm -f $(logfile)
	rm -rf mesh
	rm -rf results/
	rm -rf Normalvelocities.dat
