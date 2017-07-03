bindir=../../bin/

.PHONY: run clean

# Build and run Sphere reference test case.
run:
			$(bindir)/preProc .
			$(bindir)/solver .
			$(bindir)/postProc .

# Delete results of the computations
clean:
			rm -rf mesh
			rm -rf results/
			rm -rf Normalvelocities.dat
