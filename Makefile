all: par seq bench

clean:
	@echo ==== Cleaning package ====
	(cd par; make clean)
	(cd seq; make clean)
	(cd bench; make clean)

par:
	(cd par; make all)

seq:
	(cd seq; make all)

bench:
	(cd bench; make all)
