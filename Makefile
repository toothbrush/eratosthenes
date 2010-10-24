all: par seq bench

clean: force-look
	@echo ==== Cleaning package ====
	(cd par; make clean)
	(cd seq; make clean)
	(cd bench; make clean)

par: force-look
	(cd par; make all)

seq: force-look
	(cd seq; make all)

bench: force-look
	(cd bench; make all)

force-look:
	true
