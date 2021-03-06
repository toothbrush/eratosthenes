all: par seq bench report

clean: force-look
	@echo ==== Cleaning package ====
	(cd par; make clean)
	(cd seq; make clean)
	(cd bench; make clean)
	(cd report; make clean)

par: force-look
	(cd par; make all)

seq: force-look
	(cd seq; make all)

report: force-look
	(cd report; make all)

bench: force-look
	(cd bench; make all)

force-look:
	true
