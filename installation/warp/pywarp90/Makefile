# Possible targets: (install is the default)
#     install, debug, clean, cleandebug: for serial, python2 or 3
#     pinstall, pdebug, pclean, pcleandebug: for parallel, python2 or 3
#     installall, debugall, cleanall: for serial and parallel, python2 or 3
#     install3, debug3, clean3, clean3debug: for serial, python3 (obsolete, saved for backward compatibility)
#     pinstall3, pdebug3, pclean3, pclean3debug: for parallel, python3 (obsolete, saved for backward compatibility)

# For now, assume using python (instead of python3)
PYTHON = python
-include Makefile.local

PYTHON_VERSION := $(shell $(PYTHON) -c "import sys;print(sys.version_info.major)")

# Check if python ran successfully, whether PYTHON_VERSION was given a value.
ifeq ($(PYTHON_VERSION), )
	# Running python failed.
    ifeq ($(PYTHON), python)
	    # PYTHON was not set by the user, so try python3
        PYTHON = python3
	    # Note that the version is fetched again to make sure that python3 can be found.
        PYTHON_VERSION := $(shell $(PYTHON) -c "import sys;print(sys.version_info.major)")
    endif
endif

install:
ifeq ($(PYTHON_VERSION), 2)
	make -f Makefile.Forthon install PYTHON=$(PYTHON)
else ifeq ($(PYTHON_VERSION), 3)
	make -f Makefile.Forthon3 install PYTHON=$(PYTHON)
else
	@echo $(PYTHON) executable not found
endif

pinstall:
ifeq ($(PYTHON_VERSION), 2)
	make -f Makefile.Forthon.pympi install PYTHON=$(PYTHON)
else ifeq ($(PYTHON_VERSION), 3)
	make -f Makefile.Forthon3.pympi install PYTHON=$(PYTHON)
else
	@echo $(PYTHON) executable not found
endif

debug:
ifeq ($(PYTHON_VERSION), 2)
	make -f Makefile.Forthon install DEBUG='-g --farg "-O0"' BUILDBASEDIR=builddebug PYTHON=$(PYTHON)
else ifeq ($(PYTHON_VERSION), 3)
	make -f Makefile.Forthon3 install DEBUG='-g --farg "-O0"' BUILDBASEDIR=build3debug PYTHON=$(PYTHON)
else
	@echo $(PYTHON) executable not found
endif

pdebug:
ifeq ($(PYTHON_VERSION), 2)
	make -f Makefile.Forthon.pympi install DEBUG='-g --farg "-O0"' BUILDBASEDIR=buildparalleldebug PYTHON=$(PYTHON)
else ifeq ($(PYTHON_VERSION), 3)
	make -f Makefile.Forthon3.pympi install DEBUG='-g --farg "-O0"' BUILDBASEDIR=build3paralleldebug PYTHON=$(PYTHON)
else
	@echo $(PYTHON) executable not found
endif

clean:
ifeq ($(PYTHON_VERSION), 2)
	make -f Makefile.Forthon clean
else ifeq ($(PYTHON_VERSION), 3)
	make -f Makefile.Forthon3 clean
else
	@echo $(PYTHON) executable not found
endif

pclean:
ifeq ($(PYTHON_VERSION), 2)
	make -f Makefile.Forthon.pympi clean
else ifeq ($(PYTHON_VERSION), 3)
	make -f Makefile.Forthon3.pympi clean
else
	@echo $(PYTHON) executable not found
endif

cleandebug:
ifeq ($(PYTHON_VERSION), 2)
	make -f Makefile.Forthon clean BUILDBASEDIR=builddebug
else ifeq ($(PYTHON_VERSION), 3)
	make -f Makefile.Forthon3 clean BUILDBASEDIR=build3debug
else
	@echo $(PYTHON) executable not found
endif

pcleandebug:
ifeq ($(PYTHON_VERSION), 2)
	make -f Makefile.Forthon.pympi clean BUILDBASEDIR=buildparalleldebug
else ifeq ($(PYTHON_VERSION), 3)
	make -f Makefile.Forthon3.pympi clean BUILDBASEDIR=build3paralleldebug
else
	@echo $(PYTHON) executable not found
endif

installall: install pinstall

debugall: debug pdebug

cleanall: clean pclean cleandebug pcleandebug
	rm -rf build* ../scripts/build*

install3:
	make -f Makefile.Forthon3 install

pinstall3:
	make -f Makefile.Forthon3.pympi install

debug3:
	make -f Makefile.Forthon3 install DEBUG='-g --farg "-O0"' BUILDBASEDIR=build3debug

pdebug3:
	make -f Makefile.Forthon3.pympi install DEBUG='-g --farg "-O0"' BUILDBASEDIR=build3paralleldebug

clean3:
	make -f Makefile.Forthon3 clean

pclean3:
	make -f Makefile.Forthon3.pympi clean

clean3debug:
	make -f Makefile.Forthon3 clean BUILDBASEDIR=build3debug

pclean3debug:
	make -f Makefile.Forthon3.pympi clean BUILDBASEDIR=build3paralleldebug

help:
	@echo "Possible targets: (install is the default)"
	@echo "    install, debug, clean, cleandebug: for serial"
	@echo "    pinstall, pdebug, pclean, pcleandebug: for parallel"
	@echo "    installall, debugall, cleanall: for serial and parallel"
	@echo "    help: this message"

