EXTRA_DIST = lhapdf.cpp lhapdf.pyx clhapdf.pxd

if WITH_CYTHON
lhapdf.cpp: lhapdf.pyx clhapdf.pxd
	cython lhapdf.pyx --cplus
else
lhapdf.cpp:
	@echo "Not (re)generating lhapdf.cpp since Cython is not installed"
endif


if WITH_PYTHON

PYEXT_ENV = CC=$(CC) CXX=$(CXX)

## Always force setup.py, it's not good at guessing what needs to rebuild
all-local: lhapdf.cpp
	$(PYEXT_ENV) $(PYTHON) setup.py install --install-lib=build/ --force

install-exec-local:
	$(PYEXT_ENV) $(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix) --force

uninstall-local:
	rm -rf $(DESTDIR)$(PYTHON_PATH)/LHAPDF-$(PACKAGE_VERSION)*.egg-info

clean-local:
	$(PYEXT_ENV) $(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc
	@rm -rf $(builddir)/build
	@rm -rf dist

distclean-local:
	@rm -rf LHAPDF.egg-info

endif
