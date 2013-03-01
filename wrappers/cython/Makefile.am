if WITH_CYTHON
if WITH_PYTHON

all-local:
	$(PYTHON) setup.py install --install-lib=build/

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(DESTDIR)$(PYTHON_PATH)/LHAPDF-$(PACKAGE_VERSION)*.egg-info

clean-local:
	$(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc
	@rm -rf $(builddir)/build
	@rm -rf dist

distclean-local:
	@rm -rf LHAPDF.egg-info

endif
endif
