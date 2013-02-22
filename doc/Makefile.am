#if WITH_PDF
man:
	(cd manual && make)
#endif


## Doxygen
if WITH_DOXYGEN
doxy:
	@echo "Running Doxygen..."
	@$(DOXYGEN) || true
install-doxy:
	$(install_sh_DATA) html/* $(DESTDIR)$(pkgdatadir)/doc/html
mostlyclean-local:
	@rm -rf $(top_builddir)/html $(top_builddir)/latex
uninstall-local:
	@rm -rf $(DESTDIR)$(pkgdatadir)/doc; \
	test ! -d $(DESTDIR)$(pkgdatadir) || rmdir --ignore-fail-on-non-empty $(DESTDIR)$(pkgdatadir)
endif
