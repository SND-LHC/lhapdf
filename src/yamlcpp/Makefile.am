EXTRA_DIST = yaml-cpp-0.3.0.tar.gz

yaml-cpp/liblhapdf-yaml-cpp.a: yaml-cpp-0.3.0.tar.gz
	test -d yaml-cpp || tar xf $^; true
	grep COMPILE_FLAGS yaml-cpp/CMakeLists.txt || $(SED) -i -e $$'s/LABEL_SUFFIX}"/LABEL_SUFFIX}"\\\n\\\tCOMPILE_FLAGS "-fPIC"/' yaml-cpp/CMakeLists.txt; true
	for i in `find yaml-cpp -type f -name "*.h" -o -name "*.cpp"`; do (grep LHAPDF_YAML $$i > /dev/null) || \
      $(SED) -i -e 's/namespace\ \+YAML/namespace LHAPDF_YAML/' -e 's/YAML::/LHAPDF_YAML::/g' $$i; done; true
	cd yaml-cpp && cmake . && $(MAKE) && mv libyaml-cpp.a liblhapdf-yaml-cpp.a

all: yaml-cpp/liblhapdf-yaml-cpp.a
	@true
