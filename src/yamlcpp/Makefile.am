EXTRA_DIST = yaml-cpp-0.3.0.tar.gz

yaml-cpp/libyaml-cpp.a: yaml-cpp-0.3.0.tar.gz
	test -d yaml-cpp || tar xf $^; true
	grep COMPILE_FLAGS yaml-cpp/CMakeLists.txt || $(SED) -i -e 's/LABEL_SUFFIX}"/LABEL_SUFFIX}"\n\tCOMPILE_FLAGS "-fPIC"/' yaml-cpp/CMakeLists.txt; true
	$(SED) -i -e 's/namespace\ \+YAML/namespace LHAPDF_YAML/' -e 's/YAML::/LHAPDF_YAML::/g' `find yaml-cpp -name "*.cpp" -o -name "*.h"`; true
	cd yaml-cpp && cmake . && $(MAKE)

all: yaml-cpp/libyaml-cpp.a
	@true
