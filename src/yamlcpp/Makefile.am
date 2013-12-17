EXTRA_DIST = yaml-cpp-0.3.0.tar.gz

yaml-cpp/libyaml-cpp.a: yaml-cpp-0.3.0.tar.gz
	test -d yaml-cpp || tar xf $^; true
	grep COMPILE_FLAGS yaml-cpp/CMakeLists.txt || sed -i -e 's/LABEL_SUFFIX}"/LABEL_SUFFIX}"\n\tCOMPILE_FLAGS "-fPIC"/' yaml-cpp/CMakeLists.txt; true
	cd yaml-cpp && cmake . && $(MAKE)

all: yaml-cpp/libyaml-cpp.a
	@true
