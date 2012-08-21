INC=../../include/
SRC=../../src/

settest: settest.cpp $(SRC)*.cpp $(INC)*.h
	clang++ -I$(INC) $(SRC)*.cpp -I/Users/Martin/Desktop/lhapdfv6/deps/yaml-cpp/include/ -L/Users/Martin/Desktop/lhapdfv6/deps/yaml-cpp/lib -lyaml-cpp settest.cpp -o settest

settestall: settest.cpp $(SRC)*.cpp $(INC)*.h
	clang++ -Wall -I$(INC) $(SRC)*.cpp -I/Users/Martin/Desktop/lhapdfv6/deps/yaml-cpp/include/ -L/Users/Martin/Desktop/lhapdfv6/deps/yaml-cpp/lib -lyaml-cpp settest.cpp -o settest
