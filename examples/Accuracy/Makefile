INC=../../include/
SRC=../../src/

accuracy: accuracy.cpp $(SRC)*.cpp $(INC)*.h
	clang++ -g -I/Users/Martin/Desktop/lhapdfv6/deps/yaml-cpp/include/ -L/Users/Martin/Desktop/lhapdfv6/deps/yaml-cpp/lib/ -lyaml-cpp -I$(INC) $(SRC)*.cpp -I/usr/local/include/LHAPDF -L/usr/local/lib/ -lLHAPDF accuracy.cpp -o accuracy
