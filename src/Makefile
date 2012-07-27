INC=../include

#Dummy
dummy: libLHAPDF.so

#Applications
accuracy: accuracy.cpp libLHAPDF.a

analytic: analytic.cpp libLHAPDF.a

#Library
libLHAPDF.a: Bicubic.o Bilinear.o Nearest.o PDFGrid.o PDFSet.o
	
libLHAPDF.so: Bicubic.o Bilinear.o Nearest.o PDFGrid.o PDFSet.o

#Source
Bicubic.o: Bicubic.cpp
	$(CC) -c -I $(INC) $? -o $@ 

Bilinear.o: Bilinear.cpp
	$(CC) -c -I $(INC) $? -o $@ 

Nearest.o: Nearest.cpp
	$(CC) -c -I $(INC) $? -o $@ 

PDFGrid.o: PDFGrid.cpp
	$(CC) -c -I $(INC) $? -o $@ 

PDFSet.o: PDFSet.cpp
	$(CC) -c -I $(INC) $? -o $@