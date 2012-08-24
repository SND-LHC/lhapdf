CFLAGS=-I ../include -I /scratch/s0821167/build/include -fPIC

#Dummy
dummy: libLHAPDF.so

#Applications
accuracy: accuracy.cpp libLHAPDF.a

analytic: analytic.cpp libLHAPDF.a

#Library
libLHAPDF.a: Bicubic.o Bilinear.o Nearest.o PDFGrid.o PDFSet.o AlphaS.o
	
libLHAPDF.so: Bicubic.o Bilinear.o Nearest.o PDFGrid.o PDFSet.o AlphaS.o
	$(CC) -shared -fPIC $? -o $@

#Source
Bicubic.o: Bicubic.cpp
	$(CC) -c $(CFLAGS) $< -o $@ 

Bilinear.o: Bilinear.cpp
	$(CC) -c $(CFLAGS) $< -o $@ 

Nearest.o: Nearest.cpp
	$(CC) -c $(CFLAGS) $< -o $@ 

PDFGrid.o: PDFGrid.cpp
	$(CC) -c $(CFLAGS) $< -o $@ 

PDFSet.o: PDFSet.cpp
	$(CC) -c $(CFLAGS) $< -o $@

AlphaS.o: AlphaS.cpp
	$(CC) -c $(CFLAGS) $< -o $@
