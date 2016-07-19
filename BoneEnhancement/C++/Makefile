#chb2ab, Crispin Bernier, lab06

CXX=clang++
CXXFLAGS= -O2
OFILES=main.cpp
Exe=./a.out
Args=1original.pgm 2blurred.pgm 3differential.pgm 4ReflectionImage.pgm 5ShadowImage.pgm 6ProbabilityImage.pgm 7Thresholded.pgm >out.txt
.SUFFIXES: .o .cpp

make: $(OFILES)
	$(CXX) $(CXXFLAGS) $(OFILES)

main.o: main.cpp image.h image.cpp

run: $(Exe)
	$(Exe) $(Args)
