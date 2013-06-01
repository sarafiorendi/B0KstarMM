CC        = g++
CCFLAGS   = -g -Wall
SOURCES   = 
ROOTFLAGS = `root-config --cflags --libs` -lMinuit

all:
	$(CC) $(CCFLAGS) -o MyPixelAnalysisRatios $(ROOTFLAGS) $(SOURCES) MyPixelAnalysisRatios.cc

clean:
	rm -rf MyPixelAnalysisRatios
