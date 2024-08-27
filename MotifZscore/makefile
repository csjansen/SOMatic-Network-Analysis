CC = g++
FILES = main.cpp
FILES2 = main2.cpp
OUT_EXE = MotifZScore
OUT_EXE2 = MotifZScore2

build: $(FILES)
	$(CC) $(FILES) -O3 -fpermissive -o $(OUT_EXE) -std=c++11 -pthread
	$(CC) $(FILES2) -O3 -fpermissive -o $(OUT_EXE2) -std=c++11 -pthread

clean:
	rm -f MotifZScore
	rm -f MotifZScore2

