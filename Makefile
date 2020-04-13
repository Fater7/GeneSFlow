all: main
main: main.o GeneSFlow.o GeneBody.o GSFRange.o
	g++ -g -std=c++11 ./bin/main.o ./bin/GeneSFlow.o ./bin/GeneBody.o ./bin/GSFRange.o -o ./main 

main.o: main.cpp GeneSFlow.cpp GeneBody.cpp GSFRange.cpp
	g++ -g -std=c++11 -c main.cpp -o ./bin/main.o

GeneSFlow.o: GeneSFlow.cpp GeneBody.cpp GSFRange.cpp
	g++ -g -std=c++11 -c GeneSFlow.cpp -o ./bin/GeneSFlow.o

GeneBody.o: GeneBody.cpp
	g++ -g -std=c++11 -c GeneBody.cpp -o ./bin/GeneBody.o

GSFRange.o: GSFRange.cpp
	g++ -g -std=c++11 -c GSFRange.cpp -o ./bin/GSFRange.o

clean:
	-rm -r ./bin/*.o main
