all:
	g++ ./src/heatbath.cpp -o ./bin/heatbath.o -lm
	./bin/heatbath.o

test:
	g++ -Ofast ./src/heatbath.cpp -o ./bin/heatbath.o -lm
	./bin/heatbath.o

profile:
	g++ -Ofast ./src/heatbath.cpp -o ./bin/heatbath.o -lm -pg
	./bin/heatbath.o
	gprof heatbath

plot:
	g++ -Ofast ./src/heatbath.cpp -o ./bin/heatbath.o -lm
	./bin/heatbath.o
