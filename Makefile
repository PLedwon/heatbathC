all:
	g++ -Ofast -o ./bin/heatbath ./src/heatbath.cpp -lm
	./bin/heatbath

compile:
	g++ -Ofast  -o ./bin/heatbath ./src/heatbath.cpp -lm

profile:
	g++ -Ofast ./src/heatbath.cpp -o ./bin/heatbath -lm -pg
	./bin/heatbath
	gprof ./bin/heatbath

plot:
	python3 ./plots/plots.py
