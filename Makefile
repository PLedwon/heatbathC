all:
	g++ -Ofast  -o ./bin/heatbath ./src/heatbath.cpp -lm
	./bin/heatbath

profile:
	g++ -Ofast ./src/heatbath.cpp -o ./bin/heatbath -lm -pg
	./bin/heatbath
	gprof heatbath

plot:
	g++ -Ofast ./src/heatbath.cpp -o ./bin/heatbath -lm
	./bin/heatbath
