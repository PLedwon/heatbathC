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

clean:
	rm ../csvData/*.csv
	rm ./data/log/*.txt

slurm:
	g++ -Ofast -o ./bin/heatbath ./src/heatbath.cpp -lm
	./slurmSubmit.sh
