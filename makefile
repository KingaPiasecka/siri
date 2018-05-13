compile:
	/opt/nfs/mpich-3.2/bin/mpicxx -ggdb -std=c++1y Graph.cpp main.cpp -o main
	
expSmall:
	/opt/nfs/mpich-3.2/bin/mpiexec -n 4 ./main cases/exampleSmall

expBig:
	/opt/nfs/mpich-3.2/bin/mpiexec -n 4 ./main cases/exampleBig

clean:
	rm -f main