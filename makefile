#run:
#	source ~/srir/source_bash.sh && /opt/nfs/mpich-3.2/bin/mpicxx -ggdb -std=c++1y map.cpp dijkstra.cpp main.cpp -DDEBUG -o main && /opt/nfs/mpich-3.2/bin/mpiexec -n 2 ./main testcases/testcase0
all:
	/opt/nfs/mpich-3.2/bin/mpicxx -ggdb -std=c++1y map.cpp main.cpp -o main
	
test0:
	/opt/nfs/mpich-3.2/bin/mpiexec -n 3 ./main testcases/testcase0

test1:
	/opt/nfs/mpich-3.2/bin/mpiexec -n 3 ./main testcases/testcase1

test3:
	/opt/nfs/mpich-3.2/bin/mpiexec -n 4 ./main testcases/testcase3
