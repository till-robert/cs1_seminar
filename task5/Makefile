all:
	g++ -O3 task5_genloops.cpp -o task5_genloops
	./task5_genloops
	g++ -O3 task5.cpp task5_loops.cpp -o task5
	@echo
	./task5
clean:
	rm -f task5_genloops task5 task5_loops.cpp task5_loops.h