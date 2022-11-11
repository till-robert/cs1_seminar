all:
	mkdir -p build
	g++ -O3 task5_genloops.cpp -o build/task5_genloops
	./build/task5_genloops
	g++ -O3 task5.cpp build/task5_loops.cpp -o build/task5
	@echo
	./build/task5
clean:
	rm -rf build