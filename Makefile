all:
	cd main_paper && make
	cd supplements && make

clean:
	cd main_paper && make clean
	cd supplements && make clean
