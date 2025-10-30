all: simulation_1.rds simulation_2.rds
	Rscript figure_4-5.R

simulation_1.rds:
	Rscript simulation_1.R

simulation_2.rds:
	Rscript simulation_2.R

clean:
	rm -f figure_4.pdf
	rm -f figure_5.pdf

clean-all: clean
	rm -f simulation_1.rds
	rm -f simulation_2.rds 
