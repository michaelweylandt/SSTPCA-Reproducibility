all: simulation_computational_convergence_results.rds
	Rscript figure_3.R

simulation_computational_convergence_results.rds:
	Rscript simulation_computational_convergence_figure_3.R

clean:
	rm -f figure_3.pdf

clean-all: clean
	rm simulation_computational_convergence_results.rds
