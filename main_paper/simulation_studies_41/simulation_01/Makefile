all: simulation_methods_comparison_results.rds
	Rscript figure_2.R

simulation_methods_comparison_results.rds:
	Rscript simulation_methods_figure_2.R

clean:
	rm -f figure_2.pdf

clean-all: clean
	rm simulation_methods_comparison_results.rds
