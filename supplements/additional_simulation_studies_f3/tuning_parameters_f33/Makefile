all: select_rank_BIC_error.rds select_K_error.rds
	Rscript figure_A9.R
	Rscript figure_A10.R

clean:
	rm -f figure_A9.pdf
	rm -f figure_A10.pdf

select_rank_BIC_error.rds:
	Rscript simulation_select_rank_BIC.R

select_K_error.rds:
	Rscript simulation_select_K.R

clean-all: clean
	rm -f select_rank_BIC_error.rds
	rm -f select_K_error.rds
