all: multi_rank_error_simul.rds
	cd multifactor_sstpca_sec_f31/ && make
	cd deflation_sec_f32/ && make
	cd tuning_parameters_sec_f33/ && make

multi_rank_error_simul.rds:
	Rscript simulate_f3_data.R

clean:
	cd multifactor_sstpca_sec_f31/ && make clean
	cd deflation_sec_f32/ && make clean
	cd tuning_parameters_sec_f33/ && make clean

clean-all: clean
	rm -f multi_rank_error_simul.rds
