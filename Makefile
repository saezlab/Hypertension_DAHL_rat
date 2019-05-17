analysis:
	@echo "Data pre-processing"
	@echo "==================="
	@echo
	Rscript 1_preprocessing.R
	@echo
	@echo "Computing differential expression"
	@echo "================================="
	@echo
	Rscript 2_diff_exp.R
	@echo
	@echo "Running PHONEMeS"
	@echo "================"
	@echo
	Rscript 3_phonemes.R
	@echo
	@echo "*** Resulting network saved in results/ ***"
	@echo
	@echo "+===================+"
	@echo "|     FINISHED!     |"
	@echo "+===================+"
