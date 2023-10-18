CONFIG_FILE=/Users/leework/Documents/Research/projects/project_ace/data/raw/Misha_HT_BBN963_ELISpot_Results/ace_ready/ACE_Configuration_08062023.xlsx
SPOT_COUNT_FILE=/Users/leework/Documents/Research/projects/project_ace/data/raw/Misha_HT_BBN963_ELISpot_Results/ace_ready/ACE_Spot_Counts_08062023.csv
OUTPUT_EXCEL_FILE_EMPIRICAL=/Users/leework/Documents/Research/projects/project_ace/data/raw/Misha_HT_BBN963_ELISpot_Results/ace_outputs/ACE_Deconvolution_Empirical_08062023.xlsx
OUTPUT_EXCEL_FILE_EM=/Users/leework/Documents/Research/projects/project_ace/data/raw/Misha_HT_BBN963_ELISpot_Results/ace_outputs/ACE_Deconvolution_EM_08062023.xlsx
MIN_SPOT_COUNT=343

# Empirical
#ace deconvolve \
#	--readout-file-type pool_id \
#	--readout-files $SPOT_COUNT_FILE \
#	--assignment-excel-file $CONFIG_FILE \
#	--mode empirical \
#	--min-spot-count $MIN_SPOT_COUNT \
#	--output-excel-file $OUTPUT_EXCEL_FILE_EMPIRICAL

# Statistical
ace deconvolve \
	--readout-file-type pool_id \
	--readout-files $SPOT_COUNT_FILE \
	--assignment-excel-file $CONFIG_FILE \
	--mode em \
	--output-excel-file $OUTPUT_EXCEL_FILE_EM
