cd ..
make clean all
cd bin
# ./drug_sim -input_deck input_deck_test.txt -hill_file control/IC50_drug_control.csv
./drug_sim -input_deck input_deck_test.txt -hill_file IC50_control_multi.csv
# ./drug_sim -input_deck input_deck_test.txt -hill_file IC50_quinidine.csv
# ./drug_sim -input_deck input_deck_test.txt -hill_file IC50_quinidine_16.csv