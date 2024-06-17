#!/usr/bin/env bash


if [ $1 == "ld" ] ; then
	mkdir -p data_results
	# loading Gold standards:
	#------------------------#
	cp $SCRATCH/executions/GraphPrioritizer/report/zampieri/ranking_report/control_pos ./data_results/zampieri_pos
	cp $SCRATCH/executions/GraphPrioritizer/report/menche/ranking_report/control_pos ./data_results/buphamalai_pos
	cp $SCRATCH/executions/BackupGenes/report/ranking_report/control_pos ./data_results/compensation_pos

	# loading results
	#-----------------#
	# # load data from GraphPrioritizer embeddings data
	# cp $SCRATCH/executions/GraphPrioritizer/report/menche/kernel_report/parsed_uncomb_kernel_metrics ./data_results/embeddings_nonInt_stats
	# cp $SCRATCH/executions/GraphPrioritizer/report/menche/kernel_report/parsed_comb_kernel_metrics ./data_results/embeddings_Int_stats
	# # load buphamalai data
	# cp $SCRATCH/executions/GraphPrioritizer/report/menche/ranking_report/parsed_non_integrated_rank_pos_cov ./data_results/buphamalai_nonInt_coverage
	# cp $SCRATCH/executions/GraphPrioritizer/report/menche/ranking_report/parsed_integrated_rank_pos_cov ./data_results/buphamalai_Int_coverage
	# # load zampieri data
	# cp $SCRATCH/executions/GraphPrioritizer/report/zampieri/ranking_report/parsed_non_integrated_rank_pos_cov ./data_results/zampieri_nonInt_coverage
	# cp $SCRATCH/executions/GraphPrioritizer/report/zampieri/ranking_report/parsed_integrated_rank_pos_cov ./data_results/zampieri_Int_coverage
	# # load from compensatory data
	# cp $SCRATCH/executions/BackupGenes/report/ranking_report/parsed_non_integrated_rank_pos_cov ./data_results/compensation_nonInt_coverage
	# cp $SCRATCH/executions/BackupGenes/report/ranking_report/parsed_integrated_rank_pos_cov ./data_results/compensation_Int_coverage
fi

if [ $1 == "rp" ] ; then
	source ~soft_bio_267/initializes/init_python
	name_dir=`date +%d_%m_%Y`
    mkdir -p ./meta_reports/
	report_html -t ./metareport.py  -d `ls data_results/* | tr -s [:space:] "," | sed 's/,*$//g'` -o "$name_dir"
fi