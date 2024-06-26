<%
    import py_exp_calc.exp_calc as pxc
    import pandas as pd
    import re
    import os
    import sys 
    sys.path.append("./aux_report")
    import pyreport_helper as prh

    parse_name = {
                "string_ppi_combined": "STRING combined", 
                "string_ppi_textmining":"STRING textmining",
                "string_ppi_coexpression": "STRING coexpression",
                "string_ppi_neighborhood": "STRING neighborhood",
                "string_ppi_experimental": "STRING experiments",
                "string_ppi_cooccurence": "STRING cooccurrence",
                "phenotype": "HPO",
                "disease":"Disease",
                "pathway": "Pathway",
                "DepMap_effect_pearson": "DepMap Pearson",
                "string_ppi_database":"STRING databases",
                "DepMap_effect_spearman":"DepMap Spearman",
                "hippie_ppi": "Hippie",
                "DepMap_Kim":"DepMap Kim",
                "string_ppi_fusion":"STRING fusion",
                "gene_hgncGroup": "HGNC group",
                "integration_mean_by_presence": "IMP",
                "mean": "Mean",
                "max": "Max",
                "median": "Median",
                "el": "EL",
                "raw_sim": "GSM",
                "node2vec": "n2v",
                "rf": "RF"
                }
            
    if plotter.hash_vars["selection"][0][0] == "original_databases":
        whitelist = ["string_ppi_combined", "phenotype", "disease", "pathway",
        "DepMap_effect_pearson", "DepMap_effect_spearman","hippie_ppi","DepMap_Kim",
        "gene_hgncGroup", "biological_process","molecular_function","cellular_component"]
    elif plotter.hash_vars["selection"][0][0] == "string_channels":
        whitelist = ["string_ppi_textmining", "string_ppi_combined", "string_ppi_coexpression",
        "string_ppi_neighborhood","string_ppi_experimental","string_ppi_cooccurence","string_ppi_database","string_ppi_fusion"]

    # Selecting databases of interest
    for table in plotter.hash_vars.keys(): 
        if not re.search("_pos$",table):
            print(whitelist)
            print("\n----\n"*5)
            print(len(plotter.hash_vars[table]))
            plotter.hash_vars[table]=prh.select_columns(plotter.hash_vars[table], 1, whitelist)
            print(len(plotter.hash_vars[table]))

    # Parsing tables names and order
    for table in plotter.hash_vars.keys(): 
        if table != "selection": prh.parse_table(plotter.hash_vars[table],parse_name=parse_name)
    prh.order_columns(plotter.hash_vars["embeddings_nonInt_stats"], 2, {"GSM": 1, "n2v": 4, "EL":2, "RF": 3})
    prh.order_columns(plotter.hash_vars["zampieri_nonInt_coverage"], 2, {"GSM": 1, "n2v": 4, "EL":2, "RF": 3})
    prh.order_columns(plotter.hash_vars["buphamalai_nonInt_coverage"], 2, {"GSM": 1, "n2v": 4, "EL":2, "RF": 3})
    prh.order_columns(plotter.hash_vars["compensation_nonInt_coverage"], 2, {"GSM": 1, "n2v": 4, "EL":2, "RF": 3})

    # edge density
    plotter.hash_vars["embeddings_nonInt_stats"] = prh.add_column_values(plotter.hash_vars["embeddings_nonInt_stats"],
    [23,24,25],["Edge density","separator1","Genome"])
    edge_density_heatmap = prh.parse_heatmap_from_flat(plotter.hash_vars['embeddings_nonInt_stats'][1:],
        1,2,6,{22:"Gene coverage"},{23:"Metric",24:"separator",25:"Dataset"},100)
    # zampieri
    zampieri_positives = len(plotter.hash_vars["zampieri_pos"])
    plotter.hash_vars["zampieri_nonInt_coverage"] = prh.add_column_values(plotter.hash_vars["zampieri_nonInt_coverage"],
    [4,5,6],["GS-Coverage","separator2","Zampieri"])
    zampieri_cov_heatmap = prh.parse_heatmap_from_flat(plotter.hash_vars["zampieri_nonInt_coverage"][1:],1,2,3,None,{4:"Metric",5:"separator",6:"Dataset"},100/zampieri_positives)
    # buphamalai
    buphamalai_positives = len(plotter.hash_vars["buphamalai_pos"])
    plotter.hash_vars["buphamalai_nonInt_coverage"] = prh.add_column_values(plotter.hash_vars["buphamalai_nonInt_coverage"],
    [4,5,6],["GS-Coverage","separator2","Buphamalai"])
    buphamalai_cov_heatmap = prh.parse_heatmap_from_flat(plotter.hash_vars["buphamalai_nonInt_coverage"][1:],
        1,2,3,None,{4:"Metric",5:"separator",6:"Dataset"}, 100/buphamalai_positives)
    # compensation
    compensation_positives = len(plotter.hash_vars["compensation_pos"])
    plotter.hash_vars["compensation_nonInt_coverage"] = prh.add_column_values(plotter.hash_vars["compensation_nonInt_coverage"],
    [4,5,6],["GS-Coverage","separator2","Compensation"])
    # concatenate all heatmaps
    compensation_cov_heatmap = prh.parse_heatmap_from_flat(plotter.hash_vars["compensation_nonInt_coverage"][1:],
        1,2,3,None,{4:"Metric",5:"separator",6:"Dataset"}, 100/compensation_positives)
    plotter.hash_vars["heatmap_concat"] = prh.concat_heatmap([edge_density_heatmap,zampieri_cov_heatmap,buphamalai_cov_heatmap,compensation_cov_heatmap])
    plotter.hash_vars["heatmap_concat"].insert(0,plotter.hash_vars["heatmap_concat"][0][:])
    plotter.hash_vars["heatmap_concat"][1][0] = "eGSM type"
    for idx, value in enumerate(plotter.hash_vars["heatmap_concat"][0]):
        if idx > 1:
            plotter.hash_vars["heatmap_concat"][0][idx] = plotter.hash_vars["heatmap_concat"][0][idx] + f"_{idx}"
%>

${ plotter.heatmap(id = 'heatmap_concat', title="",header = True, row_names = True, smp_attr=[1],var_attr=[1,2,3,4],
                                        config= {"varTextRotate":45,
                                        "splitVariablesBy":"separator",
                                        "smpOverlays":["Gene coverage"],
                                        "varOverlays":["Metric","eGSM type","Dataset"],
                                        "smpOverlayProperties": {"Gene coverage":{"type":"Bar","thickness":200,"color":"#A1A1A1"}},
                                        "varOverlayProperties":{
                                            "Metric":{"position":"top","scheme":"Matlab","showLegend":"true","type":"Default"},
                                            "eGSM type":{"position":"bottom","scheme":"Matlab","showLegend":"true","type":"Default"},
                                            "Dataset":{"position":"bottom","scheme":"Matlab","showLegend":"true","type":"Default"}
                                            },
                                        "setMinX":0,
                                        "setMaxX":100, 
                                        "xAxisTitle": "Density", 
                                        "smpTextScaleFontFactor":0.65,
                                        "showLegendTitle":False,
                                        "showVariableNames":False,
                                        "samplesClustered":True,
                                        "showSmpDendrogram":True})} 
