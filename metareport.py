<%
    import py_exp_calc.exp_calc as pxc
    import pandas as pd
    import re

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

    def parse_heatmap_from_flat(data,nrow,ncol,nvalue,smp_attr,var_attr,scale_factor=1):
            pairs = {}
            sample_attributes = {}
            var_attributes = {}
            for row in data:
                    if not pairs.get(row[nrow]):
                            pairs[row[nrow]] = {}
                    pairs[row[nrow]][row[ncol]] = row[nvalue]
                    if smp_attr: sample_attributes[row[nrow]] = [row[attr] for attr in smp_attr.keys()]
                    if var_attr: var_attributes[row[ncol]] = [row[attr] for attr in var_attr.keys()]
            mat, row, col = pxc.to_wmatrix_rectangular(pairs)
            mat = scale_factor*mat
            col_attrs = []
            col_attrs = [i for i in smp_attr.values()] if smp_attr else []
            table = [["-", *col_attrs, *col]]
            # Adding var attributes
            if var_attributes:
                for i, attr_name in enumerate(var_attr.values()):
                    table.append([attr_name]+["-"]*(len(table[0])-1))
                for cidx, cid in enumerate(col):
                    cidx = cidx + len(col_attrs) + 1
                    for attr_idx, attr in enumerate(var_attributes[cid]):
                        attr_idx = attr_idx + 1
                        table[attr_idx][cidx] = attr
            if sample_attributes:
                for idx,elem in enumerate(row): table.append([elem,*sample_attributes.get(elem), *mat[idx,:].tolist()])
            else:
                for idx,elem in enumerate(row): table.append([elem, *mat[idx,:].tolist()])
            return table

    def add_column_value(table,column,value):
        for row in table:
            row.insert(column, value)
        return table

    def add_column_values(table,colums, values):
        for i, value in enumerate(values):
            table = add_column_value(table,colums[i],value)
        return table

    def concat_heatmap(heatmaps):
        heatmaps_df = []
        for i, heatmap in enumerate(heatmaps):
            if i == 0:
                table = heatmap
            else:
                table = [row[1:] for row in heatmap]
            heatmaps_df.append(pd.DataFrame(table))
        heatmap_concat = pd.concat([*heatmaps_df], axis=1).reindex(heatmaps_df[0].index)
        return heatmap_concat.values.tolist()

    def parsed_string(data, blacklist = ["sim"]):
        words = []
        for word in data.split("_"):
                for blackword in blacklist:
                        word = re.sub(blackword,"",word)
                word = word.capitalize()
                words.append(word)
        parsed_data = " ".join(words)
        return parsed_data

    def parse_data(table, blacklist = ["sim"], column = "all"):
        parsed_table = []
        for i,row in enumerate(table):
                parsed_table.append(row)
                for j,data in enumerate(row):
                        if type(data) == str and not data.startswith("HGNC:"):
                                if parse_name.get(data):
                                        parsed_table[i][j] = parse_name[data]
                                else:
                                        parsed_table[i][j] = parsed_string(data, blacklist)
                        else:
                                continue
        return parsed_table
                
    def parse_table(name, blacklist=["sim"], include_header = False):
            if not include_header:
                    tab_header = plotter.hash_vars[name].pop(0)
                    plotter.hash_vars[name] = parse_data(plotter.hash_vars[name])
                    plotter.hash_vars[name].insert(0, tab_header)
            else:
                    plotter.hash_vars[name] = parse_data(plotter.hash_vars[name])
    
    def order_columns(name, column):
        tab_header = plotter.hash_vars[name].pop(0)
        plotter.hash_vars[name].sort(key=lambda x: x[column])
        plotter.hash_vars[name].insert(0, tab_header)

    for table in plotter.hash_vars.keys(): parse_table(table)
    order_columns("embeddings_nonInt_stats", 2)
    order_columns("zampieri_nonInt_coverage", 2)
    order_columns("buphamalai_nonInt_coverage", 2)
    order_columns("compensation_nonInt_coverage", 2)

    # edge density
    plotter.hash_vars["embeddings_nonInt_stats"] = add_column_values(plotter.hash_vars["embeddings_nonInt_stats"],
    [23,24,25],["Edge density","separator1","Genome"])
    edge_density_heatmap = parse_heatmap_from_flat(plotter.hash_vars['embeddings_nonInt_stats'][1:],
        1,2,6,{22:"Gene coverage"},{23:"Metric",24:"separator",25:"Dataset"},100)
    # zampieri
    zampieri_positives = len(plotter.hash_vars["zampieri_pos"])
    plotter.hash_vars["zampieri_nonInt_coverage"] = add_column_values(plotter.hash_vars["zampieri_nonInt_coverage"],
    [4,5,6],["GS-Coverage","separator2","Zampieri"])
    zampieri_cov_heatmap = parse_heatmap_from_flat(plotter.hash_vars["zampieri_nonInt_coverage"][1:],1,2,3,None,{4:"Metric",5:"separator",6:"Dataset"},100/zampieri_positives)
    # buphamalai
    buphamalai_positives = len(plotter.hash_vars["buphamalai_pos"])
    plotter.hash_vars["buphamalai_nonInt_coverage"] = add_column_values(plotter.hash_vars["buphamalai_nonInt_coverage"],
    [4,5,6],["GS-Coverage","separator2","Buphamalai"])
    buphamalai_cov_heatmap = parse_heatmap_from_flat(plotter.hash_vars["buphamalai_nonInt_coverage"][1:],
        1,2,3,None,{4:"Metric",5:"separator",6:"Dataset"}, 100/buphamalai_positives)
    # compensation
    compensation_positives = len(plotter.hash_vars["compensation_pos"])
    plotter.hash_vars["compensation_nonInt_coverage"] = add_column_values(plotter.hash_vars["compensation_nonInt_coverage"],
    [4,5,6],["GS-Coverage","separator2","Compensation"])
    compensation_cov_heatmap = parse_heatmap_from_flat(plotter.hash_vars["compensation_nonInt_coverage"][1:],
        1,2,3,None,{4:"Metric",5:"separator",6:"Dataset"}, 100/compensation_positives)
    print("+-+-"*30)
    print(edge_density_heatmap)
    print("+-+-"*30)
    print(zampieri_cov_heatmap)
    print("+-+-"*30)
    print(buphamalai_cov_heatmap)
    print("+-+-"*30)
    print(compensation_cov_heatmap)
    plotter.hash_vars["heatmap_concat"] = concat_heatmap([edge_density_heatmap,zampieri_cov_heatmap,buphamalai_cov_heatmap,compensation_cov_heatmap])
    plotter.hash_vars["heatmap_concat"].insert(0,plotter.hash_vars["heatmap_concat"][0])
    plotter.hash_vars["heatmap_concat"][1][0] = "eGSM type"
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
                                        "showSmpDendrogram":False})} 
