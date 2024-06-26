import re
import py_exp_calc.exp_calc as pxc

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

def parse_data(table, blacklist = ["sim"], column = "all", parse_name={}):
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
            
def parse_table(original_table, blacklist=["sim"], include_header = False, parse_name={}):
        original_table 
        if not include_header:
                tab_header = original_table.pop(0)
                original_table = parse_data(table=original_table, parse_name=parse_name)
                original_table.insert(0, tab_header)
        else:
                original_table = parse_data(table=original_table, parse_name=parse_name)

def order_columns(original_table, column, custom_orden):
    tab_header = original_table.pop(0)
    if custom_orden:
        original_table.sort(key=lambda x: custom_orden[x[column]])
    else:
        original_table.sort(key=lambda x: x[column])
    original_table.insert(0, tab_header)

def select_columns(original_table, column, whitelist):
    tab_header = original_table.pop(0)
    filtered_table = []
    for row in original_table:
        if row[column] in whitelist: filtered_table.append(row)
    filtered_table.insert(0,tab_header)
    original_table = filtered_table