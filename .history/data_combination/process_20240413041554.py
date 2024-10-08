import os
import pandas as pd
from pandarallel import pandarallel

## 先去除is_edit为-1的行，只保留长度23的行和is_mis_synthesis为0的行（去除indel）

## 重新命名columns

## 筛选canonical的seq

## 计算总编辑效率 = 编辑/总 = 1 - (未编辑/总)

## 计算每一种outcome的编辑效率 = 某一种outcome的条数/相同的gRNA的“在3-10上有编辑 & canonical的条数”

# head -n 1000000 /dev/urandom > test_data.txt
# time python workflow.py test_data.txt



def process_data(file_path):
    pandarallel.initialize(nb_workers=20)
    data_list = []
    for data in pd.read_csv(file_path, header=0, sep='\t',chunksize=1e5):
        # data = pd.pivot_table(data, values='count', index=['gRNA', 'target'], columns=['is_edit'], aggfunc='sum').fillna(0)
        data = data[(data['is_edit'] != -1) & (data['tgt_len'] == 23) & (data['is_mis_synthesis'] == 0)]
        data_list.append(data)
    data = pd.concat(data_list,ignore_index=True)
    data.columns = ['gRNA', 'target','tgt_len','is_mis_synthesis','is_edit']
    data['gRNA'] = data['gRNA'].str[:-3]
    data['target'] = data['target'].str[:-3]
    
    # 筛选合法编辑        
        
    file_path = os.path.basename(file_path)

    if 'ABE' in file_path:
        file_df = data[data.parallel_apply(check_canonical_ABE, axis=1)]
    else:
        file_df = data[data.parallel_apply(check_canonical_CBE, axis=1)]
    
    #计算总编辑效率
        
    file_df['count'] = file_df.groupby(['gRNA', 'target'])['gRNA'].transform('count')
    file_df = file_df.drop_duplicates(subset=['gRNA', 'target'], keep='first') # 去除重复的行
    file_df['all_count'] = file_df.groupby(['gRNA'])['count'].transform('sum')
    file_df['unedited_count'] = file_df.groupby(['gRNA','target'])['count'].transform('sum')
    file_df['unedited_proportion'] = file_df['unedited_count'] / file_df['all_count']
    file_df['overall_efficiency'] = 1 - file_df['unedited_proportion'][file_df['is_edit'] == 0]
    
    # 赋值给每一种gRNA相同的总编辑效率
    
    # groups = pd.pivot_table(data, index=['gRNA'])
    groups = file_df.groupby(['gRNA'])
    for name, group in groups:
        equal_rows = group[group['gRNA'] == group['target']]
        if not equal_rows.empty:
            c_value = equal_rows.iloc[0]['overall_efficiency']
            file_df.loc[(file_df['gRNA'] == name) & (file_df['gRNA'] != file_df['target']), 'overall_efficiency'] = c_value
                
    
    # 计算每一种outcome的编辑效率
    
    cano_edited = file_df[file_df['gRNA'].str[2:10] != file_df['target'].str[2:10]].copy()
    cano_edited["deno"] = cano_edited.groupby('gRNA')['count'].transform('sum')
    cano_edited["proportion"] = cano_edited["count"] / cano_edited["deno"]
    unedited = (file_df['gRNA'] == file_df['target'])
    cano_unedited = file_df[unedited]
    processed_cano_seq = pd.concat([cano_edited, cano_unedited], axis=0)
    processed_cano_seq = processed_cano_seq.fillna(0)
    
    # 去除重复和重新排序
    
    processed_cano_seq = processed_cano_seq.sort_values(by=['gRNA', 'target'])
    processed_cano_seq = processed_cano_seq.reset_index(drop=True)
    processed_cano_seq = processed_cano_seq.drop(columns=['tgt_len','is_mis_synthesis','is_edit','unedited_count','unedited_proportion','deno'])
    
    processed_cano_seq.to_csv(file_path[0:-4] + '_processed.csv', header=True, encoding='utf-8', index=False)
    

def check_canonical_ABE(file_data):
    gRNA = file_data['gRNA']
    target = file_data['target']
    for i in range(len(gRNA)):
        if gRNA[i] != target[i]:
            if not (gRNA[i] == 'A' and target[i] == 'G'):
                return False
    return True

def check_canonical_CBE(file_data):
    gRNA = file_data['gRNA']
    target = file_data['target']
    for i in range(len(gRNA)):
        if gRNA[i] != target[i]:
            if not (gRNA[i] == 'C' and target[i] == 'T'):
                return False
    return True



