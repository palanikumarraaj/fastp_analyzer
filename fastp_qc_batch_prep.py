# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 17:14:28 2026

@author: User
"""

import json
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import os

def extract_fastp_metrics(json_path):
    with open(json_path) as f:
        data = json.load(f)
    
    row = {
        'sample': json_path.stem,   # or extract from filename as you prefer
        'fastp_version': data.get('summary', {}).get('fastp_version'),
        'sequencing': data.get('summary', {}).get('sequencing'),
        
        # Before filtering
        'before_total_reads': data['summary']['before_filtering']['total_reads'],
        'before_total_bases': data['summary']['before_filtering']['total_bases'],
        'before_q20_rate': data['summary']['before_filtering']['q20_rate'],
        'before_q30_rate': data['summary']['before_filtering']['q30_rate'],
        'before_gc_content': data['summary']['before_filtering']['gc_content'],
        'before_read1_mean_length': data['summary']['before_filtering'].get('read1_mean_length'),
        'before_read2_mean_length': data['summary']['before_filtering'].get('read2_mean_length'),
        
        # After filtering
        'after_total_reads': data['summary']['after_filtering']['total_reads'],
        'after_total_bases': data['summary']['after_filtering']['total_bases'],
        'after_q20_rate': data['summary']['after_filtering']['q20_rate'],
        'after_q30_rate': data['summary']['after_filtering']['q30_rate'],
        'after_gc_content': data['summary']['after_filtering']['gc_content'],
        'after_read1_mean_length': data['summary']['after_filtering'].get('read1_mean_length'),
        'after_read2_mean_length': data['summary']['after_filtering'].get('read2_mean_length'),
        
        # Filtering rates
        'percentage_reads_survived': (data['summary']['after_filtering']['total_reads'] / 
                              data['summary']['before_filtering']['total_reads']) * 100,
        'percentage_bases_survived': (data['summary']['after_filtering']['total_bases'] / 
                              data['summary']['before_filtering']['total_bases']) * 100,
        
        # Other useful sections
        'adapter_trimmed_reads': data.get('adapter_cutting', {}).get('adapter_trimmed_reads', 0),
        'duplication_rate': (data.get('duplication', {}).get('rate', None)) * 100,

        # Insert size
        'insert_size' : data.get('insert_size', {}).get('peak',None),

        # Read1_before_filtering
        'read1_before_filtering_total_reads' : data['read1_before_filtering']['total_reads'],

        # Read2_before_filtering
        'read2_before_filtering_total_reads' : data['read2_before_filtering']['total_reads'],


        # Read1_after_filtering
        'read1_after_filtering_total_reads' : data['read1_after_filtering']['total_reads'],


        # Read2_after_filtering
        'read2_after_filtering_total_reads' : data['read2_after_filtering']['total_reads'],
        
        # Filtering_results
        'passed_filter_reads_overall' : data['filtering_result']['passed_filter_reads'],
        'low_quality_reads' : data['filtering_result']['low_quality_reads'],
        'too_many_N_reads' : data['filtering_result']['too_many_N_reads'],
        'too_short_reads' : data['filtering_result']['too_short_reads'],
        'too_long_reads' : data['filtering_result']['too_long_reads'],
             
    }
    
    return row


# Main processing
folder = input("Please enter the batch location (e.g., /mnt/ngs1): ")
os.chdir(folder)
folder = Path(folder)
json_files = list(folder.glob("*.json"))

records = []
for f in tqdm(json_files):
    try:
        records.append(extract_fastp_metrics(f))
    except Exception as e:
        print(f"Error processing {f}: {e}")
        
        
df = pd.DataFrame(records)

## new varaibles

df['read1_total_reads_diff'] = df['read1_before_filtering_total_reads'] - df['read1_after_filtering_total_reads']
df['read2_total_reads_diff'] = df['read2_before_filtering_total_reads'] - df['read2_after_filtering_total_reads']
df['read1and2_afterQC'] = df.apply(lambda row: 'No_reads_issue' if row['read2_total_reads_diff'] == row['read1_total_reads_diff'] else 'Reads_anomaly', axis=1)
df['adapter_timmed_reads_percentage'] = (df['adapter_trimmed_reads']) / df['before_total_reads'] * 100
df['proper_read_length_before'] =  df.apply(lambda row: 'Proper_paried_length' if row['before_read1_mean_length'] == row['before_read2_mean_length'] else 'Reads_anamaly', axis=1)
df['read1_length_loss'] = (df['before_read1_mean_length'] - df['after_read1_mean_length']) / df['before_read1_mean_length'] * 100
df['read2_length_loss'] = (df['before_read2_mean_length'] - df['after_read2_mean_length']) / df['before_read2_mean_length'] * 100
df['total_q30_before_after_comp'] = (df['after_q30_rate'] - df['before_q30_rate']) / df['after_q30_rate'] * 100
df['total_q20_before_after_comp'] = (df['after_q20_rate'] - df['before_q20_rate']) / df['after_q20_rate'] * 100
df['total_gc_before_after_comp'] = (df['before_gc_content'] - df['after_gc_content']) / df['before_gc_content'] * 100
df['total_reads_before_after_comp'] = (df['before_total_reads'] - df['after_total_reads']) / df['before_total_reads'] * 100
df['total_bases_before_after_comp'] = (df['before_total_bases'] - df['after_total_bases']) / df['before_total_bases'] * 100
df['read1_total_reads_retention'] = (df['read1_before_filtering_total_reads'] - df['read1_total_reads_diff']) / df['read1_before_filtering_total_reads'] * 100
df['read2_total_reads_retention'] = (df['read2_before_filtering_total_reads'] - df['read2_total_reads_diff']) / df['read2_before_filtering_total_reads'] * 100


# Select Specific Columns from DataFrame
f1f = df[['sample', 'read1_total_reads_retention', 'read2_total_reads_retention', 'total_q30_before_after_comp', 'total_q20_before_after_comp', 'total_gc_before_after_comp', 'total_reads_before_after_comp', 'total_bases_before_after_comp', 'duplication_rate', 'low_quality_reads', 'too_many_N_reads', 'too_short_reads', 'read1_length_loss', 'read2_length_loss', 'percentage_reads_survived',  'percentage_bases_survived' , 'insert_size', 'adapter_timmed_reads_percentage']]


## To transform based on model_sample.csv format
f1f_transformed = pd.DataFrame(f1f).set_index('sample').T

f1f_transformed.to_csv("samples_batch_values.csv", index=True)















