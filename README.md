# fastp_analyzer
fastp batch analyzer to flag issues on QC metrics for Human WES samples

**The goal of this project is to create useful varaibles, utilize certain thresholds derived from a larger cohort, recalibrate the thresholds and flag the QC varaibles respectively related to biology or libray prep based.**

NOTE : The thresholds and distribution file is created based on a large cohort of 730 fastp json files and its QC parameters. So, it can be directly used or can also be created and calibrated.

## Case 1 - Direct Usage

If the given source files for distribution levels and threshold values are relevant to Human WES fastp results and range, then proceed directly.

**Step 1**

```
$ python3 fastp_qc_batch_prep.py
```
The script *python3 fastp_qc_batch_prep.py* can be used directly and can request only the location of your batch samples fastp json files.

The script will look for **json** file extensions and load the json file for collecting required parameters and analysis.

