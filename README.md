# fastp_analyzer
use fastp analyzer to flag issues on QC metrics for Human WES samples

**The goal of this project is to create useful variables, utilize certain thresholds derived from a larger cohort, recalibrate the thresholds and flag the QC varaibles respectively related to biology or libray preparation based.**

NOTE : The thresholds and distribution file is created based on a large cohort of 730 fastp json files and its QC parameters. So, it can be directly used or can also be created and calibrated.

The entire script, method and dataset is prepared, tested on a single fastp version format.

Hence, the usage is limited or focused on fastq file QC files. The concept of WES capture kits or assembly version is not useful.

Kindly note that Human WES is targeted in the range of 6 to 8 GB and hence the distibution is best suited for similar size. For targeted small panels or highly deviated levels like less than 4 or higher than 10 gb, please check or modify threshold file accordingly or create your own datasets and prepare threshold levels accordingly.

## Case 1 - Direct Usage

If the given source files for distribution levels and threshold values are relevant to Human WES fastp results and range, then proceed directly.

Kindly explore the distribution related provided file to check with your relevant json files or batch files if required.


**Step 1**

```
$ python3 fastp_qc_batch_prep.py
```
The script *python3 fastp_qc_batch_prep.py* can be used directly and can request only the location of your batch samples fastp json files.

The script will look for **json** file extensions and load the json file for collecting required parameters and analysis.

Kindly use the location where only fastp output json files present and avoid any other json files from different tools.



