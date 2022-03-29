# ctg-sc-multi-10x
10x GEX and cmo and/or adt


- Run one project at a time only
- Demux -> cellranger -> QC -> delivery
- Supports ADT/HTO, CMO and RNA libraries sequenced on same flowcell.
- Supports different indexing of RNA and ADT/HTO library (e.g. RNA dual 10bp and ADT/HTO single 6bp). See `Handle dual and single indexing in same sequencing run` for more info.
- TotalSeq feature IDs used for ADT/HTO in experiment can be specified in samplesheet. Pipeline creates feature reference csv from IDs. See `2. Feature reference (feature.ref.csv)` and `1. Samplesheet` sections below.
- CMO ids are specified in samplesheet, pr sample
- If antibody/cmo features are not in reference, a custom feature-ref csv can be created and specified with drivers -f tag.
- For ADT/HTO: Currently, BioLegend TotalSeq ADT A,B and C as well as TotalSeq HTO A are in references, located at `/projects/fs1/shared/references/cellranger_totalseq`.
- For CMO: CMO301-312 (default cmo reference) are supported.


1. Clone and build the Singularity container for this pipeline: https://github.com/perllb/ctg-sc-multi-10x/tree/master/container
2. If not using reference totalSeq/cmo tags/antibodies, prepare the feature ref csv. See section `Feature reference` below
3. Edit your samplesheet to match the example samplesheet. See section `SampleSheet` below
4. Edit the nextflow.config file to fit your project and system. 
5. Run pipeline 
```
nohup nextflow run pipe-sc-multi-10x.nf > log.pipe-sc-multi-10x.txt &
```
## Driver
- Driver for pipeline located in `.../shared/ctg-pipelines/ctg-sc-multi-10x/sc-multi-10x-driver`

## Input files

The following files must be in the runfolder to start pipeline successfully.

1. Samplesheet  (see `SampleSheet` section below)
2. (OPTIONAL): Feature reference csv (see `Feature Reference` section below) if not using reference features on lsens4 `/projects/fs1/shared/references/cellranger_totalseq`.

### 1. Samplesheet (CTG_SampleSheet.sc-multi-10x.csv):

- If using other names than CTG_SampleSheet.**sc-multi-10x**.csv - then specify which sheet when starting driver: e.g. `sc-multi-10x-driver -s CTG_SampleSheet.2022_102.csv`

#### Example sheet

##### With both ADT and CMO
```
[Header]
metaid,2021_067_citeseqTestmulti
antibodies,"ADT_A0574,ADT_A0052,ADT_A0394,ADT_A0161,ADT_A0063,ADT_A0576,ADT_A0054,ADT_A0048"
cmo,y
email,per.a@med.lu.se
autodeliver,y
[Data]
Lane,Sample_ID,index,Sample_Project,Sample_Species,Sample_Lib,Sample_Pair,Sample,CMO,cmotype
,pool_1,SI-TT-A6,2022_228,mouse,gex,1,Pool1,CMO301|CMO302|CMO303|CMO308|CMO311,single
,pool_2,SI-TT-B6,2022_228,mouse,gex,2,Pool2,CMO302|CMO304|CMO307|CMO310|CMO312,single
,pool_1_CP,SI-NN-A2,2022_228,mouse,cp,1,Pool1,,,
,pool_2_CP,SI-NN-B2,2022_228,mouse,cp,2,Pool2,,,
,pool_1_ADT,ACAGTG,mouse,2022_228,adt,1,Pool1,,,
,pool_2_ADT,TGACCA,mouse,2022_228,adt,2,Pool2,,,
```
#### [Header] section

- `email` : Email to customer (or ctg staff) that should retrieve email with qc and deliver info upon completion of pipeline. Note: only lu emails works (e.g. @med.lu.se or @lth.se.
- `autodeliver` : set to `y` if email should be sent automatically upon completion. Otherwise, set to `n`.
- `metaid` : optional. set to create working folders with this name. otherwise, it will be created based on runfolder name/date.
- `cmo` : Set to `y` if there are CMO libraries in the project, otherwise `n`. CMO IDs are listed in [Data] Section, pr sample (see [Data] section below)
- `antibodies` : set to 'n' if ADT/HTO are NOT in experiment. Otherwise, list of ADT/HTO IDs (must match reference) used in experiment (across all samples). The sc-multi-10x-driver will use these IDs to extract all info from cellranger_totalseq references, and create the "features" csv file needed in count analysis.
  - recommended if using antibodies that are defined in the totalSeq human cocktail csv files (/projects/fs1/shared/references/cellranger_totalseq). 
  	- Set IDs of antibodies that were used in experiment (and will be used to create the feature reference file for `count` analysis.)
  	- IMPORTANT: They MUST match the IDs of the totalSeq human cocktail references on lsens 
	- The list of antibodies should be comma-separated, and best if quoted (should also work without quote, but not tested).
  - Alternative is to
  	- create a feature.ref.csv file, and add it to runfolder. It must have the standard required cellranger format (see: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis). Note: The sc-multi-10x-driver will first look for a `feature.ref.csv` file in runfolder. If not found, it will look for `antibodies,` field in header section. 
  	- The last alternative is to run the driver with -f flag, pointing to any file that serve as feature-ref csv. e.g. `sc-multi-10x-driver -f /path/to/feature.ref.csv`

#### [Data] section

 | Lane | Sample_ID | index | Sample_Species | Sample_Project | Sample_Lib | Sample_Pair | Sample | CMO | cmotype |
 | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
 | | Sr1 | SI-GA-D9 | human | 2022_022 | gex | 1 | S1 | CMO301\|CMO302 | multi |
 | | Sr1_CP | SI-GA-H9 | human | 2022_022 | cp | 1 | S1 | | |
 | | Sr1_ADT | ATGCAG | human | 2022_022 | adt | 1 | S1 | |
 | | Sr2 | SI-GA-C9 | human | 2022_022 | gex | 2 | S2 | CMO303\|CMO304 | single |
 | | Sr2_CP | SI-GA-C9 | human | 2022_022 | cp | 2 | S2 | |
 | | Sr2_ADT | ATGTTG | human | 2022_022 | adt | 2 | S2 | |


- The nf-pipeline takes the following Columns from samplesheet to use in channels:

- `Sample_ID` : ID of sample. Sample_ID can only contain a-z, A-Z and "_".  E.g space and hyphen ("-") are not allowed! If 'Sample_Name' is present, it will be ignored. IMPORTANT: CellPlex (CP) Sample_ID MUST have corresponding name to GEX Sample_ID, with `_CP` suffix! So if GEX Sample_ID is "Sample_A", CP sample MUST be "Sample_A_CP". If GEX Sample_ID is "Cascvf", CP Sample_ID MUST be "Cascvf_CP".
- `index` : Must use index ID (10x ID) if dual index. For single index, the index sequence works too.
- `Sample_Project` : Project ID. E.g. 2021_033, 2021_192.
- `Sample_Species` : Only 'human'/'mouse'/'custom' are accepted. If species is not human or mouse, set 'custom'. This custom reference genome has to be specified in the nextflow config file. See below how to edit the config file.
- `Sample_Lib` : 'gex'/'cp'/'adt'. Specify whether sample is RNA (gex) or CellPlex (cp) or ADT/HTO (adt) library. 
- `Sample_Pair` : To match the GEX sample with the corresponding CP and ADT/HTO sample. e.g. in the example above, sample 'Sr1' is the GEX library, that should be matched with 'Sr1_CP' which is the CP library of the sample and 'Sr1_ADT' the corresponding ADT library.
- `Sample` : A common name for GEX+CellPle+ADT sample. Must be the same for the corresponing GEX, CP and ADT samples (such as in example above). So all samples with "Sample_Pair" = 1, must have same "Sample" etc.
- `CMO` : CMO IDs used for the given sample. NOTE: This is specified in the `GEX` sample only in the sheet. MUST be specified for the GEX, otherwise it will crash. Leave the CP sample empty here.
- `cmotype` : `multi` or `single`. 
  - `multi`: If there are multiple CMOs pr "sample". The same "biological sample" were split and a CMO was added to each "subsample", which were then pooled to one "sample_id" (one row) in samplesheet. This would cause the CMO-"declaration" in the library.csv to be 
  ```
  [samples]
  sample_id,cmo_ids
  sample1,CMO301|CMO302|CMO303|CMO304
  ```
  - `single`: Each CMO represent one sample. Typically, different "biological samples" were treated separately with different CMOs, and then pooled to one "samplesheet sample". So the declaration will be:
  ```
  [samples]
  sample_id,cmo_ids
  sample1,CMO301
  sample2,CMO302
  sample3,CMO303
  sample4,CMO304
  ```
 Note that for cmotype=`single`, the pipeline will generate the following nomenclature on this "CMO-declaration" in the library.csv file:
```
[samples]
  sample_id,cmo_ids
  ${sid}-301,CMO301
  ${sid}-302,CMO302
  ${sid}-303,CMO303
  ${sid}-304,CMO304
  etc..
```

In the samplesheet template **below** it will generate these two csv:
- For CellPlex1: Take the `Sample` column: `CellP1`, and it is cmotype=multi:
 ```
  [samples]
  sample_id,cmo_ids
  CellP1,CMO301|CMO302
  ```
- For CellPlex3: Take the `Sample` column: `CellP3`, and it is cmotype=single: (one CMO pr sample)
 ```
  [samples]
  sample_id,cmo_ids
  CellP3-303,CMO303
  CellP3-304,CMO304
  ```


### Samplesheet template

#### Samplesheet name : `CTG_SampleSheet.sc-multi-10x.csv`

##### With CMO and ADT
```
[Header],,,,,,
metaid,2022_112,,,,,
cmo,y,
antibodies,"ADT_A0933, ADT_A0432, ADT_A0395, ADT_A0572, ADT_A0815, ADT_A0062, ADT_A0944, ADT_A0244, ADT_A0145, ADT_A0068, ADT_A0245, ADT_A0155, ADT_A0817",,,,,
email,per.a@med.lu.se,,,,,
autodeliver,y,,,,,
[Data],,,,,,
Lane,Sample_ID,index,Sample_Species,Sample_Project,Sample_Lib,Sample_Pair,Sample,CMO,cmotype
,s17,SI-TT-B7,human,2022_112,gex,1,Sample17,CMO301|CMO302,single
,s17_CP,SI-NN-H1,human,2022_112,cp,1,Sample17,,
,s17_ADT,ACATCG,human,2022_112,adt,1,Sample17,,
```

##### With ADT, without CMO
```
[Header],,,,,,
metaid,2022_112,,,,,
cmo,n,
antibodies,"ADT_A0933, ADT_A0432, ADT_A0395, ADT_A0572, ADT_A0815, ADT_A0062, ADT_A0944, ADT_A0244, ADT_A0145, ADT_A0068, ADT_A0245, ADT_A0155, ADT_A0817",,,,,
email,per.a@med.lu.se,,,,,
autodeliver,y,,,,,
[Data],,,,,,
Lane,Sample_ID,index,Sample_Species,Sample_Project,Sample_Lib,Sample_Pair,Sample,CMO,cmotype
,s17,SI-TT-B7,human,2022_112,gex,1,Sample17,,
,s17_ADT,ACATCG,human,2022_112,adt,1,Sample17,,
```

##### With CMO, without ADT
```
[Header],,,,,,
metaid,2022_112,,,,,
cmo,y,
antibodies,n,
email,per.a@med.lu.se,,,,,
autodeliver,y,,,,,
[Data],,,,,,
Lane,Sample_ID,index,Sample_Species,Sample_Project,Sample_Lib,Sample_Pair,Sample,CMO,cmotype
,s17,SI-TT-B7,human,2022_112,gex,1,Sample17,CMO301|CMO302,single
,s17_CP,SI-NN-H1,human,2022_112,cp,1,Sample17,,
```

### 2. Feature reference (feature.ref.csv)

- Pipeline driver will create `metadata/feature.ref.csv` file based on samplesheet:
1. Add antibodies declared in samplesheet, by matching the IDs with the totalseq reference (`/projects/fs1/shared/references/cellranger_totalseq`)
2. Add default CMO reference (CMO301-312) to this csv. It is located at `/projects/fs1/shared/references/cellranger_multiplex`.
3. This resulting `feature.ref.csv` will be added to `metadata` folder in ctg-project folder.

So, with both ADT and CMOs it will look something like: 
```
id,name,read,pattern,sequence,feature_type
ADT_A0933,Annexin.A1,R2,5P(BC),CCCACTGGAGCAATT,Antibody Capture
ADT_A0432,CD230,R2,5P(BC),CAGGTCCCTTATTTC,Antibody Capture
ADT_A0395,Hu.B7.H4,R2,5P(BC),TGTATGTCTGCCTTG,Antibody Capture
ADT_A0572,Hu.C5L2,R2,5P(BC),ACAATTTGTCTGCGA,Antibody Capture
ADT_A0815,Hu.CCR10,R2,5P(BC),ATCTGTATGTCACAG,Antibody Capture
ADT_A0062,Hu.CD10,R2,5P(BC),CAGCCATTCATTAGG,Antibody Capture
ADT_A0944,Hu.CD101,R2,5P(BC),CTACTTCCCTGTCAA,Antibody Capture
CMO301,CMO301,R2,5P(BC),ATGAGGAATTCCTGC,Multiplexing Capture
CMO302,CMO302,R2,5P(BC),CATGCCAATAGAGCG,Multiplexing Capture
CMO303,CMO303,R2,5P(BC),CCGTCGTCCAAGCAT,Multiplexing Capture
CMO304,CMO304,R2,5P(BC),AACGTTAATCACTCA,Multiplexing Capture
CMO305,CMO305,R2,5P(BC),CGCGATATGGTCGGA,Multiplexing Capture
CMO306,CMO306,R2,5P(BC),AAGATGAGGTCTGTG,Multiplexing Capture
CMO307,CMO307,R2,5P(BC),AAGCTCGTTGGAAGA,Multiplexing Capture
CMO308,CMO308,R2,5P(BC),CGGATTCCACATCAT,Multiplexing Capture
CMO309,CMO309,R2,5P(BC),GTTGATCTATAACAG,Multiplexing Capture
CMO310,CMO310,R2,5P(BC),GCAGGAGGTATCAAT,Multiplexing Capture
CMO311,CMO311,R2,5P(BC),GAATCGTGATTCTTC,Multiplexing Capture
CMO312,CMO312,R2,5P(BC),ACATGGTCAACGCTG,Multiplexing Capture
```
#### NB: 
- Recommended is to use the `antibodies,` declarartion in the header: see examples above. 
- The ADT/HTO references are found in `/projects/fs1/shared/references/cellranger_totalseq`. 
- The IDs must match "id" fields in these references
- The driver will use these ids to extract the information needed for the feature.ref.csv file 
- The driver will create the feature ref file and add it to nextflow.config to use in pipeline.
- If CMOs are used in the project, the entire cmo_ref.csv will be added to the feature.ref.csv (as in example above). This cmo_ref.csv is located at `/projects/fs1/shared/references/cellranger_multiplex/cmo`  
 
#### If not using antibodies, list in header
- If for some reason you will not list the antibodies to be used in the samplesheet, you must create such a file manually.
- Csv that declares the molecule structure and unique Feature Barcode sequence of each feature present in your experiment. Format as above.  
- This should either be added in runfolder (must be named /path/to/runfolder/**feature.ref.csv**).
- Or it could be specified by driver with -f flag: sc-multi-10x-driver -f /path/to/feature.ref.csv. 


See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis for more info

Example (TotalSeq A): 
| id | name | read | pattern | sequence | feature_type | 
| --- | --- | --- | --- | --- | --- | 
| CD235a | CD235a | R2 | ^(BC) | AGAGTATGTATGGGA | Antibody Capture | 
| CD33 | CD33 | R2 | ^(BC) | TAACTCAGGGCCTAT | Antibody Capture | 
| CD71 | CD71 | R2 | ^(BC) | CCGTGTTCCTCATTA | Antibody Capture | 
| CD11b | CD11b | R2 | ^(BC) | GACAAGTGATCTGCA | Antibody Capture | 
| CD45RA | CD45RA | R2 | ^(BC) | TCAATCCTTCCGCTT | Antibody Capture | 
| CD34 | CD34 | R2 | ^(BC) | GCAGAAATCTCCCTT | Antibody Capture | 
| CD49d | CD49d | R2 | ^(BC) | CCATTCAACTTCCGG | Antibody Capture | 
| CD45 | CD45 | R2 | ^(BC) | TCCCTTGCGATTTAC | Antibody Capture | 
| CD49d | CD49d | R2 | ^(BC) | CCATTCAACTTCCGG | Antibody Capture | 
| CD45 | CD45 | R2 | ^(BC) | TCCCTTGCGATTTAC | Antibody Capture | 
| CMO222 | CMO222 | R2 | 5P(BC) | GGCGTTAATCACTCT | Multiplexing Capture | 
| CMO112 | CMO112 | R2 | 5P(BC) | CCCGTTAATCACTCT | Multiplexing Capture | 
| CMO322 | CMO322 | R2 | 5P(BC) | TTCGAAAATCACTCT | Multiplexing Capture | 
| CMO422 | CMO422 | R2 | 5P(BC) | GACGGGAATCACTCT | Multiplexing Capture | 
| CMO522 | CMO522 | R2 | 5P(BC) | AGCGAAAATCACTCT | Multiplexing Capture | 


### Feature reference template 
#### Filename : `feature.ref.csv`
```
id,name,read,pattern,sequence,feature_type
CD235a,CD235a,R2,^(BC),AGAGTATGTATGGGA,Antibody Capture
CD33,CD33,R2,^(BC),TAACTCAGGGCCTAT,Antibody Capture
CD71,CD71,R2,^(BC),CCGTGTTCCTCATTA,Antibody Capture
CD11b,CD11b,R2,^(BC),GACAAGTGATCTGCA,Antibody Capture
CD45RA,CD45RA,R2,^(BC),TCAATCCTTCCGCTT,Antibody Capture
CD34,CD34,R2,^(BC),GCAGAAATCTCCCTT,Antibody Capture
CD49d,CD49d,R2,^(BC),CCATTCAACTTCCGG,Antibody Capture
CD45,CD45,R2,^(BC),TCCCTTGCGATTTAC,Antibody Capture
CMO304,CMO304,R2,5P(BC),AACGTTAATCACTCA,Multiplexing Capture
CMO305,CMO305,R2,5P(BC),CGCGATATGGTCGGA,Multiplexing Capture
CMO306,CMO306,R2,5P(BC),AAGATGAGGTCTGTG,Multiplexing Capture
CMO307,CMO307,R2,5P(BC),AAGCTCGTTGGAAGA,Multiplexing Capture
CMO308,CMO308,R2,5P(BC),CGGATTCCACATCAT,Multiplexing Capture

``` 

## Pipeline steps:

Cellranger version: cellranger v6.1.2 

* `delivery_info`: Generates the delivery info csv file needed to create delivery email (see deliverAuto below)
* `parse samplesheets`: Creates samplesheets (one for RNA, and one for ADT/HTO) for demux based on the input samplesheet. 
* `generate library csv`: Creates library.csv file based on input samplesheet. One .csv per matched RNA and ADT/HTO sample.
* `Demultiplexing` (cellranger mkfastq): Converts raw basecalls to fastq, and demultiplex samples based on index (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq). Does this separately for RNA and ADT/HTO (since they often have different index types (dual/single)
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). MultiQC summarizes FastQC reports into one document (https://multiqc.info/).
* `Align` + `Counts` + `Feature Barcoding` (cellranger multi): Aligns fastq files to reference genome, counts genes for each cell/barcode, and quantifies ADT/HTO features per barcode - Then performs secondary analysis such as clustering and generates the cloupe files (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi).
* `multiQC`: Compile fastQC and cellranger count metrics in multiqc report
* `md5sum`: md5sum of all generated files
* `deliverAuto`: executes delivery script in bin/ that edits template html delivery email, creates user, pass and executes mutt command to send email
* `pipe_done`: marks pipeline as done (puts ctg.sc-cite-seq-10x.$projid.done in runfolder and logs)


## Handle dual and single indexing in same sequencing run

RNA and ADT/HTO libraries must often have different indexing. It is handled by the pipeline by:
- It looks up the length of `adt` index, and setting the --use-bases-mask accordingly. If adt index is found to be 6 bases, it will set --use-bases-mask=Y28n*,I6n*,N10,Y90n* during mkfastq_adt. 
- By default, it will assume that RNA sample indices are dual, and ADT indices are single. It will thus set --filter-single-index during mkfastq_adt, and --filter-dual-index during mkfastq_rna. 
	
## Container
- `cellranger-v6.1.2`
https://github.com/perllb/ctg-sc-multi-10x/tree/master/container

Build container:
NOTE: Environment.yml file has to be in current working directory
```
sudo -E singularity build cellranger_v6.1.2.sif cellranger_v6.1.2-builder
```

Add path to .sif in nextflow.config

## Output:
* ctg-PROJ_ID-output
    * `qc`: Quality control output. 
        * cellranger metrics: Main metrics summarising the count / cell output 
        * fastqc output (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * multiqc output: Summarizing FastQC output and demultiplexing (https://multiqc.info/)
    * `fastq`: Contains raw fastq files from cellranger mkfastq.
    * `count-cr`: Cellranger multi output. Here you find gene/cell count matrices, feature quantification, secondary analysis output, and more. See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi for more information on the output files.
    * `summaries`: 
        * web-summary files which provide an overview of essential metrics from the 10x run. 
        * cloupe files which can be used to explore the data interactively in the Loupe browser (https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/what-is-loupe-cell-browser)  
    * `ctg-md5.PROJ_ID.txt`: text file with md5sum recursively from output dir root    



## Custom genome 

If custom genome (not hg38 or mm10) is used

1. Set "Sample_Species" column to 'custom' in samplesheet:

Example:
 | Lane | Sample_ID | index | Sample_Species | Sample_Project | Sample_Lib | Sample_Pair | Sample | CMO | cmotype |
 | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
 | | Sr1 | SI-GA-D9 | **custom** | 2022_022 | gex | 1 | S1 | CMO301\|CMO302 | multi |
 | | Sr1_CP | SI-GA-H9 | **custom** | 2022_022 | cp | 1 | S1 | | |
 | | Sr1_ADT | ATGCAG | **custom** | 2022_022 | adt | 1 | S1 | |
 | | Sr2 | SI-GA-C9 | **custom** | 2022_022 | gex | 2 | S2 | CMO303\|CMO304 | single |
 | | Sr2_CP | SI-GA-C9 | **custom** | 2022_022 | cp | 2 | S2 | |
 | | Sr2_ADT | ATGTTG | **custom** | 2022_022 | adt | 2 | S2 | |

2. In nextflow.config, set 
 `custom_genome=/PATH/TO/CUSTOMGENOME`
 
### Add custom genes (e.g. reporters) to cellranger annotation

Use the `ctg-cellranger-add2ref` script. 

https://github.com/perllb/ctg-cellranger-add2ref

