# SampleQC
The Snakemake pipeline integrates NTSM and VerifyBamID to verify sample identity and detect potential contamination in sequencing data.

## Requirements
- [Snakemake 7+](https://snakemake.readthedocs.io/)
- Singularity (or Apptainer) for running the ntsm Docker image
- Python libraries (for plotting step):
  - `pandas`
  - `numpy`
  - `seaborn`
  - `matplotlib`

---

## Input Data
### Manifest (`manifest.tab`)
The manifest must include:
- `ID` — Unique identifier for each dataset
- `FOFN` — File-of-filenames (list of FASTQ/FASTA files)
- `TYPE` - Sequencing platform or data type
Choose one of the following:
  - `PacBio`
  - `ONT`
  - `Illumina`

Example:
```
ID	FOFN	TYPE
SampleA	fofn/SampleA.fofn	PacBio
SampleB	fofn/SampleB.fofn	ONT
SampleC	fofn/SampleC.fofn	Illumina
```


### Config (`config.yaml`)
Important keys:
- `MANIFEST`: Path to the manifest file
- `REF_SITE`: Reference sites fasta (default: `db_source/human_sites_n10.fa`, relative to the Snakefile directory)
- `EXTERNAL_COUNTS_DIR`: Directory containing external count files (optional)
- `COUNT_FILE_EXP`: File extension for the external count files(default: `count`)

---

## Running the Pipeline
1. Edit `config.yaml` and `manifest.tab` to reflect your datasets.
2. Run Snakemake:
```bash
   ln -s /net/eichler/vol28/software/pipelines/ntsm_smk/runcluster .
   ./runcluster 30
```


### References
- NTSM: https://github.com/JustinChu/ntsm  
- VerifyBamID: https://github.com/Griffan/VerifyBamID.git
