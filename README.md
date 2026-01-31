# CReSIL: Accurate Identification of Extrachromosomal Circular DNA from Long-read Sequences

CReSIL is a tool for detecting eccDNA from long-read sequencing data.

:sparkles: Version `V1.2.0+hifi` - Now supports both **Oxford Nanopore (ONT)** and **PacBio HiFi** platforms!

For a complete list of changes, please see the [full changelog](https://github.com/YxZhang-XHCY/cresil-hifi/releases).



## Installation instructions
```
## Install Conda environment
cd cresil
conda env create -f environment.yml

## Activate CReSIL environment
conda activate cresil

## Install CReSIL
pip install .
```

## Testing (Human eccDNA) - ONT Data
For this demonstration, we use the `r941_min_hac_g507` model with a default breakpoint overlapping length of 50 bp. The required human reference genome can be loaded from [the Sample data section](#sample-data).
> [!NOTE]
> The default consensus model has been updated to `r1041_e82_400bps_sup_v4.3.0` since release `V1.1.0`
```bash
## Go to an example folder
cd example

## Run trim
cresil trim -t 4 -fq exp_reads.fastq -r reference.mmi -o cresil_result

## Run eccDNA identification for enriched data
cresil identify -t 4 -fa reference.fa -fai reference.fa.fai -fq exp_reads.fastq -cm r941_min_hac_g507 -trim cresil_result/trim.txt

## Run eccDNA annotation
cresil annotate -t 4 -rp reference.rmsk.bed -cg reference.cpg.bed -gb reference.gene.bed -identify cresil_result/eccDNA_final.txt

## Run visualize eccDNA (ec1)
cresil visualize -t 4 -c ec1 -identify cresil_result/eccDNA_final.txt

## Run Circos (ec1)
cd cresil_result/for_Circos/ec1
circos -noparanoid -conf circos.conf
```

## Testing (Human eccDNA) - PacBio HiFi Data

For PacBio HiFi data, use the `-platform hifi` parameter. HiFi mode automatically:
- Uses `map-hifi` preset for minimap2
- Skips Medaka sequence correction (HiFi has high accuracy ~99.9%)

```bash
## Run trim with HiFi mode
cresil trim -t 4 -platform hifi -fq hifi_reads.fastq -r reference.mmi -o cresil_result

## Run eccDNA identification with HiFi mode (no -cm needed, Medaka is skipped)
cresil identify -t 4 -platform hifi -fa reference.fa -fai reference.fa.fai -fq hifi_reads.fastq -trim cresil_result/trim.txt

## Run eccDNA annotation (same as ONT)
cresil annotate -t 4 -rp reference.rmsk.bed -cg reference.cpg.bed -gb reference.gene.bed -identify cresil_result/eccDNA_final.txt

## Run visualize eccDNA (same as ONT)
cresil visualize -t 4 -c ec1 -identify cresil_result/eccDNA_final.txt
```

The structure of files and directories generated from CReSIL is as follows:

```
cresil_result/
├── trim.txt
├── eccDNA_final.txt
├── cresil_run/
│   ├── subGraphs.summary.txt
│   ├── tmp/
│   └── assemGraph/
│       ├── ec1/
│       ├── ec2/
│       └── ..
├── cresil_gAnnotation/
│   ├── gene.annotate.txt
│   ├── CpG.annotate.txt
│   ├── repeat.annotate.txt
│   └── variant.annotate.txt
└── for_Circos/
    └── ec1/
        ├── circos.conf
        ├── ticks.conf
        ├── ideogram.conf
        ├── karyotype.conf
        ├── data/
        └── tmp/
```
For more details about mapped reads for each eccDNA, please see the eccDNA folder in **assemGraph**.

## Additional useful commands

To run CReSIL to identify eccDNA without sequence correction and variant calling steps - skip mode (from the trimming step)

```
## Run eccDNA identification without sequence correction and variant calling steps
cresil identify -t 4 -fa reference.fa -fai reference.fa.fai -fq exp_reads.fastq -cm r941_min_hac_g507 -s -trim cresil_result/trim.txt
```

To run CReSIL to identify eccDNA in a relaxed mode (from the trimming step)

```
## Run eccDNA identification with minimizing parameters (minimum size of eccDNA, average depth, number of supported breakpoints)
cresil identify -t 4 -minrsize 40 -depth 1 -break 1 -fa reference.fa -fai reference.fa.fai -fq exp_reads.fastq -cm r941_min_hac_g507 -trim cresil_result/trim.txt
```

To use CReSIL to identify eccDNA in whole-genome long-read (WGLS) sequencing data (from the trimming step)

```bash
## Run eccDNA identification for WGLS data (ONT)
cresil identify_wgls -t 4 -r reference.mmi -fa reference.fa -fai reference.fa.fai -fq exp_reads.fastq -trim cresil_result/trim.txt

## Run eccDNA identification for WGLS data (HiFi)
cresil identify_wgls -t 4 -platform hifi -r reference.mmi -fa reference.fa -fai reference.fa.fai -fq hifi_reads.fastq -trim cresil_result/trim.txt
```

## Interface
 - `cresil trim` - find and trim potential eccDNA regions from long reads
 - `cresil identify` - identify and verify eccDNA from enriched data
 - `cresil identify_wgls` - identify and verify eccDNA from whole-genome long-read (WGLS) sequencing data
 - `cresil annotate` - annotate identified eccDNA
 - `cresil visualize` - generate Circos configuration files for specified eccDNA

### Platform Parameter

The `-platform` parameter is available for `trim`, `identify`, and `identify_wgls` commands:

| Platform | Description |
|----------|-------------|
| `ont` (default) | Oxford Nanopore data, uses `map-ont` preset, runs Medaka for sequence correction |
| `hifi` | PacBio HiFi data, uses `map-hifi` preset, skips Medaka (high accuracy data) |

## Sample data
[References and simulated data](https://app.box.com/s/5leixacmp1xx8qs7qtcuyz93lgqpzgkn) for the user to go through the example of CReSIL pipeline and used in the CReSIL benchmarks.
> [!WARNING]
>We encourage users to use a primary assembly without patches or alternative loci as a reference to reduce false positives and false negatives from predictions

## Citation
Please cite the following article if you use CReSIL in your research
> Visanu Wanchai, Piroon Jenjareonpun, Thongpun Leangapichat, Gerard Arrey, Charles M Burnham, Maria C Tümmle, Jesus Delgado-Calle, Birgitte Regenberg, Intawat Nookaew, CReSIL: Accurate Identification of Extrachromosomal Circular DNA from Long-read Sequences, *Briefings in Bioinformatics.* 2022;, bbac422, https://doi.org/10.1093/bib/bbac422.<br>

See our recent work using CReSIL `V1.2.0` at
> Charles M. Burnham, Alongkorn Kurilung, Visanu Wanchai, Birgritte Regenberg, Jesus Delgado-Calle, Alexei G. Basnakian, Intawat Nookaew, An Enhancement of Extrachromosomal Circular DNA Enrichment and Amplification to Address the Extremely Low Overlap Between Replicates, *bioRxiv* 2025.06.28.662146; doi: https://doi.org/10.1101/2025.06.28.662146.<br>

## License and Copyright

CReSIL is distributed under the terms of the MIT License
