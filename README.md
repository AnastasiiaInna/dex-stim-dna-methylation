# __DEX-stimulated DNAm arrays processing__

## **Issues**

__On cluster:__

1. Installation `minfi` : 

```R
Sys.setenv(LC_CTYPE="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

BiocManager::install("minfi")

## OR

install.packages("remotes")
remotes::install_github("hansenlab/minfi")
````
or
```sh
export LC_CTYPE="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
```
m
```sh
srun -u --pty --part=pe -c 12 --mem=50G R --vanilla
```

## __Note:__

_from https://bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html_

- For each CpG, there are two measurements: a methylated intensity (denoted by M) and an unmethylated intensity (denoted by U). These intensity values can be used to determine the proportion of methylation at each CpG locus. Methylation levels are commonly reported as either **beta** values (_β = M/(M+U)_) or **M**-values (_M-value = log2(M/U)_).

- For practical purposes, a small offset, α, can be added to the denominator of the β value equation to avoid dividing by small values, which is the default behaviour of the getBeta function in minfi. The default value for α is 100. It may also be desirable to add a small offset to the numerator and denominator when calculating M-values to avoid dividing by zero in rare cases, however the default getM function in minfi does not do this. 

- Beta values and M-values are related through a logit transformation. 

- **Beta values** are generally preferable for describing the level of methylation at a locus or for graphical presentation because percentage methylation is easily interpretable. However, due to their distributional properties, **M-values** are more appropriate for statistical testing 

## **Input data:**

- iData data : `/binder/mgp/workspace/2020_DexStim_Array_Human/methylation/`
- RData: 

`/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/`

`/home/ahryhorzhevska/mpip/datasets/methyl/rData` (starting with normalization step, the data are stored in this folder)

    * `rgSet_dex`: raw data from the dex IDAT files; organized at the probe (not CpG locus) level and has two channels (Red and Green).
    
    The following data are obtained from __rgSet_dex__:
    
    * `Mset_original`: data organized by the CpG locus level, but not mapped to a genome and has two channels, _Meth (methylated)_ and _Unmeth (unmethylated)_;
    * `RatioSet_original`: data organized by the CpG locus level, but not mapped to a genome, and has at least one of two channels, _Beta_ and/or _M (logratio of Beta)_;
    * `RawBetas`: _Beta_ values matrix, obtained from _RatioSet_;
    * `gRatioSet_original`: the same as _RatioSet_ but mapped to a genome;
    * `pd_original`: phenotype data

To get _rgSet_ and the rest of the data, run:

```sh
cd ~/01_prepare_data_formats/
chmod +x ./getRgSetFormats.sh
sbatch ./getRgSetFormats.sh
```

## **QC:**

1.  Calculate detection p-values. We can generate a detection p-value for every CpG in every sample, which is indicative of the quality of the signal. The method used by minfi to calculate detection p-values compares the total signal (_M_ + _U_) for each probe to the background signal level, which is estimated from the negative control probes. Very small p-values are indicative of a reliable signal whilst large p-values, for example _> 0.01_, generally indicate a poor quality signal.

2. Remove poor quality samples with a mean detection p-value _> 0.01_:

* from _rgSet_
* from _beta_ matrix
* from _p-values_ table
* from _targets_ table

_Result:_

**Data**

> `rgSet_qc.Rdata`

> `RawBetas_qc.Rdata`

> `detP_qc.Rdata` : _p-values_ table;

**Reports**

>  `detectionP.pdf`

> `qcReport.pdf`

3. Distribution artefacts

_Result:_

**Reports**

> `03_beta_densities_report.pdf`

**Suspicious observations:**

> From array _R01C01_;

> Partially from _R03C01_;

> Sample _200712160065_.

4. Sex mismatches

By looking at the median total intensity of the X chromosome-mapped probes, denoted _xMed_, and the median total intensity of the Y-chromosome-mapped probes, denoted _yMed_, one can observe two different clusters of points corresponding to which gender the samples belong to. To predict the gender, minfi separates the points by using a cutoff on _log2(yMed)−log2(xMed)_. The default cutoff is −2. 

_Result:_

_The result should be cheked with genotype data_

**Data**

> `sex_predicted.Rdata`

5. Remove sampleID found at the steps _3_ and _4_

_Result:_

**Data**

> `rgSet_clean.Rdata`

> `RawBetas_clean.Rdata`

> `detP_clean.Rdata`

> `pd_celan.Rdata`

> `annotated_data_clean.Rdata`

## **Normalization:**

_Result:_

**Data**

Folder: `/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/rData`:

> `BMIQ.quantileN.Rdata`

>`quantileN.Rdata`

>`quantileNBetas.Rdata`

>`quantileNMs.Rdat`

Folder: `/home/ahryhorzhevska/mpip/datasets/methyl/rData`

> `BMIQ_quantileN.Rdata`

**Reports**

> `BetaValue_Distributions_Norm_Quantile.pdf`

## **Filtering:**

_Result:_

**Data**

Folder: `/home/ahryhorzhevska/mpip/datasets/methyl/rData`

> `quantileN_filtered.Rdata`

>`BMIQ.quantileN_filtered.Rdata` 

>`Betas_quantileN_filtered.Rdata`

>`Ms_quantileN_filtered.Rdata`

**Reports**

> `BetaValue_Distributions_Norm_quantile_Filter.png`

## **Batch effects correction:**

Description will come soon

## **MixupMapper:**

Required files description can be found here:

https://github.com/molgenis/systemsgenetics/wiki/File-descriptions

**Data**

Methylation data :

> `/home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper`

SNPs:

> `/home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps/mixupmapper/`

_Initial step:_

> `imputed genotypes: /home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps/Dex_genoData_SNPs.bed`
>` batch-adjusted beta values: /binder/mgp/datasets/2020_DexStim_Array_Human/methylation/rData/Betas_combated.Rdata`

_Intermediate step:_

> `TRITYPER genotypes: /home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps/mixupmapper/`
> `Trait file: batch-adjusted beta values: /home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/beta_combated_for_mixupmapper.txt`
> `Genotype - phenotype coupling: /home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/genotypemethylationcoupling.txt`

_Results:_

**Run:**

1. Convert genotypes _.bed_ to _TriTyper_ format:

```sh
DIR_GENOTYPES_BED="/home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps/Dex_genoData_SNPs"
DIR_GENOTYPES_TRITYPER="/home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps/mixupmapper"

wget https://github.com/molgenis/systemsgenetics/releases/download/1.4.0_20-8.1/GenotypeHarmonizer-1.4.23-dist.tar.gz
tar -xvf GenotypeHarmonizer-1.4.23-dist.tar.gz
cd GenotypeHarmonizer-1.4.23

java -jar GenotypeHarmonizer.jar -i $DIR_GENOTYPES_BED -I PLINK_BED -o $DIR_GENOTYPES_TRITYPER -O TRITYPER  --update-id
```

2. Beta values normalization

```sh
screen -S mixupmapper

BETA_VALUES_FILENAME="/home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/betas_combat_veh_mixupmapper_final.txt"
MIXUPMAPPER_DATA_DIR="/home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/"

MIXUPMAPPER_DIR="/home/ahryhorzhevska/github/dex-stim-dna-methylation/06_mixupmapper/eqtl-mapping-pipeline-1.2.4E-SNAPSHOT"
java -jar $MIXUPMAPPER_DIR/eqtl-mapping-pipeline.jar --mode normalize --in $BETAS_VALUE_FILENAME  --out $MIXUPMAPPER_DATA_DIR --centerscale
```

Result : Exception in thread "main" java.lang.NoClassDefFoundError: umcg/genetica/io/Gpio


