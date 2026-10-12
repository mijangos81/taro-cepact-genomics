# taro-cepact-genomics

Analysis code for:

> Restall J, Mijangos JL, Kilian A, Campbell B. Genome-wide SNPs reveal clonal
> redundancy and connectivity in the world's largest taro collection.
> *Horticulture Research* (submitted).

`reproduce_manuscript.R` regenerates every figure, table and reported number of
the paper and its Supplementary Information from the DArTseq genotypes of the
taro (*Colocasia esculenta*) collection at the Centre for Pacific Crops and Trees
(CePaCT). It runs from the repository root in ten sections, in the order the
analyses depend on each other; the section list is at the top of the script.

## Data

Download the data deposit (Zenodo, DOI [10.5281/zenodo.23307371](https://doi.org/10.5281/zenodo.23307371)), decompress the `.gz` files and
place them as follows:

| Files | Folder |
|---|---|
| `Report_DTa25-11061_*_SNP_mapping_{4,5}.csv`, `SEQ_SNPs_counts_{0,1,2}_Target.csv`, `taro*.rds` | `input_files/` |
| `metadata_DTa23-8697.csv`, `coor_countries.csv`, `ID_targets.csv` | `info/` |
| `ploidy.csv` | `outputs/` |

The script downloads the annotated Ahi assembly (GCA_009445465.1) and Swiss-Prot
itself. Read alignments are not deposited, so the nQuack block and the allele-depth
re-test on BAM files (both in section 3) cannot be re-run; they are marked in the
script.

## Software

Developed with R 4.4.2 and dartRverse 1.2.6 (dartR.base, dartR.popgen), adegenet
2.1.12, poppr 2.9.8, LEA 3.18.0, ape 5.8-1, nQuack 1.0.2, Qploidy2 1.9.1, rjags
4-17 and vcfR 1.15.0; the full package list is at the top of the script.

- `gl.find.genes.for.loci` needs dartR.popgen from the dev branch:
  `remotes::install_github("green-striped-gecko/dartR.popgen@dev")`.
- gbs2ploidy is installed from the CRAN archive
  (https://cran.r-project.org/src/contrib/Archive/gbs2ploidy/).
- Binaries on `PATH`: plink; OHANA (`qpas`, `nemeco`, `selscan`, `convert`; put the
  OHANA folder first on `PATH`, because these names clash with the iHS selscan
  and ImageMagick); PHYLIP (`fitch`, `consense`); BLAST+ (`blastn`, `blastp`,
  `makeblastdb`); samtools; bcftools; PopCluster; curl.

## Citation

Cite the paper and this archive (Zenodo, DOI [CODE_DOI]); see `CITATION.cff`.

## Licence

Code: MIT (see `LICENSE`). Data: CC BY 4.0 (see the data deposit).
