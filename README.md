# GeneiASE_Tools
R scripts to process input/output related to GeneiASE software.

GeneiASE software can be downloaded from: https://github.com/edsgard/geneiase

---

## createInputForGeneiASE_static.R

Create the input file for the GeneiASE static analysis. The script needs the file of filtered SNPs, the vep annotation file, the depth of alignment at each position and the rejected depth of alignement for the fitering conditions (both generated using mpileupTools (https://github.com/adeschen/mpileupTools)), the minimum depth of alignment for the SNP to be used in the analysis.

