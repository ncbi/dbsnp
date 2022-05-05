# GrafPop Software Documentation

GrafPop is a software tool that can be used for inferring subject ancestry with genotypes.

A C++ executable `grafpop` and other files are included in the package `GrafPop.tar.gz` and visible as separate files after the user executes the command:
```sh
tar zxvf GrafPop.tar.gz
```

GrafPop calculates genetic distances from each subject (or sample) to several reference populations and estimates subject ancestry and ancestral proportions based on these distances. See [Jin et al., 2019](https://www.g3journal.org/content/9/8/2447.long) for the algorithm, which was first implemented in the GRAF-pop feature of the [GRAF](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/Software.cgi) software package, using the 10,000 fingerprint SNPs ([Jin et al., 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0179106)).

Four genetic distances scores, GD1, GD2, GD3, GD4, are used in ancestry inference in the current version of GrafPop. Subjects in the input datasets are clustered using these scores and plotted on scatter plots. GrafPop assumes that each subject is an admixture of three ancestries: European (E), African (F), and Asian (A), and estimates ancestral proportions _P<sub>e</sub>, P<sub>f</sub>, P<sub>a</sub>_ based on GD1 and GD2 scores using barycentric coordinates. It also assigns a population ID (PopID) to each subject using the cutoff values shown in Tables 1 and 2, by default.
 
_Table 1. Grouping subjects based on the ancestry proportions_

| PopID | Population | Cutoff standard |
| --- | --- | --- |
| 1 | European | P<sub>e</sub> ≥ 90% |
| 2 | African | P<sub>f</sub> ≥ 95% |
| 3 | East Asian | P<sub>a</sub> ≥ 95% |
| 4 | African American | 40% ≤ P<sub>f</sub> < 95% and P<sub>a</sub> < 14% |
| 5 | Latin American 1 | P<sub>f</sub> < 50% and P<sub>e</sub> < 90% and P<sub>a</sub> < 14% and GD1 < 1.48 |
| 6,7,8 | (Three populations) | Otherwise and P<sub>f</sub> < 14% |
| 9 | Other | P<sub>a</sub> ≥ 14% and P<sub>f</sub> ≥ 14% |

_Table 2. Separating Asians and Latin Americans using GD1 and GD4 scores_

| PopID | Population | Cutoff standard |
| --- | --- | --- |
| 7 | Asian-Pacific Islander | GD1 > 30 × (GD4)<sup>2</sup> + 1.58 |
| 8 | South Asian | GD4 > 5 × (GD1 - 1.524)<sup>2</sup> + 0.0575 |
| 6 | Latin American 2 | GD1 + GD4 < 1.525 and PopID is not 7 |


### Input files

GrafPop takes genotype datasets in either PLINK format (`.fam`, `.bim`, `.bed`) or VCF format (`.vcf` or `.vcf.gz`). In addition, GrafPop can read self-reported races/ethnicities from an input file and compare the populations inferred from genotypes with the self-reported ones. The input file should be a plain text file with two columns (without column header), containing subject IDs (or sample IDs, depending on the type of IDs in the genotype dataset) and the self-reported races/ethnicities, respectively.

### Running `grafpop` to infer subject ancestry

`grafpop` is a command line executable that can be run under GNU/LINUX 64 bit systems.  Brief usage is given when the program is executed without parameters:

```sh
$ grafpop

Usage: grafpop <Binary PLINK set or VCF file> <output file>

```

The following command determines population structures and saves results to the output file:

```sh
$ grafpop data/TGP_anc_geno.bed results/TGP_pop_scores.txt 
```
The SNPs in the PLINK set are entered as RS IDs. If the SNP IDs are not provided, or provided but not in RS IDs, it is acceptable to `grafpop` if GRCh 37 or GRCh 38 chromosome positions are included in the genotype dataset, e.g.,
```sh
$ grafpop data/TG_10_1000_g37.bed results/TG_g37_pops.txt
$ grafpop data/TG_10_1000_g38.bed results/TG_g38_pops.txt
```

The input dataset can also be a VCF file, e.g.,
```sh
$ grafpop data/TGP_anc_geno.bed results/TGP_pop_scores.txt 
```

`grafpop` also accepts a zipped VCF file, e.g.,
```sh
$ grafpop data/TG_2_zip_chr2.vcf.gz results/TG_2_zip_pops.txt
```
If the VCF file includes many more SNPs than those being used by GrafPop, e.g., containing whole genome sequencing data,  `grafpop` can still read the data and do ancestry inference. However, it is recommend that the Perl script `ExtractAncSnpsFromVcfGz.pl` be used to extract the genotypes before `grafpop` is run (see instructions below for usage of the Perl script).

`grafpop` and the Perl scripts included in the package can be called from other directories, e.g.,

```sh
$ cd results
$ ../grafpop ../data/TG_10_1000_g38.bed TG_g38_pops.txt
```
GrafPop C++ executable and Perl scripts need to find some information included in the `data` directory when being run. If the executable is moved away from the `data` directory, the user can set environment variable `GRAFPATH` to the directory where GrafPop `data` directory and Perl packages (`.pm` files) are located and call `grafpop` and Perl scripts from any location. 

### Running `PlotGrafPopResults.pl` to plot population results

The results generated by `graf` can be passed to `PlotGrafPopResults.pl` for further processing. The following instructions are displayed on the screen when the script is run without parameters:

```sh
$ PlotGrafPopResults.pl

Usage: PlotGrafPopResults.pl <input file> <output file> [Options]

    Note:
          Input file is the file generated by the C++ grafpop program that includes subject ancestry scores.
          Output file should be a .png file.
          Options should be entered after the two required parameters.

    Options:
        Specify the input file with self-reported subject race information
            -spf     text file with two columns (no header): subject and self-reported population

        Set window size in pixels
            -gw      graph width (500 - 2000, default 800)

        Set graph axis limits, max - min should be between 0.1 and 1.5
            -xmin    min x value
            -xmax    max x value
            -ymin    min y value
            -ymax    max y value

        Set minimum and maximum numbers of genotyped fingerprint SNPs for samples to be processed
            -minsnp  minimum number of SNPs with genotypes
            -maxsnp  maximum number of SNPs with genotypes

        Set population cutoff lines
            -ecut    proportion: cutoff European proportion dividing Europeans from other populations. Default 90%.
            -fcut    proportion: cutoff African proportion dividing Africans from other populations. Default 95%.
            -acut    proportion: cutoff East Asian proportion dividing East Asians from other populations. Default 95%.
            -ohcut   proportion: cutoff African proportion dividing Latin Americans from Other population. Default 14%.
            -fhcut   proportion: cutoff African proportion dividing Latin Americans from African Americans. Default 40%.

        Select some self-reported populations (by IDs) to be highlighted on the graph (for studies with multiple races)
            -pops    comma separated population IDs, e.g., -pops 1,3,4 -> highlight populations #1, #3 and #4

        Select self-reported populations (by IDs) to show areas including 95% dbGaP subjects with genotypes of
        at least 50,000 ancestry SNPs, to help estimate subject populations
            -areas   comma separated dbGaP self-reported population IDs, e.g., -areas 1,3
                         -> show areas that include 95% dbGaP subjects with self-reported populations #1 and #3
                            1: European                 2: African
                            3: East Asian               4: African American
                            5: Latin American 1         6: Latin American 2
                            7: Asian-Pacific Islander   8: South Asian

        Select which score to show on the y-axis
            -gd4:    show GD4 on y-axis (GD4 separates South Asians from Latin Americans and other Asians)

        Set population cutoff lines
            -cutoff: show cutoff lines

        Rotate the plot on x-axis by a certain angle
            -rotx    angle (in degrees) to rotate the GD2/GD3 vs. GD1 plot on the x-axis (0 - 360)

        Set the size (diameter) of each dot that represents each subject
            -dot     dot size in pixels (1 - 10, default 2)
```

The script takes two required parameters, which must be the first two arguments and are not preceded by flags, unlike all the optional arguments, which are preceded by a flag.  The first parameter should be the name of the file that is generated by `grafpop` and contains subject genetic distance scores. The second parameter is the output file, expected to be `.png` file.  The script processes the scores and saves the results to the output file.  The default graph is GD2 vs. GD21, e.g.,

```sh
$ grafpop data/TGP_anc_geno.bed results/TGP_pop_scores.txt
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops.png
```

When option `-gd4` is specified, the script generates a graph of GD4 vs. GD1:

```sh
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_gd4.png -gd4
```

When self-reported populations are available, the information can be passed to the script with `-spf` option so that the script can color-code the subjects using the self-reported populations, e.g.,

```sh
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_sp.png -spf data/TGP_SbjSuperPop.txt 
```

The format of the input race/ethnicity/population file is described above. In the graph generated by the script, the populations are numbered and color coded.

The cutoff lines used to partition the subjects are drawn on the graphs when option `-cutoff` is set, e.g.,

```sh
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_sp_cut.png -spf data/TGP_SbjSuperPop.txt -cutoff 
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_gd4_sp_cut.png -gd4 -spf data/TGP_SbjSuperPop.txt -cutoff
```

If multiple subjects appear at the same locations in the x-y plane, the user can use option `-pops` to bring some populations to the front, while setting some populations to the back and fade them out in the graph.  For example, the following command generates a graph with the European populations (No. 3, 4, 14, 16, 21) to the front with different colors and other populations in the back and colored yellow. The assignments of colors to populations are currently hard-coded.

```sh
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_p.png -spf data/TGP_SbjPop.txt
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_p_eur.png  -spf data/TGP_SbjPop.txt -pops 3,4,14,16,21
```

The population numbers following `-pops` should be separated by commas without spaces. All the populations are listed on the screen and numbered: 
```sh
Self-reported races/ethnicities
1: GWD (n=113)
2: YRI (n=108)
3: TSI (n=107)
4: IBS (n=107)
5: CHS (n=105)
...
```

One can also use the `-rotx` option to rotate the graph of GD2 vs. GD1 around x-axis by a certain angle specified in degrees (can be any real number).  For example, the following command generates a graph showing the subjects rotated by 90°:

```sh
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_sp_rot90.png -spf data/TGP_SbjSuperPop.txt -rotx 90
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_sp_rot120.png -spf data/TGP_SbjSuperPop.txt -rotx 120
```

Options `-gw, -xmin, -xmax, -ymin, -ymax, -dot`, can be used to adjust the graph size, specify axis limits, and set the dot size, e.g.,

```sh
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_p_eur_gw1200.png  -spf data/TGP_SbjPop.txt -pops 3,4,14,16,21 -gw 1200
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_p_eur_minmax.png  -spf data/TGP_SbjPop.txt -pops 3,4,14,16,21 -xmin 1.42 -xmax 1.54 -ymin 1.38 -ymax 1.5 
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_p_eur_minmax_dot5.png  -spf data/TGP_SbjPop.txt -pops 3,4,14,16,21 -xmin 1.42 -xmax 1.54 -ymin 1.38 -ymax 1.5 -dot 5

$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_p_eas_rng.png -spf data/TGP_SbjPop.txt -pops 5,6,8,12,20 -xmin 1.66 -ymin 1.05 -ymax 1.19 -dot 5
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_p_afr_rng.png -spf data/TGP_SbjPop.txt -pops 1,2,13,15,17,24,26, -xmax 1.4 -ymin 1.08 -ymax 1.48 -dot 5
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_p_sas_amr.png -spf data/TGP_SbjPop.txt -gd4 -pops 9,10,11,18,22,7,19,23,25
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_p_sas_amr_rng.png -spf data/TGP_SbjPop.txt -gd4 -pops 9,10,11,18,22,7,19,23,25 -xmin 1.4 -xmax 1.6 -ymin -0.04 -ymax 0.14 -dot 5
```

 One can use the option `-areas` to select populations to show the expected oval areas that include 95% of dbGaP subjects with at least 50,000 ancestry SNPs with genotypes, e.g.,

```sh
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_sp_areas.png -spf data/TGP_SbjSuperPop.txt -areas 1,4,7,8
```

The integers in the comma-delimited string represent the eight self-reported population groups in dbGaP, with most common ancestry terms in each group shown below:

```sh
    1: European                 2: African
    3: East Asian               4: African American
    5: Latin American 1         6: Latin American 2
    7: Asian-Pacific Islander   8: South Asian
```

GrafPop uses the ancestry proportions shown in Tables 1 and 2 as default cutoff values.  The user can use options `-ecut`, `-fcut`, `-acut`, `-ohcut`, `-fhcut` to set the cutoff values to different numbers, e.g.,

```sh
$ PlotGrafPopResults.pl results/TGP_pop_scores.txt results/TGP_pops_sp_popcut.png -spf data/TGP_SbjSuperPop.txt -ecut 85 -fcut 85 -fhcut 30 -ohcut 15
```

### Running `SaveSamples.pl` to save samples and their ancestry scores and population assignments to a file
Similar to `PlotGrafPopResults.pl`, `SaveSamples.pl` takes two parameters, and the first one is the file outputted by `grafpop`. The second parameter is the output file to save the samples.

```sh
$ SaveSamples.pl

Usage: SaveSamples.pl <input file> <output file> [Options]

    Note:
          Input file is the file generated by the C++ grafpop program that includes subject ancestry scores.
          Samples and ancestry scores will be saved to the output file as plain texts.

    Options:
        Set a rectangle area to retrieve subjects from graph of GD2 (y) vs. GD1 (x)
            -xcmin   min x value
            -xcmax   max x value
            -ycmin   min y value
            -ycmax   max y value
            -isByd:  retrieve subjects whose values are beyond the above rectangle

        Set minimum and maximum numbers of genotyped fingerprint SNPs for samples to be processed
            -minsnp  minimum number of SNPs with genotypes
            -maxsnp  maximum number of SNPs with genotypes

        Set population cutoff lines
            -ecut    proportion: cutoff European proportion dividing Europeans from other populations. Default 90%.
            -fcut    proportion: cutoff African proportion dividing Africans from other populations. Default 95%.
                                 Set it to -1 to combine African and African American populations
            -acut    proportion: cutoff East Asian proportion dividing East Asians from other populations. Default 95%.
                                 Set it to -1 to combine East Asian and Other Asian populations
            -ohcut   proportion: cutoff African proportion dividing Latin Americans from Other population. Default 13%.
            -fhcut   proportion: cutoff African proportion dividing Latin Americans from African Americans. Default 40%.

        The input file with self-reported subject race information
            -spf     a file with two columns: subject and self-reported population
```

The following command saves all samples into the output file:
```sh
$ SaveSamples.pl results/TGP_pop_scores.txt results/TGP_pop_smps.txt
```
Options `-xcmin`, `-xcmax`, `-ycmin`, `-ycmax`, `-isByd` can be used to specify a rectangular area and let the script retrieve samples whose x(GD1), y(GD2) scores are either within or beyond this area. For example, the following command saves all samples with 1.7 < GD1 < ∞, which are all EAS samples.

```sh
$ SaveSamples.pl results/TGP_pop_scores.txt results/TGP_pop_eas_smps.txt -xcmin 1.7
```

When option `-isByd` is set to 1, the script retrieves subjects whose values are beyond the rectangular area specified by options `-xcmin`, `-xcmax`, `-ycmin`, `-ycmax`.  For example, the following command excludes most of the 1000 Genomes Project's subjects with super populations AMR (Ad Mixed American) and SAS (South Asian):

```sh
$ SaveSamples.pl results/TGP_pop_scores.txt results/TGP_pop_efa_smps.txt -xcmin 1.3 -xcmax 1.66 -ycmax 1.4 -isByd
```

Other options are the same as those in `PlotGrafPopResults.pl`, e.g.,

```sh
$ SaveSamples.pl results/TGP_pop_scores.txt results/TGP_pop_minsnp.txt -spf data/TGP_SbjPop.txt -minsnp 99980
$ SaveSamples.pl results/TGP_pop_scores.txt results/TGP_pop_popcut.txt -spf data/TGP_SbjPop.txt -ecut 85 -ohcut 15 -fhcut 30
```

### Running `ExtractAncSnpsFromVcfGz.pl` to extract genotypes of ancestry SNPs from one or more VCF files 
When the VCF file contains many more SNPs than those used by GrafPop, one can use `ExtractAncSnpsFromVcfGz.pl` to extract genotypes of ancestry SNPs from the file and then pass the output file to `grafpop`. When the genotype data is saved in multiple VCF files, `ExtractAncSnpsFromVcfGz.pl` can also be used to extract genotypes and combine them into one VCF file, so that the file can be passed to `grafpop`. 

`ExtractAncSnpsFromVcfGz.pl` takes two required parameters. The first parameter specifies the input VCF file, and second one is the output VCF file. 

```sh
$ ExtractAncSnpsFromVcfGz.pl
Usage: ExtractAncSnpsFromVcfGz.pl <vcf_or_vcf_gz_file> <output_vcf_file> [keyword]
```

For example, one can run the following commands to extract genotypes of ancestry SNPs and save the data to the output file, then pass the file to `grafpop` for ancestry inference:
```sh
$ ExtractAncSnpsFromVcfGz.pl data/TG_2smps_chr4.vcf results/TG_2smps_chr4_anc.vcf
$ grafpop results/TG_2smps_chr4_anc.vcf results/TG_2smps_chr4_anc_pops.txt
```

If the genotypes of the same set of samples are saved in a set of multiple VCF files, e.g., one file for one chromosome, `ExtractAncSnpsFromVcfGz.pl` can find all these files and extract genotypes of ancestry SNPs and save the results into one VCF file, so that it can be used by `grafpop`. If all file names differ from one another only by an embedded integer, one can run `ExtractAncSnpsFromVcfGz.pl` with the "keyword" before the integer as the optional third parameter to let the script extract genotypes from all these VCF files, e.g., 
```sh
$ ExtractAncSnpsFromVcfGz.pl data/TG_2smps_chr4.vcf results/TG_2smps_chr4_anc.vcf chr
$ grafpop results/TG_2smps_chr4_anc.vcf results/TG_2smps_chr4_anc_pops.txt
```
The first parameter can be any one of these VCF files. The script also works with zipped VCF files, e.g.,
```sh
$ ExtractAncSnpsFromVcfGz.pl data/TG_2_zip_chr17.vcf.gz results/TG_2_zip_anc.vcf chr
```
The above command extracts genotypes of ancestry SNPs from the four files with name like "TG_2_zip_chr<`number>`.vcf.gz" and save the results to the output file "TG_2_zip_anc.vcf". 

## References

Jin Y, Schäffer AA, Sherry ST, and Feolo M (2017). [Quickly identifying identical and closely related subjects in large databases using genotype data.](https://www.ncbi.nlm.nih.gov/pubmed/?term=28609482) PLoS One. 12(6):e0179106.

Jin Y, Schäffer AA, Feolo M, Holmes JB and Kattman BL (2019). [GRAF-pop: A Fast Distance-based Method to Infer Subject Ancestry from Multiple Genotype Datasets without Principal Components Analysis.](https://www.g3journal.org/content/9/8/2447.long) G3: Genes | Genomes | Genetics. August 1, 2019  vol. 9  no. 8  2447-2461.

