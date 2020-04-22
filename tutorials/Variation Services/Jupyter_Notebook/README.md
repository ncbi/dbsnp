# Variation services Notebooks

You can launch a live Jupyter Notebook in your web browser using `Binder` (beta) below to test and experiment.
The experimental `Binder` server may take a few minutes to launch and can be slow, run out of memory, and may not always work.

## `SPDI` notebooks

* Launch `spdi_batch.ipynb` notebook on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncbi/dbsnp/master?filepath=%2Ftutorials%2FVariation%20Services%2FJupyter_Notebook%2Fspdi_batch.ipynb)

Batch rs ID annotation and conversion HGVS to SPDI, HGVS to RS ID, and retrieve RS JSON objects using Variation Services.

* Launch `navs_spdi_demo.ipynb` notebook on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncbi/dbsnp/master?filepath=tutorials%2FVariation%20Services%2FJupyter_Notebook%2Fnavs_spdi_demo.ipynb)

Use NCBI API Variation Service (navs) module to convert variant notations (HGVS, RS, and SPDI), remap variants between sequences, and normalized variants.

## The Allele Frequency Aggregator (ALFA) project's notebooks

* Launch `by_gene.ipynb` notebook on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncbi/dbsnp/master?filepath=tutorials%2FVariation%20Services%2FJupyter_Notebook%2Fby_gene.ipynb)

Retrieve frequency data by gene using eUtils and Variation Service

* Launch `by_rsid.ipynb` notebook on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncbi/dbsnp/master?filepath=tutorials%2FVariation%20Services%2FJupyter_Notebook%2Fby_rsid.ipynb)

Retrieve frequency data by rsid in JSON format

* Launch `frequencies_for_vcf.ipynb` notebook on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncbi/dbsnp/master?filepath=tutorials%2FVariation%20Services%2FJupyter_Notebook%2Ffrequencies_for_vcf.ipynb)

Compare frequencies for two populations for variations from a VCF file

* Launch `metadata_as_hash.ipynb` notebook on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncbi/dbsnp/master?filepath=tutorials%2FVariation%20Services%2FJupyter_Notebook%2Fmetadata_as_hash.ipynb)
Retrieve and transform ALFA project and population meta data

* Launch `plot.ipynb` notebook on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncbi/dbsnp/master?filepath=tutorials%2FVariation%20Services%2FJupyter_Notebook%2Fplot.ipynb)

Plot population minor allele frequencies (MAFs) for variants in a gene location

* Launch `querying_subsets_ftp.ipynb` notebook on [Google Colab](https://colab.research.google.com/github/ncbi/dbsnp/blob/master/tutorials/Variation%20Services/Jupyter_Notebook/querying_subsets_ftp.ipynb)

Retrieving subsets of VCF data from remote NCBI FTP site

    Notes:

    * A (free) Google account is required to run a notebook in Google Colab.

    * We do not use Binder for this notebook because they have a policy of not allowing FTP connections: <https://github.com/binder-examples/getting-data#large-public-files>
