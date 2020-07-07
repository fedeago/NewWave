# NewWave

A scalable R package for dimensionality reduction and batch effect removal
of single-cell RNA-seq data.

## How to install 

Clone this repository on your machine and than install it

```
install.packages("path/to/NewWave", repos = NULL, type="source")
```

## How to use

If you have your data stored in a SummarizedExperiment object with 
batch effect variable called "batch" stored in colData then:

```
newWave(data,X = "~batch")
```

A more deteiled case of use is showned in the vignette.

## Useful links

[ZINB-WaVE article](https://www.nature.com/articles/s41467-017-02554-5)

[ZINB-WaVE package](https://bioconductor.org/packages/release/bioc/html/zinbwave.html)

[SharedObject package](https://bioconductor.org/packages/release/bioc/html/SharedObject.html)