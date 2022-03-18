# NewWave

A scalable R package for dimensionality reduction and batch effect removal
of single-cell RNA-seq data.

## How to install 

NewWave is available on Bioconductor!

```
BiocManager::install("NewWave")
```

## How to use

If you have your data stored in a SummarizedExperiment object with 
batch effect variable called "batch" stored in colData then:

```
newWave(data,X = "~batch")
```

In the X matrix you can teoretically store any variable related to the cell, both quantitative and qualitative. If you a batch variable, it must be a factor.

A complete case of use is shown in the vignette.

In the link to the workshop you can find a more detailed explanation plus an 
example of the usage with DelayedArray and HDF5 files.

## Useful links

[Bioc2021 workshop](https://fedeago.github.io/SurfingNewWave/articles/vignette.html)

[ZINB-WaVE article](https://www.nature.com/articles/s41467-017-02554-5)

[ZINB-WaVE package](https://bioconductor.org/packages/release/bioc/html/zinbwave.html)

[SharedObject package](https://github.com/Jiefei-Wang/SharedObject)
