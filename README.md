# SMR

The SMR software tool was originally developed to implement the SMR & HEIDI methods to test for pleiotropic association between the expression level of a gene and a complex trait of interest using summary-level data from GWAS and expression quantitative trait loci (eQTL) studies ([Zhu et al. 2016 Nature Genetics](https://www.nature.com/articles/ng.3538/)). The SMR & HEIDI methodology can be interpreted as an analysis to test if the effect size of a SNP on the phenotype is mediated by gene expression. This tool can therefore be used to prioritize genes underlying GWAS hits for follow-up functional studies. The methods are applicable to all kinds of molecular QTL (xQTL) data, including DNA methylation QTL (mQTL) and protein abundance QTL (pQTL).


## Requirement  

* Eigen

* libz  

## Installation

### Install requirments  
Install Eigen and libz. If you want build smr staticlly, you need static library of those them. And the static version of omp is needed too.

### Build

Simply, in smr directoy,

```
make
```

If you want compile static version, type:

```
make smr_static
```

If your Eigen or libz is not installed in standard head file and library searching path, you can
specific them as following:

```
make -DEIGEN_PATH="where the Eigen's head file located" -DZLIB_INCLUDE="path of zlib head file" -DZLIB_LIB="path of zlib library"
```

To build smr by debug mode, using:

```
make -DDEBUG=ON
```

## Usage  
visit https://yanglab.westlake.edu.cn/software/smr/ for software's document.
