# GBM_SV_anaysis

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Setting up the development environment](#setting-up-the-development-environment)
- [License](#license)
- [Issues](https://github.com/neurodata/mgcpy/issues)

# Overview

Offering code&commands for Micro-C and CUT&tag analysis
``mgcpy`` aims to be a comprehensive independence testing package including commonly used independence tests and additional functionality such as two sample independence testing and a novel random forest-based independence test. These tests are not only included to benchmark MGC but to have a convenient location for users if they would prefer to utilize those tests instead. The package utilizes a simple class structure to enhance usability while also allowing easy extension of the package for developers. The package can be installed on all major platforms (e.g. BSD, GNU/Linux, OS X, Windows)from Python Package Index (PyPI) and GitHub.


# System Requirements

## Software requirements
### OS Requirements
This package is supported for *Linux*. The package has been tested on the following systems:
+ Linux: Ubuntu 16.04


### software Dependencies


```
-bwa mem (v0.7.17;https://bio-bwa.sourceforge.net/)
-runHiC (v0.8.4-r1; https://zenodo.org/badge/doi/10.5281/zenodo)
-cooler (v0.8.6;https://cooler.readthedocs.io/en/latest/index.html) 
-EagleC (v0.1.3;https://github.com/XiaoTaoWang/EagleC)
-Neoloopfinder (v0.3.0.post4;https://github.com/XiaoTaoWang/NeoLoopFinder)
-HiCrep (v0.2.3;https://github.com/TaoYang-dev/hicrep)
-cooltools (v0.3.2;https://cooltools.readthedocs.io/en/latest/#)
-peakachu (v1.2.0;https://github.com/open2c/cooltools)
-STAR (v2.6.0c;https://github.com/alexdobin/STAR)
-RSEM (v1.3.3;https://github.com/deweylab/RSEM)
-deepTools2 (v3.5.1;https://github.com/deeptools/deepTools)
-Bowtie2 (v2.3.4.1;https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
-cBioPortal (v3.5.3;https://www.cbioportal.org/)
-Picard tools (v2.20.7;https://broadinstitute.github.io/picard/)
-SEACR (v1.3;https://github.com/FredHutch/SEACR)
-CNVkit (v0.9.9;https://github.com/etal/cnvkit)
-GEPIA2 (v7.0;http://gepia2.cancer-pku.cn/)
-DisGenet Database (v7.0;https://www.disgenet.org/)
```


## Power Curves
- Recreated Figure 2 in https://arxiv.org/abs/1609.05148, with the addition of MDMR and Fast MGC
![Power Curves](https://raw.githubusercontent.com/neurodata/mgcpy/master/power_curves_dimensions.png)

# License

This project is covered under the **Apache 2.0 License**.
