# ASCA: ANOVA-Simultaneous Component Analysis in Python
[![DOI](https://zenodo.org/badge/569291256.svg)](https://zenodo.org/badge/latestdoi/569291256)


<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
* [Getting Started](#getting-started)
* [Simple Examples](#simple-examples)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [References](#references)


<!-- ABOUT THE PROJECT -->
## About The Project
ASCA is a multivariate approach to the standard ANOVA, using simultaneous component analysis to interprete the underlying factors and interaction from a design of experiment dataset. This project implements ASCA in python to support open source development and a wider application of ASCA.


<!-- GETTING STARTED -->
## Getting Started

Install this library either from the official pypi or from this Github repository:
```
pip install ASCA
```

## Install most updated version from Github

In a environment terminal or CMD:
```bat
pip install git+https://github.com/tsyet12/ASCA
```




### Simple Example
```python

    X = [[1.0000,0.6000], 
    [3.0000,0.4000],
    [2.0000,0.7000],
    [1.0000,0.8000],
    [2.0000,0.0100],
    [2.0000,0.8000],
    [4.0000,1.0000],
    [6.0000,2.0000],
    [5.0000,0.9000],
    [5.0000,1.0000],
    [6.0000,2.0000],
    [5.0000,0.7000]]
    X=np.asarray(X)

    F = [[1,     1],
     [1,     1],
     [1,     2],
     [1,     2],
     [1,     3],
     [1,     3],
     [2,     1],
     [2,     1],
     [2,     2],
     [2,     2],
     [2,     3],
     [2,     3]]
    F=np.asarray(F)
    interactions = [[0, 1]]

    ASCA=ASCA()
    ASCA.fit(X,F,interactions)
    ASCA.plot_factors()
    ASCA.plot_interactions()

```


![Figure_1](https://user-images.githubusercontent.com/19692103/205870275-df745bee-125d-4fa4-8e2a-00fa96ce9e2c.png)
![Figure_2](https://user-images.githubusercontent.com/19692103/205870291-960146ac-02f6-4852-b3d5-71c666550259.png)
![Figure_3](https://user-images.githubusercontent.com/19692103/205872428-245e778e-c805-4dfc-b5d4-0af7c890c9f2.png)


```python
# How the variables are syntaxed
ASCA.factors
ASCA.factors_explained
ASCA.factors_loadings
ASCA.factors_projected
ASCA.factors_scores
ASCA.factors_singular

ASCA.interactions
ASCA.interaction_explained
ASCA.interaction_loadings
ASCA.interaction_projected
ASCA.interaction_scores
ASCA.interaction_singular

ASCA.residuals
ASCA.effects

```


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b testbranch/prep`)
3. Commit your Changes (`git commit -m 'Improve testbranch/prep'`)
4. Push to the Branch (`git push origin testbranch/prep`)
5. Open a Pull Request


<!-- LICENSE -->
## License

Distributed under the Open Sourced BSD-2-Clause License. See [`LICENSE`](https://github.com/tsyet12/ASCA/blob/master/LICENSE) for more information.


<!-- CONTACT -->
## Contact
Main Developer:

Sin Yong Teng sinyong.teng@ru.nl or tsyet12@gmail.com
Radboud University Nijmegen

<!-- References -->
## References
Smilde, Age K., et al. "ANOVA-simultaneous component analysis (ASCA): a new tool for analyzing designed metabolomics data." Bioinformatics 21.13 (2005): 3043-3048.

Jansen, Jeroen J., et al. "ASCA: analysis of multivariate data obtained from an experimental design." Journal of Chemometrics: A Journal of the Chemometrics Society 19.9 (2005): 469-481.


## Acknowledgements
The research contribution from S.Y. Teng is supported by the European Union's Horizon Europe Research and Innovation Program, under Marie Skłodowska-Curie Actions grant agreement no. 101064585 (MoCEGS).

## Cite this software
Sin Yong Teng. (2022). tsyet12/ASCA: ASCA ANOVA-Simultaneous Component Analysis (v1.0). Zenodo. https://doi.org/10.5281/zenodo.7404343
