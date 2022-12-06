# ASCA: ANOVA-Simultaneous Component Analysis in Python

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
* [Getting Started](#getting-started)
* [Usage Examples](#usage-examples)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [References](#references)


<!-- ABOUT THE PROJECT -->
## About The Project



<!-- GETTING STARTED -->
## Getting Started

Install this library either from the official pypi or from this Github repository:


## Install most updated version from Github

In a environment terminal or CMD:
```bat
pip install git+https://github.com/tsyet12/ASCA
```




### Simple Example
'''python
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

'''


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

Distributed under the Open Sourced BSD-2-Clause License. See [`LICENSE`](https://github.com/tsyet12/Chemsy/blob/main/LICENSE) for more information.


<!-- CONTACT -->
## Contact
Main Developer:

Sin Yong Teng sinyong.teng@ru.nl or tsyet12@gmail.com
Radboud University Nijmegen

<!-- References -->
## References
