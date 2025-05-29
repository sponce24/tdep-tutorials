Carrier conductivity including anharmonic effects with TDEP and EPW
===

This tutorials cover how to compute the drift and Hall mobility of c-BN using the Boltzmann transport equation (BTE) within the [EPW code](https://epw-code.org) and including anharmonic effects computed with sTDEP.
The theory for the BTE can be found in the review [[Ponce2020]](#suggested-reading).

We will show how to compute the harmonic phonon dispersion of c-BN using density functional perturbation theory (DFPT) with [Quantum ESPRESSO](https://www.quantum-espresso.org/), the anharmonic temperature dependent phonon dispersion using stochastic TDEP and to convert the effective
harmonic interatomic force constants (IFC) into real-space Quantum ESPRESSO XML format.
We will then use these to compute transport properties with the [EPW code](https://epw-code.org).

For the description of inputs and more information please follow the link:
- [TDEP inputs](https://tdep-developers.github.io/tdep/files/)
- [Quantum ESPRESSO inputs](https://epwdoc.gitlab.io/source/doc/Inputs.html)
- [EPW inputs](https://epwdoc.gitlab.io/source/doc/Inputs.html)


## Suggested reading
- [[Ponce2020] S. Poncé, W. Li, S. Reichardt and F. Giustino, Reports on Progress in Physics **83**, 036501 (2020)](https://iopscience.iop.org/article/10.1088/1361-6633/ab6a43)
