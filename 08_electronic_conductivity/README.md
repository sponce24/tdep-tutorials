Carrier conductivity including anharmonic effects with TDEP and EPW
===

This tutorials cover how to compute the drift and Hall mobility of c-BN using the Boltzmann transport equation (BTE) within the [EPW code](https://epw-code.org) and including anharmonic effects computed with sTDEP.
The theory for the BTE can be found in the review [[Ponce2020]](#suggested-reading).

We will show how to compute the harmonic phonon dispersion of c-BN using density functional perturbation theory (DFPT) with [Quantum ESPRESSO](https://www.quantum-espresso.org/), the anharmonic temperature dependent phonon dispersion using stochastic TDEP and to convert the effective harmonic interatomic force constants (IFC) into real-space Quantum ESPRESSO XML format.
We will then use these to compute transport properties with the [EPW code](https://epw-code.org).

For the description of inputs and more information please follow the link:
- [TDEP inputs](https://tdep-developers.github.io/tdep/files/)
- [Quantum ESPRESSO inputs](https://epwdoc.gitlab.io/source/doc/Inputs.html)
- [EPW inputs](https://epwdoc.gitlab.io/source/doc/Inputs.html)

The TDEP to Quantum ESPRESSO interface is described in [[Yin2025]](#suggested-reading).

## Theory

In this example we are going to calculate the drift and Hall hole carrier mobility of c-BN.
The drift mobility is obtained with:

$$
  \mu^{\mathrm{d}}_{\alpha\beta} = \frac{-1}{V^{\mathrm{uc}}n^{\text{c}}} \sum_n \! \int \! \frac{\mathrm{d}^3 k}{\Omega^{\mathrm{BZ}}} \, v_{n\mathbf{k}\alpha} \partial_{E_{\beta}} f_{n\mathbf{k}}
$$

where the out of equilibrium occupations are obtained by solving the BTE:

$$
  \partial_{E_{\beta}} f_{n\mathbf{k}} = e v_{n\mathbf{k}\beta} \frac{\partial f^0_{n\mathbf{k}}}{\partial \varepsilon_{n\mathbf{k}}} \tau_{n\mathbf{k}} + \frac{2\pi \tau_{n\mathbf{k}}}{\hbar}
  \sum_{m\nu} \!\int\! \frac{\mathrm{d}^3 q}{\Omega_{\mathrm{BZ}}} | g_{mn\nu}(\mathbf{k},\mathbf{q})|^2 \\
 \!\!\! \times \Big[(n_{\mathbf{q}\nu}+1-f_{n\mathbf{k}}^0)\delta(\varepsilon_{n\mathbf{k}}-\varepsilon_{m\mathbf{k}+\mathbf{q}}  + \hbar \omega_{\mathbf{q}\nu} )
 +  (n_{\mathbf{q} \nu}+f_{n\mathbf{k}}^0)\delta(\varepsilon_{n\mathbf{k}}-\varepsilon_{m\mathbf{k}+\mathbf{q}}  - \hbar \omega_{\mathbf{q}\nu} ) \Big]  \partial_{E_{\beta}} f_{m\mathbf{k}+\mathbf{q}}.
$$

The scattering rate is defined as:

$$
\tau_{n\mathbf{k}}^{-1} \equiv  \frac{2\pi}{\hbar} \sum_{m\nu} \!\int\! \frac{d^3q}{\Omega_{\mathrm{BZ}}} | g_{mn\nu}(\mathbf{k,q})|^2 \big[ (n_{\mathbf{q}\nu} +1 - f_{m\mathbf{k+q}}^0 ) \\
 \times  \delta( \varepsilon_{n\mathbf{k}}-\varepsilon_{m\mathbf{k}+\mathbf{q}} - \hbar \omega_{\mathbf{q}\nu}) +  (n_{\mathbf{q}\nu} +   f_{m\mathbf{k+q}}^0 )\delta( \varepsilon_{n\mathbf{k}}-\varepsilon_{m\mathbf{k}+\mathbf{q}} +  \hbar \omega_{\mathbf{q}\nu}) \big].
$$

A common approximation to Eq.~\eqref{eq:eqBTE} is called the self-energy relaxation time approximation (SERTA) and consists in neglecting the second term in the right-hand of the equation which gives:
$$
  \mu_{\alpha\beta}^{\mathrm{SERTA}} = \frac{-e}{V^{\mathrm{uc}}n^{\text{c}}} \sum_n \! \int \! \frac{\mathrm{d}^3 k}{\Omega^{\mathrm{BZ}}} \frac{\partial f^0_{n\mathbf{k}}}{\partial \varepsilon_{n\mathbf{k}}}  v_{n\mathbf{k}\alpha} v_{n\mathbf{k}\beta} \tau_{n\mathbf{k}}.
$$

The the low-field phonon-limited carrier mobility in the presence of a small finite magnetic field \textbf{B} is given by:

$$
  \mu_{\alpha\beta}(B_\gamma) = \frac{-1}{V^{\mathrm{uc}}n^{\text{c}}} \sum_n \! \int \! \frac{\mathrm{d}^3 k}{\Omega^{\mathrm{BZ}}} \, v_{n\mathbf{k}\alpha} [\partial_{E_{\beta}} f_{n\mathbf{k}}(B_\gamma) - \partial_{E_{\beta}} f_{n\mathbf{k}}],
$$

again solving the BTE with finite (small) magnetic field:

$$
 \Big[ 1 - \frac{e}{\hbar}\tau_{n\mathbf{k}} ({\bf v}_{n\mathbf{k}} \times {\bf B}) \cdot \nabla_{\bf k}
\Big] \partial_{E_{\beta}} f_{n\mathbf{k}}(B_\gamma) = e v_{n\mathbf{k}\beta} \frac{\partial f^0_{n\mathbf{k}}}{\partial \varepsilon_{n\mathbf{k}}} \tau_{n\mathbf{k}} + \frac{2\pi \tau_{n\mathbf{k}}}{\hbar}
  \sum_{m\nu} \!\int\! \frac{\mathrm{d}^3 q}{\Omega^{\mathrm{BZ}}} | g_{mn\nu}(\mathbf{k},\mathbf{q})|^2 \\
  \!\!\!\! \times \Big[(n_{\mathbf{q}\nu}\!+\!1\!-\!f_{n\mathbf{k}}^0)\delta(\varepsilon_{n\mathbf{k}}\!-\!\varepsilon_{m\mathbf{k}\!+\!\mathbf{q}}  \!+\! \hbar \omega_{\mathbf{q}\nu} )
 +  (n_{\mathbf{q} \nu}+f_{n\mathbf{k}}^0)\delta(\varepsilon_{n\mathbf{k}}-\varepsilon_{m\mathbf{k}+\mathbf{q}}  - \hbar \omega_{\mathbf{q}\nu} ) \Big] \partial_{E_{\beta}} f_{m\mathbf{k}+\mathbf{q}}(B_\gamma).
$$

The Hall factor and Hall mobility are then obtained as:

$$
  r_{\alpha\beta}(\hat{\mathbf{B}}) &\equiv  \lim_{\mathbf{B} \rightarrow 0} \sum_{\delta\epsilon} \frac{[\mu_{\alpha\delta}^{\rm d}]^{-1} \, \mu_{\delta\epsilon}(\mathbf{B}) \, [\mu_{\epsilon\beta}^{\rm d}]^{-1}}{|\mathbf{B}|} \\
  \mu_{\alpha\beta}^{\rm Hall}(\hat{\mathbf{B}}) &= \sum_{\gamma} \mu_{\alpha\gamma}^{\rm d} r_{\gamma\beta}(\hat{\mathbf{B}}),
$$

where $\hat{\mathbf{B}}$ is the direction of the magnetic field.

## Preliminary calculations with Quantum ESPRESSO

For this tutorial, you will need to compile [Quantum ESPRESSO](https://www.quantum-espresso.org/) and we will assume that the following executables are in your path: `pw.x, ph.x, epw.x`.

First, go to the folder `example_cBN`

Then download the standard boron and nitrogen PBE pseudopotential from [PseudoDojo](https://www.pseudo-dojo.org/) in upf format and renamed them `B-PBE.upf` and `N-PBE.upf`
   ```bash
   cd 1_qe/
   mpirun -np 4 pw.x -in scf.in | tee scf.out
   ```


## Suggested reading
- [[Ponce2020] S. Poncé, W. Li, S. Reichardt and F. Giustino, Reports on Progress in Physics **83**, 036501 (2020)](https://iopscience.iop.org/article/10.1088/1361-6633/ab6a43)
- [[Yin2020]  J. Yin, O. Hellman, and S. Poncé, arXiv:2505.20092 (2025)](https://doi.org/10.48550/arXiv.2505.20092)
