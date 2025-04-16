# Tutorial for the Thermal Conductivity
 
In this tutorial, you will learn how to compute lattice thermal conductivity in solids using the TDEP method. 
Useful references on the topic:

[R.E. Peierls, Quantum Theory of Solids (1955)](https://books.google.es/books?id=WvPcBUsSJBAC&redir_esc=y)

[M. Omini and A. Sparavigna Phys. Rev. B 53, 9064 (1996)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.53.9064)

[D. A. Broido et al., Appl. Phys. Lett. 91, 231922 (2007)](https://pubs.aip.org/aip/apl/article-abstract/91/23/231922/334217/Intrinsic-lattice-thermal-conductivity-of?redirectedFrom=fulltext)

[O. Hellman and D.A. Broido, Phys. Rev. B 90, 134309 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.134309)

[A. Castellano et al., Phys. Rev. B 111, 094306 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.094306)

[J. Klarbring et al., Phys Rev Lett 125, 045701 (2020)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.045701)


# General scope
In this tutorial, you will learn how  to calculate the lattice thermal conductivity in the mode-coupling formalism, including collective and off-diagonal contributions up to fourth-order interactions.


The tutorial covers:
* the basic features included in the thermal conductivity routine of TDEP
* thermal conductivity as a function of temperature
* thermal conductivity for a given isotope distribution
* thermal conductivity vs q-point grid
* plot and analysis of the results
* thermal conductivity including the four-phonon contributions to the scattering

The tutorial **does not cover**:
* Relax the structure
* Supercell convergence
* Extract force constants

  
There are three test cases: 


1. [Aluminium](#al)
2. [MgO](#mgo)
3. [Graphene](#graphene)




The data provided includes the IFCs and the unitcell obtained from previous calculations. **This tutorial is not meant to produce converged results.**

# Before running the tutorial:
* TDEP installed
* [H5PY](https://docs.h5py.org/en/stable/) for Python installed
* Access to a plotting tool (for example  [matplotlib](https://matplotlib.org/), if you are using Python)
* Have a converged set of IFCs (you can use some previous examples or the data provided in this tutorial)

## Input files:
* [infile.ucposcar](https://tdep-developers.github.io/tdep/files/#infile.ucposcar)
* [infile.forceconstant](https://tdep-developers.github.io/tdep/program/extract_forceconstants/#outfileforceconstant)
* [infile.forceconstant_thirdorder](https://tdep-developers.github.io/tdep/program/extract_forceconstants/#outfileforceconstant_thirdorder)

## Optional input file:

* [infile.isotopes](https://tdep-developers.github.io/tdep/files/#infile.isotopes)
* [infile.forceconstant_fourthorder](https://tdep-developers.github.io/tdep/program/extract_forceconstants/#infile.forceconstant_fourthorder)


# Basic steps

## Preparation

Read the documentation for the [thermal conductivity](https://tdep-developers.github.io/tdep/program/thermal_conductivity/)

Go now into your work directory and copy the files provided.


#### Note: The tutorial files can be downloaded from the school webpage [tdep_school](https://github.com/tdep-developers/tdep-tutorials/tree/main)

 
Inspect the content of the folder: 

You can see some examples. 

# Al


Go in Examples/Al

It contains the minimum input files needed for the thermal conductivity
* [infile.ucposcar](https://tdep-developers.github.io/tdep/files/#infile.ucposcar)
* [infile.forceconstant](https://tdep-developers.github.io/tdep/program/extract_forceconstants/#outfileforceconstant)
* [infile.forceconstant_thirdorder](https://tdep-developers.github.io/tdep/program/extract_forceconstants/#outfileforceconstant_thirdorder)

## Compute the thermal conductivity

If you run it using:

```
mpirun thermal_conductivity > kappa.log
```

you will be able to get *output _thermal_conductivity* which contains the thermal conductivity tensor   $\kappa_{\alpha \beta}$ , with the decomposition from all contributions,  at a given temperature. This file can be parsed with tools such as numpy.

It looks like this

```
# Unit:               W/m/K
# Temperature:          0.300000000000E+03
# Single mode approximation
#                      kxx                      kyy                      kzz                      kxy                      kxz                      kyz

# Collective contribution
#                      kxx                      kyy                      kzz                      kxy                      kxz                      kyz

# Off diagonal (coherence) contribution
#                      kxx                      kyy                      kzz                      kxy                      kxz                      kyz

# Total thermal conductivity
#                      kxx                      kyy                      kzz                      kxy                      kxz                      kyz


```

As explained in the documentation, without specifying any optional flag, you will obtain the thermal conductivity for a natural isotope distribution, with a q-mesh of 26 26 26 (default value), at 300K.

You can check the status of the calculation by looking at the file kappa.log printed
```
Initialize calculation
 ... read unitcell poscar
 ... read second order forceconstant
 ... read third order forceconstant
 ... generating q-point mesh
 ... generating harmonic properties on the q-point mesh
 ... done in        0.291 s

 Calculating scattering events
 ... walltime() used to generate random state
 ... creating Monte-Carlo grid
 ... distributing q-point/modes on MPI ranks
 ... everything is ready, starting scattering computation
  ... computing scattering amplitude     100.0% |========================================|  105.56641s
 ... symmetrizing scattering matrix
 ... done in      120.580 s

 Thermal conductivity calculation
 ... computing kappa in the single mode approximation
 ... computing off diagonal (coherence) contribution
 ... solving iteratively the collective contribution
 iter         kxx            kyy            kzz            kxy            kxz            kyz       DeltaF/F
    0        16.6454        16.6454        16.6454         0.0000         0.0000         0.0000
    1        19.2847        19.2847        19.2847         0.0000         0.0000         0.0000   1.143E+00
    2        20.8637        20.8637        20.8637         0.0000         0.0000         0.0000   3.357E+00
    3        21.4629        21.4629        21.4629         0.0000         0.0000         0.0000   5.783E-02
    4        21.7904        21.7904        21.7904         0.0000         0.0000         0.0000   5.922E-03
    5        21.9321        21.9321        21.9321         0.0000         0.0000         0.0000   6.913E-03
    6        22.0077        22.0077        22.0077         0.0000         0.0000         0.0000   6.154E-04
    7        22.0425        22.0425        22.0425         0.0000         0.0000         0.0000   1.642E-04
    8        22.0607        22.0607        22.0607         0.0000         0.0000         0.0000   2.955E-05
    9        22.0694        22.0694        22.0694         0.0000         0.0000         0.0000   1.140E-05
   10        22.0738        22.0738        22.0738         0.0000         0.0000         0.0000   2.358E-06
   11        22.0760        22.0760        22.0760         0.0000         0.0000         0.0000   7.698E-07
   12        22.0771        22.0771        22.0771         0.0000         0.0000         0.0000   1.646E-07
 ... done in        0.645 s


... symmetrizing the thermal conductivity tensors

 Decomposition of the thermal conductivity (in W/m/K)
 Single mode approximation (SMA)
              kxx            kyy            kzz            kxy            kxz            kyz
             16.6454        16.6454        16.6454         0.0000         0.0000         0.0000
 Correction to include collective contribution via iterative procedure
              kxx            kyy            kzz            kxy            kxz            kyz
              5.4317         5.4317         5.4317         0.0000         0.0000         0.0000
 Off diagonal (coherence) contribution
              kxx            kyy            kzz            kxy            kxz            kyz
              0.0236         0.0236         0.0236         0.0000         0.0000         0.0000
 Total thermal conductivity
              kxx            kyy            kzz            kxy            kxz            kyz
             22.1007        22.1007        22.1007         0.0000         0.0000         0.0000

 ... computing cumulative kappa
 ... estimating kappa with boundary
 ... computing density of state
 ... computing spectral kappa
 ... computing angular momentum

 ... dumping auxiliary data to files

Scattering rates can be found in                                outfile.thermal_conductivity_grid.hdf5
Thermal conductivity tensor can be found in                     outfile.thermal_conductivity
Cumulative and spectral thermal conductivity can be found in    outfile.cumulative_thermal_conductivity.hdf5



```
The first step (iter 0) represents the RTA solution. The second is the converged iterative solution within the given tolerance (default value 1e-5 Tolerance for the iterative solution). The third value accounts for the coherence contribution. The last one is the total thermal conductivity given by the sum of the previous ones.
For reference, take a look at the documentation about the [thermal conductivity](https://tdep-developers.github.io/tdep/program/thermal_conductivity/).


Repeat the run for several temperatures, and plot the results. Here an example of scripts you can use to automate the process:
```
#!/bin/bash

# Define the range of temperatures
temperatures=(100 200 300 400 500)  # Replace these with the actual temperatures you want

# Create the combined output file
> combined_kappa.kappa

# Row counter
row=1

for temp in "${temperatures[@]}"; do
  echo "Running thermal_conductivity at temperature $temp..."
  mpirun thermal_conductivity -qg 8 8 8 --temperature $temp

  mv outfile.thermal_conductivity "output_${temp}.kappa"
  thermal_data=$(awk '/^# Total thermal conductivity/ {getline; print; getline; print}' "output_${temp}.kappa")

  thermal_data=$(echo "$thermal_data" | sed '/^#/d' | sed 's/^[ \t]*//;s/[ \t]*$//')

  echo "$row     $temp   $thermal_data" >> combined_kappa.kappa

  ((row++))
done

echo "All runs completed. Relevant outputs with temperatures appended to combined_kappa.kappa."
```

plot the results with gnuplot
```
gnuplot
p 'combined_kappa.kappa' u 2:3 w l
```
and study the plot. 

![here you can see how the plot should look like](https://github.com/RobertaFarris93/tdep-tutorials/blob/thermal_conductivity/04_thermal_conductivity/Plots/Al_kappa_T.png)

## Understand the results

* Does the trend look reasonable?
* Are the values comparable with the literature ([thermal conductivity of pure aluminum is estimated to be in the range between 220 and 250 (W/mK)](https://thermtest.com/thermal-resources/materials-database))? Why?

Get familiar with the optional flags available for thermal conductivity.


## Study the convergence 

### q-points

* Perform the same calculation using different grids of q-points
* Plot the thermal conductivity as a function of 1/q

The calculation of thermal conductivity, being an integrated quantity, requires its evaluation under the assumption of an infinitely refined q-point grid. Unfortunately, this is impossible from a computational point of view, but using progressively finer grids, the behavior of thermal conductivity should scale linearly with q. Thus, in order to converge the thermal conductivity value, we could perform the calculation for a set of q-grids and then study the convergence by plotting the thermal conductivity against 1/q and extrapolating the value for 1/q at 0.  The point of intersection on the y-axis resulting from this regression corresponds to the thermal conductivity within the hypothetical context of an infinitely dense q-point grid. 
For more details, see [Esfarjani, K. et. al., Phys. Rev. B 84, 085204 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.085204).


![Here the convergence test for thermal conductivity of aluminum](https://github.com/RobertaFarris93/tdep-tutorials/blob/thermal_conductivity/04_thermal_conductivity/Plots/Al_convergence_q-points.png)

This is a python script you can customize and use for that purpose:
```
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import Polynomial

q = []
kappa = []

# Open and read the file line by line
with open("combined_kappa.q-points", "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 3:
            continue  # skip lines that are too short
        try:
            q_val = float(parts[1])
            kappa_val = float(parts[2])
            q.append(q_val)
            kappa.append(kappa_val)
        except ValueError:
            continue  # skip lines with non-numeric data

q = np.array(q)
kappa = np.array(kappa)
inv_q = 1.0 / q

# Linear fit
fit = Polynomial.fit(inv_q, kappa, deg=1)
kappa_inf = fit(0)

print(f"Extrapolated thermal conductivity at 1/q = 0 (infinite q): {kappa_inf:.6f} W/m·K")

# Plot
x_range = np.linspace(0, max(inv_q)*1.1, 100)
y_fit = fit(x_range)

plt.plot(inv_q, kappa, 'o', label='Data')
plt.plot(x_range, y_fit, '-', label='Linear Fit')
plt.axvline(0, color='gray', linestyle='--')
plt.axhline(kappa_inf, color='red', linestyle='--', label=f'Extrapolated: {kappa_inf:.3f}')
plt.xlabel('1 / q-points')
plt.ylabel('Thermal Conductivity (W/m·K)')
plt.title('Extrapolation to Infinite q-points')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("Al_convergence_q-points.png")
plt.show()

```

# MgO

### Thermal conductivity of MgO

Now, you should be familiar with the thermal_conductivity routine. So far, you’ve learned how to run calculations for different grids of q-points and various temperature ranges.
Let’s now explore the available output information you can extract from TDEP.

We will analyze the thermal conductivity of a well-known semiconductor.
Go to the Examples/MgO directory.

Create a new folder and copy all the input files there.
Perform the calculation using the following command:

```
mpirun thermal_conductivity -qg 12 12 12 

```

By default, TDEP uses the isotope natural distribution.(tabulated in the code, taken from the symbol in infile.ucposcar). In case you want to specify some other distribution, you can write in the same directory the [inpute.isotopes file](https://tdep-developers.github.io/tdep/files/#infileisotopes) following the example linked here. 
```
1         # number of isotopes for first atom in infile.ucposcar
1 28.0855 # concentration, mass, one line per isotope
2         # number of isotopes for second atom
0.5 12.0  # concentration, mass
0.5 13.0  # concentration, mass
...
```

Now repeat the calculation in a different folder for the case of pure MgO by using:

```
mpirun thermal_conductivity -qg 12 12 12 --noisotope

```

Compare the outputs. 
The isotope scattering is known to decrease the thermal conductivity of MgO by 30%-40% at 300K. Did you observe that? For reference, see [Florian Knoop et.al, PRB 107, 224304 (2023)](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.107.224304). 


![MgO: comparison betweeen natural isotope distribution and pure cases](https://github.com/RobertaFarris93/tdep-tutorials/blob/thermal_conductivity/04_thermal_conductivity/Plots/MgO_T.png)


#### Extrapolation for an infinite grid of q-points

Compute the thermal conductivity using different grids, as in the previous example. 
![Fit the k_xx against 1/qx and extrapolate the value for qx=0](https://github.com/RobertaFarris93/tdep-tutorials/blob/thermal_conductivity/04_thermal_conductivity/Plots/MgO_convergence_q-points.png). 

**In order to reach convergence, you may need access to a cluster/HPC.**


#### qg3ph MC grid: speed-up the calculations

You may have noticed that TDEP offers the possibility to select Density of q-point mesh for Brillouin zone integration and the dimension of the grid for the threephonon integration through the flag **--qpoint_grid3ph value#1 value#2 value#3**, **-qg3ph value#1 value#2 value#3**.
For more details, have a look at the manual: [Monte Carlo integration for the scattering rates](https://tdep-developers.github.io/tdep/program/thermal_conductivity/#monte-carlo-integration-for-the-scattering-rates)

The idea is, we generate a full grid, on which the thermal conductivity will be integrated. A subset of this full grid can then be selected to perform the scattering integration. In order to improve the convergence, these points are not selected entirely at random but using a stratified approach to sample more uniformly the Brillouin zone.

Converging the grids is an important step to ensure accurate results. Since the convergence of the Monte-Carlo grids is not related to the convergence of the full grid, their determination can be done independently. To reduce the computational cost of the convergence, an approach is to fix the full grid to a moderately large density, and converge the third-order grid. Once the Monte-Carlo grid densities are known, then the full grid density can be determined.

![What is a good value for qg3ph?](https://github.com/RobertaFarris93/tdep-tutorials/blob/thermal_conductivity/04_thermal_conductivity/Plots/MgO_convergence_qg3ph.png). 

# Post-processing options
So far we have seen how to extract the thermal conductivity tensor using TDEP routine.  Let's have a look at the other output files you may find in the working directory. In case your calculations are not finished yet, you can use the output file provided in the MgO directory, using a q-grid of 28x28x28 points and using a temperature of 300K. 

*outfile.cumulative_thermal_conductivity.hdf5*

The file contains the information on the cumulative thermal conductivity plots described in the manual, per each computed temperature. You can inspect the information contained in the file using the following minimal Python script:
```
import h5py
fn = h5py.File('outfile.cumulative_thermal_conductivity.hdf5', 'r')
print("what's inside:", list(fn.keys()))
fn.close()
```

```
['boundary_scattering_kappa', 'boundary_scattering_lengths', 'cumulative_kappa_vs_mean_free_path', 'cumulative_kappa_vs_mean_free_path_per_atom', 'cumulative_kappa_vs_mean_free_path_per_mode', 'dos', 'frequencies', 'generating_angular_momentum_tensor', 'generating_angular_momentum_tensor_vs_frequency', 'mean_free_path_axis', 'spectral_kappa_vs_frequency', 'spectral_kappa_vs_frequency_per_atom', 'spectral_kappa_vs_frequency_per_mode']
```

##### 1. Cumulative thermal conductivity vs mean free path
The cumulative thermal conductivity can then be computed as a sum of the fraction of heat that is carried by phonons with mean free paths smaller than l, where l is:
```math
l_{\lambda} = \left| v_{\lambda} \right| \tau_{\lambda} \,,
```

```math
\kappa_{\alpha\beta}^{\textrm{acc}}(l)= \frac{1}{V} \sum_{\lambda} C_{\lambda} v^{\alpha}_{\lambda} v^{\beta}_{\lambda} \tau_{\lambda} \Theta(l- l_{\lambda} ) \,,
```
##### 2.  The spectral thermal conductivity 

which is a measure which frequencies contribute most to thermal transport.
```math
\kappa_{\alpha\beta}(\omega)=\frac{1}{V} \sum_{\lambda}C_{\lambda} v^{\alpha}_{\lambda} v^{\beta}_{\lambda} \tau_{\lambda} \delta(\omega- \omega_{\lambda})
```

Note that all quantities do not include the off-diagonal coherence contribution.

#### Read and analyze the output
Let's analyze the output file with the following script

```
import numpy
import matplotlib.pyplot as plt
import h5py

params = {'legend.fontsize': 20,
          'figure.figsize': (15, 5),
         'axes.labelsize': 30,
         'axes.titlesize':30,
         'xtick.labelsize':20,
         'ytick.labelsize':20}
plt.rcParams.update(params)

# Open the HDF5 file in read mode
fn = h5py.File('outfile.cumulative_thermal_conductivity.hdf5', 'r')


#plot one component of the spectral thermal conductivity as a function of frequency
plt.plot(fn['frequencies'][:], fn['spectral_kappa_vs_frequency'][0, 0, :], lw = 4)
#add labels, titles etc.
plt.title('Spectral thermal conductivity of MgO')
plt.xlabel('Frequency [THz]')
plt.ylabel(r'$\kappa$ [W/K/m/THz]')
plt.tight_layout()
plt.savefig('Spectral_thermal_conductivity_MgO.png')
plt.show()
#plot the cumulative thermal conductivity as a function of the total mean free path
plt.semilogx(fn['mean_free_path_axis'][:], 
             fn['cumulative_kappa_vs_mean_free_path'][0, 0, :], lw = 4)
#add labels, titles etc.
plt.title('Cumulative kappa vs mean free path of MgO')
plt.xlabel('Mean Free Path [m]')
plt.ylabel(r'$\kappa$ [W/K/m]')
plt.xlim(1E-9,1E-5)
plt.tight_layout()
plt.savefig('thermal_conductivity_vs_mfp_MgO.png')
```
The script should give you the plots reported below. 
![Spectral thermal conductivity of MgO](https://github.com/RobertaFarris93/tdep-tutorials/blob/thermal_conductivity/04_thermal_conductivity/Plots/Spectral_thermal_conductivity_MgO.png)


![Cumulative kappa vs mean free path of MgO](https://github.com/RobertaFarris93/tdep-tutorials/blob/thermal_conductivity/04_thermal_conductivity/Plots/thermal_conductivity_vs_mfp_MgO.png)


# Convergence of thermal conductivity 

As for the other physical quantities, the thermal conductivity needs to be tested against all the parameters used in the calculations. In particular:

* supercell size
* number of configurations used
* range of the forces cutoffs (number of neighbors used included in the integral)
* number of iterations in the self-consistent loop 

## Supercell size convergence

you can test the size of the supercells used by following the steps explained in [Tutorial 02 Sampling](https://github.com/tdep-developers/tdep-tutorials/blob/main/02_sampling/sTDEP/01_MgO/README.md).


## Self-consistent loop

In order to converge the thermal conductivity, we should test it against the sampling in a iterative way. For doing that, we should repeat the steps explained in the  [Tutorial 02 Sampling](https://github.com/tdep-developers/tdep-tutorials/blob/main/02_sampling/sTDEP/01_MgO/README.md) and test the goodness of our fit for the desired property, in this case the thermal conductivity (sTDEP scheme). 
[Example of MgO for sampling](https://github.com/tdep-developers/tdep-tutorials/blob/main/02_sampling/sTDEP/01_MgO/README.md)

Remember to extract the 3-order IFCs, needed for calculating the thermal conductivity, at each step:
```
mpirun extract_forceconstants  -rc2 8 -rc3 4
```
 
**How many steps did you need to reach convergence?**

## Cutoffs of IFcs

In the data provided, you will find a directory named ```convergence_tests/input_MgO```.

Inspect the folder. It contains the necessary TDEP input files required for fitting the force constants and calculating related properties.

In Tutorial 01, you learned how to increase the second-order force constant cutoff (**rc2**) and observe how the phonon dispersion changes.

Now, repeat the same steps, but this time keep rc2 fixed and change rc3, the cutoff for third-order force constants.

To get started, use the following command:

```
mpirun extract_forceconstants  -rc2 8 -rc3 xx
```

Explore how changing rc3 influences the results and assess the convergence of the thermal conductivity or other properties of interest.


**How the thermal conductivity changes with the 3rd order IFCs cutoff?**

# Next steps
 

# Use your material of interest

You can now use your own structure to calculate the thermal conductivity.

Copy your primitive structure and the force constants into your working folder, and repeat the previous steps.

Alternatively, you can use the provided input files for silicon (Si), which will be used as an example in the next tutorial.

# Graphene

So far, we have focused on bulk materials and explored the capabilities of the thermal conductivity routine in **TDEP**.  
For **two-dimensional materials**, such as **graphene**, it is essential to include the **four-phonon scattering contribution** to obtain accurate thermal conductivity values.

We provide an example for graphene, including the necessary **interatomic force constants (IFCs)**.

Perform a calculation including the four-phonon scattering contribution using the `--fourthorder` flag, and compare the results to a calculation **without** this contribution:

```bash
mpirun thermal_conductivity -qg 8 8 1 --fourthorder
```

*Are the results similar?*

**Note**: Graphene is a two-dimensional material, so the thermal conductivity must be corrected to account for its effective thickness.  
 The results should be multiplied by a factor of **15.0 / 3.335**.


Below what you should have obtained, as a function of T

![Here an example of the thermal conductivity convergence.](https://github.com/RobertaFarris93/tdep-tutorials/blob/thermal_conductivity/04_thermal_conductivity/Plots/Graphene_kappa_vs_temp_4ph.png)
![Here an example of the thermal conductivity convergence.](https://github.com/RobertaFarris93/tdep-tutorials/blob/thermal_conductivity/04_thermal_conductivity/Plots/Graphene_kappa_vs_temp_3ph.png)

