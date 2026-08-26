# NetBiF
Network Bifurcation Finder: MATLAB Code to implement Structural Bifurcation Analysis

NetBiF is a MATLAB-based computational tool to predict and classify steady-state bifurcation in biochemical systems. The software uses the method of **Structural Bifurcation Analysis** to structurally identify bifurcation conditions.

## Overview
NetBiF takes a biochemical reaction network as input and returns a list of bifurcating species and bifurcation parameters, along with the minimal labelled buffering structures localizing the bifurcation behavior.

The main function is ``bifurcation_finder(model).``

## Specific Features
In order to generate the bifurcation summary, NetBiF proceeds to both numerical and symbolic techniques. Currently, NetBiF provides the following:

- Model inspection
  - Generates a list of all species and reactions (including reverse reactions) in the network.
- Stoichiometric and kinetic matrices construction (uses methods in [COMPILES](https://github.com/pvnlubenia/COMPILES/))
  - Constructs the stoichiometric matrix, reactant and product complex matrices, and kinetic-order matrix.
- Flux-augmented matrix construction
  - Builds the flux-augmented matrix (A) from the network's stoichiometry.
  - Computes kernel and cokernel matrices of A
- Labelled buffering structure (LBS) enumeration
  - Generates candidate subnetworks for LBS
  - Filters candidates using output-completeness and index criteria
- Flux-augmented matrix symbolic representation
  - Converts nonzero entries of the flux-augmented matrix into symbolic variables for analysis
- Determinant structure analysis
  - Organizes LBS into nested determinant structures for bifurcation analysis
- Bifurcation candidate detection
  - Computes and factors symbolic determinants of determinant structures
  - Identifies factors that may change sign under assumption of positive rate parameters
  - Reports corresponding bifurcating species, bifurcation parameters, and minimal LBS containing them.

NetBiF conveniently provides these intermediate results, along with a final bifurcation summary directly in MATLAB command window.

## Input Model Format
The input for NetBiF uses the syntax from [COMPILES](https://github.com/pvnlubenia/COMPILES/). To add reactions to the network, NetBiF uses the addReaction function whose output is ``model.`` We refer the reader to the documentation of [COMPILES](https://github.com/pvnlubenia/COMPILES/).
