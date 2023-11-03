# A Decomposition Approach to Multi-Agent Systems with Bernoulli Packet Loss

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5997060.svg)](https://doi.org/10.5281/zenodo.5997060)

## General

This repository contains an implementation of the algorithms described in the paper

> C. Hespe, H. Saadabadi, A. Datar, H. Werner and Y. Tang, "A Decomposition Approach to Multi-Agent Systems with Bernoulli Packet Loss," in *IEEE Transactions on Control of Network Systems*, doi: [10.1109/TCNS.2023.3275917](https://doi.org/10.1109/TCNS.2023.3275917).

It may be used to recreate and validate the figures from the paper.
To do so, run either of the three main entry points in the repository, the scripts `conservatism.m`, `scaling_large.m` and `scaling_small.m`.
Be advised that each of these scripts has a runtime of at least one hour.
The raw data used in the figures in the paper is available in the subdirectory `figures`.

## Prerequisites

To run the scripts in this repository, you will need a working copy of [*Yalmip*](https://yalmip.github.io/) together with a suitable SDP solver in your *Matlab* path.

The code in this repository was tested in the following environment:

* *Windows 10* Version 20H2
* *Matlab* 2021
* *Yalmip* 16-January-2020

The *Matlab* [`parfor`](https://de.mathworks.com/help/parallel-computing/parfor.html) feature from the *Parallel Computing Toolbox* is used to speed up the calculations.
*Matlab* should automatically detect if that toolbox is not available and run the iterations sequentially in that case.
However, this will drastically prolong the runtime of the scripts to up to a day or more!
You may want to reduce the number of sampling points for the figures or run the calculations for smaller networks.
