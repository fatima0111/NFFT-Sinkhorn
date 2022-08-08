# NFFT-Sinkhorn

Acceleration of the Sinkhorn algorithm for multi-marginal optimal transport with tree and circle structured cost functions. The new algorithm is faster than the usual Sinkhorn algorithm. Numerical examples such as Euler flows and the tree-shaped Barycenter problem illustrate the accuracy of the new algorithm.

## Prerequisites
To run the code the fastsum package has be installed. For more details refer to https://www-user.tu-chemnitz.de/~potts/nfft/ and https://github.com/NFFT/nfft. 

#### Make sure that the fastsum package is in the same parent directory the NFFT-Sinkhorn directories.

## Overview

## Reference

When you are using this code, please cite the paper

[Accelerating the Sinkhorn algorithm for sparse multi-marginal optimal transport by fast Fourier transforms](https://arxiv.org/abs/2208.03120), Fatima Antarou Ba, Michael Quellmalz, 2022. 

This paper also explains the algorithms in more detail.

## Directory structure

| File/Folder   | Purpose                                                                            |
| ------------- |:----------------------------------------------------------------------------------:| 
| Libs          | Sinkhorn Algorithm from Section 4 and NFFT-Sinkhorn algorithm from Section 5.      |
| images        | Test images for Wasserstein barycenters with general tree                          |
| Output_Circle | Output of the numerical examples for MOT problem with tree-structured cost function|
| Output_Tree   | Output of the numerical examples for MOT problem with tree-structured cost function|
| Test_functions| Implementation of numerical examples in Section 6. of the paper                    |
| Utils         |                                                                                    | 


## Legal Information & Credits

Copyright (c) 2022 [Fatima Antarou Ba](https://www.math.tu-berlin.de/fachgebiete_ag_modnumdiff/angewandte_mathematik/v_menue/team/fatima_antarou_ba/v_menue/homepage/) and [Michael Quellmalz](https://page.math.tu-berlin.de/~quellm/index.php)

This software was written by Fatima Antarou Ba and Michael Quellmalz. It was developed at the Institute of Mathematics, TU Berlin. The first mentioned author acknowledges support by the German Research Foundation within the Bundesministerium für Bildung und Forshung within the Sale project.

NFFT-Sinkhorn is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. If not stated otherwise, this applies to all files contained in this package and its sub-directories.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
