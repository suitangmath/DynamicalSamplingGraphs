# Randomized dynamical sampling for bandlimited graph signals

This repository is for our work:
Longxiu Huang, Deanna Needell, and Sui Tang. <a href=https://arxiv.org/abs/2109.14079> Robust recovery of bandlimited graph signals via randomized dynamical sampling </a>  Information and Inference: A Journal of the IMA, 2024.

###### To display math symbols properly, one may have to install a MathJax plugin. For example, [MathJax Plugin for Github](https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima?hl=en).

## Introduction
construct_sampling.m file contains the code for generating random samples according to the three random space-time sampling regimes introduced in our study.

construct_opt_sampling_distribution.m file contains the code for generating the optimal sampling distribution for the three random space-time sampling regimes introduced in our study.

## Syntex
```
samplingInfo = construct_sampling(sysInfo,obsInfo,samplingInfo)
```
Input:
'sysInfo' includes the name of the graph example, graph structures (e.g., number of nodes, graph Laplacian operator, bandwidth), and the operator driving the dynamic graph signal.

'obsInfo' details the total evolution time, the noise variance on the samples, and the cumulative number of samples across all time points.

'samplingInfo' specifies the sampling disbution: optimal distribution or uniform distribution. By default, samplingInfo contains the uniform distributions for the three proposed sampling regimes.

Output:
`samplingInfo': add $S_j(\pi_{A,T}),S_j, P_{\Omega_j}, W_j$ (see Table 4.1 of our work) for the three proposed sampling regimes to samplingInfo

```
samplingInfo =  construct_opt_sampling_distribution(sysInfo,obsInfo,samplingInfo)
```

Input:
'sysInfo' includes the name of the graph example, graph structures (e.g., number of nodes, graph Laplacian operator, bandwidth), and the operator driving the dynamic graph signal.

'obsInfo' details the total evolution time, the noise variance on the samples, and the cumulative number of samples across all time points.

'samplingInfo' contains the uniform distributions for the three proposed sampling regimes.

Output:
samplingInfo: the optimal sampling distributions for the three proposed sampling regimes.


## Demo

Run `RandomDS_Reconstruction.m` for a test demo.
