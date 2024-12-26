# Parametrizing a ionizable lipid to Martini3
This tutorial aims to showcase the parametrization of ionizable lipids to Martini3.

Adapted from the [Parametrization of a new small molecule](https://github.com/ricalessandri/Martini3-small-molecules/blob/main/tutorials/M3tutorials--parameterizing-a-new-small-molecule.md) Martini 3 Tutorial.

Files required and worked examples for this tutorial can be downloaded here.

## Table of contents

Introduction

[1) Generate atomistic reference data](#introduction)

2) Atom-to-bead mapping

3) Generate the CG mapped trajectory from the atomistic simulation

4) Create the initial CG itp and tpr files

5) Generate target CG distributions from the CG mapped trajectory

6) Create the CG simulation

7) Optimize CG bonded parameters

8) Comparison to experimental results, further refinements, and final considerations

References and notes

## Introduction

In this tutorial we will learn how to build a Martini 3 topology for an ionizable lipid. The aim is to have a pragmatic description of the Martini 3 coarse-graining (CGing) principles [1,2], which follow the main ideas outlined in the seminal Martini 2 work [3].

We will use as an example the molecule SM-102 (Fig. 1), and make use of Gromacs versions 2024.x or later. Required files and worked examples can be downloaded here.

<p align="center">
<img src="SM-102.png" width="500" alt="SM-102">
</p>

## Generate atomistic reference data

