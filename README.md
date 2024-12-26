# Parametrizing a ionizable lipid to Martini3
This tutorial aims to showcase the parametrization of ionizable lipids to Martini3.

Adapted from the [Parametrization of a new small molecule](https://github.com/ricalessandri/Martini3-small-molecules/blob/main/tutorials/M3tutorials--parameterizing-a-new-small-molecule.md) Martini 3 Tutorial.

Files required and worked examples for this tutorial can be downloaded here.

## Table of contents

[Introduction](#introduction)

[1) Generate atomistic reference data](#generate-atomistic-reference-data)

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

We will need atomistic reference data to extract bonded parameters for the CG model. Note that we will need all the hydrogen atoms to extract bond lengths further down this tutorial, so make sure that your atomistic structure contains all the hydrogens.

Here, we will use the LigParGen server [6] as a way to obtain an atomistic (or all-atom, AA) structure and molecule topology based on the OPLS-AA force field, but of course feel free to use your favorite atomistic force field. Other web-based services such as the automated topology builder (ATB) or CHARMM-GUI can also be used to obtain reference topologies based on other AA force fields. Another important option is to look in the literature for atomistic studies of the molecule you want to parametrize: if you are lucky, somebody might have already published a validated atomistic force field for the molecule, which you can then use to create reference atomistic simulations.

Start by feeding the SMILES string for 1-ethylnaphthalene (namely, CCc1cccc2ccccc21) to the LigParGen server, and pick the "1.14*CM1A-LBCC (Neutral molecules)" charge model (nothing special about this choice of charge model). After submitting the molecule, the server will generate input parameters for several molecular dynamics (MD) packages. Download the structure file (PDB) as well as the OPLS-AA topology in the GROMACS format (TOP) and rename them ENAP_LigParGen.pdb and ENAP_LigParGen.itp, respectively. You can now unzip the zip archive provided:

```
unzip  files-m3-newsm.zip
```

which contains a folder called ENAP-in-water that contains some template folders and useful scripts. We will assume that you will be carrying out the tutorial using this folder structure and scripts. Note that the archive contains also a folder called ENAP-worked where you will find a worked version of the tutorial (trajectories not included). This might be useful to use as reference to compare your files (e.g., to compare the ENAP_LigParGen.itp you obtained with the one you find in ENAP-worked/1_AA-reference).

We can now move to the first subfolder, 1_AA-reference, and copy over the files you just obtained from the LigParGen server:

cd     ENAP-in-water/1_AA-reference
[move here the obtained ENAP_LigParGen.pdb and ENAP_LigParGen.itp files]

Input files obtained from LigParGen may come with unknown residue names. Before launching the AA MD simulation, we will substitute the UNK residue name by ENAP. To do so, open the ENAP_LigParGen.pdb with your text editor of choice and replace the UNK entries on the 4th column of the ATOM records section. This column defines the residue name on a PDB file. Now open the ENAP_LigParGen.itp file and replace the UNK entries under the [ moleculetype ] directive and on the 4th column of the [ atoms ] directive. These define the residue name in a GROMACS topology file. (A lengthier discussion on GROMACS topology files will be given in section 4).) Alternatively, the following command - that relies on the Unix utility sed - will replace any UNK occurrence with ENAP (note the extra space after UNK which is important to keep the formatting of the pdb file!):

sed -i 's/UNK /ENAP/' ENAP_LigParGen.pdb
sed -i 's/UNK /ENAP/' ENAP_LigParGen.itp

Now launch the AA MD simulation:

bash prepare_1mol_AA_system.sh  ENAP_LigParGen.pdb  spc216.gro  SOL  3

The last command will run an energy-minimization, followed by an NPT equilibration of 250 ps, and by an MD run of 10 ns (inspect the script and the various mdp files to know more). Note that 10 ns is a rather short simulation time, selected for speeding up the current tutorial. You should rather use at least 50 ns, or an even longer running time in case of more complex molecules (you can try to experiment with the simulation time yourself!). In this case, the solvent used is water; however, the script can be adapted to run with any other solvent, provided that you input also an equilibrated solvent box. You should choose a solvent that represents the environment where the molecule will spend most of its time.
