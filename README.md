# CosSq

<img src ="https://github.com/chen3262/CosSq/blob/master/pic.png" width="500">

This is a C++ program to analyze orientations of water molecules from multiple-frame [GROMACS](http://www.gromacs.org/) structure files (*.gro). Fined-grained [OpenMp](http://www.openmp.org) is implented to fulfill parallel computing in multiple-processors computers.

## Requirements
[GROMACS](http://www.gromacs.org/)

## Preparing multi-frame .gro file

Use the **trjconv** tool installed with GROMACS to convert a trajectory file (**.trr**) into a multi-frame **.gro** file. Type

```bash
gmx trjconv -f *.trr -s &.tpr -n *.ndx -pbc while -o sample.gro
```

A sample output **"sample.gro"** is provided for testing **CosSq**.

## Running CosSq

We use **make** to compile and execuate **CosSq** in once. 

```bash
make Makefile_CosSq_OpenMP
```

After the job finish, three output files will be generated:

```bash
sample.cos
sample.dip
sample.log
```

The **sample.cos** file prints out distant-dependent <cos^n> of water (n=1,2,3,4). The **sample.dip** file gives the z-component of polarization density of water in along z-direction.

## License

Copyright (C) 2017 Si-Han Chen chen.3262@osu.edu
