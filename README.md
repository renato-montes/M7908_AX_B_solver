Linear Solver
=============

A solver of the linear AX = B matrix equation, given as an assignment for
BCIT's COMP 7908 (Linear Algebra).

The program takes the size of a square matrix as an argument, and the input
is a .csv with the values, where every row's last element is part of the
"B" constant (so a 3x3 would have an input .csv of 4x3).

On a Windows 7 machine with 4GB RAM and a 2.6 GHz CPU, this single-threaded
program solves the sample 100x100 matrix equation in roughly 4 seconds. It uses
the algorithm taught in class, which is not the most computationally efficient
either.

Usage: 
```
./solver.exe 3 < three.csv
./solver.exe 100 < hundred.csv
```
