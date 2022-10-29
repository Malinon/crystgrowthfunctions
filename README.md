# Introduction

This is Python (SageMath) library dedicated to growth functions of periodic tessellations.
The package finds growth functions and visualize how crystallographic growth functions depend on frame's anchor point.
All results are exact. Numerical errors are eliminated by using symbolic computations.

# Implemented features

1. Finding topological growth functions
2. Finding crystallographic growth functions
    * For frames defined by translation vectors
    * For frames with scaled sides
3. Finding and visualizing regions within which growth functions are identical

# Examples

You can find example usages of the library in example/examples.ipynb.
You can use some of the functionality of the library in terminal using located in the same folder script.
```console
sage ./examples/example_app
```

# Installation
Go to main directory of the repository and type in terminal following commands

```console
python3 -m build
sage --pip install ./dist/crystgrowthpoly-1.0.0.tar.gz
```

# Input file format
The library read tessellation's description from file

Input File Format:
    The first line shows the number of cells in each dimension sorted in ascending order with respect to the cell dimension.
    The following lines encode the coordinates of the tessellation vertices.
    Coordinates can be given in crystallographic or Cartesian coordinates.
    When describing vertices, the SageMath syntax can be used.
    Higher dimensional cells are described after the vertices.
    Each line encodes one cell.
    Higher dimensional cells are encoded with vertex indices.
    The last two lines contain information about translation vectors

Example input file:
    6 6 1
        0     0
        1/3  -1/3
        1     0
        4/3   2/3
        1     1
        1/3   2/3
    1 2
    2 3
    3 4
    4 5
    5 6
    6 1
    1 2 3 4 5 6
    1 0
    -1/2 sqrt(3)/2

# Contact

If you have any questions regarding the package, comments or found error, please contact me.   
E-mail: jakmal7@st.amu.edu.pl


