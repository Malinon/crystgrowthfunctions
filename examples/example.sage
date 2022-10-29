import crystgrowthpoly as cryst
import sys

WELCOME_MESSAGE = """ Welcome in crystgrowthpoly example application!
For and format of input file please run this application with argument "--help"
"""
HELP_MESSAGE = """ This application presents basic features of crystgrowthpoly library dedicated to growth functions of periodic tessellation.

Arguments:
--help                Print this message
--dim <dimmension>    Dimmension of tessellation. By default dim=2
--no_diagrams         Do not create orphic diagrams
--no_description      Textual description of orphic diagrams is not printed
--one_variable        One variable polynomials are computed. The frame grows identically in each direction
--cartesian           Use Cartesian coordinates to describe vertices
--full_plot           Create 2 diagrams. First for 2-d domains and second for 1/0-domains 

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
    sqrt(3)*2 0
    -sqrt(3) 3
"""

dim = 2
description = True
symmetric_growth = False
create_diagrams = True
crystallographic_coordinates = True
full_plot = False


def handle_argument(index):
    if sys.argv[index] == "--help":
        print(HELP_MESSAGE)
    elif sys.argv[index] == "--dim":
        global dim
        print("aaa")
        dim = int(sys.argv[index + 1])
        return index + 2
    elif sys.argv[index] == "--no_diagrams":
        global create_diagrams
        create_diagrams = False
    elif sys.argv[index] == "--one_variable":
        global symmetric_growth
        symmetric_growth = True
    elif sys.argv[index] == "--cartesian":
        global crystallographic_coordinates
        crystallographic_coordinates = False
    elif sys.argv[index] == "--full_plot":
        global full_plot
        full_plot = True
    return index + 1

index = 1
while index < len(sys.argv):
    index = handle_argument(index)

print("Input file:")
file_path = input()
while True:
    print("Does the file contain translation vectors y/N")
    answer = input()
    cartesian_vectors_included = (answer == "y")
    if answer == "N" or answer == "y":
        break
# Read tessellation
print("dim", dim)
tes = cryst.read_tessellation_from_file(file_path=file_path, cartesian_vectors_included=cartesian_vectors_included,
                                        dim=dim, crystallographic_coordinates=crystallographic_coordinates)
# Present topological growth functions
top_funcs = tes.get_growth_polynomials(symmetric_growth)
print("Topological growth functions")
for f in top_funcs:
    f.show()

if create_diagrams:
    print("")
    print("Generating orphic diagrams")
    plots = tes.plot_domains(symmetric_growth = symmetric_growth, full_plot=full_plot, description = description,
                             export_format="png")