############################################
#
# High level driver code for Monte carlo calculations
#
# Daniel Ben-Zion
# 2019
#
############################################


#goal:
# what determines a model at the highest level?
# its VMC so theres no hamiltonian or temperature

# level zero:
# lattice
#   - lattice needs to know: position of each site
#   - each site should know who its neighbors are
#   - thats it?
# degrees of freedom
#
# wavefunction
#


# what determines monte carlo?
# number of steps, method (ie metropolis)



dims = (4); # dimensions of the lattice, how many unit cells in each direction
lattice = ?;
model = ?;
