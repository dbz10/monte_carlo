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

# level zero:
# lattice
#   - lattice needs to know: position of each site
#   - each site should know who its neighbors are
#   - thats it?
# hamiltonian
#   - or something like that, more generally boltzmann weight
# degrees of freedom
#   - ising spin, rotor, etc?

# what determines monte carlo?
# temperature, number of steps, method (ie metropolis)



dims = (4); # dimensions of the lattice
lattice = ?;
model = ?;
