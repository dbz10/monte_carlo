""" Misc math functions used in FreeFermionGutzwiller.jl"""

function get_Determinant(r,states)
    """ Computes slater determinant |ϕ_j(r_i)|"""
    mat = states[r,:]
    return det(mat)
end

function get_Inverse_Matrix(r,states)
    """ Compute the inverse of the matrix ϕ_j(r_i) """
    mat = states[r,:]
    return inv(mat)
end

function get_update_vectors(state::Gutzwillerstate,move::SwapNeighborMove)
    """ Gets vectors u,v which perform rank 1 update to the slater determinant
    specifically, u_i = ϕ_i(r') - ϕ_i(r) and v is zero except at column i """
    # need the yellowbook here to know which electron is being moved.

end

function det_ratio_factor(A_inv,u,v)
    """ Compute ratio of determinants using Matrix Determinant Lemma"""
end
