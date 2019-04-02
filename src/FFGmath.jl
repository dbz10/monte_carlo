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

function det_ratio_factor(A_inv,u,v)
    """ Compute ratio of determinants using Matrix Determinant Lemma"""
    return 1+ v'*A_inv*u
end

function sherman_morrison_inverse_update(A_inv,u,v)
    """ Compute (A')^-1 where A' = A + u v' """
    return A_inv - (A_inv*u*v'*A_inv)/(1+v'*A_inv*u)
end
