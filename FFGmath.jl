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

function get_update_vectors(state::GutzwillerState,move::SwapNeighborMove)
    """ Gets vectors which perform rank 1 row update to the slater determinant
    specifically, u is zero except at column i and v_i = ϕ_i(r') - ϕ_i(r)
    u1 and v1 are for up spin, and u2 and v2 are for down spin"""
    # get electron indices from business directory
    bd = state.business_directory

    site_A = move.sites[1]
    site_B = move.sites[2]

    electron_A = bd[site_A]
    electron_B = bd[site_B]

    u1 = zeros(length(state.r_up))
    v1 = zeros(length(state.r_up))
    u2 = zeros(length(state.r_down))
    v2 = zeros(length(state.r_down))

    if electron_A < electron_B
        # then A is up (1) and B is down (2)
        # up electron is moving from site A to site B
        electron_A = get_electron_index(electron_A,state)
        u1[electron_A] = 1
        v1 = state.wavefunctions[site_B,:] - state.wavefunctions[site_A,:]

        # down electron is moving from site B to site A
        electron_B = get_electron_index(electron_B,state)
        u2[electron_B] = 1
        v2 = state.wavefunctions[site_A,:] - state.wavefunctions[site_B,:]

    else
        # A is down (2) and B is up (1)
        electron_A = get_electron_index(electron_A,state)
        u2[electron_A] = 1
        v2 = state.wavefunctions[site_B,:] - state.wavefunctions[site_A,:]

        electron_B = get_electron_index(electron_B,state)
        u1[electron_B] = 1
        v1 = state.wavefunctions[site_A,:] - state.wavefunctions[site_B,:]
    end

    return u1,v1,u2,v2
end

function det_ratio_factor(A_inv,u,v)
    """ Compute ratio of determinants using Matrix Determinant Lemma"""
    return 1+ v'*A_inv*u
end

function sherman_morrison_inverse_update(A_inv,u,v)
    """ Compute (A')^-1 where A' = A + u v' """
    return A_inv - (A_inv*u*v'*A_inv)/(1+v'*A_inv*u)
end
