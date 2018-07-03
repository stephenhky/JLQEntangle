# define the bipartite state
tensor = Array(Complex128, 2, 2)
tensor[1, 2] = sqrt(0.6)im
tensor[2, 1] = sqrt(0.4)
tensor[1, 1] = tensor[2, 2] = 0.0

# function for finding reduced density matrix
function reduced_density_matrix(bipartite_tensor, keep)
    tensor_shape = size(tensor)
    reddenmat = Array(Complex128, tensor_shape[keep], tensor_shape[keep])
    for i = 1:tensor_shape[keep]
        for ip = 1:tensor_shape[keep]
            reddenmat[i, ip] = 0
            if keep==1
                for j = 1:tensor_shape[2]
                    reddenmat[i, ip] += bipartite_tensor[i, j]*conj(bipartite_tensor[ip, j])
                end
            elseif keep==2
                for j = 1:size(tensor)[1]
                    reddenmat[i, ip] += bipartite_tensor[j, i]*conj(bipartite_tensor[j, ip])
                end
            else
                reddenmat[i, ip] /= 0
            end
        end
    end

    return reddenmat
end

# Schmidt decomposition
function schmidt_decomposition(bipartite_tensor)
    tensor_shape = size(bipartite_tensor)
    mindim = minimum(tensor_shape)

    rho1 = reduced_density_matrix(bipartite_tensor, 1)
    eigvals, eigvecs = eig(rho1)

    coefmat0 = bipartite_tensor * eigvecs
    modes = [Dict("weights"=>eigvals[i], "modeA"=>coefmat0[:, i]/norm(coefmat0[:, i]), "modeB"=>eigvecs[:, i])
        for i in 1:mindim]
    return modes
end

# testing
reduced_density_matrix(tensor, 1)
reduced_density_matrix(tensor, 1)
