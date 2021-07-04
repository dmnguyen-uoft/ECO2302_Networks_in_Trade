using LinearAlgebra
###################################################
# Required Functions to Generate all possible matrices
# Function constructing network matrix with 0 diagonal from N*(N-1) array
function to_matrix(x,N)
    matrix = zeros.(N,N)
    index = 1
    for i in 1:N
        for j in 1:N
            if i != j
                matrix[i,j] = x[index]
                index = index + 1
            end
        end
    end
    return matrix
end

# Recursive Function to generate all permutations of (0,1) for given length
function recurse(len, index, value, current_array, all_array)
    append!(current_array, value)
    # Base case
    if index == len
        push!(all_array, deepcopy(current_array))
    else
        recurse(len, index + 1, 0, current_array, all_array)
        recurse(len, index + 1, 1, current_array, all_array)
    end
    pop!(current_array)
end

# Generate list of all possible matrices
function matrix_generate(N)
    current_array = []
    all_array = []

    # Generate n*(n-1) array
    recurse(N*(N-1), 1, 0, current_array, all_array)
    recurse(N*(N-1), 1, 1, current_array, all_array)

    # Turn into matrix
    matrix_list = []
    for i in 1:length(all_array) 
        push!(matrix_list, to_matrix(all_array[i], N))
    end
    return matrix_list
end

###################################################
# Simple model of endogenous networks
@time begin     # Measure running time
# Set Parameters here
N = 2
# A = [1,2]
A = ones.(N)
α = 0.1
f = 0.1
B = 1
Iden = Diagonal(ones.(N))      # Identity matrix

# Initialize Network structure
E = zeros.(N,N)         

# Generate List of all possible Networks
network_list = matrix_generate(N)

stable_NE = ones.(length(network_list))      # stable = 1 meaning E is stable NE
list_stable_matrix = []

for (index,E) in enumerate(network_list)
    # Solve for vector T: A = (I-αE)T
    T = (Iden - α*E)\A
    
    # Calculate profit
    Π = zeros.(N)
    for k in 1:N
        Π[k] = B*T[k] - f*sum(E[k,:])
    end

    for i in 1:N
    # Check if the profit for firm i is maximized
        for j in 1:N
            E_hat = copy(E)
            # Calculate profit Π_hat by changing relationship E[i,j]
            if j != i
                qw = abs(E[i,j]-1)
                E_hat[i,j] = qw
                T_hat = (Iden - α*E_hat)\A
                Π_hat = B*T_hat[i] - f*sum(E_hat[i,:])
                if Π_hat > Π[i]
                    stable_NE[index] = 0
                end
            end
        end
    end
    if stable_NE[index] != 0
       push!(list_stable_matrix, E) 
    end
end

end

#######
# Show results
# This vector shows which matrix is stable
stable_NE

# Number of equilibria
sum(stable_NE)

# Equilibria
list_stable_matrix[2]

network_list
# Results: N = 5 takes 65.184888 seconds
