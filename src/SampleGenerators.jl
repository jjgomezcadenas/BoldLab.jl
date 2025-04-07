using DataFrames

"""Sample_3D: theta and phi are randomized."""
function Sample_3D(N, dimensions)
    # dimensions is assumed to be a tuple or vector like (dim_x, dim_y)
    x = rand(N) .* dimensions[1]
    y = rand(N) .* dimensions[2]
    theta = acos.(1 .- 2 .* rand(N))
    phi = 2 * π .* rand(N)
    return DataFrame(x = x, y = y, theta = theta, phi = phi)
end

"""Sample_2D: theta is given and constant, phi is randomized. """
function Sample_2D(N, dimensions, theta)
    x = rand(N) .* dimensions[1]
    y = rand(N) .* dimensions[2]
    theta_vec = fill(theta, N)
    phi = 2 * π .* rand(N)
    return DataFrame(x = x, y = y, theta = theta_vec, phi = phi)
end

"""Sample_1D: both theta and phi are given and constant."""
function Sample_1D(N, dimensions, theta, phi)
    x = rand(N) .* dimensions[1]
    y = rand(N) .* dimensions[2]
    theta_vec = fill(theta, N)
    phi_vec = fill(phi, N)
    return DataFrame(x = x, y = y, theta = theta_vec, phi = phi_vec)
end