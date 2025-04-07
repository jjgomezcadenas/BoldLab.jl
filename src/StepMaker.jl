# StepMaker.jl
# Translated from MATLAB StepMaker.m
# Generates synthetic stepwise signals with optional noise and flat/gaussian step styles

module StepMaker

using Random, Statistics

export generate_steps

"""
    generate_steps(n_steps::Int; 
                   minstep::Float64=1.0, 
                   maxstep::Float64=4.0, 
                   flat::Bool=true, 
                   gauss_noise::Bool=true,
                   dwell::Int=100) -> Vector{Float64}

Generate a stepwise trace with `n_steps` transitions, with random or flat step heights,
optional Gaussian noise, and constant dwell time.
"""
function generate_steps(n_steps::Int; 
                        minstep::Float64=1.0, 
                        maxstep::Float64=4.0, 
                        flat::Bool=true, 
                        gauss_noise::Bool=true,
                        dwell::Int=100)
    total_length = n_steps * dwell
    trace = Float64[]
    current_level = 0.0

    for i in 1:n_steps
        step_size = flat ? minstep : (minstep + rand() * (maxstep - minstep))
        direction = rand(Bool) ? 1 : -1
        current_level += direction * step_size
        segment = fill(current_level, dwell)
        if gauss_noise
            segment .+= randn(dwell)
        end
        append!(trace, segment)
    end

    return trace
end

end # module
