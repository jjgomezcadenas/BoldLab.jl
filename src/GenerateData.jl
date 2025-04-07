using DataFrames, Distributions, Random

"""Compute the radiative rate"""
function Rad_rate(molecule, exc, row)
    B = molecule.Cross_section(exc.wl) * molecule.QY
    I = exc.Intensity(row[:x], row[:y])
    return B * sin(row[:theta])^2 * sin(row[:phi])^2 * I
end

# Compute the MO efficiency
function MOeff(NA, row)
    delta = asin(NA)
    tetha = row[:theta]
    sincosm = sin(tetha - delta) * cos(tetha - delta)
    sincosp = sin(tetha + delta) * cos(tetha + delta)
    sincos0 = sin(delta) * cos(delta)  # computed but not used further
    return (2 * delta - sincosp + sincosm) / π / 2
end

# Compute a history of time events based on exponential waiting times
function History(molecule, r)
    t = 0.0
    hist = Float64[]
    dpb = r / molecule.M
    ddark = r / molecule.Darks[1]  # Julia indexing: first element of Darks corresponds to Python's Darks[0]
    
    tdark = rand(Exponential(ddark))
    tpb = rand(Exponential(dpb))
    
    while tdark < tpb
        t += tdark
        push!(hist, t)
        # For tground, Python uses lambd=1/molecule.Darks[1] (i.e. second element in Python)
        # so here we assume molecule.Darks[2] corresponds to that value.
        tground = rand(Exponential(1 / molecule.Darks[2]))
        t += tground
        push!(hist, t)
        tdark = rand(Exponential(ddark))
        tpb = rand(Exponential(dpb))
    end
    t += tpb
    push!(hist, t)
    return hist
end

# Generate data based on the molecule, excitation, and a DataFrame sample.
function Generate_Data(molecule, exc, sample::DataFrame, NA, wl0, wl1)
    data = deepcopy(sample)
    rs = Float64[]
    moeffs = Float64[]
    hists = Vector{Vector{Float64}}()
    names = String[]
    
    for row in eachrow(sample)
        rate = Rad_rate(molecule, exc, row)
        push!(rs, rate)
        push!(moeffs, MOeff(NA, row))
        push!(hists, History(molecule, rate))
        push!(names, molecule.name)
    end
    
    data.R = rs
    data.MOeff = moeffs
    data.History = hists
    data.Name = names
    # Assuming Spectral_fraction returns a scalar; fill a column with that constant.
    data.Spectral_fraction = fill(molecule.Spectral_fraction(wl0, wl1), nrow(data))
    
    return data
end