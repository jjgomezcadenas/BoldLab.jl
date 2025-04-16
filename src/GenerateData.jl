using DataFrames, Distributions, Random
using .SimpleLogger
"""Compute the radiative rate"""
function rad_rate(fm::FBMolecule, exc::LaserExcitation, row)
    B = cross_section(fm, exc.λ)
    I = intensity(exc, row[:x], row[:y])
    B * sin(row[:theta])^2 * sin(row[:phi])^2 * I
end

"""Compute the MO efficiency"""
function mo_eff(NA, row)
    delta = asin(NA)
    tetha = row[:theta]
    sincosm = sin(tetha - delta) * cos(tetha - delta)
    sincosp = sin(tetha + delta) * cos(tetha + delta)
    #sincos0 = sin(delta) * cos(delta)  # computed but not used further
    return (2 * delta - sincosp + sincosm) / π / 2
end

""" Compute a time trajector based on exponential waiting times"""
function trajectory(molecule::FBMolecule, r::typeof(1.0Hz))
    t = 0.0s
    hist = []
    dpb = r / molecule.M             # rate/PB-rate 

    debug("typeof(molecule.M) = $(typeof(molecule.M))")
    ddark = r / molecule.Darks[1]    # rate/PB-rate

    debug("typeof(ddark) = $(typeof(ddark))")
    
    argtdark = (1.0/ddark*Hz)

    debug("typeof(argtdark) = $(typeof(argtdark))")

    argtpb = (1.0/dpb*Hz)
    argtgr = (molecule.Darks[2])/1.0s

    debug("dpb = $(dpb), ddark = $(ddark)")
    tdark = rand(Exponential(argtdark)) * 1.0s
    tpb = rand(Exponential(argtpb)) * 1.0s

    debug("tdark = $(tdark), tpb = $(tpb)")
    debug("typof (dark) = $(typeof(tdark)), tpb = $(tpb)")
    while tdark < tpb
        t += tdark
        push!(hist, t)
        debug("time in dark = $(t)")
       
        tground = rand(Exponential(argtgr)) * 1.0s
        debug("tground = $(tground) ")
        t += tground
        push!(hist, t)
        debug("time ground = $(t) ")
        tdark = rand(Exponential(argtdark)) * 1.0s
        tpb = rand(Exponential(argtpb)) * 1.0s

        debug("new dice: tdark = $(tdark), tpb = $(tpb)")
    end
    t += tpb
    push!(hist, t)
    debug("time to photobleach = $(t)")
    hist
end

function spectral_fraction(fb::FBMolecule, λmin::typeof(1.0nm), λmax::typeof(1.0nm))
	quadgk(fb.fm.pdf, λmin/nm, λmax/nm)[1]
end


"""Generate data based on the molecule, excitation, and a DataFrame sample"""
function generate_data(molecule::FBMolecule, exc::LaserExcitation,
					   sample::DataFrame, NA::Float64, 
					   λmin::typeof(1.0nm), λmax::typeof(1.0nm))
    data = deepcopy(sample)
    rs = []
    moeffs = Float64[]
	spx  = Float64[]
    trajs = []
    names = String[]
	#molecule = fmb.fm
    
    for row in eachrow(sample)
		rate = rad_rate(molecule, exc, row)
        push!(rs, rate)
		push!(spx, spectral_fraction(molecule, λmin, λmax))
        push!(moeffs, mo_eff(NA, row))
        push!(trajs, trajectory(molecule, rate))
        push!(names, molecule.name)
    end
    
    data.R = rs
    data.MOeff = moeffs
    data.Trajectories = trajs
    data.Name = names
	data.Spx = spx
    
    return data
end