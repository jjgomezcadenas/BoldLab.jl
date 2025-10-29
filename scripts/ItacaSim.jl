#!/usr/bin/env julia
# ItacaSim.jl
# Simulation of ITACA detector focusing system
# Converted from Pluto notebook to command-line script

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using BoldLab
using CSV
using DataFrames
using Plots
using Printf
using Statistics
using StatsBase
using Images, FileIO
using Interpolations
using QuadGK
using Distributions, Random
using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018

import Unitful:
    nm, μm, mm, cm, m,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    V, kV,
    μJ, mJ, J,
    mW, W

# ----------------------------
# Simple named-argument parser (--key=value or --key value)
# ----------------------------
function parse_args(defaults::Dict{String,Any})
    args = copy(ARGS)
    conf = copy(defaults)
    i = 1
    while i <= length(args)
        tok = args[i]
        if startswith(tok, "--")
            kv = split(tok[3:end], "=", limit=2)
            key = kv[1]
            if length(kv) == 2
                val = kv[2]
                conf[key] = val
                i += 1
            else
                if i == length(args)
                    error("Missing value for --$key")
                end
                conf[key] = args[i+1]
                i += 2
            end
        else
            i += 1
        end
    end
    return conf
end

function parse_float!(conf::Dict{String,Any}, key::String)
    conf[key] = parse(Float64, string(conf[key]))
end

function parse_int!(conf::Dict{String,Any}, key::String)
    conf[key] = parse(Int, string(conf[key]))
end

function parse_bool!(conf::Dict{String,Any}, key::String)
    val = lowercase(string(conf[key]))
    conf[key] = val in ["true", "1", "yes"]
end

# ----------------------------
# Main simulation function
# ----------------------------
function main()
    # Default parameters
    defaults = Dict{String,Any}(
        "N" => 1000000,              # Number of molecules (integer)
        "NM" => 1.0e6,               # Number of molecules (float, for display)
        "with_noise" => false,       # Add noise to simulation
        "laser_power" => 1.0,        # Laser power in W
        "laser_lambda" => 450.0,     # Laser wavelength in nm
        "texp_sec" => 0.1,           # Laser wavelength in nm
        "dim_x" => 36.0,             # Sample x dimension in μm
        "dim_y" => 36.0,             # Sample y dimension in μm
        "sigma" => 36.0 * 2.5,             # Gaussian beam sigma in μm
        "center_x" => 0.0,           # Beam center x in μm
        "center_y" => 0.0,           # Beam center y in μm
        "pbcycles" => 1e4,           # Photobleaching cycles
        "dkcycles" => 5e5,           # Dark cycles
        "tdark" => 2.0,              # Time in dark states in s
        "theta" => pi/2.0,           # Theta angle in radians (default π/2)
        "phi" => pi/2.0,             # Phi angle in radians (default π/4)
        "model" => "1D",
        "output_file" => "itaca_1e6_noise_1W_50mu.csv",
        "data_dir" => joinpath(dirname(@__FILE__), "..", "data")
    )

    conf = parse_args(defaults)

    # Parse numeric values
    for k in ["laser_power", "laser_lambda", "dim_x", "dim_y", "sigma",
              "center_x", "center_y", "pbcycles", "dkcycles", "tdark", "NM",
              "theta", "phi", "texp_sec"]
        parse_float!(conf, k)
    end
    parse_int!(conf, "N")
    parse_bool!(conf, "with_noise")

    # Extract parameters with proper units
    N = conf["N"]      # take N as maximum 
    NM = conf["NM"]
    with_noise = conf["with_noise"]
    laser_power = conf["laser_power"] * W
    laser_λ = conf["laser_lambda"] * nm
    dimensions = [conf["dim_x"] * μm, conf["dim_y"] * μm]
    sigma = conf["sigma"] * μm
    center = (conf["center_x"] * μm, conf["center_y"] * μm)
    pbcycles = conf["pbcycles"]
    dkcycles = conf["dkcycles"]
    tdark = conf["tdark"] * s
    theta = conf["theta"]
    phi = conf["phi"]
    output_file = string(conf["output_file"])
    model = string(conf["model"])
    data_dir = string(conf["data_dir"])
    texp_s =conf["texp_sec"]

    # Print simulation parameters
    println("\n" * "="^70)
    println("ITACA DETECTOR SIMULATION")
    println("="^70)
    println("\n### Track and Detector Parameters ###")

    # Define barn unit
    barn = 1E-24 * cm^2

    # Setup laser
    println("### Laser Setup ###")
    lsr = Laser(laser_λ, laser_power)
    I = uconvert(Hz*cm^-2, photon_density(laser_λ, laser_power, sigma, sigma))
    println("  - Wavelength: $(lsr.λ)")
    println("  - Power: $(lsr.P)")
    println("  - Photon energy: $(photon_energy(lsr.λ))")
    println("  - Intensity: $(I)")

    # Setup objective and filter
    println("\n### Objective and Filter ###")
    obj = Objective("highNA", 0.95, 100.0)
    bf = BandFilter(450.0nm, 800.0nm)
    println("  - NA: $(obj.NA)")
    println("  - Magnification: $(obj.M)")
    println("  - Filter λmin: $(bf.λmin)")
    println("  - Filter λmax: $(bf.λmax)")

    # Setup CMOS
    println("\n### CMOS Camera ###")
    of = CMOS("ORCA_Flash4", 6.5μm, 2049, 1.0, 0.06Hz, 60, oflash4_eff)
    println("  - Name: $(of.name)")
    println("  - Pixel size: $(of.pixelsize)")
    println("  - Readout noise: $(of.readout_noise) p.e.")
    println("  - Dark current: $(of.dark_current) p.e./pixel/s")

    # Load molecule data
    println("\n### Molecule Setup ###")
    naph3f, naph3c = fm_from_csv(data_dir,
                                  "naph3_emission.csv",
                                  "naph3_abs.csv")

    naph3fB = FBMolecule(naph3c, "naph3-ba", pbcycles, (dkcycles, tdark))
    println("  - Name: $(naph3fB.name)")
    println("  - QY: $(naph3fB.QY)")
    println("  - M: $(naph3fB.M)")
    println("  - Darks: $(naph3fB.Darks)")
    println("  - Cross section @ $(lsr.λ): $(cross_section(naph3fB, lsr.λ)/barn) barn")

    # Define laser excitation
    lexc = LaserExcitation(lsr, sigma, center)

    # Generate sample
    println("\n### Generating Sample ###")
    println("Creating sample of $N molecules...")
    println("  - Dimensions: $(dimensions[1]) × $(dimensions[2])")

    println("  - Using model => $(model)")
    if model == "1D"
        sample = Sample_1D(N, dimensions, theta, phi)
    elseif model == "2D"
        sample = Sample_2D(N, dimensions)
    elseif model == "3D"
        sample = Sample_3D(N, dimensions)
    else
        error("Unknown model: $(model). Choose '1D', '2D', or '3D'")
    end
    println("✓ Sample generated")

    # Generate data with progress reporting
    println("\n### Generating Integrated Signal Data ###")
    println("Processing $N molecules...")
    println("This may take a while for large samples...")

    # Start timing
    t_start = time()

    # Call generate_data_integrated_signal with progress reporting
    # We'll wrap the function to add progress
    data = generate_data_integrated_signal_with_progress(
        naph3fB, lexc, sample, obj.NA, bf.λmin, bf.λmax, texp_s
    )

    t_elapsed = time() - t_start

    println("✓ Data generation complete!")
    println(@sprintf("  Time elapsed: %.2f seconds", t_elapsed))
    println(@sprintf("  Processing rate: %.0f molecules/second", N/t_elapsed))

    # Save results (strip units from columns)
    println("\n### Saving Results ###")
    println("Writing data to: $output_file")

    # Create a copy of data without units
    data_no_units = DataFrame(
        x_μm = ustrip.(μm, data.x),
        y_μm = ustrip.(μm, data.y),
        theta = data.theta,
        phi = data.phi,
        R_Hz = ustrip.(Hz, data.R),
        Tb_s = ustrip.(s, data.Tb),
        MOeff = data.MOeff,
        Name = data.Name,
        Spx = data.Spx,
        Nphot = data.Nphot
    )

    CSV.write(output_file, data_no_units)
    println("✓ Data saved successfully (units stripped for CSV compatibility)")

    # Summary statistics
    println("\n### Summary Statistics ###")
    println(@sprintf("  Total molecules: %d", nrow(data)))
    println(@sprintf("  Mean emission rate: %.2e Hz", mean(data.R)))
    println(@sprintf("  Mean MO efficiency: %.4f", mean(data.MOeff)))
    println(@sprintf("  Mean spectral fraction: %.4f", mean(data.Spx)))

    println("\n" * "="^70)
    println("SIMULATION COMPLETE")
    println("="^70 * "\n")

    return data
end

# Wrapper function to add progress reporting to generate_data_integrated_signal
function generate_data_integrated_signal_with_progress(
    molecule::FBMolecule, exc::LaserExcitation,
    sample::DataFrame, NA::Float64,
    λmin::typeof(1.0nm), λmax::typeof(1.0nm), texp_s::Float64
)
    data = deepcopy(sample)
    rs = []
    pbt = []
    nevt = []
    moeffs = Float64[]
    spx = Float64[]
    names = String[]
    Nphot = Float64[]

    n_total = nrow(sample)
    report_interval = max(1, n_total ÷ 20)  # Report ~20 times

    texp = texp_s*1.0s
    for (idx, row) in enumerate(eachrow(sample))
        rate = BoldLab.rad_rate(molecule, exc, row)
        push!(rs, rate)
        pb = BoldLab.photobleaching_time(molecule, rate)
        push!(pbt, pb)
        sf = BoldLab.spectral_fraction(molecule, λmin, λmax)
        push!(spx, sf)
        moeff = BoldLab.mo_eff(NA, row)
        push!(moeffs, moeff)
        push!(names, molecule.name)
        nevent = rate * texp * sf * moeff
        if pb < texp
            nevent = rate * pb * sf * moeff
        end
        push!(Nphot, ustrip.(nevent))

        # Progress reporting
        if idx % report_interval == 0 || idx == n_total
            pct = round(100 * idx / n_total; digits=1)
            @printf("  Progress: %6d/%d molecules (%.1f%%)\n", idx, n_total, pct)
        end
    end

    data.R = rs
    data.Tb = pbt
    data.MOeff = moeffs
    data.Name = names
    data.Spx = spx

    # Calculate total number of photons: R * Tb * MOeff * Spx
    # R (Hz) * Tb (s) = dimensionless count
    data.Nphot = Nphot

    return data
end

# Run main function
main()
