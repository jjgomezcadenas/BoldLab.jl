using Revise
using CSV
using DataFrames
	
	
using Statistics
using StatsBase
using DelimitedFiles
using Images, FileIO, ImageIO
using SparseArrays

using BoldLab
using SimpleLogger
using JStepFinder
using StepAnalysis
using LabStepAnalysis
using histos
using NPZ
using Unitful

import Unitful:
    nm, μm, mm, cm, m,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W


list_subfolders(dir) = [basename(e) for e in readdir(dir; join=true) if isdir(e) && !startswith(basename(e), ".")]

BASEDIRS = Dict("Data"=>"/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign", 
				"MonteCarlo" => "/Users/jjgomezcadenas/Projects/BoldLab/pluto/npy")

base_directory(label) = BASEDIRS[label]


function scan_level(dir::AbstractString)
    entries = readdir(dir; join=true, sort=true)
    vis     = filter(e -> !startswith(basename(e), "."), entries)

    subdirs = filter(isdir, vis)
    npys    = filter(e -> endswith(e, ".npy"), vis)
    tiffs   = filter(e -> endswith(lowercase(e), ".tif")  ||
                           endswith(lowercase(e), ".tiff"), vis)

    return basename.(subdirs), basename.(npys), basename.(tiffs)
end


function update_nstep!(df::DataFrame)
    @assert "nmol" in names(df) "Missing column 'nmol'"
    @assert "nstep" in names(df) "Missing column 'nstep'"

    gdf = groupby(df, :nmol)
    for g in gdf
        g[!, :nstep] .= nrow(g)  # assign group size to all rows in group
    end
    return df
end


function cut_length(df; maxl=150)

	function is_bad_molecule(gdf)
    	any(gdf.stepLength .> maxl)
	end

	grouped = groupby(df, :nmol)
	
	# Keep only groups that do not satisfy the bad condition
	good_groups = filter(gdf -> !is_bad_molecule(gdf), grouped)
	
	# Reassemble the cleaned DataFrame
	vcat(good_groups...)
end


function get_vals_from_sparse(sm)
	rows = Int[]
	cols = Int[]
	vals = Float64[]

	for j in 1:size(sm, 2)
	    for idx in sm.colptr[j]:(sm.colptr[j+1]-1)
	        i = sm.rowval[idx]
	        v = sm.nzval[idx]
	        push!(rows, i)
	        push!(cols, j)
	        push!(vals, v)
	    end
	end
	return rows, cols, vals
end


function reduced_chi_squared_per_element(DX::SparseMatrixCSC{<:AbstractVector, Int},
                                          FX::SparseMatrixCSC{<:AbstractVector, Int})

    @assert size(DX) == size(FX) "DX and FX must have the same size"

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    for col in 1:size(DX, 2)
        for idx in DX.colptr[col]:(DX.colptr[col+1] - 1)
            row = DX.rowval[idx]
            dvec = DX[row, col]
            fvec = FX[row, col]

            @assert length(dvec) == length(fvec) "Vectors at ($row,$col) must be same length"

            chi2 = 0.0
            dof = 0

            for k in eachindex(dvec)
                d = dvec[k]
                f = fvec[k]
                if d > 0
                    chi2 += (d - f)^2 / d
					#chi2 += (d - f)^2 
                    dof += 1
                end
            end

            if dof > 0
                push!(rows, row)
                push!(cols, col)
                push!(vals, chi2 / dof)
            end
        end
    end

    return sparse(rows, cols, vals, size(DX, 1), size(DX, 2))
end


function trz_anal(imst::AbstractArray{T,3}; nsigma::Real=10,
                  niter::Real=3, thr::Real=0.1, sel="base",
                  maxl::Real=20) where {T<:Real}

    xmc, ymc, tmc = size(imst)
    TRZ = build_traces(imst)

    println("Reading MC TRZ:  size = $(size(TRZ))")

    df, npixel, ngfit, nfailed, _, _, DX, FX, _, _ = fit_data3(TRZ; ped=0.0, niter=niter, thr=thr, sel=sel)

    println("""
    - niter = $(niter), thr=$(thr) 
	- Fitted $(npixel) pixels 
	- failed fit fraction =$(nfailed/npixel)
	- Found $(length(nonzeros(FX))) good fits
	""")

    println("keeping molecules where niter= $(niter)")
    sdf = filter(row -> !(row.nstep <niter), df)
    update_nstep!(sdf)	

    nsdf = length(unique(sdf.nmol))
    println("good fits ater cut = $(nsdf)")

    mstl = maximum(sdf.stepLength)
    println("maximum step length in df = $(mstl)")
    mx = maximum([100, mstl - maxl])

    println("Cutting molecules where at least one step is larger than $(mx)")
    sdf2 = cut_length(sdf; maxl=mx)

    nsdf2 = length(unique(sdf2.nmol))
    println("good fits ater cut = $(nsdf2)")


    bchi2_matrix = reduced_chi_squared_per_element(DX, FX)
	I, J, Vc = get_vals_from_sparse(bchi2_matrix)

	println("""
	### bFIT
	- mean chi2 = $(mean(Vc))
	- std chi2 = $(std(Vc))
	- max chi2 = $(maximum(Vc))
	""")
end

niter = 5
thr=0.1

mc_file = "/Users/jjgomezcadenas/Projects/BoldLab/pluto/npy/n3d_mc_1e3_noise.npy"
root_dir = "/Users/jjgomezcadenas/BoldLab/BOLD/Alfonso/2025_03_Blank_measurements_campaign/"
blank_dir = "/Week2/Data/0_31_Monday/Blank1/"
evap_dir = "/Week3/Data/10_Thursday/W3_SVE0_2/"
free_dir ="/Week3/Data/10_Thursday/W3_SVN0_2/"


immc = npzread(mc_file)
println("Reading MC file: Image size = $(size(immc))")

for thr in 0.1:0.1:0.5  
    trz_anal(immc; niter=niter, thr=thr, sel="base")
end



