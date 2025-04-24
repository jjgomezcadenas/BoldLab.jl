# """
# Module dffunctions.jl provides tools to read data frames expressing
# data series relevant for LaserLab software
# (e.g, molecular cross sections, filters), and describe them when
# relevant as functions.
# """

#module dffunctions
#using Revise
#export spG, enG, enG2
#export load_df_from_csv, dftof, qpdf, ftopdf
#export FMolecule, fm_from_csv
using DataFrames
using CSV
using Interpolations
using QuadGK
using Unitful
using Printf

import Unitful:
    nm, μm, mm, cm, m,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, W

"""
    generate_filename(N, laser_power, pbcycles, dkcycles, tdark; with_noise=true)

Generates a descriptive filename like:
"n3d_mc_1e1_5mW_pb_1e5_dk_1e6_noise.npy"

# Arguments
- `N`: Number of molecules.
- `laser_power`: Laser power (e.g., `5mW` from Unitful).
- `pbcycles`: Photobleaching cycles.
- `dkcycles`: Dark cycles.
- `tdark`: Time in dark state (not included in filename).
- `with_noise`: Boolean flag to include `_noise` in the filename.

# Returns
- A filename string.
"""
function generate_filename(N, laser_power, pbcycles, dkcycles, tdark; with_noise=true)
    # Convert numbers to compact scientific notation
    function sci_str(x)
        str = @sprintf("%.0e", x)
        replace(str, r"e\+?0*" => "e")  # compact form
    end

    # Extract laser power value and unit
    power_val = ustrip(mW, laser_power)
    laser_str = "$(Int(round(power_val)))mW"

    # Build base string
    fname = "n3d_mc_$(sci_str(N))_$(laser_str)_pb_$(sci_str(pbcycles))_dk_$(sci_str(dkcycles))"
    if with_noise
        fname *= "_noise"
    end
    return fname * ".npy"
end


#-----------

"""
	CsvG

A struct which describes a simple grammar needed to read CSV files
written with spanish/english characters

# Fields
- `delim::Char`  : Delimiter used (e.g., ',' in english, ';' in spanish)
- `decimal::Char`: Symbol to represent decinal point (',' '.')
"""
struct CsvG  # stands for CSV Grammar
	delim::Char
	decimal::Char
end

spG = CsvG(';',',')  # spanish: decimals represented with ',' delimited with ';'
enG = CsvG(',','.')  # english
enG2 = CsvG('\t','.')  # english with tabs


"""
	Gf

A struct representing a generalized function:

# Fields
- `N::Number`: Normalization constant: pdf(x) = f(x)/N.
- `f::Any`: Any integrable function.
- `pdf::Any`: Normalized to area 1 (PDF).
"""
struct Gf          # general function
	N  ::Number
	f  ::Any
	pdf::Any
end

"""
	load_df_from_csv(path::String, fname::String, csvg::CsvG)

Load a dataframe from a csv file.

# Arguments
- `path::String`: a path to the data.
- `fname::String`: name of the (csv) file expressing the df.
- `csvg=spG`: CVS coding (spanish, default) or english.
"""
function load_df_from_csv(path::String, fname::String, csvg::CsvG; header=1)
    name   = string(path, "/", fname)
    csvf   = CSV.File(name; header=header, delim=csvg.delim, decimal=csvg.decimal)
    return DataFrame(csvf)
end


"""
	dftof(wl, df::DataFrame, cname::String,
			   bkgnd::Float64=0.0)

Return an interpolated function, valid in range wl.

# Arguments
- `wl`::interpolation range.
- `df::DataFrame`: data frame holding the data.
- `cname::String`: name of the column holding data to be interplated.
- `overflow::Real`: value of the data outside interpolation range.
- `scale::Real`: value to re-scale the data 
"""
function dftof(wl, df::DataFrame, cname::String, overflow::Real=0.0, scale::Real=1.0)
    li = LinearInterpolation(wl, df[!,cname], extrapolation_bc=Line())
	return gfpdf_(li, wl[1], wl[end], overflow, scale)
end


function gfpdf_(fi, xmin::Real, xmax::Real, bkgnd::Real=0.0, scale::Real=1.0)
	function fn(x)
	   	if x < xmin || x > xmax
			return bkgnd
		else
			return fi(x) * scale
		end
	end
	return fn
end


"""
	qpdf(f, λmin::Real, λmax::Real)

Return the integral of f in the interval (λmin, λmax)
(Syntactic sugar for quadgk)

# Arguments
- `f::Function`: Function to be integrated.
- `λmin::Number`: lower bound of range.
- `λmax::Number`: upper bound of range.
"""
function qpdf(f, λmin::Real, λmax::Real)
	return quadgk(f, λmin, λmax)[1]
end


"""
	ftopdf(wl, f)

Compute the PDF of function f in range wl.

# Arguments
- `wl`: Range of application.
- `f::Function`: Input function.
"""
function ftopdf(wl, f)
	function pdf(x)
		return f(x) / N
	end
    N = qpdf(f, wl[1], wl[end])
	return N, pdf
end


#end #module



