### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 80a4e1cc-4d59-4c8a-9d65-1fe0d2d77bf2
using Pkg; Pkg.activate(ENV["JBoldLab"])

# ╔═╡ 23b4e644-8331-4aa6-8369-47ce3ff0e143
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Images
	using ImageBinarization
	using Colors
	using Plots
	using Printf
	using Interpolations
	using QuadGK
	using Markdown
	using InteractiveUtils
	using LsqFit
	using Statistics
	using StatsBase
	using Unitful 
	using UnitfulEquivalences 
	using PhysicalConstants
	using Peaks
	using FFTW
	using DSP
	using Clustering
	import Glob
end

# ╔═╡ 3a7182c2-8ffd-4d6a-975f-42f9ff1b0744
using Test

# ╔═╡ 9f71bc31-df8f-48c5-815a-c8db9e46100f
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 19035c0d-93a7-4bb9-be05-d2a9b9ac4619
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ f4e379c4-a2b2-4703-bcbc-f6d7c996354a
lbl = ingredients("../src/BoldLab.jl")

# ╔═╡ d88ba94f-dd9e-4c46-adac-4c3f7e556cad
lbl.BoldLab.print_greeting()

# ╔═╡ 51710bb8-9c4a-4abc-9fe9-02e87bd4e4c5
PlutoUI.TableOfContents(title="Laser Lab CMOS analysis", indent=true)

# ╔═╡ d1ace15e-fe1a-4552-9144-a0824ae8ae0f
md"""
# Contributors

- J.J. Gómez Cadenas (jjgomezcadenas@gmail.com)
- M. Elorza
- P. Herrero
"""

# ╔═╡ 0937a6fc-6936-47d0-80de-8a38bb9a6a37
md"""
## Dropbox path
Here you must enter de path of your LaserLab Dropbox. With this path, the program can access to every file in the Dropbox.
"""

# ╔═╡ 6d04e8fe-c174-4be5-bffc-40945e8074e5
LaserLabp="C:\\Users\\Mikel\\LaserLab Dropbox\\"

# ╔═╡ 8b4fa5c3-0efe-4aa7-98cf-d37c1e12bc74
md"""
# Scheme
- Not Band Pass Filter (NBPF): Two wheels with filters, from 2 to 11. No band pass filter applied. 
- Band Pass Filter: First wheel includes a BPF and filters 3 to 6. Second wheel filters 7 to 11. When taking data with second wheel, first wheel includes BPF.
- New Filter Scheme (NFS): Two wheels with 10 filters, from 3 to 11. Band pass filters before the wheels. Filter set has been changed. 
"""

# ╔═╡ e6de7c36-135f-43f6-a0e2-3eb65fd6cf57
schema = ["NFS", "BPF", "NBPF"]

# ╔═╡ 7789a654-af1b-408a-86fc-26cc8ff221da
md""" Select scheme : $(@bind sch Select(schema))"""

# ╔═╡ a2acebd8-c197-4418-855f-cb08ee24021b
md"""
# G2SL Characteristics
#"""

# ╔═╡ 5e02e22e-7c05-4af2-b7d5-f5c61f2cee11
g2slp=joinpath(LaserLabp,"Proyectos\\FLUORI\\G2\\G2SL\\BOLD_104_SIL_GB_onquartz")

# ╔═╡ 025d7b0e-241a-4cf5-bf1a-47bc32a7e955
begin
g2df = lbl.BoldLab.load_df_from_csv(g2slp, 
"Emisi_espec_G2SIL_quartz.csv", lbl.BoldLab.enG)
	g2df = select!(g2df, [:L, :I])
	plot(g2df.L,g2df.I, lw=2, label="G2SL")
	xlabel!("λ (nm)")
	ylabel!("arbitrary units")
end

# ╔═╡ a01f191e-5f4a-4a56-9a39-b0494f02a0cd
md"""
- In the spectrum shown above, signal below 500 nm is most likely an artifact.
"""

# ╔═╡ d81cc134-1a4e-45a6-8657-f82372f8b66b
md"""
## Filters
"""

# ╔═╡ 4c92771a-8684-46ab-856a-627ac859a1fc
md"""
### Filter central values in nm
"""

# ╔═╡ f855c5a4-bf43-4247-b45f-c69957575ce0
if sch == "NFS"
	xfnm = [400.0, 438.0, 503.0, 549.0, 575.0, 600.0, 630.0, 676.0, 732.0, 810.0]
	wfnm = [40.0, 24.0,   40.0,   17.0,  15.0,  14.0,  38.0,  29.0,  68.0,  10.0]
elseif sch == "BPF"
	xfnm = [438.0,465.0,503.0,550.0,600.0,650.0,692.0,732.0]
	wfnm = [24.0, 30.0, 40.0, 49.0, 52.0, 60.0, 40.0, 68.0]
else
	xfnm = [420.0,438.0,465.0,503.0,550.0,600.0,650.0,692.0,732.0,810.0]
	wfnm = [10.0, 24.0, 30.0, 40.0, 49.0, 52.0, 60.0, 40.0, 68.0, 10.0]
end

# ╔═╡ 76383212-4bb1-425f-b82c-6dca2f6db236
begin
	filtnm = (center=xfnm,
	width  = wfnm,
	left = xfnm .- 0.5*wfnm,
	right = xfnm .+ 0.5*wfnm)

	println("Filter central values (nm) = ", filtnm.center, " width (nm) =", filtnm.width)
end

# ╔═╡ 755a1675-906d-417e-a974-a6a9520f879d
md"""
### Effect of filters in spectrum
"""

# ╔═╡ 9f30ed77-3d1a-474a-9659-7683ee429b03
begin
	wr=396.0:1.0:850.0
	fg2 = lbl.BoldLab.dftof(wr, g2df, "I")
	qs = [lbl.BoldLab.qpdf(fg2, filtnm.left[l], filtnm.right[l])/filtnm.width[l] for l in 1:length(xfnm)]
	pqyd = plot(g2df.L,g2df.I, lw=2, label="G2SL")	
	pfy = plot!(collect(wr), fg2.(wr), label="")
	pqs = scatter!(xfnm, qs, label="Filters")
	plot!(xfnm, qs, lw=2, label="")
	xlabel!("Filter (nm)")
	ylabel!("STD DC counts")
end

# ╔═╡ a921505a-33c6-4d13-9ac9-f2fb7c1a7b02
md"""
- The plot shows the expected discretized distribution of G2SL passed by the filters of the laser setup.  
"""

# ╔═╡ 16d5b19e-fdcb-4c44-9a92-d7fe8c0390df
begin
	qwl = lbl.BoldLab.qpdf(fg2, 0.0, 850.0)
	qflt = [lbl.BoldLab.qpdf(fg2, filtnm.left[l], filtnm.right[l]) for l in 1:length(xfnm)]
	qx = qflt ./qwl
	scatter(xfnm, qx, label="Fraction of total charge per filter", legend=:topleft)
	plot!(xfnm, qx, label="", legend=:topleft)
	xlabel!("Filter (nm)")
	ylabel!("STD DC counts")
end

# ╔═╡ 4a08a8db-f86d-49c4-961e-bcce3ec41653
md"""
- Plot shows the fraction of charge expected in each filter bin.
"""

# ╔═╡ 0122a81d-ee51-4355-8ddb-6429f15e8ea1
md"""
# Code
"""

# ╔═╡ c45b108e-62f9-473a-9975-eec4736d5df1
md""" 
## Select run
"""

# ╔═╡ da7dcbd8-75c3-42d5-9958-5ccbcc9624f1
begin
bodir=joinpath(LaserLabp,"Proyectos\\pdata\\")
bcmdir=joinpath(LaserLabp,"Proyectos\\data\\")
end

# ╔═╡ b9361344-165b-460c-85c6-0945f40612ef
begin
cmds = filter(name->occursin("CMOS",name), readdir(bcmdir))
md""" Select CMOS dir : $(@bind scmos Select(cmds))"""
end

# ╔═╡ a25b5dfe-e076-48c3-83c0-258b04add965
cmdir=joinpath(bcmdir,scmos)

# ╔═╡ 853b8b2e-66ed-4723-8884-213e5fd4a0e7
md"""
# Tests
"""

# ╔═╡ 20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
md"""
# Functions
"""

# ╔═╡ b9970588-422f-461f-addb-5169d2e6043e
function select(dir)
	readdir(dir)
end

# ╔═╡ c9e8c0f2-2776-43b5-87c6-c6c85e264924
let
subs = select(cmdir)
md""" Select substrate : $(@bind ssub Select(subs))"""
end

# ╔═╡ 5bd48d34-2b3c-47a5-8500-242e6a1f8f6f
subsp=joinpath(cmdir,ssub)

# ╔═╡ d5c5fbca-ae04-4ad7-a530-5f6e42f3e436
let
exps = select(subsp)
md""" Select experiment : $(@bind sexp Select(exps))"""
end

# ╔═╡ 71662eae-76e8-449a-93ee-abd1c3000b34
expp=joinpath(subsp,sexp)

# ╔═╡ f90a27b0-d729-4a9a-afe6-a713c090f467
let
runs = select(expp)
md""" Select run : $(@bind srun Select(runs))"""
end

# ╔═╡ b26174ef-ebb4-4fe3-a93d-d308d49488aa
runp=joinpath(expp,srun)

# ╔═╡ efc2a88d-309d-410c-b635-3c841e695c08
function load_df_from_csv(fp::String, csvg::CsvG; header=1)
    csvf   = CSV.File(fp; header=header, delim=csvg.delim, decimal=csvg.decimal)
    return DataFrame(csvf)
end

# ╔═╡ 0bf5c4a7-b369-4438-a987-253f242dd88f
function dummy(x,y)
	x+y
end

# ╔═╡ 71ed6fbd-0342-4cf5-9d8c-aa8f791d85f1
@test dummy(3,4) == 7

# ╔═╡ Cell order:
# ╠═80a4e1cc-4d59-4c8a-9d65-1fe0d2d77bf2
# ╠═23b4e644-8331-4aa6-8369-47ce3ff0e143
# ╠═3a7182c2-8ffd-4d6a-975f-42f9ff1b0744
# ╠═9f71bc31-df8f-48c5-815a-c8db9e46100f
# ╠═19035c0d-93a7-4bb9-be05-d2a9b9ac4619
# ╠═f4e379c4-a2b2-4703-bcbc-f6d7c996354a
# ╠═d88ba94f-dd9e-4c46-adac-4c3f7e556cad
# ╠═51710bb8-9c4a-4abc-9fe9-02e87bd4e4c5
# ╠═d1ace15e-fe1a-4552-9144-a0824ae8ae0f
# ╠═0937a6fc-6936-47d0-80de-8a38bb9a6a37
# ╠═6d04e8fe-c174-4be5-bffc-40945e8074e5
# ╠═8b4fa5c3-0efe-4aa7-98cf-d37c1e12bc74
# ╠═e6de7c36-135f-43f6-a0e2-3eb65fd6cf57
# ╠═7789a654-af1b-408a-86fc-26cc8ff221da
# ╠═a2acebd8-c197-4418-855f-cb08ee24021b
# ╠═5e02e22e-7c05-4af2-b7d5-f5c61f2cee11
# ╠═025d7b0e-241a-4cf5-bf1a-47bc32a7e955
# ╠═a01f191e-5f4a-4a56-9a39-b0494f02a0cd
# ╠═d81cc134-1a4e-45a6-8657-f82372f8b66b
# ╠═4c92771a-8684-46ab-856a-627ac859a1fc
# ╠═f855c5a4-bf43-4247-b45f-c69957575ce0
# ╠═76383212-4bb1-425f-b82c-6dca2f6db236
# ╠═755a1675-906d-417e-a974-a6a9520f879d
# ╠═9f30ed77-3d1a-474a-9659-7683ee429b03
# ╠═a921505a-33c6-4d13-9ac9-f2fb7c1a7b02
# ╠═16d5b19e-fdcb-4c44-9a92-d7fe8c0390df
# ╠═4a08a8db-f86d-49c4-961e-bcce3ec41653
# ╠═0122a81d-ee51-4355-8ddb-6429f15e8ea1
# ╠═c45b108e-62f9-473a-9975-eec4736d5df1
# ╠═da7dcbd8-75c3-42d5-9958-5ccbcc9624f1
# ╠═b9361344-165b-460c-85c6-0945f40612ef
# ╠═a25b5dfe-e076-48c3-83c0-258b04add965
# ╠═c9e8c0f2-2776-43b5-87c6-c6c85e264924
# ╠═5bd48d34-2b3c-47a5-8500-242e6a1f8f6f
# ╠═d5c5fbca-ae04-4ad7-a530-5f6e42f3e436
# ╠═71662eae-76e8-449a-93ee-abd1c3000b34
# ╠═f90a27b0-d729-4a9a-afe6-a713c090f467
# ╠═b26174ef-ebb4-4fe3-a93d-d308d49488aa
# ╠═853b8b2e-66ed-4723-8884-213e5fd4a0e7
# ╠═71ed6fbd-0342-4cf5-9d8c-aa8f791d85f1
# ╠═20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
# ╠═b9970588-422f-461f-addb-5169d2e6043e
# ╠═efc2a88d-309d-410c-b635-3c841e695c08
# ╠═0bf5c4a7-b369-4438-a987-253f242dd88f
