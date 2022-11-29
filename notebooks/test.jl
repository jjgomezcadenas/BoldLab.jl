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
# Root path
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

# ╔═╡ 7565e604-d434-4ba2-9594-f8d328dcca0f
md"""
# Parameters of the CCD
"""

# ╔═╡ adf99e37-bb92-4d6d-8a15-8dd159af31da
begin
	texp = 10.0
	adctopes = 0.48
	pixelsize = 16.0
	pedestal = 100.0 * pixelsize
	dcperpixel = 0.06 * pixelsize * texp
	dcperpixelc = 0.06 * pixelsize * texp / adctopes
	enoise = 1.0 * pixelsize
	enoisec = enoise / adctopes
	
	md"""
- Time of exposition (in sec) = $texp
- adc to photoelect = $adctopes
- pixelsize = $pixelsize
- pedestal = $pedestal
- dark current per pixel in e- = $dcperpixel
- dark current + noise per pixel in e- = $(dcperpixel + enoise)
- dark current in counts $dcperpixelc 
- dark current + noise per pixel in e- = $(dcperpixelc + enoisec)
	
"""
end

# ╔═╡ c45b108e-62f9-473a-9975-eec4736d5df1
md""" 
# Select run
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

# ╔═╡ fd7b867f-2d49-4036-9b6f-0bb6966c32e6
md"""
runp and runop are the input and output data folders
"""

# ╔═╡ 0122a81d-ee51-4355-8ddb-6429f15e8ea1
md"""
# Analysis
"""

# ╔═╡ 54faa211-0796-4899-b56a-90d0ea73ce4a
md"""
## Position of points 
- Positions folder contains measurements of the poistions of each point. 
"""

# ╔═╡ b896d27f-1a5f-4c25-86fa-4717bf94770b
md"""
## Algorithm
- Find the position of the measured points
- Measure the dark current (filter and point independent)
- Select the spot (filter independent)
- For each point, measure the response for all the filters
"""

# ╔═╡ 838eafb7-5c5b-4df5-98b1-dfbc94e58d1e
md""" 
## Analysis of single point
"""

# ╔═╡ 9646801d-da05-4101-b453-8d12ed570475
md"""
### Dark current 
- Dark folder contains measurements of the dark current for each filter. 
"""

# ╔═╡ ff92b592-9cdd-4b8c-a256-f09d166febc6
md""" Select nsigma for DC suppression: $(@bind spsigma NumberField(1.0:100.0, default=3.0))"""

# ╔═╡ 853b8b2e-66ed-4723-8884-213e5fd4a0e7
md"""
# Tests
"""

# ╔═╡ 68eaa46d-94c0-43d9-b9e3-9bd84db26b63


# ╔═╡ 20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
md"""
# Functions
"""

# ╔═╡ b9970588-422f-461f-addb-5169d2e6043e
function select(dir)
	readdir(dir)
end

# ╔═╡ c9e8c0f2-2776-43b5-87c6-c6c85e264924
begin
	cmdir=joinpath(bcmdir,scmos)
	odir=joinpath(bodir,scmos)
	subs = select(cmdir)
	md""" Select substrate : $(@bind ssub Select(subs))"""
end

# ╔═╡ d5c5fbca-ae04-4ad7-a530-5f6e42f3e436
begin
	subsp=joinpath(cmdir,ssub)
	subsop=joinpath(odir,ssub)
	exps = select(subsp)
	md""" Select experiment : $(@bind sexp Select(exps))"""
end

# ╔═╡ f90a27b0-d729-4a9a-afe6-a713c090f467
begin
	expp=joinpath(subsp,sexp)
	expop=joinpath(subsop,sexp)
	runs = select(expp)
	md""" Select run : $(@bind srun Select(runs))"""
end

# ╔═╡ b26174ef-ebb4-4fe3-a93d-d308d49488aa
begin
runp=joinpath(expp,srun)
runop=joinpath(expop,srun)
end

# ╔═╡ b0bf7aed-3234-44ca-ac05-023545ba0987
begin
	ppath = joinpath(runp, "Positions")
	pfb = readdir(ppath)[1]
	pdf = lbl.BoldLab.load_df_from_csv(ppath, pfb, lbl.BoldLab.enG; header=0)
	zp = pdf.Column1
	xp = pdf.Column2
	yp = pdf.Column3
	psxyz = scatter(xp, xp, zp, xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)")
	xysp =scatter(xp, yp, label="")
	xlabel!("x (mm)")
	ylabel!("y (mm)")
	xzp1 = scatter(xp, zp, label="")
	xlabel!("x (mm)")
	ylabel!("z (mm)")
	yzp1 = scatter(yp, zp, label="")
	xlabel!("y (mm)")
	ylabel!("z (mm)")
	plot(size=(750,750), psxyz, xysp, xzp1,yzp1, layout=(2,2))
end

# ╔═╡ 6d108da6-353f-47dc-8684-3bea555e3921
"""
Given a absolute path of a directory returns a list of files (of type dfiles) found
in the directory.

"""
function return_files(path::String, 
	                  dfiles="*.csv")
	xfiles=readdir(path,join=true)
	nxfiles=readdir(path)
	xfiles, nxfiles
	
end

# ╔═╡ a2ce4e48-73d0-491d-9de0-86a4409d358a
@test typeof(return_files(runp,"Dark"))==Tuple{Vector{String}, Vector{String}}

# ╔═╡ 411b89b3-178f-4e4a-bc83-af42e489b3eb
"""
Given the absolute path of the run returns the filter names by searching in the "Dark" folder.
"""
function flt_names(runp::String)
	p=joinpath(runp,"Dark")
	xnames=return_files(p)[2]
	fxnb = [split(name, "_")[2] for name in xnames]
	fxint = sort([parse(Int64, ff) for ff in fxnb])
	[string("Filter_", string(i)) for i in fxint]
end

# ╔═╡ f68c5806-6c44-4a13-8811-a193d59e45ba
[@test occursin("Filter",name) for name in flt_names(runp)]

# ╔═╡ 100fa551-fe67-4cba-b9e5-77676c2c2c6f
"""
Given the absolute path of the run returns the filter names by searching in the "Filter1" folder.
"""
function point_names(runp::String)
	p=joinpath(runp,"Filter1")
	xnames=return_files(p)[2]
	ns=[String(split(pd, "_")[1]) for pd in xnames]
end

# ╔═╡ 4d478749-935d-40a5-9431-0565ffa19e11
"""
"""
function get_image_path(runp::String, point_name::String, flt_name::String)
	ppath=joinpath(runp,point_name)
	fs=return_files(ppath)[2]
	fname=fs[findall([occursin(flt_name,name) for name in fs])]
	path=joinpath(ppath,fname[1])
end

# ╔═╡ 09dc9013-eec5-42ea-8085-87ccc71a54aa
"""
"""
function get_image(runp::String, point_name::String, flt_name::String)
	path=get_image_path(runp,point_name,flt_name)
	imgdf = DataFrame(CSV.File(path, header=false,delim="\t"))
	df1=Matrix(imgdf)
end
	

# ╔═╡ 28885012-f552-4ded-bfe4-c2507a685146
"""
"""
function get_nfimage(runp::String, point_name::String)
	path=get_image_path(runp,"Filter1",point_name)
	imgdf = DataFrame(CSV.File(path, header=false,delim="\t"))
	df1=Matrix(imgdf)
end

# ╔═╡ 1a0727c6-ea80-4026-a24a-9682737b90ec
"""
"""
function get_dark(runp::String, flt_name::String,darkfolder::String="Dark")
	path=get_image_path(runp,darkfolder,flt_name)
	imgdf = DataFrame(CSV.File(path, header=false,delim="\t"))
	df1=Matrix(imgdf)
end

# ╔═╡ e5d1f8a3-bf9f-49c1-8b0d-5e745f6c5f34
"""
Given a ROI (eg, a matrix representing a section of an image) it returns an histogram of the signal in the ROI
"""
function histo_signal(iroi::Matrix{Float64}, nbin::Int=100)
	vroi = vec(iroi)
	mxvroi = maximum(vroi)
	mnvroi = minimum(vroi)
	lbl.BoldLab.hist1d(vroi, "signal in roi", nbin, mnvroi,mxvroi)
end


# ╔═╡ d53e52dc-bd7e-454e-8426-9c0dcfb5ac42
begin
	fltns = flt_names(runp)
	DRK = [get_dark(runp,fltn) for fltn in fltns]
	DRK0 = mean(DRK)  # mean pixel by pixel 
	AVG = [mean(Float64.(drk)) for drk in DRK]
	STD = [std(drk) for drk in DRK]
	dkavg = mean(AVG)
	HDRK = [histo_signal(drk) for drk in DRK]
	PDRK = [HDRK[i][2] for i in 1:length(HDRK)]

	dkstd = mean(STD)
	drk0 = mean([drk .- dkavg for drk in DRK])
end

# ╔═╡ 15fa1119-f331-4cd8-825b-a878d4210ac7
heatmap((DRK[7]))

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
# ╠═d81cc134-1a4e-45a6-8657-f82372f8b66b
# ╠═4c92771a-8684-46ab-856a-627ac859a1fc
# ╠═f855c5a4-bf43-4247-b45f-c69957575ce0
# ╠═a2acebd8-c197-4418-855f-cb08ee24021b
# ╠═5e02e22e-7c05-4af2-b7d5-f5c61f2cee11
# ╠═025d7b0e-241a-4cf5-bf1a-47bc32a7e955
# ╠═a01f191e-5f4a-4a56-9a39-b0494f02a0cd
# ╠═76383212-4bb1-425f-b82c-6dca2f6db236
# ╠═755a1675-906d-417e-a974-a6a9520f879d
# ╠═9f30ed77-3d1a-474a-9659-7683ee429b03
# ╠═a921505a-33c6-4d13-9ac9-f2fb7c1a7b02
# ╠═16d5b19e-fdcb-4c44-9a92-d7fe8c0390df
# ╠═4a08a8db-f86d-49c4-961e-bcce3ec41653
# ╠═7565e604-d434-4ba2-9594-f8d328dcca0f
# ╠═adf99e37-bb92-4d6d-8a15-8dd159af31da
# ╠═c45b108e-62f9-473a-9975-eec4736d5df1
# ╠═da7dcbd8-75c3-42d5-9958-5ccbcc9624f1
# ╠═b9361344-165b-460c-85c6-0945f40612ef
# ╠═c9e8c0f2-2776-43b5-87c6-c6c85e264924
# ╠═d5c5fbca-ae04-4ad7-a530-5f6e42f3e436
# ╠═f90a27b0-d729-4a9a-afe6-a713c090f467
# ╠═fd7b867f-2d49-4036-9b6f-0bb6966c32e6
# ╠═b26174ef-ebb4-4fe3-a93d-d308d49488aa
# ╠═0122a81d-ee51-4355-8ddb-6429f15e8ea1
# ╠═54faa211-0796-4899-b56a-90d0ea73ce4a
# ╠═b0bf7aed-3234-44ca-ac05-023545ba0987
# ╠═b896d27f-1a5f-4c25-86fa-4717bf94770b
# ╠═838eafb7-5c5b-4df5-98b1-dfbc94e58d1e
# ╠═9646801d-da05-4101-b453-8d12ed570475
# ╠═ff92b592-9cdd-4b8c-a256-f09d166febc6
# ╠═d53e52dc-bd7e-454e-8426-9c0dcfb5ac42
# ╠═15fa1119-f331-4cd8-825b-a878d4210ac7
# ╠═853b8b2e-66ed-4723-8884-213e5fd4a0e7
# ╠═71ed6fbd-0342-4cf5-9d8c-aa8f791d85f1
# ╠═a2ce4e48-73d0-491d-9de0-86a4409d358a
# ╠═f68c5806-6c44-4a13-8811-a193d59e45ba
# ╠═68eaa46d-94c0-43d9-b9e3-9bd84db26b63
# ╠═20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
# ╠═b9970588-422f-461f-addb-5169d2e6043e
# ╠═6d108da6-353f-47dc-8684-3bea555e3921
# ╠═411b89b3-178f-4e4a-bc83-af42e489b3eb
# ╠═100fa551-fe67-4cba-b9e5-77676c2c2c6f
# ╠═4d478749-935d-40a5-9431-0565ffa19e11
# ╠═09dc9013-eec5-42ea-8085-87ccc71a54aa
# ╠═28885012-f552-4ded-bfe4-c2507a685146
# ╠═1a0727c6-ea80-4026-a24a-9682737b90ec
# ╠═e5d1f8a3-bf9f-49c1-8b0d-5e745f6c5f34
# ╠═0bf5c4a7-b369-4438-a987-253f242dd88f
