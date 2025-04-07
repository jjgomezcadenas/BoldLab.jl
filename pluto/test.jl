### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
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
LaserLabp=ENV["JLaserData"]

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
g2slp=joinpath(LaserLabp,"FLUORI/G2/G2SL/BOLD_104_SIL_GB_onquartz")

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
bodir=joinpath(LaserLabp,"pdata")
bcmdir=joinpath(LaserLabp,"data")
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

# ╔═╡ 54faa211-0796-4899-b56a-90d0ea73ce4a
md"""
# Position of points 
- Positions folder contains measurements of the poistions of each point. 
"""

# ╔═╡ a021d53f-c408-49d6-bb3c-6f8758c5e183
md"""
# Image visualization
"""

# ╔═╡ f3669d36-cfd5-40b7-a09e-d73e8159b031
md"""
## Filter1 
"""

# ╔═╡ 71eab91d-df5c-405f-b8af-cea1de77a7c5
md"""
## Filtered image
"""

# ╔═╡ c36851fc-1bd2-4e5d-baad-d97555642850
md"""
### Dark
"""

# ╔═╡ b6a43048-7f11-43c2-9f6a-3a462e5f1299
md"""
### Image
"""

# ╔═╡ 49d52be1-0803-49f7-8a67-6e46f847aa83
md"""
# Clusterization
"""

# ╔═╡ 853b8b2e-66ed-4723-8884-213e5fd4a0e7
md"""
# Tests
"""

# ╔═╡ 7d759e14-e9a9-440c-ac8a-de92ff312526
@test ("data" in readdir(ENV["JLaserData"]))

# ╔═╡ 2c1ca8f4-ed1a-46d1-9a34-1e76c9005a87
@test ("pdata" in readdir(ENV["JLaserData"]))

# ╔═╡ 905ca07b-5778-4d8a-815d-7a8bdd4b73d4
@test ("FLUORI" in readdir(ENV["JLaserData"]))

# ╔═╡ a2ce4e48-73d0-491d-9de0-86a4409d358a
@test typeof(return_files(runp,"Dark"))==Tuple{Vector{String}, Vector{String}}

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

# ╔═╡ 4dd6f397-0cd3-4c48-a8a2-20984418fb6f
begin
	pnts=point_names(runp)
	md""" Select point : $(@bind pnt Select(pnts))"""
end

# ╔═╡ 34c9c12c-2f39-4d4c-bc99-74b93304d724
begin
	nfim=get_nfimage(runp,pnt)
	heatmap(nfim)
end

# ╔═╡ 966c7a16-ef71-49ff-a16a-9d2a9c22731e
begin
	flts=flt_names(runp)
	md""" Select filter : $(@bind flt Select(flts))"""
end

# ╔═╡ b04549e8-2a06-4bf3-a8d7-9c48e352a1a7
begin
	drk=get_dark(runp,flt)
	heatmap(drk)
end

# ╔═╡ 6f1141c9-908d-40dc-ab2b-e3bbfdcd74cb
begin
	im=get_image(runp,pnt,flt)
	heatmap(im)
end

# ╔═╡ 9084f5bc-6b3a-412a-8d90-0dbf4df74ee3
begin
	imn=Float64.(im./maximum(im))
	img_edge=Float64.(sujoy(imn,four_connectivity=true))
	img_edgeb=binarize(img_edge,Otsu())
	iedge = Tuple.(findall(x->x==1,img_edgeb))  #indexed of the edge
	medge, xedge, yedge = permutedims(hcat(first.(iedge), last.(iedge))),first.(iedge), last.(iedge)  # edge expressed as matrix, x,y
	heatmap(img_edgeb)
end

# ╔═╡ f68c5806-6c44-4a13-8811-a193d59e45ba
[@test occursin("Filter",name) for name in flt_names(runp)]

# ╔═╡ abca5b48-0ee3-4085-b233-3154f405dd90
[@test occursin("Point",name) for name in point_names(runp)]

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

# ╔═╡ 2e19ae1a-8bff-4d3b-a406-cc710370029f
"""
"""
function sujoy(img; four_connectivity=true)
    img_channel = Gray.(img)

    min_val = minimum(img_channel)
    img_channel = img_channel .- min_val
    max_val = maximum(img_channel)

    if max_val == 0
        return img
    end

    img_channel = img_channel./max_val

    if four_connectivity
        krnl_h = centered(Gray{Float32}[0 -1 -1 -1 0; 0 -1 -1 -1 0; 0 0 0 0 0; 0 1 1 1 0; 0 1 1 1 0]./12)
        krnl_v = centered(Gray{Float32}[0 0 0 0 0; -1 -1 0 1 1;-1 -1 0 1 1;-1 -1 0 1 1;0 0 0 0 0 ]./12)
    else
        krnl_h = centered(Gray{Float32}[0 0 -1 0 0; 0 -1 -1 -1 0; 0 0 0 0 0; 0 1 1 1 0; 0 0 1 0 0]./8)
        krnl_v = centered(Gray{Float32}[0 0 0 0 0;  0 -1 0 1 0; -1 -1 0 1 1;0 -1 0 1 0; 0 0 0 0 0 ]./8)
    end

    grad_h = imfilter(img_channel, krnl_h')
    grad_v = imfilter(img_channel, krnl_v')

    grad = (grad_h.^2) .+ (grad_v.^2)

    return grad
end

# ╔═╡ Cell order:
# ╠═80a4e1cc-4d59-4c8a-9d65-1fe0d2d77bf2
# ╠═23b4e644-8331-4aa6-8369-47ce3ff0e143
# ╠═3a7182c2-8ffd-4d6a-975f-42f9ff1b0744
# ╠═9f71bc31-df8f-48c5-815a-c8db9e46100f
# ╠═19035c0d-93a7-4bb9-be05-d2a9b9ac4619
# ╠═f4e379c4-a2b2-4703-bcbc-f6d7c996354a
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
# ╠═54faa211-0796-4899-b56a-90d0ea73ce4a
# ╠═b0bf7aed-3234-44ca-ac05-023545ba0987
# ╠═a021d53f-c408-49d6-bb3c-6f8758c5e183
# ╠═4dd6f397-0cd3-4c48-a8a2-20984418fb6f
# ╠═f3669d36-cfd5-40b7-a09e-d73e8159b031
# ╠═34c9c12c-2f39-4d4c-bc99-74b93304d724
# ╠═71eab91d-df5c-405f-b8af-cea1de77a7c5
# ╠═966c7a16-ef71-49ff-a16a-9d2a9c22731e
# ╠═c36851fc-1bd2-4e5d-baad-d97555642850
# ╠═b04549e8-2a06-4bf3-a8d7-9c48e352a1a7
# ╠═b6a43048-7f11-43c2-9f6a-3a462e5f1299
# ╠═6f1141c9-908d-40dc-ab2b-e3bbfdcd74cb
# ╠═49d52be1-0803-49f7-8a67-6e46f847aa83
# ╠═9084f5bc-6b3a-412a-8d90-0dbf4df74ee3
# ╠═853b8b2e-66ed-4723-8884-213e5fd4a0e7
# ╠═7d759e14-e9a9-440c-ac8a-de92ff312526
# ╠═2c1ca8f4-ed1a-46d1-9a34-1e76c9005a87
# ╠═905ca07b-5778-4d8a-815d-7a8bdd4b73d4
# ╠═a2ce4e48-73d0-491d-9de0-86a4409d358a
# ╠═f68c5806-6c44-4a13-8811-a193d59e45ba
# ╠═abca5b48-0ee3-4085-b233-3154f405dd90
# ╠═20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
# ╠═b9970588-422f-461f-addb-5169d2e6043e
# ╠═6d108da6-353f-47dc-8684-3bea555e3921
# ╠═411b89b3-178f-4e4a-bc83-af42e489b3eb
# ╠═100fa551-fe67-4cba-b9e5-77676c2c2c6f
# ╠═4d478749-935d-40a5-9431-0565ffa19e11
# ╠═09dc9013-eec5-42ea-8085-87ccc71a54aa
# ╠═28885012-f552-4ded-bfe4-c2507a685146
# ╠═1a0727c6-ea80-4026-a24a-9682737b90ec
# ╠═2e19ae1a-8bff-4d3b-a406-cc710370029f
