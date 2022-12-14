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
cmds = filter(name->occursin("CMOS",name), lbl.BoldLab.return_files(bcmdir))
md""" Select CMOS dir : $(@bind scmos Select(cmds))"""
end

# ╔═╡ c9e8c0f2-2776-43b5-87c6-c6c85e264924
begin
	cmdir=joinpath(bcmdir,scmos)
	odir=joinpath(bodir,scmos)
	subs = lbl.BoldLab.return_files(cmdir)
	md""" Select substrate : $(@bind ssub Select(subs))"""
end

# ╔═╡ d5c5fbca-ae04-4ad7-a530-5f6e42f3e436
begin
	subsp=joinpath(cmdir,ssub)
	subsop=joinpath(odir,ssub)
	exps = lbl.BoldLab.return_files(subsp)
	md""" Select experiment : $(@bind sexp Select(exps))"""
end

# ╔═╡ f90a27b0-d729-4a9a-afe6-a713c090f467
begin
	expp=joinpath(subsp,sexp)
	expop=joinpath(subsop,sexp)
	runs = lbl.BoldLab.return_files(expp)
	md""" Select run : $(@bind srun Select(runs))"""
end

# ╔═╡ fd7b867f-2d49-4036-9b6f-0bb6966c32e6
md"""
runp and runop are the input and output data folders
"""

# ╔═╡ b26174ef-ebb4-4fe3-a93d-d308d49488aa
begin
runp=joinpath(expp,srun)
runop=joinpath(expop,srun)
end

# ╔═╡ 54faa211-0796-4899-b56a-90d0ea73ce4a
md"""
# Position of points 
- Positions folder contains measurements of the poistions of each point. 
"""

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

# ╔═╡ 3c22d436-e9b3-45a1-a87d-2d2b08580561
md"""
# Dark analysis
"""

# ╔═╡ 38080120-e90e-4d09-8328-e54e9dd7e9ff
md"""
## Mean of darks
"""

# ╔═╡ 69c529a7-7582-4830-9432-ff6e117ac0b6
md"""
Histograms of individual darks
"""

# ╔═╡ 0988bac4-9faa-40aa-ba73-e2538ad3b6e3
md"""
Histograms of the darks-mean_dark
"""

# ╔═╡ f07f9cf5-39e3-4464-98fd-ec3c7f58e235
md"""
## Stability of darks
"""

# ╔═╡ a021d53f-c408-49d6-bb3c-6f8758c5e183
md"""
# Single point
"""

# ╔═╡ 4dd6f397-0cd3-4c48-a8a2-20984418fb6f
begin
	pnts=lbl.BoldLab.point_names(runp)
	md""" Select point : $(@bind pnt Select(pnts))"""
end

# ╔═╡ f3669d36-cfd5-40b7-a09e-d73e8159b031
md"""
## Non-filtered image 
"""

# ╔═╡ 34c9c12c-2f39-4d4c-bc99-74b93304d724
begin
	nfim=lbl.BoldLab.get_nfimage(runp,pnt)
	pp0=heatmap(nfim)
end

# ╔═╡ 49d52be1-0803-49f7-8a67-6e46f847aa83
md"""
### Define ROI
"""

# ╔═╡ 9084f5bc-6b3a-412a-8d90-0dbf4df74ee3
begin
img_edge,ROI=lbl.BoldLab.Image_edge(nfim)
pp1=heatmap(img_edge)
end

# ╔═╡ dfeb71a9-d39e-4ac7-a3b2-fc954e0c4860
plot(pp0,pp1,size=(1000,400))

# ╔═╡ 71eab91d-df5c-405f-b8af-cea1de77a7c5
md"""
## Filtered image
"""

# ╔═╡ 966c7a16-ef71-49ff-a16a-9d2a9c22731e
begin
	flts=lbl.BoldLab.flt_names(runp)
	md""" Select filter : $(@bind flt Select(flts))"""
end

# ╔═╡ c36851fc-1bd2-4e5d-baad-d97555642850
md"""
### Dark
"""

# ╔═╡ 6767a952-97ba-49f2-abcc-a812406ea702
md"""
- Left: Dark image
- Right: Histogram of dark image
"""

# ╔═╡ b04549e8-2a06-4bf3-a8d7-9c48e352a1a7
begin
	drk=lbl.BoldLab.get_dark(runp,flt,"Dark")
	vlsdrk=[vl for vl in vec(drk)]
	histdrk=stephist(vlsdrk,bins=:100,yaxis=:log)
	heatdrk=heatmap(drk)
	plot(heatdrk,histdrk,size=(1000,400))
end

# ╔═╡ e06c0eee-beee-4137-acae-fb3ce6201268
md"""
Histogram of dark in ROI
"""


# ╔═╡ 447f648f-63ec-49bf-9f30-2624e86d5ff6
begin
	sdrk=[drk[i,j] for (i,j) in ROI]
	histdrkROI=stephist(sdrk,bins=:100,yaxis=:log,size=(1000,400))
end

# ╔═╡ b6a43048-7f11-43c2-9f6a-3a462e5f1299
md"""
### Image
"""

# ╔═╡ a8603145-c510-4de5-bb3d-ba3409c3f8ac
md"""
- Left: Image
- Right: Histogram of image
"""

# ╔═╡ 6f1141c9-908d-40dc-ab2b-e3bbfdcd74cb
begin
	im=lbl.BoldLab.get_image(runp,pnt,flt)
	vls=[vl for vl in vec(im)]
	histim=stephist(vls,bins=:100,yaxis=:log)
	heatim=heatmap(im)
	plot(heatim,histim,size=(1000,400))
end

# ╔═╡ 5644b297-d344-45ec-a801-1aa8d1c428a1
md"""
Histogram of image in ROI
"""

# ╔═╡ 95bea89c-720e-47a1-aa7e-a11300d585f5
begin
	sim=[im[i,j] for (i,j) in ROI]
	histimROI=stephist(sim,bins=:100,yaxis=:log,size=(1000,400))
end

# ╔═╡ 704a83dd-3914-4803-a08d-6528892263b6
md"""
### Image-Dark
"""

# ╔═╡ 019ab3d7-a98d-4277-87a2-1f8a67e90b80
md"""
- Left: Image-Dark
- Right: Histogram of image-dark
"""

# ╔═╡ d64ff72f-7e99-48d8-8dcc-3c69c5f496b3
begin
	im_drk=im-drk
	vlsim_drk=[vl for vl in vec(im_drk)]
	histim_drk=stephist(vlsim_drk,bins=:100,yaxis=:log)
	heatim_drk=heatmap(im_drk)
	plot(heatim_drk,histim_drk,size=(1000,400))
end

# ╔═╡ bac91759-91d2-44aa-ad6c-7281a4f4d8be
begin
	sim_drk=[im_drk[i,j] for (i,j) in ROI]
	histim_drkROI=stephist(sim_drk,bins=:100,yaxis=:log,size=(1000,400))
end

# ╔═╡ 967e407f-9919-4497-822e-a1f428c1e61f
md"""
## Spectrum
"""

# ╔═╡ 1b401135-decc-4ca3-85ae-139a55c9ad64
md"""
# Spectrum for all points
"""

# ╔═╡ 27be6cae-4c5d-4f0e-bb0e-259adb4ade59
# ╠═╡ disabled = true
#=╠═╡
md"""
## Total charge in each filter
"""
  ╠═╡ =#

# ╔═╡ f6329740-963a-41b7-8a3a-d08ead791aba
md"""
## Total charge in each filter / filter width
"""

# ╔═╡ 853b8b2e-66ed-4723-8884-213e5fd4a0e7
md"""
# Tests
"""

# ╔═╡ 20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
md"""
# Functions
"""

# ╔═╡ b45e7a24-4906-40ac-8468-6d40e10663b8
function sum_signal(runp,pnt,flt)
	drk=lbl.BoldLab.get_dark(runp,flt)
	im=lbl.BoldLab.get_image(runp,pnt,flt)
	im_drk=im-drk
	img_edge,ROI=lbl.BoldLab.Image_edge(nfim)
	s=[im_drk[i,j] for (i,j) in ROI]
	sumation=sum(s)
end

# ╔═╡ ac4670b5-a160-4207-8842-c029f38d4470
@test typeof(sum_signal(runp,pnts[1],flts[1]))==Int64

# ╔═╡ a8739f1f-6cd9-4be4-ae17-037fa3fc048f
function signal_flts(runp,pnt)
	nfim=lbl.BoldLab.get_nfimage(runp,pnt)
	flts=lbl.BoldLab.flt_names(runp)
	store=Vector{Float64}()
	for flt in flts
		s=sum_signal(runp,pnt,flt)
		push!(store,s)
	end
	store
end

# ╔═╡ 1021f5b4-fb1a-48ea-a560-3488697ec045
begin
p0=plot(xfnm,signal_flts(runp,pnts[1]),title="Total charge",xlabel="wavelength (nm)")
scatter!(xfnm,signal_flts(runp,pnts[1]))
p1=plot(xfnm,signal_flts(runp,pnts[1])./wfnm,title="Total charge/filter width",xlabel="wavelength (nm)")
scatter!(xfnm,signal_flts(runp,pnts[1])./wfnm)
plot(p0,p1,size=(1000,500))
end

# ╔═╡ 8889f58d-fb86-4b55-a36e-d6a59ac94aa1
# ╠═╡ disabled = true
#=╠═╡
begin
plots0=[]
plots1=[]
for pnt in pnts
	s=signal_flts(runp,pnt)
	p0=plot(xfnm,s,xlabel="wavelength (nm)")
	scatter!(xfnm,s)
	push!(plots0,p0)
	p1=plot(xfnm,s./wfnm,xlabel="wavelength (nm)")
	scatter!(xfnm,s./wfnm)
	push!(plots1,p1)
end
end
  ╠═╡ =#

# ╔═╡ e767fdda-bc58-423d-ab32-3dc64f80f681
# ╠═╡ disabled = true
#=╠═╡
begin
plot(plots0...,size=(1000,1000))
title!("Total charge")
end
  ╠═╡ =#

# ╔═╡ a0ab96ff-d512-436a-9510-4c9eb61aecf0
# ╠═╡ disabled = true
#=╠═╡
begin
plot(plots1...,size=(1000,1000))
title!("Total charge")
end
  ╠═╡ =#

# ╔═╡ b90b15df-cf24-45f7-81a9-38598ca93fa2
[@test size(signal_flts(runp,pnt))==size(flts) for pnt in pnts]

# ╔═╡ 4b175426-4d2a-4e23-bdc2-32ca92be3bcc
function dark_mean(runp,darkfolder::String="Dark")
	flts=lbl.BoldLab.flt_names(runp)
	drk0=lbl.BoldLab.get_dark(runp,flts[1])
	sum=zero(drk0)
	for flt in flts
		sum+=lbl.BoldLab.get_dark(runp,flt,darkfolder)
	end
	sum/size(flts)[1]
end
		
	

# ╔═╡ 06885865-c1bf-46a0-89eb-63e73f62fef5
begin
drkm=dark_mean(runp)
drkmheat=heatmap(drkm)
drkmhist=stephist(vec(drkm),bins=:1000,yaxis=:log)
plot(drkmheat,drkmhist,size=(1000,400))
end

# ╔═╡ 40f3dbff-47ed-434f-b3a3-7f98d9f1d7e7
begin
drkhists0=[]
drkhists1=[]
drkms=[]
drkstds=[]
for flt in flts
	drk=lbl.BoldLab.get_dark(runp,flt)
	drkms=push!(drkms,mean(drk))
	drkstds=push!(drkstds,std(drk))
	drkhist0=stephist(vec(drk),bins=:1000,yaxis=:log,title=flt)
	push!(drkhists0,drkhist0)
	drkhist1=stephist(vec(drk-drkm),bins=:1000,yaxis=:log,title=flt)
	push!(drkhists1,drkhist1)
end
end

# ╔═╡ 59978ac0-b955-437a-9378-54944119898b
begin
plot(drkhists0...,size=(1000,1000))
end

# ╔═╡ 510d5320-0e3f-4fd0-a096-4b25ffca84a5
begin
plot(drkhists1...,size=(1000,1000))
xlims!(-200,200)
end

# ╔═╡ 95b7485e-3b8e-469f-8697-55b43847f663
begin
	drkmsp=plot(drkms)
	scatter!(drkms)
	drkstdsp=plot(drkstds)
	scatter!(drkstds)
	plot(drkmsp,drkstdsp)
end

# ╔═╡ 18aac698-2fe7-4434-86ed-eee1bceb16e8
@test typeof(dark_mean(runp))==Matrix{Float64}

# ╔═╡ f78784d2-4fe5-4dbc-80a7-33d44718c5a7
@test 1500<mean(dark_mean(runp))<1700

# ╔═╡ Cell order:
# ╠═80a4e1cc-4d59-4c8a-9d65-1fe0d2d77bf2
# ╠═23b4e644-8331-4aa6-8369-47ce3ff0e143
# ╠═3a7182c2-8ffd-4d6a-975f-42f9ff1b0744
# ╠═9f71bc31-df8f-48c5-815a-c8db9e46100f
# ╠═19035c0d-93a7-4bb9-be05-d2a9b9ac4619
# ╠═f4e379c4-a2b2-4703-bcbc-f6d7c996354a
# ╠═51710bb8-9c4a-4abc-9fe9-02e87bd4e4c5
# ╟─d1ace15e-fe1a-4552-9144-a0824ae8ae0f
# ╟─0937a6fc-6936-47d0-80de-8a38bb9a6a37
# ╠═6d04e8fe-c174-4be5-bffc-40945e8074e5
# ╟─8b4fa5c3-0efe-4aa7-98cf-d37c1e12bc74
# ╟─e6de7c36-135f-43f6-a0e2-3eb65fd6cf57
# ╠═7789a654-af1b-408a-86fc-26cc8ff221da
# ╠═d81cc134-1a4e-45a6-8657-f82372f8b66b
# ╠═4c92771a-8684-46ab-856a-627ac859a1fc
# ╠═f855c5a4-bf43-4247-b45f-c69957575ce0
# ╟─a2acebd8-c197-4418-855f-cb08ee24021b
# ╠═5e02e22e-7c05-4af2-b7d5-f5c61f2cee11
# ╠═025d7b0e-241a-4cf5-bf1a-47bc32a7e955
# ╠═a01f191e-5f4a-4a56-9a39-b0494f02a0cd
# ╠═76383212-4bb1-425f-b82c-6dca2f6db236
# ╠═755a1675-906d-417e-a974-a6a9520f879d
# ╠═9f30ed77-3d1a-474a-9659-7683ee429b03
# ╠═a921505a-33c6-4d13-9ac9-f2fb7c1a7b02
# ╠═16d5b19e-fdcb-4c44-9a92-d7fe8c0390df
# ╟─4a08a8db-f86d-49c4-961e-bcce3ec41653
# ╟─7565e604-d434-4ba2-9594-f8d328dcca0f
# ╠═adf99e37-bb92-4d6d-8a15-8dd159af31da
# ╟─c45b108e-62f9-473a-9975-eec4736d5df1
# ╟─da7dcbd8-75c3-42d5-9958-5ccbcc9624f1
# ╟─b9361344-165b-460c-85c6-0945f40612ef
# ╟─c9e8c0f2-2776-43b5-87c6-c6c85e264924
# ╟─d5c5fbca-ae04-4ad7-a530-5f6e42f3e436
# ╟─f90a27b0-d729-4a9a-afe6-a713c090f467
# ╟─fd7b867f-2d49-4036-9b6f-0bb6966c32e6
# ╟─b26174ef-ebb4-4fe3-a93d-d308d49488aa
# ╠═54faa211-0796-4899-b56a-90d0ea73ce4a
# ╠═b0bf7aed-3234-44ca-ac05-023545ba0987
# ╠═3c22d436-e9b3-45a1-a87d-2d2b08580561
# ╠═38080120-e90e-4d09-8328-e54e9dd7e9ff
# ╠═06885865-c1bf-46a0-89eb-63e73f62fef5
# ╠═40f3dbff-47ed-434f-b3a3-7f98d9f1d7e7
# ╟─69c529a7-7582-4830-9432-ff6e117ac0b6
# ╠═59978ac0-b955-437a-9378-54944119898b
# ╟─0988bac4-9faa-40aa-ba73-e2538ad3b6e3
# ╠═510d5320-0e3f-4fd0-a096-4b25ffca84a5
# ╠═f07f9cf5-39e3-4464-98fd-ec3c7f58e235
# ╠═95b7485e-3b8e-469f-8697-55b43847f663
# ╠═a021d53f-c408-49d6-bb3c-6f8758c5e183
# ╟─4dd6f397-0cd3-4c48-a8a2-20984418fb6f
# ╟─f3669d36-cfd5-40b7-a09e-d73e8159b031
# ╠═34c9c12c-2f39-4d4c-bc99-74b93304d724
# ╟─49d52be1-0803-49f7-8a67-6e46f847aa83
# ╠═9084f5bc-6b3a-412a-8d90-0dbf4df74ee3
# ╠═dfeb71a9-d39e-4ac7-a3b2-fc954e0c4860
# ╟─71eab91d-df5c-405f-b8af-cea1de77a7c5
# ╟─966c7a16-ef71-49ff-a16a-9d2a9c22731e
# ╟─c36851fc-1bd2-4e5d-baad-d97555642850
# ╠═6767a952-97ba-49f2-abcc-a812406ea702
# ╠═b04549e8-2a06-4bf3-a8d7-9c48e352a1a7
# ╠═e06c0eee-beee-4137-acae-fb3ce6201268
# ╠═447f648f-63ec-49bf-9f30-2624e86d5ff6
# ╟─b6a43048-7f11-43c2-9f6a-3a462e5f1299
# ╠═a8603145-c510-4de5-bb3d-ba3409c3f8ac
# ╠═6f1141c9-908d-40dc-ab2b-e3bbfdcd74cb
# ╟─5644b297-d344-45ec-a801-1aa8d1c428a1
# ╠═95bea89c-720e-47a1-aa7e-a11300d585f5
# ╠═704a83dd-3914-4803-a08d-6528892263b6
# ╠═019ab3d7-a98d-4277-87a2-1f8a67e90b80
# ╠═d64ff72f-7e99-48d8-8dcc-3c69c5f496b3
# ╠═bac91759-91d2-44aa-ad6c-7281a4f4d8be
# ╟─967e407f-9919-4497-822e-a1f428c1e61f
# ╠═1021f5b4-fb1a-48ea-a560-3488697ec045
# ╠═1b401135-decc-4ca3-85ae-139a55c9ad64
# ╠═8889f58d-fb86-4b55-a36e-d6a59ac94aa1
# ╠═27be6cae-4c5d-4f0e-bb0e-259adb4ade59
# ╠═e767fdda-bc58-423d-ab32-3dc64f80f681
# ╠═f6329740-963a-41b7-8a3a-d08ead791aba
# ╠═a0ab96ff-d512-436a-9510-4c9eb61aecf0
# ╠═853b8b2e-66ed-4723-8884-213e5fd4a0e7
# ╠═ac4670b5-a160-4207-8842-c029f38d4470
# ╠═b90b15df-cf24-45f7-81a9-38598ca93fa2
# ╠═18aac698-2fe7-4434-86ed-eee1bceb16e8
# ╠═f78784d2-4fe5-4dbc-80a7-33d44718c5a7
# ╠═20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
# ╠═b45e7a24-4906-40ac-8468-6d40e10663b8
# ╠═a8739f1f-6cd9-4be4-ae17-037fa3fc048f
# ╠═4b175426-4d2a-4e23-bdc2-32ca92be3bcc
