### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 3a4634d0-063c-11f0-2c98-ed94f1170f79
using Pkg; Pkg.activate(ENV["JBoldLab"])

# ╔═╡ fb1cc299-6238-445e-99be-9c167823417c
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Plots 
	using Measures
	using Printf
	using InteractiveUtils
	using LsqFit
	using Statistics
	using StatsBase
	using DelimitedFiles
end

# ╔═╡ 4cbeae67-5ee0-4e33-98dd-4abb98874786
begin
	rdir = "histogramsZFFL138/"
	state = ["Free", "Chelated"]
	samples = [string(x) for x in range(1,8)]
	blanks = ["Blank1", "Blank2"]
	tsamples =vcat(samples, blanks)
	fields = ["Field1", "Field2", "Field3"]
end

# ╔═╡ 7693af9f-b4ca-4127-93c0-771f3aad54c1
tsamples

# ╔═╡ c185fa05-d717-42be-81df-9b6f0e30379b
function getfiles(rdir, state, sample, field)
	path = joinpath(rdir, state, sample, field)
	#fp1 = "histogramsZFFL138/Free/1/Field1/Points.csv"
	fss = joinpath(path,"Stepsizes.csv")
	fst = joinpath(path,"Steptimes.csv")
	fss, fst
end

# ╔═╡ 56d4de31-efc1-4b34-9990-dd0b4ff54073
function read_poinst(fname)
	# Read the file, filtering out comment lines (lines starting with '#')
	lines = filter(line -> !startswith(line, "#"), 
				   readlines(fname))

	# Replace multiple whitespace characters with a single space in each line
	cleaned_lines = map(line -> replace(line, r"\s+" => " "), lines)

	# Join the cleaned lines into a single string separated by newlines
	cleaned_text = join(cleaned_lines, "\n")
	
	# Use readdlm with a single space as the delimiter
	data = readdlm(IOBuffer(cleaned_text), ' ')
	
	# Convert the resulting array into a DataFrame (columns are auto-named)
	df = DataFrame(data, :auto)
	df = select(df, Not(:x1))
	# Display the DataFrame
	df
end


# ╔═╡ be9c2108-5b23-48f2-8ff3-3edaf899568e
function readss(fname)
	# Read all lines from the file
	lines = readlines(fname)
	
	# Filter out the comment line(s) and parse each remaining line to a Float64
	ss = [parse(Float64, strip(line)) for line in lines if !startswith(strip(line), "#")]
	
	# Display the resulting array of numbers
	ss
end

# ╔═╡ f7bb460d-0cfc-40bb-9dc6-56b4636530a4
function getplot_ss_st(rdir, state, sample, fields; 
					   smin=-400.0, smax=400.0, nbinss=30, nbinst=30)
	
	fstep(X) = filter(x -> smin <= x <= smax, X)

	
	SS = Float64[]
	ST = Float64[]

	for i in range(1,3)
		fss, fst = getfiles(rdir, state, sample, fields[i])
		ss =readss(fss)
		st =readss(fst) 
		append!(SS, fstep(ss))
		append!(ST, st)
	end
		
	lgnd = string(state, "_", sample)
	pss = histogram(SS, bins=nbinss, title="Step Sizes", 
				   xlabel="Step size", ylabel="# steps", label=lgnd)
	pst = histogram(ST, bins=nbinst, title="Step Times", xlabel="# steps", 
				   ylabel="# events", label=lgnd)
	hss = fit(Histogram, SS, nbins=nbinss)
	hst = fit(Histogram, ST, nbins=nbinst)
	pss, pst, hss, hst
end

# ╔═╡ 909c89fa-5383-4db7-ac2f-9d4e2590b602
function plot_ss_st(rdir, state, sample, fields; 
					smin=-400.0, smax=400.0, nbinss=30, nbinst=30)
	
	pss, pst, hss, hst = getplot_ss_st(rdir, state, sample, fields; 
									   smin, smax, nbinss, nbinst)
	plot(pss, pst, layout=(1,2), margin=2mm, spacing=2mm)
end

# ╔═╡ 906abfe8-24fe-40f7-abfc-9ce58a717430
function histos_ss_st(rdir, state, sample, fields; smin=-400.0, smax=400.0,
					  nbinss=30, nbinst=30)
	pss, pst, hss, hst = getplot_ss_st(rdir, state, sample, fields; smin, smax,
									   nbinss, nbinst)
	hss, hst
end

# ╔═╡ fd3a641a-78eb-4cb9-a495-4d487b02b834
function plot_samples(rdir, state, samples, fields; 
					  smin=-400.0, smax=400.0, nbinss=30, nbinst=30)
	PSS=[]
	PST=[]
	for i in range(1,length(samples))
		if i == 2
			continue
		end
		pss, pst, hss, hst = getplot_ss_st(rdir, state, samples[i], fields; 
										   smin, smax, nbinss, nbinst)
		push!(PSS, pss)
		push!(PST, pst)
	end

	PSS, PST
		
end

# ╔═╡ dc92f1d9-93a4-430a-bf72-334da4cbe2ba
function get_histos(rdir, state, samples, fields; smin=-400.0, smax=400.0,
				    nbinss=30, nbinst=30)
	PSS=[]
	PST=[]
	for i in range(1,length(samples))
		if i==2
			continue
		end
		hss, hst =histos_ss_st(rdir, state, samples[i], fields; smin, smax,
							   nbinss, nbinst)
		push!(PSS, hss)
		push!(PST, hst)
	end

	PSS, PST
		
end

# ╔═╡ 8933c28d-be0e-482b-b72c-b83798178fff
function get_ss_st(rdir, state, blanks, fields; smin=-400.0, smax=400.0,
				   nbinss=30, nbinst=30)
	hsb, hst =histos_ss_st(rdir, state, blanks, fields; smin, smax, nbinss, nbinst)
	mxss, iss = findmax(hsb.weights)
	eds = collect(hsb.edges[1])
	mxst, ist = findmax(hst.weights)
	edt = collect(hst.edges[1])
	mxss, eds[iss], mxst, edt[ist], hst.weights[1]
end

# ╔═╡ ebb1bdf6-9eac-41cd-8513-caae11875c2a
function get_ss_st_stats(rdir, state, samples, fields; smin=-400.0, smax=400.0,
						 nbinss=30, nbinst=30)
	HSB, HST = get_histos(rdir, state, samples, fields; smin, smax, nbinss, nbinst)

	MXSS = []
	EDS = []
	MXST = []
	EDT = []
	VT0 = []
	for hsb in HSB
		mxss, iss = findmax(hsb.weights)
		eds = collect(hsb.edges[1])
		push!(MXSS, mxss)
		push!(EDS, eds[iss])
	end
	for hst in HST
		mxst, ist = findmax(hst.weights)
		edt = collect(hst.edges[1])
		push!(MXST, mxst)
		push!(EDT, edt[ist])
		push!(VT0, hst.weights[1])
	end
	MXSS, EDS, MXST, EDT, VT0
	
end

# ╔═╡ 5c74b873-5716-4a8a-848d-ff9da84d8250
function subtract_blank(mxsf, mxsc, mxtf, mxtc, vt0f, vt0c; iend=7)
	xmxsf = [m-mxsf[8] for m in mxsf[1:iend]]
	xmxsc = [m-mxsc[8] for m in mxsc[1:iend]]
	xmxtf = [m-mxtf[8] for m in mxtf[1:iend]]
	xmxtc = [m-mxtc[8] for m in mxtc[1:iend]]
	xvt0f = [m-vt0f[8] for m in vt0f[1:iend]]
	xvt0c = [m-vt0c[8] for m in vt0c[1:iend]]
	return xmxsf, xmxsc, xmxtf, xmxtc, xvt0f, xvt0c
end

# ╔═╡ 524a2bde-bc60-4463-91ce-ddbd20b6a574
md"""
## Analysis
"""

# ╔═╡ 865b5b90-15d4-4c67-94c5-6f30b650bff0
md"""
### Free: Blank 1
"""

# ╔═╡ e99aca88-68b6-4ee7-8a71-99d0a9261450
plot_ss_st(rdir, state[1], blanks[1], fields; nbinss=20, nbinst=20)

# ╔═╡ 61f5f861-60c8-4bf4-93fe-34f8233a79f0
let
	mxss, vss, mxst, vst, st0 = get_ss_st(rdir, state[1], blanks[1], fields; 
										 nbinss=20, nbinst=20)
	md"""
	- peak number of steps =$(mxss)
	- step length at peak =$(vss)
	- peak number of steps =$(mxst)
	- step time at peak =$(vst)
	- number of steps at t=0:  $(st0)
	"""
end

# ╔═╡ 1b1c5630-1e22-4fe3-be28-2c0591c3a063
md"""
### Free: Blank 2
"""

# ╔═╡ bee8c3c3-4512-403e-af20-b76ae3534a10
plot_ss_st(rdir, state[1], blanks[2], fields; nbinss=20, nbinst=20)

# ╔═╡ 5cc2980f-3181-4840-ba93-83cf5745f566
let
	mxss, vss, mxst, vst, st0 = get_ss_st(rdir, state[1], blanks[2], fields;
										 nbinss=20, nbinst=20)
	md"""
	- peak number of steps =$(mxss)
	- step length at peak =$(vss)
	- peak number of steps =$(mxst)
	- step time at peak =$(vst)
	- number of steps at t=0:  $(st0)
	"""
end

# ╔═╡ 403736b5-3474-45bb-b71d-ad7639c6dd4d
md"""
### Chelated: Blank 1
"""

# ╔═╡ 32cc64c3-9214-4b24-9d3b-23d42cafe574
plot_ss_st(rdir, state[2], blanks[1], fields; nbinss=20, nbinst=20)

# ╔═╡ 205364db-aaf3-4af6-b3d8-114ba05faf16
let
	mxss, vss, mxst, vst, st0 = get_ss_st(rdir, state[2], blanks[1], fields;
										 nbinss=20, nbinst=20)
	md"""
	- peak number of steps =$(mxss)
	- step length at peak =$(vss)
	- peak number of steps =$(mxst)
	- step time at peak =$(vst)
	- number of steps at t=0:  $(st0)
	"""
end

# ╔═╡ b55a32e3-8015-4f50-b34f-cdf2fed5559e
md"""
### Chelated: Blank 2
"""

# ╔═╡ f16ba4df-3398-4a34-b12f-bcf848aad609
plot_ss_st(rdir, state[2], blanks[2], fields; nbinss=20, nbinst=20)

# ╔═╡ 5d43c3bb-be52-4a0c-a638-0a2cf64424e7
let
	mxss, vss, mxst, vst, st0 = get_ss_st(rdir, state[2], blanks[2], fields;
										 nbinss=20, nbinst=20)
	md"""
	- peak number of steps =$(mxss)
	- step length at peak =$(vss)
	- peak number of steps =$(mxst)
	- step time at peak =$(vst)
	- number of steps at t=0:  $(st0)
	"""
end

# ╔═╡ 2c4e2d24-fc1e-4d24-a672-7c1c3b1306dc
plot_ss_st(rdir, state[1], samples[1], fields, nbinss=20, nbinst=20)

# ╔═╡ ce0bca19-1679-4431-9b49-b74171745bb6
let
	mxss, vss, mxst, vst, st0 = get_ss_st(rdir, state[1], samples[1], fields, nbinss=20, nbinst=20)
	md"""
	- peak number of steps =$(mxss)
	- step length at peak =$(vss)
	- peak number of steps =$(mxst)
	- step time at peak =$(vst)
	- number of steps at t=0:  $(st0)
	"""
end

# ╔═╡ ca368f79-bb07-49d4-9c53-83c46e69e9b7
MXSS, EDS, MXST, EDT, VT0 = get_ss_st_stats(rdir, state[1], tsamples, 
											fields; smin=-400.0, smax=400.0, nbinss=30, nbinst=30)

# ╔═╡ 2e7b7742-38dc-4dad-8dae-0fab9bcc486f
VT0

# ╔═╡ 47f153ad-3943-4314-9e38-22f4af02f017
begin
	PSS, PST = plot_samples(rdir, state[1], tsamples, fields; smin=-300.0, smax=300.0, nbinss=30, nbinst=30)
	nothing 
end

# ╔═╡ 7595358d-2264-4318-ad48-c7bbc8871300
plot(PSS..., layout=(3,3), margin=2mm, spacing=2mm, size=(800,600))

# ╔═╡ 4de31e96-5c89-46cf-81eb-6eaa1dc0712c
plot(PST..., layout=(3,3), margin=2mm, spacing=2mm, size=(800,600))

# ╔═╡ ab77a9ac-b60c-4f35-b63e-c5edd2fb1464
let
	s1 = scatter(MXSS,  xlabel="Sample number", ylabel="#steps at sx")
	s1 = plot!(s1, MXSS)
	s2 = scatter(MXST,  xlabel="Sample number", ylabel="#steps at t=max")
	s2 = plot!(s2, MXST)
	s3 = scatter(VT0,  xlabel="Sample number", ylabel="#steps at t=0")
	s3 = plot!(s3, VT0,)
	plot(s1,s2, s3, layout=(3,1), margin=2mm, spacing=2mm, size=(800,600))
end


# ╔═╡ d76192c0-ef94-4b44-8b7c-5e6468a54928
begin
	PSS2, PST2 = plot_samples(rdir, state[2], tsamples, fields; smin=-300.0, smax=300.0, nbinss=40, nbinst=40)
	nothing 
end

# ╔═╡ e236f801-bb66-4658-bada-4e5d3ef9fb2b
plot(PSS2..., layout=(3,3), margin=2mm, spacing=2mm, size=(800,600))

# ╔═╡ 6832eeaf-3ffc-46c3-a2fc-c0c1d016c8ab
plot(PST2..., layout=(3,3), margin=2mm, spacing=2mm, size=(800,600))

# ╔═╡ 6071de63-a274-431f-8361-a905e4294f27
MXSS2, EDS2, MXST2, EDT2, VT02 = get_ss_st_stats(rdir, state[2], tsamples, fields;
												 smin=-400.0, smax=400.0,
												 nbinss=30, nbinst=30)

# ╔═╡ f7e92472-160a-49b1-8e78-8cae4ee886d4
let
	s1 = scatter(MXSS2,  xlabel="Sample number", ylabel="#steps at sx",
				ylims=(0,500), label="chelated")
	s1 = scatter!(s1, MXSS, ylims=(0,500), label="free")
	s1 = plot!(s1, MXSS2, label="chelated")
	s1 = plot!(s1, MXSS, label="free")
	
	
	s2 = scatter(MXST2,  xlabel="Sample number", ylabel="#steps at t=max",
				ylims=(0,500), label="chelated")
	s2 = scatter!(s2, MXST,  ylims=(0,500), label="free")
	s2 = plot!(s2, MXST2, label="chelated")
	s2 = plot!(s2, MXST, label="free")
	
	s3 = scatter(VT02,  xlabel="Sample number", ylabel="#steps at t=0",
				ylims=(0,500), label="chelated")
	s3 = scatter!(s3, VT0, label="free")
	s3 = plot!(s3, VT02, label="chelated")
	s3 = plot!(s3, VT0, label="free")
	plot(s1,s2, s3, layout=(3,1), margin=2mm, spacing=2mm, size=(800,600))
end

# ╔═╡ de6afef4-2255-4107-b9ca-4aca85a99dc5
mxsf, mxsc, mxtf, mxtc, vt0f, vt0c = subtract_blank(MXSS, MXSS2, MXST, MXST2, VT0, VT02; iend=6)

# ╔═╡ 0dcb0d14-062d-467d-8431-c431db4f8d39
let
	s1 = scatter(mxsc,  xlabel="Sample number", ylabel="#steps at sx",
				ylims=(0,350), label="chelated-blank", legend=:topleft)
	s1 = scatter!(s1, mxsf, ylims=(0,350), label="free-blank")
	s1 = plot!(s1, mxsc, label="chelated")
	s1 = plot!(s1, mxsf, label="free")
	
	
	s2 = scatter(mxtc,  xlabel="Sample number", ylabel="#steps at t=max",
				ylims=(0,400), label="chelated", legend=:topleft)
	s2 = scatter!(s2, mxtf,  ylims=(0,400), label="free")
	s2 = plot!(s2, mxtc, label="chelated")
	s2 = plot!(s2, mxtf, label="free")
	
	s3 = scatter(vt0c,  xlabel="Sample number", ylabel="#steps at t=0",
				ylims=(0,400), label="chelated", legend=:topleft)
	s3 = scatter!(s3, vt0f, label="free")
	s3 = plot!(s3, vt0c, label="chelated")
	s3 = plot!(s3, vt0f, label="free")
	plot(s1,s2, s3, layout=(3,1), margin=2mm, spacing=2mm, size=(800,600))
end

# ╔═╡ d0834b37-2339-47de-af17-37813dcd6b70
mean(mxsc)

# ╔═╡ 6d841a1a-8e17-4362-839f-8a44e6d72a20
mean(mxsf)

# ╔═╡ 60b8ee0e-a426-4d75-9454-63738bc58537
std(mxsf)

# ╔═╡ e3460ff2-cc1a-4f2e-a581-75dce7ac9476
std(mxsc)

# ╔═╡ Cell order:
# ╠═3a4634d0-063c-11f0-2c98-ed94f1170f79
# ╠═fb1cc299-6238-445e-99be-9c167823417c
# ╠═4cbeae67-5ee0-4e33-98dd-4abb98874786
# ╠═7693af9f-b4ca-4127-93c0-771f3aad54c1
# ╠═c185fa05-d717-42be-81df-9b6f0e30379b
# ╠═56d4de31-efc1-4b34-9990-dd0b4ff54073
# ╠═be9c2108-5b23-48f2-8ff3-3edaf899568e
# ╠═f7bb460d-0cfc-40bb-9dc6-56b4636530a4
# ╠═909c89fa-5383-4db7-ac2f-9d4e2590b602
# ╠═906abfe8-24fe-40f7-abfc-9ce58a717430
# ╠═fd3a641a-78eb-4cb9-a495-4d487b02b834
# ╠═dc92f1d9-93a4-430a-bf72-334da4cbe2ba
# ╠═8933c28d-be0e-482b-b72c-b83798178fff
# ╠═ebb1bdf6-9eac-41cd-8513-caae11875c2a
# ╠═5c74b873-5716-4a8a-848d-ff9da84d8250
# ╠═524a2bde-bc60-4463-91ce-ddbd20b6a574
# ╠═865b5b90-15d4-4c67-94c5-6f30b650bff0
# ╠═e99aca88-68b6-4ee7-8a71-99d0a9261450
# ╟─61f5f861-60c8-4bf4-93fe-34f8233a79f0
# ╠═1b1c5630-1e22-4fe3-be28-2c0591c3a063
# ╠═bee8c3c3-4512-403e-af20-b76ae3534a10
# ╟─5cc2980f-3181-4840-ba93-83cf5745f566
# ╠═403736b5-3474-45bb-b71d-ad7639c6dd4d
# ╠═32cc64c3-9214-4b24-9d3b-23d42cafe574
# ╠═205364db-aaf3-4af6-b3d8-114ba05faf16
# ╠═b55a32e3-8015-4f50-b34f-cdf2fed5559e
# ╠═f16ba4df-3398-4a34-b12f-bcf848aad609
# ╠═5d43c3bb-be52-4a0c-a638-0a2cf64424e7
# ╠═2c4e2d24-fc1e-4d24-a672-7c1c3b1306dc
# ╠═ce0bca19-1679-4431-9b49-b74171745bb6
# ╠═ca368f79-bb07-49d4-9c53-83c46e69e9b7
# ╠═2e7b7742-38dc-4dad-8dae-0fab9bcc486f
# ╠═47f153ad-3943-4314-9e38-22f4af02f017
# ╠═7595358d-2264-4318-ad48-c7bbc8871300
# ╠═4de31e96-5c89-46cf-81eb-6eaa1dc0712c
# ╠═ab77a9ac-b60c-4f35-b63e-c5edd2fb1464
# ╠═d76192c0-ef94-4b44-8b7c-5e6468a54928
# ╠═e236f801-bb66-4658-bada-4e5d3ef9fb2b
# ╠═6832eeaf-3ffc-46c3-a2fc-c0c1d016c8ab
# ╠═6071de63-a274-431f-8361-a905e4294f27
# ╠═f7e92472-160a-49b1-8e78-8cae4ee886d4
# ╠═de6afef4-2255-4107-b9ca-4aca85a99dc5
# ╠═0dcb0d14-062d-467d-8431-c431db4f8d39
# ╠═d0834b37-2339-47de-af17-37813dcd6b70
# ╠═6d841a1a-8e17-4362-839f-8a44e6d72a20
# ╠═60b8ee0e-a426-4d75-9454-63738bc58537
# ╠═e3460ff2-cc1a-4f2e-a581-75dce7ac9476
