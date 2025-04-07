### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 34233420-116a-11f0-0e9c-75681c7ddb3f
using Pkg; Pkg.activate(ENV["JBoldLab"])

# ╔═╡ 3738ca4c-c6e2-4a4d-90f2-1ee3284b7c89
push!(LOAD_PATH, ENV["JBoldLab"] * "/src")

# ╔═╡ bba8e102-398c-4191-9d88-5cbc7879d6cf
begin
	using Revise
	using BoldLab
	using SimpleLogger
	using dffunctions
end

# ╔═╡ bf3dd6ae-7cf8-45b4-b407-f480cb9e7f30
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
	using Images, FileIO
	using Interpolations
	using QuadGK
end

# ╔═╡ 934b4ebc-a4f4-4900-9d19-f8b04726a635
ENV["JBoldLab"]

# ╔═╡ 21345bd9-5cd7-4775-805c-e1407072a616
names(BoldLab)

# ╔═╡ a650449d-6f90-4655-ae83-64fe4eb67b22
names(SimpleLogger)

# ╔═╡ 41053e81-3b76-41f3-9f54-076518a64e98
names(dffunctions)

# ╔═╡ b37387b8-79f9-487e-906b-5e654efa4d9f


# ╔═╡ be94d574-ecb6-4d5e-aea4-e0f7feaf421a
"""Simple trapezoidal integral""" 
function trapz(x::Vector{Float64}, y::Vector{Float64})
	s = 0.0
	for i in 1:length(x)-1
		s += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
	end
	return s
end
	


# ╔═╡ d5662217-0a05-4206-99b1-6bfa8fe4f00b
"""A sample normalization/interpolation function"""
    function norm_and_interp(x::Vector{Float64}, y::Vector{Float64})
        # Simple trapezoidal integration for area
        area = sum(0.5 .* (y[1:end-1] .+ y[2:end]) .* diff(x))
        y_norm = y ./ area
        return LinearInterpolation(x, y_norm, extrapolation_bc=Flat())
    end


# ╔═╡ c4c28c20-b071-46dd-addc-b4d0025cc170
function interplate_emss(df)
		
	pdf_free = norm_and_interp(df.lambda, df.NAPH3_emission) 
	pdf_ba2p = norm_and_interp(df.lambda, df."NAPH3-Ba_emission") 
	pdf_free, pdf_ba2p
end


# ╔═╡ 951b11b0-cf34-465b-90b6-3b146940b8d2
function interplate_abs(df, ; xsec=false)
	eps_to_xsec = 3.82E+3 # eps in M*1 cm^-1, xsec in barn (10^-24 cm2)
	
	if xsec == true
		Yfree = eps_to_xsec * df."NAPH3_UV-vis"
		YBa2p = eps_to_xsec * df."NAPH3-Ba_Uvvis"
	else
		Yfree = df."NAPH3_UV-vis"
		YBa2p = df."NAPH3-Ba_Uvvis"
	end
	
	pdf_free = LinearInterpolation(df.lambda, Yfree, extrapolation_bc=Line()) 
	pdf_ba2p = LinearInterpolation(df.lambda, YBa2p, extrapolation_bc=Line()) 
	pdf_free, pdf_ba2p
end

# ╔═╡ ed5b380c-c08b-46ae-8339-74c8654982b4
function plot_emss(df)
	
	p1 = plot(df.lambda, df.NAPH3_emission,
	     label = "NAPH3 Emission",
	     xlabel = "Wavelength (nm)",
	     ylabel = "Emission (au)",
		 title = "NAPH3 Emission Spectra",
	     legend = :topright)
	
	p1 = plot!(p1, df.lambda, df."NAPH3-Ba_emission",
	      label = "NAPH3-Ba Emission")
	p1
end

# ╔═╡ bd1551b6-510c-42d1-a711-349e295c95de
function plot_abs(df; xsec=false)
	eps_to_xsec = 3.82E+3 # eps in M*1 cm^-1, xsec in barn (10^-24 cm2)

	if xsec == false
		p1 = plot(df.lambda, df."NAPH3_UV-vis",
		     label = "NAPH3 Absorption",
		     xlabel = "Wavelength (nm)",
		     ylabel = "\$\\sigma (M^{-1} cm^{-1})\$",
			 title ="NAPH3 Absorption Spectra",
		     legend = :topright)

		p1 = plot!(p1, df.lambda, df."NAPH3-Ba_Uvvis",
		      label = "NAPH3-Ba Absorption")
	else
		xs = eps_to_xsec * df."NAPH3_UV-vis"
		p1 = plot(df.lambda, xs,
		     label = "NAPH3 Absorption",
		     xlabel = "Wavelength (nm)",
		     ylabel = "\$\\sigma (barn)\$",
			 title ="NAPH3 Absorption Spectra",
		     legend = :topright)

		xs = eps_to_xsec * df."NAPH3-Ba_Uvvis"
		p1 = plot!(p1, df.lambda, xs,
		      label = "NAPH3-Ba Absorption")
	end
	p1
end

# ╔═╡ b74a28fe-e6d8-4156-b7b1-bd2276755b40
function plot_spectra(dfa, dfe; xsec=false)
	p1 = plot_abs(dfa; xsec)
	p2 = plot_emss(dfe)
	plot(p1,p2, layout=(1,2), size=(700, 500))
end

# ╔═╡ 88cbfd09-94ae-44d2-b64a-602578962ce1
struct Molecule{T}
    QY::Float64  # Quantum Yield
    em_spectrum::Tuple{Vector{Float64}, Vector{Float64}}  # emission spectrum pdf
    abs_spectrum::Tuple{Vector{Float64}, Vector{Float64}}   # absorption spectrum
    xs::T        # Interpolated absorption cross section function
    pdf::T       # Interpolated emission PDF function
end

# ╔═╡ 0ce2e4f2-8d13-44aa-8753-b29409ba2306
"""
Creates the FM structs from the DF. Assumes that the DF columns are:
for dfe (emission data frame):
c1: lambda, c2; free emission; c3 chelated emissions.
for dfa (absorbtion data frame):
c1: lambda, c2; free abs; c3 chelated abs.
"""
function fm_from_df(dfe, dfa)
	eps_to_xsec = 3.82E+3 # eps in M*1 cm^-1, xsec in barn (10^-24 cm2)
	
	function itp(λ, ϕ)
		area = trapz(λ, ϕ)
		y     = ϕ ./area
		LinearInterpolation(λ, y, extrapolation_bc=Line())
	end
		
	ndfe = names(dfe)
	le = Float64.(dfe[!,ndfe[1]])
	emf = Float64.(dfe[!,ndfe[2]])
	emc = Float64.(dfe[!,ndfe[3]])
	
	iemf = itp(le, emf)
	iemc = itp(le, emc)

	ndfa = names(dfa)
	la = Float64.(dfa[!,ndfa[1]])
	ema = Float64.(dfa[!,ndfa[2]]) 
	emca = Float64.(dfa[!,ndfa[3]]) 
	
	iema = LinearInterpolation(la, ema * eps_to_xsec, extrapolation_bc=Line())
	iemca = LinearInterpolation(la, emca * eps_to_xsec, extrapolation_bc=Line())
	
	Molecule(1.0, (le,emf), (la,ema), iemf, iema), Molecule(1.0, (le,emc), (la,emca), iemc, iemca)
end

# ╔═╡ 4badd96d-2461-4507-b771-1b4882cd239d
function plot_molecule(n3f, n3c)
	p1 = plot(n3f.em_spectrum[1], n3f.em_spectrum[2],
			  label = "NAPH3",
		      xlabel = "λ (nm)",
		      ylabel = "a.u.",
			  title = "NAPH3 measured emission spectrum",
		      legend = :topright)
	p1 = plot!(p1, n3c.em_spectrum[1], n3c.em_spectrum[2],
			  label = "NAPH3-Ba")
	
	ls = 400:10:800
	p2 = plot(ls, n3f.xs.(ls),
			  label = "NAPH3",
		      xlabel = "λ (nm)",
		      ylabel = "PDF",
			  title = "NAPH3 emission PDF",
		      legend = :topright)
			  
	p2 = plot!(p2, ls, n3c.xs.(ls),
			  label = "NAPH3-Ba")

	p3 = plot(n3f.abs_spectrum[1], n3f.abs_spectrum[2],
			  label = "NAPH3",
		      xlabel = "λ (nm)",
		      ylabel = "\$\\epsilon (M^{-1} cm^{-1})\$",
			  title = "NAPH3 \$\\epsilon \$",
		      legend = :topright)
	p3 = plot!(p3, n3c.abs_spectrum[1], n3c.abs_spectrum[2],
			  label = "NAPH3-Ba")
	
	ls = 200:1:500
	p4 = plot(ls, n3f.pdf.(ls),
			  label = "NAPH3",
		      xlabel = "λ (nm)",
		      ylabel = "\$\\sigma\$ (barn)",
			  title = "NAPH3 \$\\sigma\$ (barn)",
		      legend = :topright)
	p4 = plot!(p4, ls, n3c.pdf.(ls),
			  label = "NAPH3-Ba")
	plot(p1,p2, p3, p4, layout=(2,2), size=(700, 500))
end

# ╔═╡ 9288a859-e358-4c78-8947-38244f9019f7
begin
	jbl = ENV["JBoldLab"]
	fn3e = joinpath(jbl, "data", "naph3_emission.csv")
	fn3a = joinpath(jbl, "data", "naph3_abs.csv")
	
	dfn3e = CSV.read(fn3e, DataFrame)
	dfn3a = CSV.read(fn3a, DataFrame)
	plot_spectra(dfn3a, dfn3e, xsec=false)
end

# ╔═╡ 8353ef43-a2c9-4d45-bd76-df81872ed900
n3f, n3c = fm_from_df(dfn3e, dfn3a)

# ╔═╡ f65ea54a-c375-4ba3-864c-820bb98f67d9
plot_molecule(n3f, n3c)

# ╔═╡ 1d59335e-7c76-49c2-88f6-aef6e4312941
begin
	naph3f, naph3c = fm_from_csv(joinpath(jbl, "data"), 
									  "naph3_emission.csv",
									  "naph3_abs.csv")
end

# ╔═╡ a2076949-49f3-4964-a2f8-6b3d6c6303fb
plot_molecule(naph3f, naph3c)

# ╔═╡ Cell order:
# ╠═34233420-116a-11f0-0e9c-75681c7ddb3f
# ╠═bf3dd6ae-7cf8-45b4-b407-f480cb9e7f30
# ╠═3738ca4c-c6e2-4a4d-90f2-1ee3284b7c89
# ╠═934b4ebc-a4f4-4900-9d19-f8b04726a635
# ╠═bba8e102-398c-4191-9d88-5cbc7879d6cf
# ╠═21345bd9-5cd7-4775-805c-e1407072a616
# ╠═a650449d-6f90-4655-ae83-64fe4eb67b22
# ╠═41053e81-3b76-41f3-9f54-076518a64e98
# ╠═b37387b8-79f9-487e-906b-5e654efa4d9f
# ╠═be94d574-ecb6-4d5e-aea4-e0f7feaf421a
# ╠═d5662217-0a05-4206-99b1-6bfa8fe4f00b
# ╠═c4c28c20-b071-46dd-addc-b4d0025cc170
# ╠═951b11b0-cf34-465b-90b6-3b146940b8d2
# ╠═ed5b380c-c08b-46ae-8339-74c8654982b4
# ╠═bd1551b6-510c-42d1-a711-349e295c95de
# ╠═b74a28fe-e6d8-4156-b7b1-bd2276755b40
# ╠═88cbfd09-94ae-44d2-b64a-602578962ce1
# ╠═0ce2e4f2-8d13-44aa-8753-b29409ba2306
# ╠═4badd96d-2461-4507-b771-1b4882cd239d
# ╠═9288a859-e358-4c78-8947-38244f9019f7
# ╠═8353ef43-a2c9-4d45-bd76-df81872ed900
# ╠═f65ea54a-c375-4ba3-864c-820bb98f67d9
# ╠═1d59335e-7c76-49c2-88f6-aef6e4312941
# ╠═a2076949-49f3-4964-a2f8-6b3d6c6303fb
