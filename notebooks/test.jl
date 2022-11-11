### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

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

# ╔═╡ Cell order:
# ╠═80a4e1cc-4d59-4c8a-9d65-1fe0d2d77bf2
# ╠═23b4e644-8331-4aa6-8369-47ce3ff0e143
# ╠═9f71bc31-df8f-48c5-815a-c8db9e46100f
# ╠═19035c0d-93a7-4bb9-be05-d2a9b9ac4619
# ╠═f4e379c4-a2b2-4703-bcbc-f6d7c996354a
# ╠═d88ba94f-dd9e-4c46-adac-4c3f7e556cad
# ╠═51710bb8-9c4a-4abc-9fe9-02e87bd4e4c5
# ╠═d1ace15e-fe1a-4552-9144-a0824ae8ae0f
