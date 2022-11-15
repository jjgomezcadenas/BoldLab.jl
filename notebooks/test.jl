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

# ╔═╡ 0122a81d-ee51-4355-8ddb-6429f15e8ea1
md"""
# Code
"""

# ╔═╡ 0937a6fc-6936-47d0-80de-8a38bb9a6a37
md"""
## Dropbox path
Here you must enter de path of your LaserLab Dropbox. With this path, the program can access to every file in the Dropbox.
"""

# ╔═╡ 6d04e8fe-c174-4be5-bffc-40945e8074e5
LaserLabp=ENV["JLaserData"]

# ╔═╡ c45b108e-62f9-473a-9975-eec4736d5df1
md""" 
## Select run
"""

# ╔═╡ da7dcbd8-75c3-42d5-9958-5ccbcc9624f1
begin
bodir=joinpath(LaserLabp,"pdata\\")
bcmdir=joinpath(LaserLabp,"data\\")
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

# ╔═╡ 7d759e14-e9a9-440c-ac8a-de92ff312526
@test ("data" in readdir(ENV["JLaserData"]))

# ╔═╡ 2c1ca8f4-ed1a-46d1-9a34-1e76c9005a87
@test ("pdata" in readdir(ENV["JLaserData"]))

# ╔═╡ 905ca07b-5778-4d8a-815d-7a8bdd4b73d4
@test ("FLUORI" in readdir(ENV["JLaserData"]))

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
# ╠═0122a81d-ee51-4355-8ddb-6429f15e8ea1
# ╠═0937a6fc-6936-47d0-80de-8a38bb9a6a37
# ╠═6d04e8fe-c174-4be5-bffc-40945e8074e5
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
# ╠═7d759e14-e9a9-440c-ac8a-de92ff312526
# ╠═2c1ca8f4-ed1a-46d1-9a34-1e76c9005a87
# ╠═905ca07b-5778-4d8a-815d-7a8bdd4b73d4
# ╠═71ed6fbd-0342-4cf5-9d8c-aa8f791d85f1
# ╠═20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
# ╠═b9970588-422f-461f-addb-5169d2e6043e
# ╠═0bf5c4a7-b369-4438-a987-253f242dd88f
