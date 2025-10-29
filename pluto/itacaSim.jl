### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 6eee4840-a8e3-11f0-9db8-b93ffc85840e
begin
	using Pkg; 
	Pkg.activate(ENV["JBoldLab"])
end

# ╔═╡ e8f2e87c-30cf-424f-a562-0f3f61ac80f1
begin
	using BoldLab

	using CSV
	using DataFrames
	using Plots 
	using Printf
	using HDF5
	using Statistics
	using StatsBase
	using Images, FileIO
	using Interpolations
	using QuadGK
	using Distributions, Random
	using Unitful
	using UnitfulEquivalences
	using PhysicalConstants.CODATA2018
	using NearestNeighbors
	import Measures 
	import TiffImages
end

# ╔═╡ 88cdb157-b57a-40b6-ab3c-3ab8bbf1564d
begin
	cmdir=joinpath(ENV["DATA"], "HD5t/itaca")
	pdir =joinpath(ENV["PROJECTS"], "BoldLab")
end

# ╔═╡ 45fc0354-c4e7-45ac-8467-609b28cf7a61
import Unitful:
    nm, μm, mm, cm, m,
    ns, μs, ms, s, minute, hr, d, yr, Hz,
    eV,
    V, kV,
    μJ, mJ, J,
	mW, W

# ╔═╡ d75adcd8-2b38-4c02-b14f-4dfc4280a855
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

# ╔═╡ 16e942a5-20ed-4937-a5e6-987236a1e76d
jn = ingredients(string(pdir,"/src/BoldLab.jl"))

# ╔═╡ eb956cef-2a99-41fd-ad54-bb5272f8d7b7
Diffusion = ingredients("diffusion.jl")

# ╔═╡ 35dab5a3-75ec-4df6-8443-3ed9537e1dd7
md"""
- NAPH cross sections
"""

# ╔═╡ d892bb34-e552-4da2-ae1b-c096ab429709
begin
jbl = ENV["JBoldLab"]
	naph3f, naph3c = fm_from_csv(joinpath(jbl, "data"), 
									  "naph3_emission.csv",
									  "naph3_abs.csv")
	plot_molecule(naph3f, naph3c)
end

# ╔═╡ 795cdf62-31a8-497b-9cee-f1ea4077eaf9
begin
	mipa = load("mipa.png")
	mipa_small = imresize(mipa, ratio=0.5)
end

# ╔═╡ 92754e5b-ae80-4699-b1e6-3ec665f08248
begin
#P = 10.0 #in bar

	diamx(dx) = sqrt(2*dx*dx)
	trk_l = 100.0mm 
	sigma_t = 1mm
	trk_a = trk_l * 4*sigma_t
	
	n_i = 1e+5
	dsk_p = 0.5mm
	n_dsk_mm2 = 1.0/(dsk_p/mm)^2
	n_d = Int(trk_a/dsk_p^2)
	n_id = n_i/n_d
	dsk_x = 36μm
	dsk_d = diamx(dsk_x) 
	dskm = dsk_d/μm
	dsk_a = π * (dsk_d/2)^2
	
	n_mol = 1e+4/μm^2 # number of molecules per micron
	n_mol_mu2 = 1e+4
	eff_PET = 0.01
	n_dsk = n_mol * dsk_a
	
	println("- track length = $(trk_l/mm) mm")
	println("- track area = $(trk_a/mm^2) mm2")
	println("- disk pitch = $(dsk_p/mm) mm")
	println("- number of disks/mm2 = $(n_dsk_mm2)")
	println("- number of disks = $(n_d)")
	println("- average number of ions/disk = $(n_id)")
	
	println("- disk diameter = $(@sprintf("%.1f", dsk_d/μm)) μm")
	println("""- disk area = $(@sprintf("%.1f", dsk_a/μm^2)) μm2""")
	println("""- number of molecules/micron = $(@sprintf("%.1e", n_mol))""")
	println("""- number of molecules/disk = $(@sprintf("%.1e", n_dsk))""")
end

# ╔═╡ e35b6e57-f943-4488-9743-d2231095ba89
begin
	
end

# ╔═╡ 77b0fe9d-a457-4d36-8ebd-7c2b8ed4389e
md"""
## Description of problem
- Assume that the track has a length of $(trk_l/mm) mm 10 and diffusion $(sigma_t/mm) mm. 
- We then want to scan an area of (track length) x 4 mm ($\pm 2 \sigma$ of diffusion). 
- Therefore $(trk_a/mm^2) mm2 
- We focus the track using a microstructure similar to the Micro Pin Array Detector (MIPA), shown above.
- Replace the current cathode by a gate grid held at some positive voltage. In normal operations, the cathode is held at the same voltage, which we can take as ground (V=0 volts). When a $\beta\beta$ event triggers the detector, the cathode is ramped down to a more negative voltage, so that the gate is "opened".
- The gate is an hexagonal grid with hexagons of 2.5 mm diameter (like that of NEXT). The MIPA is designed to focus each one of these hexagons in disks of $(@sprintf("%.1f",dskm)) μm. 
- The disks are held at a more negative voltage than the grid and above the cathode voltage. 
- To sample properly the track across the diffusion we place the disks at a pitch of 500 μm. 
- Thus we get $(n_dsk_mm2) disks per mm2, and a maximum of $(n_d) disks to scan.  
- The number of ions in the track is $(@sprintf("%.1e", n_i)) and in average $(n_id) ions per disk.

"""

# ╔═╡ 88a1e362-2ff8-4fa0-8b28-1057c311a668
md"""
### Read data from file 
"""

# ╔═╡ e0347b2c-8f8c-46dd-a857-01a036b70073
md"""
### Load tracks
"""

# ╔═╡ 8f9d90d8-d20b-42b6-a784-ccb2516f2158
md"""
## Appendix
"""

# ╔═╡ 9b5eec23-6302-4919-9f7f-07e9565fea96


  function select_neighbors_kdtree(edf::DataFrame, idf::DataFrame; cutoff=5.0)
      # Build KDTree for efficient spatial queries (idf points)
      idf_points = hcat(idf.x, idf.y)'  # 2 x N matrix
      kdtree = KDTree(idf_points)

      # Query points from edf
      edf_points = hcat(edf.x, edf.y)'  # 2 x M matrix

      # Find all neighbors within cutoff for each edf point
      indices = inrange(kdtree, edf_points, cutoff)

      # Flatten and get unique indices from idf
      selected_indices = unique(vcat(indices...))

      # Return filtered idf
      return idf[selected_indices, :]
  end

  

# ╔═╡ 3c65af6a-ad59-42c0-9fb8-be4b2830d7b5
function parse_filename(filename::String)
      # Extract n (number of molecules)
      n_match = match(r"_n([\d.e+\-]+)", filename)
      n = n_match !== nothing ? parse(Float64, n_match.captures[1]) : nothing

      # Extract power in mW
      p_match = match(r"_([\d.]+)mW", filename)
      p = p_match !== nothing ? parse(Float64, p_match.captures[1]) : nothing

      # Extract pb (photobleaching) - handle the + sign
      pb_match = match(r"_pb([\d.e+\-]+)", filename)
      pb = pb_match !== nothing ? parse(Float64, pb_match.captures[1]) : nothing

      # Extract texp (exposure time) - handle the + sign and .csv extension
      te_match = match(r"_texp([\d.e+\-]+?)(?:\.csv)?$", filename)
      te = te_match !== nothing ? parse(Float64, te_match.captures[1]) : nothing

      return (n=n, p=p, pb=pb, te=te)
  end
 
  

# ╔═╡ 7f0e2846-74ef-4a62-82e9-070e4252514c
begin
	filem = "itaca_molecules_n1e5_5mW_pb1e+4_texp1e-1.csv"
	filem = "itaca_molecules_n1e5_10mW_pb1e+4_texp1e-1.csv"
  #filem = "itaca_molecules_n1e6_10mW_pb1e+4.csv"
#filem = "itaca_molecules_n1e6_100mW_pb1e+4.csv"
  result = parse_filename(filem)
  #println("n = $(result.n)")    # n = 1.0e6
  #println("p = $(result.p) mW")  # p = 100.0 mW
  #println("pb = $(result.pb)")   # pb = 10000.0
end

# ╔═╡ 0b28053d-2403-4160-9e58-4d74932700ec
md"""
- number of molecules in file = $(result.n)
- laser power = $(result.p) mW
- photobleaching cycles = $(result.pb)
- exposition time = $(result.te)
"""

# ╔═╡ 11b4e71e-51c0-4104-a821-3c6c9178a9ae
im1e6 = CSV.read(filem, DataFrame) # 10^6 molecules

# ╔═╡ 527a8d40-0d63-42e7-9647-3d30730b1654
begin
	#nphot = collect(im1e6.Nphot)
	nphot = collect(im1e6.R_Hz)
	rphot = collect(im1e6.R_Hz)
	tpb = collect(im1e6.Tb_s)
	xmu = collect(im1e6.x_μm)
	ymu = collect(im1e6.y_μm)
	
	hxmu, pxmu2 = step_hist(ymu;
                 nbins = 100,
                 xlim   = (0.0, 15.0),
                 xlabel = " diameter (μm)",
                 ylabel = "Frequency",
                 title=" distribution of molecules ")

	hnph, pnph = step_hist(nphot;
                 nbins = 500,
                 #xlim   = (0.0, 5e+4),
                 xlabel = " N photons",
                 ylabel = "Frequency",
                 title=" Photons emitted ")
	hrph, prph = step_hist(rphot;
                 nbins = 100,
                 #xlim   = (0.0, 5e+4),
                 xlabel = " Rate (Hz)",
                 ylabel = "Frequency",
                 title=" Emission rate ")
	htpb, ptpb = step_hist(tpb;
                 nbins = 100,
                 xlim   = (0.0, 15.0),
                 xlabel = " time (second)",
                 ylabel = "Frequency",
                 title=" Phtobleaching time ")
	plot(pxmu2, pnph, prph, ptpb,
       layout=(2, 2))
end

# ╔═╡ 61da2c80-a9a7-4ea0-a855-356c0c300021
begin
	ng_mol = mean(nphot)
	nc = n_id 
	ngc = nc * ng_mol
	nb = n_dsk * eff_PET
	ngb = nb * ng_mol 
	stn = ngc/sqrt(ngb)
	println("""- average number of photons/molecule = $(@sprintf("%.1e", ng_mol))""")
	println("""- number chelated/pin = $(@sprintf("%.1e", nc))""")
	println("""- number gammas/chelated = $(@sprintf("%.1e", ngc))""")
	println("""- number bkgnd/pin = $(@sprintf("%.1e", nb))""")
	println("""- number gammas/background = $(@sprintf("%.1e", ngb))""")
	println("""- S/N = $(@sprintf("%.1e", stn))""")
end

# ╔═╡ e465533f-8416-4b4a-9db5-5d76eb38520c
begin
	plot(pnph,  ptpb,
       layout=(1, 2), size=(600,300))
end

# ╔═╡ 825d578e-6662-456a-96aa-afd6f2107475
begin
"""
sigma_t (mm) = dtmm /sqrt(Pbar) x sqrt(Lcm)
"""
function sigma_t_mm(Lcm, Pbar; dtmm=3.5) 
    dtmm /sqrt(Pbar) * sqrt(Lcm)
end

"""
sigma_l (mm) = dlmm /sqrt(Pbar) x sqrt(Lcm)
"""
function sigma_l_mm(Lcm, Pbar; dlmm=0.9) 
    dlmm /sqrt(Pbar) * sqrt(Lcm)
end


"""
sigma_t (mm) = sqrt(2Kb/e x T x L/E) * 10
where:
- 2Kb/e in V/K
- T in K
- L in cm
- E in V/cm
"""
function sigma_t_ion_mm(Tk, Lcm, Evcm; kB_ovr_e = 8.6173E-5)  
    # T -> K, L -> cm, E -> v/cm
    return sqrt(2*kB_ovr_e*Tk*(Lcm/Evcm)) * 10 #-> mm
end
end

# ╔═╡ 22ac54a3-e73e-4a19-a0ac-31f5d24d9beb
begin 
	sigmat_xe = 3.5 # mm
	sigmat_he = 1.5 # mm
	sigmal_xe = 0.9 # mm
	sigmal_he = 0.75 # mm
	kB_ovr_e_lxe = 8.6173E-5 # V/K
	kB_ovr_e_gxe = 2.59E-2 # V/K
	sigmat_nh4_lxe_50 = sigma_t_ion_mm(170.0, 50.0, 400.0) # 50 cm drift, 400 V/cm
	sigmat_nh4_lxe_10 = sigma_t_ion_mm(170.0, 10.0, 400.0) # 50 cm drift, 400 V/cm
	
	sigmat_nh4_gxe_100 = sigma_t_ion_mm(298.0, 100.0, 400.0) # 100 cm drift, 400 V/cm
	sigmat_nh4_gxe_50 = sigma_t_ion_mm(298.0, 50.0, 400.0) # 100 cm drift, 400 V/cm
	sigmat_nh4_gxe_10 = sigma_t_ion_mm(298.0, 10.0, 400.0) # 100 cm drift, 400 V/cm
	
	sigmat_e_xe_100 = sigma_t_mm(100.0, 10.0; dtmm=3.5)
	sigmat_e_he_100 = sigma_t_mm(100.0, 10.0; dtmm=1.6)
	sigmal_e_xe_100 = sigma_l_mm(100.0, 10.0; dlmm=0.9)
	sigmal_e_he_100 = sigma_l_mm(100.0, 10.0; dlmm=0.75)
	
	sigmat_e_xe_50 = sigma_t_mm(50.0, 10.0; dtmm=3.5)
	sigmat_e_he_50 = sigma_t_mm(50.0, 10.0; dtmm=1.6)
	sigmal_e_xe_50 = sigma_l_mm(50.0, 10.0; dlmm=0.9)
	sigmal_e_he_50 = sigma_l_mm(50.0, 10.0; dlmm=0.75)
	
	sigmat_e_xe_10 = sigma_t_mm(10.0, 10.0; dtmm=3.5)
	sigmat_e_he_10 = sigma_t_mm(10.0, 10.0; dtmm=1.6)
	sigmal_e_xe_10 = sigma_l_mm(10.0, 10.0; dlmm=0.9)
	sigmal_e_he_10 = sigma_l_mm(10.0, 10.0; dlmm=0.75)
	
	println("| Parameter                  | 100 mm          | 50 mm           | 10 mm           |")
	println("|----------------------------|-----------------|-----------------|-----------------|")
	println("| Ions (σₜ)                  | $(@sprintf("%.2f", sigmat_nh4_gxe_100)) mm | $(@sprintf("%.2f", sigmat_nh4_gxe_50)) mm | $(@sprintf("%.2f", sigmat_nh4_gxe_10)) mm |")
	println("| Electrons in Xe (σₜ)       | $(@sprintf("%.2f", sigmat_e_xe_100)) mm | $(@sprintf("%.2f", sigmat_e_xe_50)) mm | $(@sprintf("%.2f", sigmat_e_xe_10)) mm |")
	println("| Electrons in Xe/He (σₜ)    | $(@sprintf("%.2f", sigmat_e_he_100)) mm | $(@sprintf("%.2f", sigmat_e_he_50)) mm | $(@sprintf("%.2f", sigmat_e_he_10)) mm |")
	println("| Electrons in Xe (σₗ)       | $(@sprintf("%.2f", sigmal_e_xe_100)) mm | $(@sprintf("%.2f", sigmal_e_xe_50)) mm | $(@sprintf("%.2f", sigmal_e_xe_10)) mm |")
	println("| Electrons in Xe/He (σₗ)    | $(@sprintf("%.2f", sigmal_e_he_100)) mm | $(@sprintf("%.2f", sigmal_e_he_50)) mm | $(@sprintf("%.2f", sigmal_e_he_10)) mm |")
end

# ╔═╡ 3ee80fd1-d3d6-4fca-a415-e37561001a2a
function get_event(hitsdf::DataFrame, event_id::Int)

    ghdf = groupby(hitsdf, :event_id)
    
    event_data = nothing
    for group in ghdf
        if first(group.event_id) == event_id
            event_data = group
            break
        end
    end

    # Check if event exists
    if event_data === nothing
        error("Event ID $event_id not found in grouped dataframe")
    end

	DataFrame(event_data)
end

# ╔═╡ 2a789937-bd3b-45dc-83ed-7d62abcfbfc7
function get_dataset_dfs(filename::String)
    h5open(filename, "r") do fid
        group = fid["MC"]
        dfs = Dict{String, DataFrame}()

        for name in keys(group)
            data = read(group[name])
            # Try to make a DataFrame if it's an array of structs or tuples
            try
                dfs[name] = DataFrame(data)
            catch
                dfs[name] = DataFrame((value=data,))
            end
        end

        return dfs  # Dict of DataFrames
    end
end

# ╔═╡ cc400dae-2b04-4941-8fea-300426618bed
begin
	xfile = "/Users/jjgomezcadenas/Data/HD5t/itaca/0nubb_100um.next.h5"
	filename = basename(xfile)
	#hitdf = CSV.read(xfile, DataFrame)
	dfs = get_dataset_dfs(xfile)
	hitdf = dfs["hits"]
	println("Loaded xfile = $(filename)")
end

# ╔═╡ 2ea5696e-6aa0-4f16-9ca3-c9ee4c1f421f
begin
	nevent =4
	evtdf = get_event(hitdf, nevent)
	evtdf.electrons = Int.(round.(evtdf.energy .* 1e5))
end

# ╔═╡ 83b8a99b-c8a8-4f70-a630-fe95491f3223
begin
	dxy_he_50 =Diffusion.diffuse_xy_image_mc(evtdf; sigma_mm=sigmat_e_he_50, 
								   nbins=100, nsigma=3.0)
	pxy_he_50 =Diffusion.plot_diffused_xy(dxy_he_50; 
					 title="""e- Xe/He: Ld=50 cm: σt =$(@sprintf("%.2f",sigmat_e_he_50)) mm""", colorbar=false)
	dxy_nh4_50 =Diffusion.diffuse_xy_image_mc(evtdf; sigma_mm=sigmat_nh4_gxe_50, nbins=100, nsigma=3.0)
	pxe_50_nh4 =Diffusion.plot_diffused_xy(dxy_nh4_50; title="""Nh4: Ld=50 cm: σt =$(@sprintf("%.2f", sigmat_nh4_gxe_50)) mm""", colorbar=false)
	plot(pxy_he_50, pxe_50_nh4,
       layout=(1, 2),
       size=(1200, 800),
       plot_title="Diffusion in Xe/He: electrons & ions",
       link=:both)
end

# ╔═╡ 34bcea1a-c135-4289-a45e-c51319f664d3
begin
	grid, ions_df = Diffusion.diffuse_xy_ions_mc(evtdf; sigma_mm=sigmat_nh4_gxe_50, pitch_mm=0.5)
	md"""
- ion grid extends in x (mm) from $(grid.x_max) to $(grid.x_min)
- ion grid extends in y (mm) from $(grid.y_max) to $(grid.y_min)
- pitch (in mm) = $(grid.dx)
- number of bins = $(grid.nbins_x * grid.nbins_y)
"""
end

# ╔═╡ 83f79fa7-ce89-4b18-b27e-fa0d2f2a1cfc
begin
	hions, pions = step_hist(ions_df.n_ions[ions_df.n_ions .>0];
                 nbins = 100,
                 #xlim   = (0.0, 5500.0),
                 xlabel = " N ions/pin",
                 ylabel = "Frequency",
                 title=" N ion/pixel")
	pions
end

# ╔═╡ 86deb2ca-f847-41fe-a30f-67469a14b93d
idf = filter(:n_ions => >=(1), ions_df)

# ╔═╡ 482911a0-da45-4851-b73b-62faea930087
Diffusion.plot_diffused_xy(idf; title="Photons/ions", intensity_col=:photons, zero_as_nan=true)

# ╔═╡ e9cd96bf-b38d-4b0a-b1d6-5fee14e9a2be
pions_df =Diffusion.plot_diffused_xy(ions_df; title="""Nh4: Ld=50 cm: σt =$(@sprintf("%.2f", sigmat_nh4_gxe_50)) mm""", intensity_col=:n_ions, colorbar=false, zero_as_nan=true)

# ╔═╡ 93ab09d6-a439-45b7-8d5f-0eb23702d3cb
begin
	Diffusion.plot_diffused_xy(ions_df; title="Photons/ions", intensity_col=:npot_bkg)
end

# ╔═╡ 0c020931-0ff4-4ccc-8440-61fd7c6470f2
begin
	pion_sgn = Diffusion.plot_diffused_xy(ions_df; title="Photons/chelated", intensity_col=:npot_sgn, colorbar=false)
end

# ╔═╡ bf476c92-549c-4d9e-a301-a13aacd8dfd8
begin
	pion_photons = Diffusion.plot_diffused_xy(ions_df; title="Photons/all", intensity_col=:photons, colorbar=false)
end

# ╔═╡ 6061084b-a47a-4621-bfe3-96bab683f685
plot(pxy_he_50, pxe_50_nh4, pion_sgn, pion_photons,
       layout=(2, 2),
       size=(1200, 800),
       plot_title="Diffusion in Xe/He: electrons & ions",
       link=:both)

# ╔═╡ b5de052a-149a-416f-8f59-e42f97f3d5eb
begin
	xu, yu, imatrix = Diffusion.get_intensity_matrix(ions_df, :photons)
	ibkgnd = Diffusion.compute_background_from_matrix(imatrix; σ=5.0)
	ionimg = Diffusion.subtract_background_from_matrix(imatrix, ibkgnd)
	Diffusion.plot_imatrix(ionimg, xu, yu; 
                      title="intensity matrix", zero_as_nan=false,
                      titlefontsize=10, colorbar=false)
end                               

# ╔═╡ 84812b46-ee00-40d2-8d0d-b59066ef7db4
begin
	ionrb, xrb, yrb = Diffusion.rebin_matrix(ionimg, xu, yu, 2)
	Diffusion.plot_imatrix(ionrb, xrb, yrb; 
                      title="intensity matrix", zero_as_nan=false,
                      titlefontsize=10, colorbar=false)
end

# ╔═╡ d9ffc411-2c3e-4317-9571-2a7b0039f362
function sample_from_histogram(h::Histo1d, n::Int)
      # Create probability weights
      pw = ProbabilityWeights(h.weights)

      # Sample bin indices
      bin_indices = StatsBase.sample(1:length(h.weights), pw, n)

      # For each bin, sample uniformly within the bin
      samples = Float64[]
      for idx in bin_indices
          # Sample uniformly between edges[idx] and edges[idx+1]
          lower = h.edges[idx]
          upper = h.edges[idx+1]
          push!(samples, lower + rand() * (upper - lower))
      end

      return samples
  end


# ╔═╡ 9f1a3620-2ed3-44f7-a9c8-de69fea363e9
begin
	nbin_grid = grid.nbins_x * grid.nbins_y
	nphot_bkg_dsk_ion = sample_from_histogram(hnph , nbin_grid) * n_mol_mu2 * eff_PET 
	nphot_sgn_dsk_ion = sample_from_histogram(hnph , nbin_grid) 
	nsgn = nphot_sgn_dsk_ion
	ntot_mean =nphot_sgn_dsk_ion + nphot_bkg_dsk_ion 
	ntot =rand.(Poisson.(ntot_mean))
	hbkg, pbkg = step_hist(nphot_bkg_dsk_ion;
                 nbins = 100,
                 #xlim   = (0.0, 5500.0),
                 xlabel = " N photons",
                 ylabel = "Frequency",
                 title=" N photons/disk/ion: background ")
	hsgn, psgn = step_hist(nphot_sgn_dsk_ion;
                 nbins = 100,
                 #xlim   = (0.0, 5500.0),
                 xlabel = " N photons",
                 ylabel = "Frequency",
                 title=" N photons/disk/ion: signal")
	hntm, pntm = step_hist(ntot_mean;
                 nbins = 100,
                 #xlim   = (0.0, 5500.0),
                 xlabel = " N photons",
                 ylabel = "Frequency",
                 title=" N photons/disk/ion: mean ")
	hntp, pntp = step_hist(ntot;
                 nbins = 100,
                 #xlim   = (0.0, 5500.0),
                 xlabel = "N photons",
                 ylabel = "Frequency",
                 title="  N photons/disk/ion: measured ")
	plot(pbkg, psgn, pntm, pntp, layout=(2,2), size=(800, 600))
end	

# ╔═╡ 2520366e-7a26-49fb-a1bf-f2c431499de6
begin
	n_bkng_in_disk = Int(floor(n_mol_mu2  * dsk_a/μm^2 * eff_PET))

	NGB = Float64[]
	NGS = Float64[]
	for row in eachrow(ions_df)
		n_gamma_bkng_in_disk = sum(sample_from_histogram(hnph , n_bkng_in_disk ))
	  	n_gamma_signal_in_disk = sum(sample_from_histogram(hnph , row.n_ions))
	  	append!(NGB, n_gamma_bkng_in_disk)
		append!(NGS, n_gamma_signal_in_disk)

	end
end

# ╔═╡ 17662c19-6bff-46ca-9533-1f3651a4062c
 std(NGB)

# ╔═╡ 23dc9540-589f-4d17-8f07-ab4f8004daad
mean(NGB)

# ╔═╡ e2e96b6c-9b2f-492c-9313-5f1f92c24dbb
mean(NGS)

# ╔═╡ 692345c0-ea48-455c-98c6-90ef935054ea
begin
	ions_df[!, "npot_bkg"] = NGB
	ions_df[!, "npot_sgn"] = NGS
	ions_df[!, "photons"] = NGS .+ NGB
end

# ╔═╡ 6fd72e99-7a10-47a6-866f-0f75d373bab0
function plot_hits(df::DataFrame; nbins::Int=100)
    
    x = Float64.(df.x)
    y = Float64.(df.y)
    z = Float64.(df.z)
    e = Float64.(df.energy)

    # Compute padded limits (1.3x range)
    xmid, xrange = mean((minimum(x), maximum(x))), maximum(x) - minimum(x)
    ymid, yrange = mean((minimum(y), maximum(y))), maximum(y) - minimum(y)
    zmid, zrange = mean((minimum(z), maximum(z))), maximum(z) - minimum(z)

    xlim = (xmid - 0.65 * xrange, xmid + 0.65 * xrange)
    ylim = (ymid - 0.65 * yrange, ymid + 0.65 * yrange)
    zlim = (zmid - 0.65 * zrange, zmid + 0.65 * zrange)

    cmap = cgrad(:viridis, alpha=1.0)

    # === XY ===
    h_xy = fit(Histogram, (x, y), nbins=nbins)
    wxy = h_xy.weights
    wxy_masked = map(v -> v == 0.0 ? NaN : v, wxy)
    xcenters_xy = diff(h_xy.edges[1]) ./ 2 .+ h_xy.edges[1][1:end-1]
    ycenters_xy = diff(h_xy.edges[2]) ./ 2 .+ h_xy.edges[2][1:end-1]
    p1 = heatmap(xcenters_xy, ycenters_xy, wxy_masked';
        xlabel="x", ylabel="y", title="XY Heatmap",
        xlims=xlim, ylims=ylim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === XZ ===
    h_xz = fit(Histogram, (x, z), nbins=nbins)
    wxz_masked = map(v -> v == 0.0 ? NaN : v, h_xz.weights)
    xcenters_xz = diff(h_xz.edges[1]) ./ 2 .+ h_xz.edges[1][1:end-1]
    zcenters_xz = diff(h_xz.edges[2]) ./ 2 .+ h_xz.edges[2][1:end-1]
    p2 = heatmap(xcenters_xz, zcenters_xz, wxz_masked';
        xlabel="x", ylabel="z", title="XZ Heatmap",
        xlims=xlim, ylims=zlim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === YZ ===
    h_yz = fit(Histogram, (y, z), nbins=nbins)
    wyz_masked = map(v -> v == 0.0 ? NaN : v, h_yz.weights)
    ycenters_yz = diff(h_yz.edges[1]) ./ 2 .+ h_yz.edges[1][1:end-1]
    zcenters_yz = diff(h_yz.edges[2]) ./ 2 .+ h_yz.edges[2][1:end-1]
    p3 = heatmap(ycenters_yz, zcenters_yz, wyz_masked';
        xlabel="y", ylabel="z", title="YZ Heatmap",
        xlims=ylim, ylims=zlim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === 3D scatter ===
    p4 = scatter(x, y, z, marker_z=e, ms=2,
        xlabel="x", ylabel="y", zlabel="z", title="3D Scatter",
        xlims=xlim, ylims=ylim, zlims=zlim,
        colorbar_title="Energy", legend=false, cgrad=cmap)

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
end

# ╔═╡ 5d3197e2-1e5c-4f70-aa7a-db88a4c57c2a
plot_hits(evtdf; nbins=100)

# ╔═╡ Cell order:
# ╠═6eee4840-a8e3-11f0-9db8-b93ffc85840e
# ╠═88cdb157-b57a-40b6-ab3c-3ab8bbf1564d
# ╠═e8f2e87c-30cf-424f-a562-0f3f61ac80f1
# ╠═45fc0354-c4e7-45ac-8467-609b28cf7a61
# ╠═d75adcd8-2b38-4c02-b14f-4dfc4280a855
# ╠═16e942a5-20ed-4937-a5e6-987236a1e76d
# ╠═eb956cef-2a99-41fd-ad54-bb5272f8d7b7
# ╠═35dab5a3-75ec-4df6-8443-3ed9537e1dd7
# ╠═d892bb34-e552-4da2-ae1b-c096ab429709
# ╠═795cdf62-31a8-497b-9cee-f1ea4077eaf9
# ╠═92754e5b-ae80-4699-b1e6-3ec665f08248
# ╠═61da2c80-a9a7-4ea0-a855-356c0c300021
# ╠═e35b6e57-f943-4488-9743-d2231095ba89
# ╠═77b0fe9d-a457-4d36-8ebd-7c2b8ed4389e
# ╠═88a1e362-2ff8-4fa0-8b28-1057c311a668
# ╠═7f0e2846-74ef-4a62-82e9-070e4252514c
# ╠═0b28053d-2403-4160-9e58-4d74932700ec
# ╠═11b4e71e-51c0-4104-a821-3c6c9178a9ae
# ╠═527a8d40-0d63-42e7-9647-3d30730b1654
# ╠═e465533f-8416-4b4a-9db5-5d76eb38520c
# ╠═e0347b2c-8f8c-46dd-a857-01a036b70073
# ╠═cc400dae-2b04-4941-8fea-300426618bed
# ╠═2ea5696e-6aa0-4f16-9ca3-c9ee4c1f421f
# ╠═5d3197e2-1e5c-4f70-aa7a-db88a4c57c2a
# ╠═22ac54a3-e73e-4a19-a0ac-31f5d24d9beb
# ╠═83b8a99b-c8a8-4f70-a630-fe95491f3223
# ╠═34bcea1a-c135-4289-a45e-c51319f664d3
# ╠═83f79fa7-ce89-4b18-b27e-fa0d2f2a1cfc
# ╠═9f1a3620-2ed3-44f7-a9c8-de69fea363e9
# ╠═2520366e-7a26-49fb-a1bf-f2c431499de6
# ╠═17662c19-6bff-46ca-9533-1f3651a4062c
# ╠═23dc9540-589f-4d17-8f07-ab4f8004daad
# ╠═e2e96b6c-9b2f-492c-9313-5f1f92c24dbb
# ╠═692345c0-ea48-455c-98c6-90ef935054ea
# ╠═86deb2ca-f847-41fe-a30f-67469a14b93d
# ╠═e9cd96bf-b38d-4b0a-b1d6-5fee14e9a2be
# ╠═93ab09d6-a439-45b7-8d5f-0eb23702d3cb
# ╠═0c020931-0ff4-4ccc-8440-61fd7c6470f2
# ╠═bf476c92-549c-4d9e-a301-a13aacd8dfd8
# ╠═482911a0-da45-4851-b73b-62faea930087
# ╠═b5de052a-149a-416f-8f59-e42f97f3d5eb
# ╠═84812b46-ee00-40d2-8d0d-b59066ef7db4
# ╠═6061084b-a47a-4621-bfe3-96bab683f685
# ╠═8f9d90d8-d20b-42b6-a784-ccb2516f2158
# ╠═9b5eec23-6302-4919-9f7f-07e9565fea96
# ╠═3c65af6a-ad59-42c0-9fb8-be4b2830d7b5
# ╠═825d578e-6662-456a-96aa-afd6f2107475
# ╠═3ee80fd1-d3d6-4fca-a415-e37561001a2a
# ╠═2a789937-bd3b-45dc-83ed-7d62abcfbfc7
# ╠═d9ffc411-2c3e-4317-9571-2a7b0039f362
# ╠═6fd72e99-7a10-47a6-866f-0f75d373bab0
