using Plots
import Measures

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
    p2 = plot(ls, n3f.pdf.(ls),
            label = "NAPH3",
            xlabel = "λ (nm)",
            ylabel = "PDF",
            title = "NAPH3 emission PDF",
            legend = :topright)

			  
	p2 = plot!(p2, ls, n3c.pdf.(ls),
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
	p4 = plot(ls, n3f.xs.(ls),
			  label = "NAPH3",
		      xlabel = "λ (nm)",
		      ylabel = "\$\\sigma\$ (barn)",
			  title = "NAPH3 \$\\sigma\$ (barn)",
		      legend = :topright)
	p4 = plot!(p4, ls, n3c.xs.(ls),
			  label = "NAPH3-Ba")
	plot(p1,p2, p3, p4, layout=(2,2), size=(700, 500))
end



"""
    plot_molecule_emission(n3f::FMolecule, n3c::FMolecule, filters::Dict{String, Vector{Float64}})

Plot emission spectra for two FMolecule objects (e.g., NAPH3 and NAPH3-Ba), and overlay
vertical lines representing the start (red, dashed) and end (green, dashed) of each filter range.

# Arguments
- `n3f`: First molecule with emission spectrum
- `n3c`: Second molecule with emission spectrum
- `filters`: Dictionary mapping filter names to `[λ_start, λ_end]` wavelength ranges

# Returns
Combined plot with both spectra and filter markers.
"""
function plot_molecule_emission(n3f::FMolecule, n3c::FMolecule, filters::Dict{String, Vector{Float64}})
    # Plot 1: NAPH3
    p1 = plot(n3f.em_spectrum[1], n3f.em_spectrum[2],
        label = "NAPH3",
        xlabel = "λ (nm)",
        ylabel = "a.u.",
        title = "NAPH3 measured emission spectrum",
        legend = :topright,
        lw = 1.5)

    # Plot 2: NAPH3-BA
    p2 = plot(n3c.em_spectrum[1], n3c.em_spectrum[2],
        label = "NAPH3-Ba",
        xlabel = "λ (nm)",
        ylabel = "a.u.",
        title = "NAPH3-BA measured emission spectrum",
        legend = :topright,
        lw = 1.5)

    # Add vertical filter markers to both plots
    for (_name, range) in filters
        λ_start, λ_end = range

        vline!(p1, [λ_start], label = "", color = :red, lw = 1.0, linestyle = :dash)
        vline!(p1, [λ_end],   label = "", color = :green, lw = 1.0, linestyle = :dash)
        vline!(p2, [λ_start], label = "", color = :red, lw = 1.0, linestyle = :dash)
        vline!(p2, [λ_end],   label = "", color = :green, lw = 1.0, linestyle = :dash)
    end

    return plot(p1, p2, layout = (1, 2), size = (800, 500))
end


function plot_frames(imst::AbstractArray{T,3}; nscale::Int=20) where {T<:Real}

	FF = []

	
    for i in 1:9
        fn = (i-1) * nscale + i
		if fn > size(imst)[3]
			warn("requested frame = $(fn) is to large, set smaller nscale")
			fn = size(imst)[3] -i
			warn("set fn = $(fn)")
		end
        push!(FF, heatmap(imst[:, :, fn],
            colorbar=false,  # Optional: removes extra space
            title="Frame $fn",
            titlefontsize=7,
            tickfontsize=6,
            guidefontsize=6,
            titlelocation=:left,
            aspect_ratio=:equal))
    end

    plot(FF...;
        layout=(3, 3),
        size=(900, 900),
        margin=1.0*Measures.mm,
        top_margin=1.0*Measures.mm,
        bottom_margin=1.0*Measures.mm,
        left_margin=1.0*Measures.mm,
        right_margin=1.0*Measures.mm,
        plot_titlefontsize=7,
        legendfontsize=6)
end


function plot_frames_with_centroids(imst::AbstractArray{T,3}, region_stats::Vector{Vector{Dict{Symbol, Any}}}; 
	nscale::Int=20) where {T<:Real}

	FF = []

	
    for i in 1:9
        fn = (i-1) * nscale + i
		if fn > size(imst)[3]
			warn("requested frame = $(fn) is to large, set smaller nscale")
			fn = size(imst)[3] -i
			warn("set fn = $(fn)")
		end
		frame = imst[:, :, fn]
		centroids = [r[:centroid] for r in region_stats[fn]]
	
		p = heatmap(frame, color=:grays, title="Frame $fn",
		titlefontsize=7,
		tickfontsize=6,
		guidefontsize=6,
		titlelocation=:left,
		aspect_ratio=:equal)
		p = scatter!(p, [c[1] for c in centroids], [c[2] for c in centroids], markersize=4, color=:red)
        push!(FF, pp)
    end

    plot(FF...;
        layout=(3, 3),
        size=(900, 900),
        margin=1.0*Measures.mm,
        top_margin=1.0*Measures.mm,
        bottom_margin=1.0*Measures.mm,
        left_margin=1.0*Measures.mm,
        right_margin=1.0*Measures.mm,
        plot_titlefontsize=7,
        legendfontsize=6)
end


function plot_traces(TRZS, peaks; ftrz=1, ltrz=9,  figsize=(1500,1500))
    function select_trace(TRZS, df::DataFrame, row::Int) 
        i = df.i[row]
        j = df.j[row]
        trace = TRZS[i, j]  # safe, only defined where needed
        i,j, trace
    end
	function pltd(tz, i, j)
		plot(1:length(tz), tz, 
		label="Trace =($(i),$(i))", color=:gray, lw=1,
		xlabel="time steps", ylabel="Intensity", title="", 
		legend=:topright, grid=true)
	end
	PP =[]
	ntrz = ltrz - ftrz + 1
	
	if ntrz <= 9 
		ly = (3,3)
	elseif ntrz <= 16 
		ly = (4,4)
	elseif  ntrz <= 25 
		ly = (5,5)
	end

	for it in 1:ntrz
		i,j, tz = select_trace(TRZS, peaks, it)
		
		push!(PP, pltd(tz, i, j))
	end
	plot(PP..., layout=ly, size=figsize)
end


function plot_traces_stats(xsum, xmean, xstd; meanmx=25.0, summx=1e+3, stdmx=50.0)
	xmin = 0.0
	xmax = meanmx
	data = vec(xmean)
	filtered = data[(data .>= xmin) .& (data .<= xmax)]
	edges = range(xmin, xmax; length=101)  # 50 bins
	p1 = histogram(filtered, bins=collect(edges), title="Mean of traces", xlabel="trace index", ylabel="Frequency")

	xmin = 0.0
	xmax = summx
	data = vec(xsum)
	filtered = data[(data .>= xmin) .& (data .<= xmax)]
	edges = range(xmin, xmax; length=101)  # 50 bins
	p2 = histogram(filtered, bins=collect(edges), title="Sum of traces", xlabel="trace index", ylabel="Frequency")

	xmin = 0.0
	xmax = stdmx
	data = vec(xstd)
	filtered = data[(data .>= xmin) .& (data .<= xmax)]
	edges = range(xmin, xmax; length=101)  # 50 bins
	p3 = histogram(filtered, bins=collect(edges), title="std of traces", xlabel="trace index", ylabel="Frequency")

	plot(p1, p2, p3;
     layout=(3, 1),
     size=(700, 900),
     margin=5 * Measures.mm,
     bottom_margin=5 * Measures.mm,
     top_margin=5 * Measures.mm,
     plot_titlefontsize=10,
     guidefontsize=9,
     tickfontsize=8,
     legendfontsize=8)
end

function plot_stats(total_intensity, mean_intensity, std_intensity)
	# Plot mean and std over time
	n_frames = length(total_intensity)
	# Create a 3-row, 1-column plot layout
	plot_layout = @layout [a; b; c]

	p1 = plot(1:n_frames, total_intensity,
          label="Total Intensity", xlabel="Frame", ylabel="Total",
          title="Total Intensity Over Time", lw=1)

	p2 = plot(1:n_frames, mean_intensity,
          label="Mean Intensity", xlabel="Frame", ylabel="Mean",
          title="Mean Intensity Over Time", lw=1)

	p3 = plot(1:n_frames, std_intensity,
          label="Std Dev", xlabel="Frame", ylabel="Std",
          title="Standard Deviation Over Time", lw=1)

	# Combine them in a single figure
	#plot(p1, p2, p3, layout=plot_layout, size=(700, 900))
	plot(p1, p2, p3;
     layout=(3, 1),
     size=(700, 900),
     margin=5 * Measures.mm,
     bottom_margin=5 * Measures.mm,
     top_margin=5 * Measures.mm,
     plot_titlefontsize=10,
     guidefontsize=9,
     tickfontsize=8,
     legendfontsize=8)
end


function plot_traces_h2d(xsum, xmean, xstd; bins=50)
	function h2d(x1, x2, x1l, x2l, xt)
		# Ensure matrices have the same shape
		@assert size(x1) == size(x2)
	
		# Flatten both to vectors
		x = vec(x1)
		y = vec(x2)
	
		# 2D histogram
		histogram2d(x, y;
		       bins=bins,
		       xlabel=x1l,
		       ylabel=x2l,
		       title=xt)
	end
	p1 = h2d(xsum, xstd, "xsum", "xstd", "Sum vs Std")
	p2 = h2d(xmean, xstd, "mean", "xstd", "Mean vs Std")
	plot(p1, p2;
     layout=(1, 2),
     margin=1 * Measures.mm,
     bottom_margin=5 * Measures.mm,
     top_margin=5 * Measures.mm,
     plot_titlefontsize=8,
     guidefontsize=9,
     tickfontsize=6,
     legendfontsize=8)
end

