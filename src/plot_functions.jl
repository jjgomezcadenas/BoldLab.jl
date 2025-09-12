using Plots
import Measures
using StatsBase: percentile

# Fit types are available from fit_functions.jl which is included before this file
# ExponentialFit, PeltFit, and AsfFit are defined in the same module scope

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


# Internal helper function to determine color limits for better contrast
function _determine_color_limits(frame::AbstractArray{T,2}, clims, use_percentile, percentile_range) where {T<:Real}
    if clims == :auto
        if use_percentile
            # Use percentile range to exclude outliers
            vmin = percentile(vec(frame), percentile_range[1])
            vmax = percentile(vec(frame), percentile_range[2])
            return (vmin, vmax)
        else
            # Use min/max but exclude zeros for minimum if present
            non_zero = frame[frame .> 0]
            if !isempty(non_zero)
                return (minimum(non_zero), maximum(frame))
            else
                return (minimum(frame), maximum(frame))
            end
        end
    else
        return clims
    end
end

"""
    plot_frames(imst::AbstractArray{T,3}, frame_number::Int; kwargs...)

Plot a single frame from a 3D image stack.

# Arguments
- `imst`: 3D array of images (height × width × frames)
- `frame_number`: The frame number to plot (1-indexed)
- `clims`: Color limits - :auto (default), or tuple (min, max)
- `percentile_range`: Tuple for percentile clipping (default: (1, 99))
- `use_percentile`: Use percentile range for color limits (default: false)
- `colormap`: Colormap to use (default: :hot)
- `colorbar`: Show colorbar (default: true for single frame)
- `size`: Figure size tuple (default: (600, 500))
"""
function plot_frames(imst::AbstractArray{T,3}, frame_number::Int; 
                    clims=:auto,
                    percentile_range=(1, 99),
                    use_percentile=false,
                    colormap=:hot,
                    colorbar=true,
                    size=(600, 500)) where {T<:Real}
    
    # Check frame number is valid
    if frame_number < 1 || frame_number > size(imst, 3)
        error("Frame number $frame_number is out of range. Valid range: 1 to $(size(imst, 3))")
    end
    
    frame = imst[:, :, frame_number]
    color_limits = _determine_color_limits(frame, clims, use_percentile, percentile_range)
    
    plt = heatmap(frame,
        colorbar=colorbar,
        title="Frame $frame_number",
        xlabel="X",
        ylabel="Y",
        aspect_ratio=:equal,
        clims=color_limits,
        color=colormap,
        size=size)
    
    # Display plot when not in Pluto environment
    if !isdefined(Main, :PlutoRunner)
        display(plt)
    end
    
    return plt
end

"""
    plot_frames(imst::AbstractArray{T,3}; kwargs...)

Plot multiple frames from a 3D image stack in a 3×3 grid.

# Arguments
- `imst`: 3D array of images (height × width × frames)
- `nscale`: Interval between frames to plot (default: 20)
- `clims`: Color limits - :auto (default), or tuple (min, max)
- `percentile_range`: Tuple for percentile clipping (default: (1, 99))
- `use_percentile`: Use percentile range for color limits (default: false)
- `colormap`: Colormap to use (default: :grays)
"""
function plot_frames(imst::AbstractArray{T,3}; 
                    nscale::Int=20, 
                    clims=:auto,
                    percentile_range=(1, 99),
                    use_percentile=false,
                    colormap=:hot) where {T<:Real}

	FF = []

	
    for i in 1:9
        fn = (i-1) * nscale + i
		if fn > size(imst)[3]
			warn("requested frame = $(fn) is to large, set smaller nscale")
			fn = size(imst)[3] -i
			warn("set fn = $(fn)")
		end
        
        frame = imst[:, :, fn]
        color_limits = _determine_color_limits(frame, clims, use_percentile, percentile_range)
        
        push!(FF, heatmap(frame,
            colorbar=false,  # Optional: removes extra space
            title="Frame $fn",
            titlefontsize=7,
            tickfontsize=6,
            guidefontsize=6,
            titlelocation=:left,
            aspect_ratio=:equal,
            clims=color_limits,
            color=colormap))
    end

    plt = plot(FF...;
        layout=(3, 3),
        size=(900, 900),
        margin=1.0*Measures.mm,
        top_margin=1.0*Measures.mm,
        bottom_margin=1.0*Measures.mm,
        left_margin=1.0*Measures.mm,
        right_margin=1.0*Measures.mm,
        plot_titlefontsize=7,
        legendfontsize=6)
    
    # Display plot when not in Pluto environment
    if !isdefined(Main, :PlutoRunner)
        display(plt)
    end
    
    return plt
end


function plot_frames_with_centroids(imst::AbstractArray{T,3}, region_stats::Vector{Vector{Dict{Symbol, Any}}}; 
	nscale::Int=20,
	clims=:auto,
	percentile_range=(1, 99),
	use_percentile=false,
	colormap=:grays) where {T<:Real}

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
		
		color_limits = _determine_color_limits(frame, clims, use_percentile, percentile_range)
	
		p = heatmap(frame, color=colormap, title="Frame $fn",
		titlefontsize=7,
		tickfontsize=6,
		guidefontsize=6,
		titlelocation=:left,
		aspect_ratio=:equal,
		clims=color_limits)
		p = scatter!(p, [c[1] for c in centroids], [c[2] for c in centroids], markersize=4, color=:red)
        push!(FF, p)
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

"""
    plot_exponential_fit(trace::Vector, fit_result::ExponentialFit; kwargs...)

Plot original trace data and exponential fit.

# Arguments
- `trace`: Original data vector  
- `fit_result`: ExponentialFit from fit_traces_exponential

# Optional keyword arguments
- `data_color`: Color for data line (default: :gray)
- `data_linewidth`: Width of data line (default: 2)
- `data_marker`: Marker style for residuals plot (default: :circle)
- `data_markersize`: Size of residual markers (default: 4)
- `fit_color`: Color for fit line (default: :blue)
- `fit_linewidth`: Width of fit line (default: 2)
- `xlabel`: X-axis label (default: "Time step")
- `ylabel`: Y-axis label (default: "Intensity")
- `title`: Plot title (default: "Exponential Fit")
- `legend`: Show legend (default: true)
- `legend_position`: Legend position (default: :topright)
- `size`: Figure size (default: (600, 400))
- `show_residuals`: Show residuals subplot (default: false)
- `grid`: Show grid (default: true)

# Returns
- Plot object

# Example
```julia
trace = [100.0, 80.0, 64.0, 51.0, 41.0, 33.0, 26.0, 21.0]
fit_result = fit_traces_exponential(trace)
plot_exponential_fit(trace, fit_result)
```
"""
function plot_exponential_fit(trace::Vector, fit_result::ExponentialFit;
                             data_color=:gray,
                             data_linewidth=1,
                             data_marker=:circle,
                             data_markersize=4,
                             fit_color=:blue,
                             fit_linewidth=2,
                             xlabel="Time step",
                             ylabel="Intensity",
                             title="Exponential Fit",
                             legend=true,
                             legend_position=:topright,
                             size=(600, 400),
                             show_residuals=false,
                             grid=true)
    
    # Time points
    t = 1:length(trace)
    
    # Create main plot with data as continuous line
    p1 = plot(t, trace,
              seriestype=:line,
              color=data_color,
              linewidth=data_linewidth,
              label="Data",
              xlabel=xlabel,
              ylabel=ylabel,
              title=title,
              legend=legend ? legend_position : false,
              grid=grid,
              framestyle=:box)
    
    # Add fit line
    plot!(p1, t, fit_result.fitX,
          color=fit_color,
          linewidth=fit_linewidth,
          label="Fit: A₀=$(round(fit_result.A0, digits=2)), τ=$(round(fit_result.tau, digits=2)), C₀=$(round(fit_result.C0, digits=2))")
    
    # Note: ExponentialFit doesn't have converged field, so skip convergence check
    
    # Optional residuals subplot
    if show_residuals
        residuals = trace .- fit_result.fitX
        p2 = plot(t, residuals,
                  seriestype=:scatter,
                  color=data_color,
                  marker=data_marker,
                  markersize=data_markersize-1,
                  label="Residuals",
                  xlabel=xlabel,
                  ylabel="Residuals",
                  legend=false,
                  grid=grid,
                  framestyle=:box)
        
        # Add zero line
        hline!(p2, [0], color=:black, linestyle=:dash, linewidth=1)
        
        # Combine plots
        p = plot(p1, p2, layout=(2,1), size=size)
        return p
    else
        plot!(size=size)
        return p1
    end
end

"""
    plot_exponential_fit(trace::Vector, fit_result::Nothing; kwargs...)

Handle the case when fitting fails and returns nothing.
"""
function plot_exponential_fit(trace::Vector, fit_result::Nothing;
                             data_color=:gray,
                             data_marker=:circle,
                             data_markersize=4,
                             xlabel="Time step",
                             ylabel="Intensity",
                             title="Exponential Fit (Failed)",
                             legend=true,
                             legend_position=:topright,
                             size=(600, 400),
                             grid=true,
                             kwargs...)
    
    t = 1:length(trace)
    
    p = plot(t, trace,
             seriestype=:scatter,
             color=data_color,
             marker=data_marker,
             markersize=data_markersize,
             label="Data (Fit Failed)",
             xlabel=xlabel,
             ylabel=ylabel,
             title=title,
             legend=legend ? legend_position : false,
             grid=grid,
             framestyle=:box,
             size=size)
    
    return p
end

"""
    plot_pelt_fit(tz::Vector, step_fit::PeltFit; kwargs...)

Plot original trace data and PELT step fit results.

# Arguments
- `tz`: Original trace vector
- `step_fit`: PeltFit result from fit_pelt

# Optional keyword arguments
- `data_color`: Color for data line (default: :gray)
- `data_linewidth`: Width of data line (default: 2)
- `data_marker`: Marker style for residuals plot (default: :circle)
- `data_markersize`: Size of residual markers (default: 4)
- `fit_color`: Color for step fit line (default: :red)
- `fit_linewidth`: Width of step fit line (default: 1)
- `step_markers`: Show markers at step change points (default: true)
- `step_marker_color`: Color for step markers (default: :orange)
- `step_marker_size`: Size of step markers (default: 6)
- `xlabel`: X-axis label (default: "Time step")
- `ylabel`: Y-axis label (default: "Intensity")
- `title`: Plot title (default: "PELT Step Fit")
- `legend`: Show legend (default: true)
- `legend_position`: Legend position (default: :topright)
- `size`: Figure size (default: (600, 400))
- `show_residuals`: Show residuals subplot (default: false)
- `grid`: Show grid (default: true)

# Returns
- Plot object

# Example
```julia
trace = [10.0, 10.1, 9.9, 5.2, 5.1, 4.9, 2.1, 2.0, 1.9]
step_result = fit_pelt(trace)
plot_pelt_fit(trace, step_result; show_residuals=true)
```
"""
function plot_pelt_fit(tz::Vector, step_fit::PeltFit;
                       data_color=:gray,
                       data_linewidth=1,
                       data_marker=:circle,
                       data_markersize=4,
                       fit_color=:red,
                       fit_linewidth=2,
                       step_markers=true,
                       step_marker_color=:orange,
                       step_marker_size=6,
                       xlabel="Time step",
                       ylabel="Intensity", 
                       title="PELT Step Fit",
                       legend=true,
                       legend_position=:topright,
                       size=(600, 400),
                       show_residuals=false,
                       grid=true)
    
    # Time points
    t = 1:length(tz)
    
    # Create main plot with data as continuous line
    p1 = plot(t, tz,
              seriestype=:line,
              color=data_color,
              linewidth=data_linewidth,
              label="Data",
              xlabel=xlabel,
              ylabel=ylabel,
              title=title,
              legend=legend ? legend_position : false,
              grid=grid,
              framestyle=:box)
    
    # Create step function from StepFit results
    step_function = zeros(Float64, length(tz))
    
    # Build step function using stepTime and stepHeight
    if step_fit.nstep > 0 && length(step_fit.stepTime) > 0
        current_pos = 1
        
        for i in 1:step_fit.nstep
            step_end = min(step_fit.stepTime[i], length(tz))
            step_function[current_pos:step_end] .= step_fit.stepHeight[i]
            current_pos = step_end + 1
        end
        
        # Fill any remaining points with the last step height
        if current_pos <= length(tz)
            last_height = step_fit.stepHeight[end]
            step_function[current_pos:end] .= last_height
        end
    end
    
    # Add step fit line
    plot!(p1, t, step_function,
          color=fit_color,
          linewidth=fit_linewidth,
          label="Steps: $(step_fit.nstep) steps (noise region: $(step_fit.si):$(step_fit.sf))")
    
    # Add step change markers
    if step_markers && step_fit.nstep > 0
        step_times = step_fit.stepTime[step_fit.stepTime .<= length(tz)]
        step_intensities = [step_function[min(st, length(tz))] for st in step_times]
        
        scatter!(p1, step_times, step_intensities,
                color=step_marker_color,
                markersize=step_marker_size,
                label="Step changes",
                markerstrokewidth=1,
                markerstrokecolor=:black)
    end
    
    # Optional residuals subplot
    if show_residuals
        residuals = tz .- step_function
        p2 = plot(t, residuals,
                  seriestype=:scatter,
                  color=data_color,
                  marker=data_marker,
                  markersize=data_markersize-1,
                  label="Residuals",
                  xlabel=xlabel,
                  ylabel="Residuals",
                  legend=false,
                  grid=grid,
                  framestyle=:box)
        
        # Add zero line
        hline!(p2, [0], color=:black, linestyle=:dash, linewidth=1)
        
        # Combine plots
        p = plot(p1, p2, layout=(2,1), size=size)
        return p
    else
        plot!(size=size)
        return p1
    end
end

"""
    plot_pelt_fit(tz::Vector, step_fit::Nothing; kwargs...)

Handle the case when PELT fitting fails and returns nothing.
"""
function plot_pelt_fit(tz::Vector, step_fit::Nothing;
                       data_color=:gray,
                       data_marker=:circle,
                       data_markersize=4,
                       xlabel="Time step",
                       ylabel="Intensity",
                       title="PELT Step Fit (Failed)",
                       legend=true,
                       legend_position=:topright,
                       size=(600, 400),
                       grid=true,
                       kwargs...)
    
    t = 1:length(tz)
    
    p = plot(t, tz,
             seriestype=:line,
             color=data_color,
             linewidth=2,
             label="Data (No Steps Found)",
             xlabel=xlabel,
             ylabel=ylabel,
             title=title,
             legend=legend ? legend_position : false,
             grid=grid,
             framestyle=:box,
             size=size)
    
    return p
end

"""
    plot_asf_fit(tz::Vector, asf_fit::AsfFit; kwargs...)

Plot original trace data and ASF (Automatic Step Finder) fit results.

# Arguments
- `tz`: Original trace vector
- `asf_fit`: AsfFit result from fit_asf

# Optional keyword arguments
- `data_color`: Color for data line (default: :gray)
- `data_linewidth`: Width of data line (default: 2)
- `data_marker`: Marker style for residuals plot (default: :circle)
- `data_markersize`: Size of residual markers (default: 4)
- `fit_color`: Color for step fit line (default: :green)
- `fit_linewidth`: Width of step fit line (default: 1)
- `step_markers`: Show markers at step change points (default: true)
- `step_marker_color`: Color for step markers (default: :orange)
- `step_marker_size`: Size of step markers (default: 6)
- `xlabel`: X-axis label (default: "Time step")
- `ylabel`: Y-axis label (default: "Intensity")
- `title`: Plot title (default: "ASF Step Fit")
- `legend`: Show legend (default: true)
- `legend_position`: Legend position (default: :topright)
- `size`: Figure size (default: (600, 400))
- `show_residuals`: Show residuals subplot (default: false)
- `grid`: Show grid (default: true)

# Returns
- Plot object

# Example
```julia
trace = [10.0, 10.1, 9.9, 5.2, 5.1, 4.9, 2.1, 2.0, 1.9]
asf_result = fit_asf(trace)
plot_asf_fit(trace, asf_result; show_residuals=true)
```
"""
function plot_asf_fit(tz::Vector, asf_fit::AsfFit;
                      data_color=:gray,
                      data_linewidth=1,
                      data_marker=:circle,
                      data_markersize=4,
                      fit_color=:green,
                      fit_linewidth=2,
                      step_markers=true,
                      step_marker_color=:orange,
                      step_marker_size=6,
                      xlabel="Time step",
                      ylabel="Intensity", 
                      title="ASF Step Fit",
                      legend=true,
                      legend_position=:topright,
                      size=(600, 400),
                      show_residuals=false,
                      grid=true)
    
    # Time points
    t = 1:length(tz)
    
    # Create main plot with data as continuous line
    p1 = plot(t, tz,
              seriestype=:line,
              color=data_color,
              linewidth=data_linewidth,
              label="Data",
              xlabel=xlabel,
              ylabel=ylabel,
              title=title,
              legend=legend ? legend_position : false,
              grid=grid,
              framestyle=:box)
    
    # Add ASF fit line (use the fitted data from the AsfFit struct)
    plot!(p1, t, asf_fit.fitX,
          color=fit_color,
          linewidth=fit_linewidth,
          label="ASF: $(length(asf_fit.stepHeights)) steps (best_shot=$(asf_fit.best_shot))")
    
    # Add step change markers
    if step_markers && length(asf_fit.stepTimes) > 0
        step_times = asf_fit.stepTimes[asf_fit.stepTimes .<= length(tz)]
        step_intensities = [asf_fit.fitX[min(st, length(tz))] for st in step_times]
        
        scatter!(p1, step_times, step_intensities,
                color=step_marker_color,
                markersize=step_marker_size,
                label="Step changes",
                markerstrokewidth=1,
                markerstrokecolor=:black)
    end
    
    # Optional residuals subplot
    if show_residuals
        residuals = tz .- asf_fit.fitX
        p2 = plot(t, residuals,
                  seriestype=:scatter,
                  color=data_color,
                  marker=data_marker,
                  markersize=data_markersize-1,
                  label="Residuals",
                  xlabel=xlabel,
                  ylabel="Residuals",
                  legend=false,
                  grid=grid,
                  framestyle=:box)
        
        # Add zero line
        hline!(p2, [0], color=:black, linestyle=:dash, linewidth=1)
        
        # Combine plots
        p = plot(p1, p2, layout=(2,1), size=size)
        return p
    else
        plot!(size=size)
        return p1
    end
end

"""
    plot_asf_fit(tz::Vector, asf_fit::Nothing; kwargs...)

Handle the case when ASF fitting fails and returns nothing.
"""
function plot_asf_fit(tz::Vector, asf_fit::Nothing;
                      data_color=:gray,
                      data_marker=:circle,
                      data_markersize=4,
                      xlabel="Time step",
                      ylabel="Intensity",
                      title="ASF Step Fit (Failed)",
                      legend=true,
                      legend_position=:topright,
                      size=(600, 400),
                      grid=true,
                      kwargs...)
    
    t = 1:length(tz)
    
    p = plot(t, tz,
             seriestype=:line,
             color=data_color,
             linewidth=2,
             label="Data (No Steps Found)",
             xlabel=xlabel,
             ylabel=ylabel,
             title=title,
             legend=legend ? legend_position : false,
             grid=grid,
             framestyle=:box,
             size=size)
    
    return p
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

"""
    plot_crops_results(cro::Dict; title="CROPS Results")

Plot penalty and constrained cost as a function of number of changepoints from CROPS output.

# Arguments
- `cro::Dict`: CROPS output dictionary with keys "penalty", "number", "constrained", "changepoints"
- `title::String`: Plot title (default: "CROPS Results")

# Returns
- Plots.Plot object showing both penalty and constrained cost vs number of changepoints
"""
function plot_crops_results(cro::Dict; title="CROPS Results")
    # Extract data from CROPS dictionary
    numbers = cro["number"]
    penalties = cro["penalty"] 
    constrained = cro["constrained"]
    
    # Create the main plot with penalty on left y-axis
    p1 = plot(numbers, penalties, 
             label="Penalty", 
             color=:blue,
             linewidth=2,
             marker=:circle,
             markersize=4,
             xlabel="Number of Changepoints",
             ylabel="Penalty", 
             title=title,
             legend=:topright)
    
    # Create twin axis for constrained cost without legend
    p2 = twinx(p1)
    plot!(p2, numbers, constrained,
          color=:red, 
          linewidth=2,
          marker=:square,
          markersize=4,
          ylabel="Constrained Cost",
          legend=false,
          grid=false)  # Disable grid on twin axis to avoid overlap
    
    # Add the constrained cost to the main plot's legend manually
    plot!(p1, NaN*numbers[1:1], NaN*penalties[1:1], 
          label="Constrained Cost", 
          color=:red, 
          linewidth=2, 
          marker=:square,
          markersize=4)
    
    return p1
end

