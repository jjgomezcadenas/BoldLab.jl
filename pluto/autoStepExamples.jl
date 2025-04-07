# ╔═╡ new-cell-load-trace
data_path = "/Users/jjgomezcadenas/Projects/BoldLab/data/trace_sample.txt"  # Adjust filename
dataX = vec(readdlm(data_path, Float64))  # load as Float64 vector

# ╔═╡ new-cell-run-stepfinder
FitX, Fits, S_curves, best_shots, steptable = auto_stepfinder(dataX; demo=0.0)

# ╔═╡ new-cell-plot-fit
plot(dataX, label="Original Trace", lw=1, legend=:topright)
plot!(Fits[end, :], label="Step Fit", lw=2)

# ╔═╡ new-cell-show-table
steptable
