
function create_dir!(dir)
	if isdir(dir) == false
		mkdir(dir)
	end
end


function output_dirs!(odir, sexp)
	create_dir!(joinpath(odir, sexp))
	csvdir = joinpath(odir, sexp,"csv")
	pngdir = joinpath(odir, sexp, "png")
	create_dir!(csvdir)
	create_dir!(pngdir)
	csvdir, pngdir
end
