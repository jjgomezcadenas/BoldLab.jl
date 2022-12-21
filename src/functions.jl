using Colors
using Images
using ImageBinarization
"""
Given a absolute path of a directory returns a list of files (of type dfiles) found
in the directory.

"""
function return_files(path::String, 
	                  dfiles="*.csv")
	nxfiles=filter(fname->fname[1]!=only("."),readdir(path))
	nxfiles
	
end

"""
Given the absolute path of the run returns the filter names by searching in the "Dark" folder.
"""
function flt_names(runp::String, dflt::String="Dark")
	p=joinpath(runp,dflt)
	xnames=return_files(p)
	fxnb = [split(name, "_")[2] for name in xnames]
	fxint = sort([parse(Int64, ff) for ff in fxnb])
	[string("Filter_", string(i)) for i in fxint]
end

"""
Given the absolute path of the run returns the filter names by searching in the "Filter1" folder.
"""
function point_names(runp::String,dflt::String="Filter1")
	p=joinpath(runp, dflt)
	xnames=return_files(p)
	ns=[String(split(pd, "_")[1]) for pd in xnames]
end

"""
Returns the full path of the image given the run path, point name and filter name.
"""
function get_image_path(runp::String, point_name::String, flt_name::String)
	ppath=joinpath(runp,point_name)
	fs=return_files(ppath)
	fname=fs[findall([occursin(flt_name,name) for name in fs])]
	path=joinpath(ppath,fname[1])
end

"""
Returns the image that corresponds to the selected point and filter in the run.
"""
function get_image(runp::String, point_name::String, flt_name::String)
	path=get_image_path(runp,point_name,flt_name)
	imgdf = DataFrame(CSV.File(path, header=false,delim="\t"))
	df1=Matrix(imgdf)
end

"""
Returns the nonfiltered image given a point and a run.
"""
function get_nfimage(runp::String, point_name::String, dflt::String="Filter1")
	path=get_image_path(runp,dflt,point_name)
	imgdf = DataFrame(CSV.File(path, header=false,delim="\t"))
	df1=Matrix(imgdf)
end

"""
Returns the dark image given the filter and run. 'darkfolder' allows to select the dark set. 
When darkfolder="Dark", the darks are obtained before the measurements, and if "Dark_after", at the end of the measurements 
"""
function get_dark(runp::String, flt_name::String, darkfolder::String="Dark")
	path=get_image_path(runp,darkfolder,flt_name)
	imgdf = DataFrame(CSV.File(path, header=false,delim="\t"))
	df1=Matrix(imgdf)
end
"""
Computes the edge of the image using sujoy algorithm, and binarizes it using Otsu threshold.
"""
function Image_edge(nfim::Matrix)
	nfimn=Float64.(nfim./maximum(nfim))
	img_edge=Float64.(sujoy(nfimn,four_connectivity=true))
	img_edgeb=binarize(img_edge,Otsu())
	iedge = Tuple.(findall(x->x==1,img_edgeb))  #indexed of the edge
	return img_edgeb, iedge
end

"""
Coomputes the edge of the image.
"""
function sujoy(img; four_connectivity=true)
    img_channel = Gray.(img)

    min_val = minimum(img_channel)
    img_channel = img_channel .- min_val
    max_val = maximum(img_channel)

    if max_val == 0
        return img
    end

    img_channel = img_channel./max_val

    if four_connectivity
        krnl_h = centered(Gray{Float32}[0 -1 -1 -1 0; 0 -1 -1 -1 0; 0 0 0 0 0; 0 1 1 1 0; 0 1 1 1 0]./12)
        krnl_v = centered(Gray{Float32}[0 0 0 0 0; -1 -1 0 1 1;-1 -1 0 1 1;-1 -1 0 1 1;0 0 0 0 0 ]./12)
    else
        krnl_h = centered(Gray{Float32}[0 0 -1 0 0; 0 -1 -1 -1 0; 0 0 0 0 0; 0 1 1 1 0; 0 0 1 0 0]./8)
        krnl_v = centered(Gray{Float32}[0 0 0 0 0;  0 -1 0 1 0; -1 -1 0 1 1;0 -1 0 1 0; 0 0 0 0 0 ]./8)
    end

    grad_h = imfilter(img_channel, krnl_h')
    grad_v = imfilter(img_channel, krnl_v')

    grad = (grad_h.^2) .+ (grad_v.^2)

    return grad
end
	
