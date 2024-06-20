using BoldLab
using Test


@testset "BoldLab.jl" begin
p1=joinpath([ENV["JLaserData"],"data","5CMOS","Quartz"])
p2=joinpath(p1,BoldLab.readdir(p1)[1])
runp=joinpath(p2,BoldLab.readdir(p2)[1])

#Check that enviroment variable "JLaserData" points the root folder of the project.
@test ("data" in readdir(ENV["JLaserData"]))
@test ("pdata" in readdir(ENV["JLaserData"]))
@test ("FLUORI" in readdir(ENV["JLaserData"]))
@test typeof(BoldLab.return_files(runp,"Dark"))==Vector{String}
[@test occursin("Filter",name) for name in BoldLab.flt_names(runp)]
[@test occursin("Point",name) for name in BoldLab.point_names(runp)]
@test maximum(BoldLab.Image_edge(BoldLab.get_nfimage(runp,BoldLab.point_names(runp)[1]))[1])==1
end
