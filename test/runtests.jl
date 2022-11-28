using BoldLab
using Test

@testset "BoldLab.jl" begin
    
#Check that enviroment variable "JLaserData" points the root folder of the project.
@test ("data" in readdir(ENV["JLaserData"]))
@test ("pdata" in readdir(ENV["JLaserData"]))
@test ("FLUORI" in readdir(ENV["JLaserData"]))

end
