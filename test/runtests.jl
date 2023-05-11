using EmeraldStudio
using Test


@testset verbose = true "EmeraldStudio" begin
    @testset verbose = true "EmeraldStudio NameSpace" begin
        include("namespace/config.jl");
        include("namespace/states.jl");
    end;
end;
