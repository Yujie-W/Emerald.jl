using Emerald
using Test


@testset verbose = true "Emerald" begin
    @testset verbose = true "Emerald NameSpace" begin
        include("namespace/config.jl");
        include("namespace/states.jl");
    end;
end;
