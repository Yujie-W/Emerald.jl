using Emerald.EmeraldLand.NameSpace

@testset verbose = true "Configuration Structs" begin
    for FT in [Float32, Float64]
        consts = NameSpace.UniversalConstants{FT}();
        @test isbits(consts);
        config = NameSpace.HyperspectralAbsorption{FT}(NameSpace.LAND_2017);
        @test isbits(config);
        config = NameSpace.HyperspectralAbsorption{FT}(NameSpace.LAND_2021);
        @test isbits(config);
    end;
end
