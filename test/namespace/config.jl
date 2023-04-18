using Emerald.EmeraldLand.NameSpace

@testset verbose = true "Configuration Structs" begin
    for FT in [Float32, Float64]
        for var in [NameSpace.CanopyStructure{FT}(),
                    NameSpace.GSVSoilAlbedo{FT}(NameSpace.LAND_2021),
                    NameSpace.HyperspectralAbsorption{FT}(NameSpace.LAND_2017),
                    NameSpace.HyperspectralAbsorption{FT}(NameSpace.LAND_2021),
                    NameSpace.HyperspectralRadiation{FT}(NameSpace.LAND_2021),
                    NameSpace.UniversalConstants{FT}(),
                    NameSpace.WaveLengthSet{FT}(NameSpace.LAND_2021),
                    NameSpace.EmeraldConfiguration{FT}(NameSpace.LAND_2021),]
            @test isbits(var);
        end;
    end;
end
