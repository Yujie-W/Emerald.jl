using EmeraldStudio.EmeraldLand.NameSpace

@testset verbose = true "Configuration Structs" begin
    for FT in [Float32, Float64]
        for static in [true, false]
            NameSpace.use_static_arrays!(static);
            for var in [NameSpace.CanopyStructure{FT}(),
                        NameSpace.GSVSoilAlbedo{FT}(),
                        NameSpace.HyperspectralAbsorption{FT}(),
                        NameSpace.HyperspectralAbsorption{FT}(),
                        NameSpace.HyperspectralRadiation{FT}(),
                        NameSpace.UniversalConstants{FT}(),
                        NameSpace.WaveLengthSet{FT}(),
                        NameSpace.EmeraldConfiguration{FT}(),
                        NameSpace.EmeraldConfiguration{FT}(NameSpace.LAND_2017),]
                static ? (@test isbits(var)) : (@test true);
            end;
        end;
    end;
end
