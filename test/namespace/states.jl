using Emerald.EmeraldLand.NameSpace

@testset verbose = true "State Variables" begin
    for FT in [Float32, Float64]
        state = NameSpace.MultipleLayerSPACState{FT,10}();
        @test isbits(state);
    end;
end
