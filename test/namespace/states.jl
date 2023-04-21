using Emerald.EmeraldLand.NameSpace

@testset verbose = true "State Variables" begin
    for FT in [Float32, Float64]
        for static in [true, false]
            NameSpace.use_static_arrays!(static);

            state = NameSpace.MultipleLayerSPACState{FT}();

            static ? (@test isbits(state)) : (@test true);
        end;
    end;
end
