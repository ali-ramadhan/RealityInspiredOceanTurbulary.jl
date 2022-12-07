using Oceananigans: @compute

function mixed_layer_depth(model)
    b = BuoyancyField(model)
    @compute ∂b̄∂z = Field(Average(∂z(b), dims=(1, 2)))

    _, I_ml = findmax(abs.(interior(∂b̄∂z)))

    mixed_layer_depth = - znode(Face(), I_ml.I[3], model.grid)

    return mixed_layer_depth
end
