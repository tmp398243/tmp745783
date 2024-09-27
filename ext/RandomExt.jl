"This module extends Lorenz63Filter with functionality from Random."
module RandomExt

using Lorenz63Filter: Lorenz63Filter
using Random

"""
    greeting()

Call [`Lorenz63Filter.greeting`](@ref) with a random name.


# Examples

```jldoctest
julia> @test true;

```

"""
Lorenz63Filter.greeting() = Lorenz63Filter.greeting(rand(5))

end
