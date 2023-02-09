module MVU

using LinearAlgebra, SparseArrays, Convex, SCS, WGLMakie

function runMVU(k, n=1000)
    x, theta = swissroll(n)
    y = mvu(x, k)
    fig = Figure(resolution=(1000,1000))
    ax = LScene(fig)
    scatter!(Point3f0.(eachcol(y)), color = theta, markersize=200)
    fig[1, 1] = ax
    display(fig)
    return y
end

function swissroll(n)
    theta = 1.5 .* pi .+ 3 .* pi .* rand(Float64, n);
    z = 10 * rand(Float64, n);
    return [theta' .* [cos.(theta'); sin.(theta')]; z'], theta;
end

### PROBLEM 4 - YOUR CODE HERE
function mvu(x, k)
    n = size(x, 2)
    y = x

    return y
end
### END HOMEWORK PROBLEM

end