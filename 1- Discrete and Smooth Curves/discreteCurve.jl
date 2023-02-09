module DiscreteCurve

using LinearAlgebra
using WGLMakie

# Sample Lissajous curve
a = 4
b = 2
delta = pi / 3
n = 100

t = range(0, stop=2*pi, length=n)
x = sin.(a .* t .* delta)
y = sin.(b .* t)
xy = [x y]

function problem2cd()
    ## Problem 2(c)
    n = length(x)
    u = zeros(n - 2)
    v = zeros(n - 2)
    ### YOUR CODE HERE TO COMPUTE GRADIENT ###
    ### END HOMEWORK PROBLEM ###

    # Plot curve and gradient vectors
    fig1, ax1, plt1 = lines(x, y, linewidth=3, color=:red, figure=(resolution=(1000,1000),))
    arrows!(x[2:end-1], y[2:end-1], u, v, arrowsize=0.05)
    display(fig1)
    
    ## Problem 2(d)
    kappa = zeros(n-2)
    ### YOUR CODE HERE TO COMPUTE KAPPA ###
    ## END HOMEWORK PROBLEM ###

    curvcolor = [kappa[1]; kappa; kappa[end]]
    fig2, ax2, plt2 = lines(x, y, linewidth=4, color=curvcolor, figure=(resolution=(1000,1000),))
    cbar = Colorbar(fig2, plt2)
    cbar.width = 30
    fig2[1, 1] = ax2
    fig2[1, 2] = cbar
    display(fig2)
end

function problem2e()
    ## Problem 2(e)
    t0 = 0
    t1 = pi * 1.25
    nsamples = 100
    # Modify nsteps if your method does not converge
    nsteps = 200

    # We provide a few examples of curves to try
    # curveFunction(t) = [cos(t)-cos(3*t).^3 sin(t)-sin(3*t).^3]
    # curveFunction(t) = [cos(t) sin(t)]
    curveFunction(t) = [t (t.-t0).*(t1.-t)]
    curve = vcat(curveFunction.(range(t0, stop=t1, length=nsamples))...)

    dispcurve = Node(Point2f0.(eachrow(curve)))
    fig, ax, plt = lines(dispcurve, linewidth=3, figure=(resolution=(1000,1000),))
    display(fig)
    for i = 1:nsteps
        ### YOUR CODE HERE TO PERFORM GRADIENT DESCENT ###
        ### END HOMEWORK PROBLEM ###
        sleep(1/30)
        dispcurve[] = Point2f0.(eachrow(curve))
    end
end

end
