module DiscreteCurve

using LinearAlgebra
using GLMakie

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
    for i = 2:n-1
        # the gradient direction is the normal vector at the point
        # the magnitude is 2 * sin(theta/2) where theta is the angle between the two tangent vectors
        # compute differences between points
        x_i = x[i] - x[i-1]
        x_i_1 = x[i+1] - x[i]
        y_i = y[i] - y[i-1]
        y_i_1 = y[i+1] - y[i]

        length_i = sqrt(x_i^2 + y_i^2)
        length_i_1 = sqrt(x_i_1^2 + y_i_1^2)

        # theta sign is important
        theta = acos((x_i * x_i_1 + y_i * y_i_1) / (length_i * length_i_1))
        theta = theta * sign(x_i * y_i_1 - x_i_1 * y_i)

        # compute gradient
        u[i-1] = (y[i+1] - y[i-1]) / (length_i + length_i_1) * 2 * sin(theta / 2)
        v[i-1] = -(x[i+1] - x[i-1]) / (length_i + length_i_1) * 2 * sin(theta / 2)
    end
    ### END HOMEWORK PROBLEM ###

    # Plot curve and gradient vectors
    fig1, ax1, plt1 = lines(x, y, linewidth=3, color=:red, figure=(resolution=(1000,1000),))
    arrows!(x[2:end-1], y[2:end-1], u, v, arrowsize=0.05)
    display(fig1)

    ## Problem 2(d)
    kappa = zeros(n-2)
    ### YOUR CODE HERE TO COMPUTE KAPPA ###
    for i = 2:n-1
        # the curvature is the magnitude of the gradient
        kappa[i-1] = sqrt(u[i-1]^2 + v[i-1]^2)
    end
    ## END HOMEWORK PROBLEM ###

    curvcolor = [kappa[1]; kappa; kappa[end]]
    fig2, ax2, plt2 = lines(x, y, linewidth=4, color=curvcolor, figure=(resolution=(1000,1000),))
    cbar = Colorbar(fig2, plt2)
    cbar.width = 30
    fig2[1, 1] = ax2
    fig2[1, 2] = cbar
    display(fig2)
end

#problem2cd()

function problem2e()
    ## Problem 2(e)
    t0 = 0
    t1 = pi * 1.25
    nsamples = 100
    # Modify nsteps if your method does not converge
    nsteps = 20000

    # We provide a few examples of curves to try
     curveFunction(t) = [cos(t)-cos(3*t).^3 sin(t)-sin(3*t).^3]
    # curveFunction(t) = [cos(t) sin(t)]
    # curveFunction(t) = [t (t.-t0).*(t1.-t)]
    curve = vcat(curveFunction.(range(t0, stop=t1, length=nsamples))...)

    dispcurve = Node(Point2f0.(eachrow(curve)))
    fig, ax, plt = lines(dispcurve, linewidth=3, figure=(resolution=(1000,1000),))
    display(fig)
    u_prev = zeros(nsamples - 2)
    v_prev = zeros(nsamples - 2)
    for i = 1:nsteps
        ### YOUR CODE HERE TO PERFORM GRADIENT DESCENT ###
        # compute gradient
        u = zeros(nsamples - 2)
        v = zeros(nsamples - 2)
        for j = 2:nsamples-1
            # the gradient direction is the normal vector at the point
            # the magnitude is 2 * sin(theta/2) where theta is the angle between the two tangent vectors
            # compute differences between points
            x_i = curve[j, 1] - curve[j-1, 1]
            x_i_1 = curve[j+1, 1] - curve[j, 1]
            y_i = curve[j, 2] - curve[j-1, 2]
            y_i_1 = curve[j+1, 2] - curve[j, 2]

            length_i = sqrt(x_i^2 + y_i^2)
            length_i_1 = sqrt(x_i_1^2 + y_i_1^2)

            # theta sign is important
            theta = acos((x_i * x_i_1 + y_i * y_i_1) / (length_i * length_i_1))
            theta = theta * sign(x_i * y_i_1 - x_i_1 * y_i)

            # compute gradient
            u[j-1] = (y_i_1 - y_i) / (length_i + length_i_1) * 2 * sin(theta / 2)
            v[j-1] = -(x_i_1 - x_i) / (length_i + length_i_1) * 2 * sin(theta / 2)
        end

        # store velocity to use in momentum
        if i == 1
            u_prev = u
            v_prev = v
        end
        # update velocity
        u = 0.995 * u_prev + 0.005 * u
        v = 0.995 * v_prev + 0.005 * v

        # Update curve position with gradient descent with momentum
        for j = 2:nsamples-1
            curve[j, 1] = curve[j, 1] - 0.5 * v[j-1] * i^1.1/5000
            curve[j, 2] = curve[j, 2] - 0.5 * u[j-1] * i^1.1/5000

        end

        u_prev = u
        v_prev = v

        ### END HOMEWORK PROBLEM ###
        sleep(1/300000000)
        dispcurve[] = Point2f0.(eachrow(curve))
    end
end

problem2e()

end

