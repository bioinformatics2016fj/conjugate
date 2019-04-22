function ConjugateGradientFSO(f, g, h, x0;
        ϵx=0.01, # precision for step size
        ϵf=0.01,
        ϵg=0.01,
        debug=false)
    #check arguments
    n, m=size(h)
    if n≠m
        error("ERROR: Matrix H is not square!")
    end
    steps=[]
    maxIterations = n
    xk = x0
    fk = f(xk...)
    gk = g(xk...)
    dk = -gk
    dh = dk'*h*dk
    α =  -dk'*gk/dh
    δ =  α .*dk
    #xn = xk .+ δ
    #fn = f(xn...)
    #gn = g(xn...)
    push!(steps, xk)
    if (norm(gk)<=ϵg)
        println("Convergence is reached after 1 iteration.")
        return xk, steps
    end
    for i in 1:maxIterations
        # iteration
        xn = xk .+ δ
        push!(steps, xn)
        fn = f(xn...)
        gn = g(xn...)
        βn = dk'*h*gn/dh
        dn = -gn .+ βn.*dk
        dh = dn'*h*dn
        α  = -dn'*gn/dh
        δ  = α.*dn
        # convegence?
        if (norm(gn)<=ϵg)
            println("Convergence is reached after ", i, " iterations.")
            return xn, steps
        end
        xk = xn
        fk = fn
        dk = dn
        if debug
            println("i=",i," x=", xn, " α=", α, " β=", βn, " gn=", gn, " d=", dn, " δ= ",δ)
        end
    end
    #println("WARN:", maxIterations, " iterations have been exceeded!")
end
(x,ConjugateGradientFSO_steps)=ConjugateGradientFSO(
    (x1,x2)->x1^2+(34/2)*x2^2+3*x1*x2-x1-30*x2,
    (x1,x2)->[2*x1+3*x2-1,34*x2+3*x1-30],
    [2 3;3 34],
    [1.,1.],
    debug=false
)


function steepest_gradient(f,g,x0;
        ϵx=0.01,
        ϵf=0.01,
        ϵg=0.01,
        maxIterations=128,
        debug=false)
    xk=x0
    fk=f(xk...)
    steepest_steps=[]
    for i in 1: maxIterations
        d=-g(xk...)
        α=search_for_alpha(f,xk,fk,d,-d)
        δ=α*d
        xn=xk.+δ
        fn=f(xn...)
        push!(steepest_steps,xk)
        if (norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg)
             println("Convergence is reached after",i,"iterations.")
        return steepest_steps
        end
         if debug
             println("i=",i,"xk=",xk,"d=",d,"δ=",δ)
        end
    xk=xn
    fk=fn
    end
    println("WARN",maxIterations,"iterations have been exceeded")
end
function search_for_alpha(f,xk,fk,d,g;
        α0=100,
        ϵ=0.5,
        t=0.5,)
    α=α0
    φ0=d'*g
    while f((xk.+α*d)...)>fk+ϵ*α*φ0
       α=t*α
    end
    return α
end
steepest_steps=steepest_gradient(
    (x1,x2)->x1^2+(34/2)*x2^2+3*x1*x2-x1-30*x2,
    (x1,x2)->[2*x1+3*x2-1,34*x2+3*x1-30],
    maxIterations = 100000,
    [1.,1.],
    debug=false
)


function Newton(f,g,h,x0;
        ϵx=0.01,
        ϵf=0.01,
        ϵg=0.01,
        maxIterations=128,
        debug=false)
    xk=x0
    fk=f(xk...)
    Newton_steps=[]
    for i in 1: maxIterations
        gk= g(xk...)
        d=-inv(h(xk...))*gk
        α = search_for_alpha(f,xk,fk,d,-d)
        δ=α*d
        xn=xk.+δ
        fn=f(xn...)
        push!(Newton_steps,xn)
        if (norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg)
              println("Convergence is reached after",i,"iterations.")
        return Newton_steps
        end
         if debug
        println("i=",i,"xk=",xk,"xn",xn,"d=",d,"δ=",δ)
        println("fk=", fk, "\tfn=", f)
        end
        xk=xn
        fk=fn
    end
    println("WARN",maxIterations,"iterations have been exceeded")
end
function search_for_alpha(f,xk,fk,d,g;
        α0=100,
        ϵ=0.5,
        t=0.5,)
    α=α0
    φ0=d'*g
    while f((xk.+α*d)...)>fk+ϵ*α*φ0
       α=t*α
    end
    return α
end
Newton_steps=Newton(
    (x1,x2)->x1^2+(34/2)*x2^2+3*x1*x2-x1-30*x2,
    (x1,x2)->[2*x1+3*x2-1,34*x2+3*x1-30],
    (x1,x2)->[2 3;3 34],
    [1.,1.],
    maxIterations = 10000000,
    debug=false
)


using Gadfly
fun=layer((x,y)->x^2+(34/2)*y^2+3*x*y-x-30*y, -2, 2, -2, 2);
conjugateGradientFSO=layer(
    x=[ConjugateGradientFSO_steps[i][1] for i in 1:length(ConjugateGradientFSO_steps)],
    y=[ConjugateGradientFSO_steps[i][2] for i in 1:length(ConjugateGradientFSO_steps)],
    label=[string(i) for i in 1:length(ConjugateGradientFSO_steps)],
    Geom.point, Geom.line, Geom.label,Theme(default_color=colorant"red"));
steepest=layer(
    x=[steepest_steps[i][1] for i in 1:length(steepest_steps)],
    y=[steepest_steps[i][2] for i in 1:length(steepest_steps)],
    label=[string(i) for i in 1:length(steepest_steps)],
    Geom.point, Geom.line, Geom.label,Theme(default_color=colorant"blue"));
newton=layer(
    x=[Newton_steps[i][1] for i in 1:length(Newton_steps)],
    y=[Newton_steps[i][2] for i in 1:length(Newton_steps)],
    label=[string(i) for i in 1:length(Newton_steps)],
    Geom.point, Geom.line, Geom.label,Theme(default_color=colorant"green"));
plot(fun,newton,steepest,conjugateGradientFSO)
