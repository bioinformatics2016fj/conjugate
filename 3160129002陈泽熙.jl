using Gadfly


norm(x)=sqrt(x'*x)

function search_for_alpha(f, xk, fk, d, g; α0=100, ϵ=0.5, τ=0.5)
    α = α0
    ϕ0= d'*g
    while f((xk .+ α*d)...) > fk + ϵ*α*ϕ0
        α = τ*α
    end
    return α
end

#牛顿法
function newton(f, g, h, x0;
        ϵx=0.01, # precision for step size
        ϵf=0.01,
        ϵg=0.01,
        maxIterations=128,
        debug=false)

    xk = x0
    fk = f(xk...)

    x=[]
    y=[]

    for i in 1:maxIterations
        # iteration
        #d =-g(xk...)
        #α = d'*d/(d'*h*d)
        #δ = α*d
        xn = xk .- h^(-1)*g(xk...)
        fn = f(xn...)
        push!(x,xk[1])
        push!(y,xk[2])
        # convegence?
        if (abs(fn-fk)<=ϵf)
            println("Convergence is reached after ", i, " iterations.")
            println("the_final_x=",xn," the_final_y=",fn)
            x=convert(Array{Float64,1},x)
            y=convert(Array{Float64,1},y)


            return x,y
        end
        if debug
            println("i=",i," xk=", xk,)
        end
        xk = xn
        fk = fn
    end
    println("WARN:", maxIterations, " iterations have been exceeded!")
end

newton_x,newton_y=newton(
    (x,y)->x^2+3x*y+3y^2-x-2y,
    (x,y)->[2x+3y-1,3x+6y-2],
    [2 3;3 6],
    [1,1],
    debug=false
    )

#最速下降法

function steepest_descent(f, g, x0;
        ϵx=0.001, # precision for step size
        ϵf=0.001,
        ϵg=0.001,
        maxIterations=128,
        debug=false)
    xk = x0
    fk = f(xk...)
    x=[]
    y=[]

    for i in 1:maxIterations
        # iteration
        d =-g(xk...)
        α = search_for_alpha(f, xk, fk, d, -d)
        δ = α*d
        xn = xk .+ δ
        fn = f(xn...)
        push!(x,xk[1])
        push!(y,xk[2])

        # convegence?
        if (norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg)
            println("Convergence is reached after ", i, " iterations.")
            println("the_final_x=",xn," the_final_y=",fn," d=",d," δ=",δ)
            x=convert(Array{Float64,1},x)
            y=convert(Array{Float64,1},y)



            return  x,y
        end
        if debug
            println("i=",i, " α=", α, " xk=", xk, " d=", d, " δ= ",δ)
        end
        xk = xn
        fk = fn
    end
    println("WARN:", maxIterations, " iterations have been exceeded!")
end
sd_x,sd_y=steepest_descent(
    (x,y)->x^2+3x*y+3y^2-x-2y,
    (x,y)->[2x+3y-1,3x+6y-2],
    [1,1] ,
    debug=false
    )
#共轭梯度法
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
    x=[]
    y=[]
    push!(x,xk[1])
    push!(y,xk[2])
    if (norm(gk)<=ϵg)
        println("Convergence is reached after 1 iteration.")
        println("the_final_x=",xk," the_final_y=",fk)

        x=convert(Array{Float64,1},x)
        y=convert(Array{Float64,1},y)
        return x,y
    end

    for i in 1:maxIterations
        # iteration
        xn = xk .+ δ
        fn = f(xn...)
        gn = g(xn...)
        βn = dk'*h*gn/dh
        dn = -gn .+ βn.*dk
        dh = dn'*h*dn
        α  = -dn'*gn/dh
        δ  = α.*dn
        push!(x,xn[1])
        push!(y,xn[2])
        # convegence?

        if (norm(gn)<=ϵg)
            println("Convergence is reached after ", i, " iterations.")
            println("the_final_x=",xn," the_final_y=",fn)
            x=convert(Array{Float64,1},x)
            y=convert(Array{Float64,1},y)
            return x,y
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
CG_x,CG_y=ConjugateGradientFSO(
    (x,y)->x^2+3x*y+3y^2-x-2y,
    (x,y)->[2x+3y-1,3x+6y-2],
    [2 3;3 6],
    [1,1],
    debug=false
    )
#画图
function final_figure(f,newton_x,newton_y,sd_x,sd_y,CG_x,CG_y)
    newton=layer(x=newton_x,y=newton_y, label=[ string(i) for i in 1:length(newton_x)],
                Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"red"))
    SD=layer(x=sd_x,y=sd_y, label=[ string(i) for i in 1:length(sd_x)],
            Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"blue"))
    GC=layer(x=CG_x,y=CG_y, label=[ string(i) for i in 1:length(CG_x)],
            Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"yellow"))
    y=layer(f,0,1,0,1)
    return plot(newton,SD,GC,y,Guide.manual_color_key("Method",["newton", "steepest_descent","ConjugateGradientFSO"],
            ["red","blue","yellow"]))
end

final_figure(
    (x,y)->x^2+3x*y+3y^2-x-2y,
    newton_x,
    newton_y,
    sd_x,
    sd_y,
    CG_x,
    CG_y
    )
