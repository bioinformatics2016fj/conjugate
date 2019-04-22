using Gadfly
using LinearAlgebra
norm(x) = sqrt(x'*x)
##
function search_for_alpha(f, xk, fk, d, g; α0=100, ϵ=0.5, τ=0.5)
    α = α0
    ϕ0= d'*g
    while f((xk .+ α*d)...) > fk + ϵ*α*ϕ0
        α = τ*α
    end
    return α
end
##
function steepest_descent(f,g,x0;
        ϵx=0.01,
        ϵf=0.01,
        ϵg=0.01,
        maxIterations=128,
        debug=false
    )
    xk=x0
    fk=f(xk...)
    ###赋予空
    x0_iteration=[]
    for i in 1:maxIterations
        # iteration
        ###将参数按入数组
        push!(x0_iteration,xk)
        ##
        d=-g(xk...)
        α=search_for_alpha(f,xk,fk,d,-d)
        δ=α*d
        xn=xk.+δ  #######.+实行并行运算，增快运行速率
        fn=f(xn...)
        # convegence?
        if (norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg)
            println("Convergence is reached after ",i," interations.")
            return x0_iteration
        end
        if debug
            println("i=",i,"xk=",xk,"d=",d,"δ=",δ)
        end
        xk=xn
        fk=fn
    end
    println("WARN:",maxIterations,"iterations have been exceeded!")
end
##
function Resizable_Newton(f, g, h, x0;
        ϵx=0.01, # precision for step size
        ϵf=0.01,
        ϵg=0.01,
        maxIterations=128,
        debug=false)
    xk = x0
    fk = f(xk...)
    ##赋予空
    x0_iteration=[]
    for i in 1:maxIterations
        ###将参数按入数组
        push!(x0_iteration,xk)
        # iteration
        gk= g(xk...)
        d =-inv(h(xk...))*gk
        α = search_for_alpha(f, xk, fk, d, gk, α0=1)
        δ = α*d
        xn = xk .+ δ
        fn = f(xn...)
        # convegence?
        if (norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg)
            println("Convergence is reached after ", i, " iterations.")
            return x0_iteration
        end
        if debug
            println("i=",i, " α=", α, " xk=", xk, " xn=", xn, " d=", d, " δ= ",δ)
            println("fk=", fk, "\tfn=", fn)
        end
        xk = xn
        fk = fn
    end
    println("WARN:", maxIterations, " iterations have been exceeded!")
end
##
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
    ###赋予空，将参数按入数组
    x0_iteration=[]
    push!(x0_iteration,xk)
    ####
    if (norm(gk)<=ϵg)
        println("Convergence is reached after 1 iteration.")
        return x0_iteration
        #return xk, fk, gk, steps
    end
    for i in 1:maxIterations
        # iteration
        ###
        ###
        xn = xk .+ δ
        ###将参数按入数组
        push!(x0_iteration,xn)
        ###
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
            return x0_iteration
            ##return xn, fn, gn, steps
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
##
function eggplot(all_iteration_point,color)###color为你想要的颜色
    col1=[all_iteration_point[i][1] for i in 1:length(all_iteration_point)]##提取第一行值
    col2=[all_iteration_point[i][2] for i in 1:length(all_iteration_point)]##提取第二行值
    first_step=[col1[1:end] col2[1:end]]##这个为迭代的前一个点
    next_step=[[col1[2:(end)];col1[end]] [col2[2:(end)];col2[end]]] ##迭代后的点，这是为了后面向量箭头有指向
    xk=[first_step next_step]
    layer1=layer(xk,x=xk[:,1], y=xk[:,2], xend=xk[:,3], yend=xk[:,4],
        Geom.point, Geom.line,Geom.vector,##使用向量表示箭头
        label=[string(i) for i in 1:length(col1)],Geom.label,
        Theme(discrete_highlight_color=x->"darkred", default_color=color)
    )
    return layer1###分别画出每个函数的layer
end
##
function egg_plot_all_method(f,g,h,x0,color1,color2,color3)###color为你想要的颜色
    method_s=steepest_descent(
        f,
        g,
        x0,
        maxIterations=2500,
        debug=true
    )

    method_r=Resizable_Newton(
        f,
        g,
        (x1,x2)->h,
        x0,
        debug=true
    )

    method_c=ConjugateGradientFSO(
        f,
        g,
        h,
        x0,
        debug=false
    )
    xsc  = Scale.x_continuous(minvalue=-1, maxvalue=1)
    ysc  = Scale.y_continuous(minvalue=0.5, maxvalue=1)
    layer1=eggplot(method_s,color1)
    layer2=eggplot(method_r,color2)
    layer3=eggplot(method_c,color3)
    layer4=layer(f,-1.,1.,0.5,1.)
    return plot(layer1,layer2,layer3,layer4,xsc,ysc,
        Guide.manual_color_key("Method", ["steepest_descent", "Resizable_Newton","ConjugateGradientFSO"],
        [color1,color2,color3])
    )###最终绘图函数
end
##
egg_plot_all_method(
    (x1,x2)->x1^2+7.5*x2^2+3x1*x2-x1-11*x2,
    (x1,x2)->[2*x1+3*x2-1; 15*x2+3*x1-11],
    [2 3;3 15],
    [1.,1.],
    "deepskyblue", "pink","green"
)
