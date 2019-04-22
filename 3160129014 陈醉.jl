norm(x)=sqrt(x'*x)
using Gadfly


function search_for_alpha(f,xk,fk,d,g;α0=100,ϵ=0.5,τ=0.5)
   α=α0
   ϕ0=d'*g
   while f((xk.+α*d)...)>fk+ϵ*α*ϕ0
       α=τ*α
   end
   return α
end

##把三个layer和原函数等高线画在一张图上，最速下降法趋势为红色，牛顿法趋势为粉色，共轭梯度法的趋势为蓝色
function Combind(f,g,h,x0;
    ϵx=0.01,
    ϵf=0.01,
    ϵg=0.01,
    maxIterations=128,
    debug=false)
    point=Steepest_Descent(f,g,x0)              ##最速下降法趋势的layer，红色线
    point1=Newton(f,g,h,x0)                     ##牛顿法趋势的layer，粉色线
    point2=ConjugateGradientFSO(f,g,h,x0)       ##共轭梯度法趋势的layer，蓝色线
    paint1=layer(f,-1,1,0,2)                    ##原函数“等高线”
    return(plot(paint1,point,point1,point2))    ##结果图
end


##输出结果图，最速下降法趋势为红色，牛顿法趋势为粉色，共轭梯度法的趋势为蓝色
Combind(
    (x,y)->x^2+8.5*y^2+3*x*y-x-14y,
    (x,y)->[2x+3y-1,17y+3x-14],
    [2 3;3 17],
    [1.,1.],
    debug=false)


    ##最速下降法的函数，return最速下降法趋势的layer
    function Steepest_Descent(f,g,x0;
        ϵx=0.01,
        ϵf=0.01,
        ϵg=0.01,
        maxIterations=128,
        debug=false)
        xk=x0
        fk=f(xk...)
        my_xy=[]
        my_xy=push!(my_xy,xk)
        for i in 1:maxIterations
            # iteration
            d=-g(xk...)
            println("d",d)
            α=search_for_alpha(f,xk,fk,d,-d;α0=100,ϵ=0.5,τ=0.5)
            println("α",α)
            δ=α*d
            xn=xk.+δ
            my_xy=push!(my_xy,xn)
            fn=f(xn...)
            # convegence?
            if ((norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg))||i==maxIterations
                println("Convergence is reached after ",i," interations.")
                my_xy=vcat(my_xy'...)
                println("vcated",my_xy)
                point=layer(x=my_xy[1:end,1],y=my_xy[1:end,2], label=[ string(i) for i in 1:length(my_xy[1:end,1])],
                            Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"red"))
                return (point)
           end
           if debug
               println("i=",i,"xk=",xk,"d=",d,"δ=",δ,"xn=",xn)
           end
           xk=xn
           fk=fn
       end
       println("WARN: ",maxIterations," iterations have been exceeded!")
    end


    #牛顿法的函数，return牛顿法趋势的layer
    function Newton(f,g,h,x0;
            ϵx=0.01,
            ϵf=0.01,
            ϵg=0.01,
            maxIterations=128,
            debug=false)
        xk=x0
        my_xy=[]
        my_xy=push!(my_xy,xk)
        fk=f(xk...)
        for i in 1:maxIterations
            # iteration
            d=-inv(h)*g(xk...)  ##-F(x)T*g(x)
            α=1
            δ=α*d
            xn=xk.+δ
            my_xy=push!(my_xy,xn)
            fn=f(xn...)
            # convegence?
            if (norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg)||i==maxIterations
                println("Convergence is reached after ",i," interations.")
                my_xy=vcat(my_xy'...)
                println("vcated",my_xy)
                point1=layer(x=my_xy[1:end,1],y=my_xy[1:end,2], label=[ string(i) for i in 1:length(my_xy[1:end,1])],
                            Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"pink"))
                return (point1)
            end
            if debug
                println("i=",i,"xk=",xk,"d=",d,"δ=",δ,"xn=",xn)
            end
            xk=xn
            fk=fn
        end
        println("WARN: ",maxIterations," iterations have been exceeded!")
    end


    #共轭梯度法的函数，return共轭梯度法趋势的layer
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
        push!(steps, xk)
        if (norm(gk)<=ϵg)
            println("Convergence is reached after 1 iteration.")
            steps=vcat(steps'...)
            point2=layer(x=steps[1:end,1],y=steps[1:end,2], label=[ string(i) for i in 1:length(steps[1:end,1])],
                    Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"blue"))
            return (point2)
        end
        for i in 1:maxIterations    # iteration
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
                steps=vcat(steps'...)
                point2=layer(x=steps[1:end,1],y=steps[1:end,2], label=[string(i) for i in 1:length(steps[1:end,1])],
                            Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"blue"))
                return (point2)
            end
            xk = xn
            fk = fn
            dk = dn
            if debug
                println("i=",i," x=", xn, " α=", α, " β=", βn, " gn=", gn, " d=", dn, " δ= ",δ)
            end
        end
        println("WARN:", maxIterations, " iterations have been exceeded!")
    end
