using Pkg
using Gadfly

##3种方法（最速下降法、牛顿法、共轭梯度法）的收敛情况图
function Three_Trends(
    f,g,h,x0;
            ϵx=0.01,
            ϵf=0.01,
            ϵg=0.01,
            maxIterations=128,
            debug=false
    )
    paint1=steepest_descent(f,g,x0)##最速下降法的趋势图 绿色的线
    paint2=newton(f,g,h,x0)##牛顿法的趋势图 红色的线
    paint3=ConjugateGradientFSO(f,g,h,x0)##共轭梯度法的趋势图 黑色的线
    paint=layer(f,-1,1,0,2)##函数图
    return(plot(paint,paint1,paint2,paint3))##把三个方法的收敛情况画出
end

Three_Trends(
        (x,y)->x^2+13.5*y^2+3*x*y-x-23y,
        (x,y)->[2x+3y-1,27y+3x-23],
        [2 3;3 27],
        [1.,1.],
        debug=false
        )


function search_for_alpha(f,xk,fk,d,g;α0=100,ϵ=0.5,τ=0.5)
   α=α0
   ϕ0=d'*g
   while f((xk.+α*d)...)>fk+ϵ*α*ϕ0
       α=τ*α
   end
   return α
end
norm(x)=sqrt(x'*x)
function steepest_descent(f,g,x0;
       ϵx=0.01,
       ϵf=0.01,
       ϵg=0.01,
       maxIterations=128,
       debug=false
   )
   println("f=",f)
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
       xn=xk.+δ  #######.+实行并行运算，增快运行速率
       my_xy=push!(my_xy,xn)##存xk值的array
       fn=f(xn...)
       # convegence?
       if ((norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg))||i==maxIterations
           println("Convergence is reached after ",i," interations.")
           my_xy=vcat(my_xy'...)
           println("vcated",my_xy)
          point=layer(x=my_xy[1:end,1],y=my_xy[1:end,2], label=[ string(i) for i in 1:length(my_xy[1:end,1])],
                    Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"green"))
           return (point)##最速梯度法的点趋势的图

       end
       if debug
           println("i=",i,"xk=",xk,"d=",d,"δ=",δ,"xn=",xn)
       end
       xk=xn
       fk=fn
   end
   println("WARN: ",maxIterations," iterations have been exceeded!")
end



function newton(f,g,h,x0;
            ϵx=0.01,
            ϵf=0.01,
            ϵg=0.01,
            maxIterations=128,
            debug=false
        )
        xk=x0
        my_xy=[]
        my_xy=push!(my_xy,xk)
        fk=f(xk...)

        for i in 1:maxIterations
            # iteration
            d=-inv(h)*g(xk...)  ##-F(x)T*g(x)
            α=1
            δ=α*d
            xn=xk.+δ  #######  新的迭代点
            my_xy=push!(my_xy,xn)
            fn=f(xn...)
            # convegence?
            if (norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg)||i==maxIterations
                println("Convergence is reached after ",i," interations.")
                my_xy=vcat(my_xy'...)
                println("vcated",my_xy)
                point1=layer(x=my_xy[1:end,1],y=my_xy[1:end,2], label=[ string(i) for i in 1:length(my_xy[1:end,1])],
                          Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"red"))
                return (point1)##牛顿法的点趋势的图

            end
            if debug
                println("i=",i,"xk=",xk,"d=",d,"δ=",δ,"xn=",xn)
            end
            xk=xn
            fk=fn
        end
        println("WARN: ",maxIterations," iterations have been exceeded!")
    end

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
            push!(steps, xk)##steps是存点的array
            if (norm(gk)<=ϵg)
                println("Convergence is reached after 1 iteration.")
                steps=vcat(steps'...)
                point2=layer(x=steps[1:end,1],y=steps[1:end,2], label=[ string(i) for i in 1:length(steps[1:end,1])],
                          Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"black"))
                return (point2)##共轭梯度法的点趋势的图
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
                    steps=vcat(steps'...)
                    point2=layer(x=steps[1:end,1],y=steps[1:end,2], label=[ string(i) for i in 1:length(steps[1:end,1])],
                              Geom.point,Geom.line,Geom.label,Theme(default_color=colorant"black"))
                    return (point2)##共轭梯度法的点趋势的图
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
