#共轭法
norm(x) = sqrt(x'*x)
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
α = -dk'*gk/dh
δ = α .*dk
#xn = xk .+ δ
#fn = f(xn...)
#gn = g(xn...)
push!(steps, xk)
if (norm(gk)<=ϵg)
println("Convergence is reached after 1 iteration.")
return xk, fk, gk, steps
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
α = -dn'*gn/dh
δ = α.*dn
# convegence?
if (norm(gn)<=ϵg)
println("Convergence is reached after ", i, " iterations.")
return xn, fn, gn, steps
end
xk = xn
fk = fn
dk = dn
if debug
println("i=",i," x=", xn, " α=", α, " β=", βn, " gn=", gn, " d=", dn)
end
end
#println("WARN:", maxIterations, " iterations have been exceeded!")
end

#画图准备
xn, fn, gn, steps=
ConjugateGradientFSO(
(x,y)->x^2+13/2y^2+3x*y-x-9y,
(x,y)->[2x+3y-1, 3x+13y-9],
[2 3; 3 13],
[1.,1.],
debug=false
)

#画图
using Gadfly
fun=layer((x,y)->x^2+13/2y^2+3x*y-x-9y, -1.2, 1.2, -1.2, 1.2);
proc=layer(
x=[steps[i][1] for i in 1:length(steps)],
y=[steps[i][2] for i in 1:length(steps)],
label=[string(i) for i in 1:length(steps)],
Geom.point, Geom.line, Geom.label);
plot(fun, proc)

function Resizable_Newton(f, g, h, x0;
ϵx=0.01, # precision for step size
ϵf=0.01,
ϵg=0.01,
maxIterations=128,
debug=false)
xk = x0
fk = f(xk...)
for i in 1:maxIterations
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
return xk, fk, d, δ
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

#画图准备
xn, fn, gn, steps=
Resizable_Newton(
(x,y)->x^2+13/2y^2+3x*y-x-9y,
(x,y)->[2x+3y-1, 3x+13y-9],
[2 3; 3 13],
[1.,1.],
maxIterations = 10000000,
debug=false
)


#画图
using Gadfly
fun=layer((x,y)->x^2+13/2y^2+3x*y-x-9y, -1.2, 1.2, -1.2, 1.2);
proc=layer(
x=[steps[i][1] for i in 1:length(steps)],
y=[steps[i][2] for i in 1:length(steps)],
label=[string(i) for i in 1:length(steps)],
Geom.point, Geom.line, Geom.label);
plot(fun, proc)

norm(x) = sqrt(x'*x)
function steepest_descent(f, g, x0;
ϵx=0.01, # precision for step size
ϵf=0.01,
ϵg=0.01,
maxIterations=128,
debug=false)
xk = x0
fk = f(xk...)
for i in 1:maxIterations
# iteration
d =-g(xk...)
α = search_for_alpha(f, xk, fk, d, -d)
δ = α*d
xn = xk .+ δ
fn = f(xn...)
# convegence?
if (norm(δ)<=ϵx)&&(abs(fn-fk)<=ϵf)&&(norm(d)<=ϵg)
println("Convergence is reached after ", i, " iterations.")
return xk, fk, d, δ
end
if debug
println("i=",i, " α=", α, " xk=", xk, " d=", d, " δ= ",δ)
end
xk = xn
fk = fn
end
println("WARN:", maxIterations, " iterations have been exceeded!")
end

function search_for_alpha(f, xk, fk, d, g; α0=100, ϵ=0.5, τ=0.5)
α = α0
ϕ0= d'*g
while f((xk .+ α*d)...) > fk + ϵ*α*ϕ0
α = τ*α
end
return α
end

#画图准备
xn, fn, gn, steps=
steepest_descent(
(x,y)->x^2+13/2y^2+3x*y-x-9y,
(x,y)->[2x+3y-1, 3x+13y-9],
[2.0 3.0; 3.0 13.0],
maxIterations = 10000,
debug=false
)

#画图
using Gadfly
fun=layer((x,y)->x^2+13/2y^2+3x*y-x-9y, -1.2, 1.2, -1.2, 1.2);
proc=layer(
x=[steps[i][1] for i in 1:length(steps)],
y=[steps[i][2] for i in 1:length(steps)],
label=[string(i) for i in 1:length(steps)],
Geom.point, Geom.line, Geom.label);
plot(fun, proc)
