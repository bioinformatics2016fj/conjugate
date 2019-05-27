#p，q，是枢轴元素
function PivotalTransform(A,p,q)
    if A[p,q]==0;
        error("ERROR:Pivotal element is 0.")
        return A
    end
        pv=A[p,:]./A[p,q]
        B=pv*(A[:,q])'
        B=A .- B'
        B[p,:]= pv
        return B
end
#select the entering basis
#r:reduced cost
function SelectQ(r)
    q=argmin(r)
    if r[q]>0
        q>0
    end
    return q
end

#select the leaving basis
function SelectP(bv,qv)
    p=0
    min=Inf
    for i=1:length(bv)
        if qv[i]>0
            temp=bv[i]/qv[i]
            if temp<min
                p=i
                min=temp
            end
        end
    end
    return p
end

#16.2
A=[3 1  0  1 4;
   6 2  1  1 5;
   2 -1 -1 0 0];
q1 = SelectQ(A[end,2:3])
println("q1=",q1)
p1=SelectP(A[1:2,end],A[1:2,q1])
println("p1=",p1)
A1=PivotalTransform(A, p1, q1)

a=(2:5)
q2=SelectQ(A1[end,a])
println("q2=",q2)
p2=SelectP(A1[1:2,end],A1[1:2,q2])
println("p2=",p2)
A2=PivotalTransform(A1, p2, q2)

b=(2,3,5)
q3=SelectQ(A2[end,b])
println("q3=",q3)
p3=SelectP(A2[1:2,end],A2[1:2,q3])
println("p3=",p3)
A3=PivotalTransform(A2, p3, q3)

#16.3
A=[1  0  1 1;
   0  1  1 2;
   -1 -1 3 0];
q1 = SelectQ(A[end,1:2])
println("q1=",q1)
p1=SelectP(A[1:2,end],A[1:2,q1])
println("p1=",p1)
A1=PivotalTransform(A, p1, q1)

a=2
q2=SelectQ(A1[end,a])
println("q2=",q2)
p2=SelectP(A1[1:2,end],A1[1:2,q2])
println("p2=",p2)
A2=PivotalTransform(A1, p2, q2)

b=2
q3=SelectQ(A2[end,b])
println("q3=",q3)
p3=SelectP(A2[1:2,end],A2[1:2,q3])
println("p3=",p3)
A3=PivotalTransform(A2, p3, q3)

#16.4
A=[1  1  0 0 5;
   0  1  1 0 7;
   1  1  0 1 9;
   -2 -1 0 0 0];
q1=SelectQ(A[end,1:2])
println("q1=",q1)
p1=SelectP(A[1:3,end],A[1:3,q1])
println("p1=",p1)
A1=PivotalTransform(A, p1, q1)
