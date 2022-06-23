with(convex):
with(LinearAlgebra):
with(ComplexityOnePackage):
with(CodeTools):

# A list of irredundant KK^*-surfaces, all isomorphic to each other
Xs1 := [
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;-1,-2,2,2,1>)),
    ComplexityOneVariety(PMatrix(1, <-3,3,0,0,0;-3,0,1,1,0;-3,0,0,0,2;2,2,-1,-2,1>)),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;-2,-1,2,2,1>)),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;-3,-4,2,2,5>)),
    ComplexityOneVariety(PMatrix(1, <-2,3,0,0,0;-2,0,3,0,0;-2,0,0,1,1;1,2,2,-2,-1>)),
    ComplexityOneVariety(PMatrix(1, <-2,3,0,0,0;-2,0,3,0,0;-2,0,0,1,1;-1,-2,-2,2,1>)),
    ComplexityOneVariety(PMatrix(1, <-2,3,0,0,0;-2,0,3,0,0;-2,0,0,1,1;-7,1,-2,4,3>)),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;2,1,-1,-1,-1>)) 
];

for i from 1 to nops(Xs1) do
    for j from 1 to nops(Xs1) do
        Test(areIsomorphic(Xs1[i], Xs1[j]), true, label = cat("areIsomorphicTest-1-", i, "-", j));
        as := areIsomorphic(Xs1[i], Xs1[j], 'operations');
        for a in as do
            Test(Equal(applyAdmissibleOperation(Xs1[i], a):-P:-mat, Xs1[j]:-P:-mat), true, label = cat("areIsomorphicOperationTest-1-", i, "-", j));
        end do;
    end do;
end do;

# A list of very similar KK^*-surfaces, but all non-isomorphic to each other
Xs2 := [
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;-1,-2,2,2,1>)),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;0,-2,2,2,1>)),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;3,-2,2,2,1>)),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,2,0;-1,-1,0,0,3;-1,-2,2,1,1>)),
    ComplexityOneVariety(PMatrix(1, <-1,-1,4,0,0;-1,-1,0,3,0;-1,-1,0,0,2;-1,-2,3,2,1>))
];

for i from 1 to nops(Xs2) do
    for j from 1 to nops(Xs2) do
        if i <> j then
            Test(areIsomorphic(Xs2[i], Xs2[j]), false, label = cat("areIsomorphicTest-2-", i, "-", j));
        end if;
    end do;
end do;

# A list of KK^*-surfaces with coefficient matrices, all isomorphic to each other
Xs3 := [
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;-1,-2,2,2,1>), <-I, 4/7, 2+9*I, 2; 2/3, 9/5, -1, 1/31>),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;-1,-2,2,2,1>), <-1100 - 1200*I, -722+40*I, -3800+600*I, 2410+620*I;620003-5700*I, 391360+2/5*I, -17102-306209*I,42780+31/5*I>),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;2,1,-1,-1,-1>), <-I, 2+9*I, 4/7, 2; 2/3, -1, 9/5, 1/31>),
    ComplexityOneVariety(PMatrix(1, <-2,3,0,0,0;-2,0,3,0,0;-2,0,0,1,1;-7,1,-2,4,3>), <-15/31+30*I, -8+2*I, -9/5+4/7*I, 1/3; 30/31+15/31*I, -2-I, 18/5+9/5*I, 4/3+2/3*I>)
];

for i from 1 to nops(Xs3) do
    for j from 1 to nops(Xs3) do
        Test(areIsomorphic(Xs3[i], Xs3[j]), true, label = cat("areIsomorphicTest-3-", i, "-", j));
        as := areIsomorphic(Xs3[i], Xs3[j], 'operations');
        for a in as do
            Test(Equal(applyAdmissibleOperation(Xs3[i], a):-P:-mat, Xs3[j]:-P:-mat), true, label = cat("areIsomorphicOperationTest-3-", i, "-", j));
        end do;
    end do; 
end do;

# A list of KK^*-surfaces with coefficient matrix, but all non-isomorphic to each other
Xs4 := [
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;-1,-2,2,2,1>), <-I, 4/7, 2+9*I, 2; 2/3, 9/5, -1, 1/31>),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;-1,-2,2,2,1>), <-1101 - 1200*I, -722+40*I, -3800+600*I, 2410+620*I;620003-5700*I, 391360+2/5*I, -17102-306209*I,42780+31/5*I>),
    ComplexityOneVariety(PMatrix(1, <-1,-1,3,0,0;-1,-1,0,3,0;-1,-1,0,0,2;2,1,-1,-1,-1>), <-I, 4/7, 2, 2+9*I; 2/3, 9/5, 1/31, -1>),
    ComplexityOneVariety(PMatrix(1, <-2,3,0,0,0;-2,0,3,0,0;-2,0,0,1,1;-7,-2,1,4,3>), <-15/31+30*I, -8+2*I, -9/5+4/7*I, 1/4; 30/31+15/31*I, -2-I, 18/5+9/5*I, 4/3+2/3*I>)
];

for i from 1 to nops(Xs4) do
    for j from 1 to nops(Xs4) do
        if i <> j then
            Test(areIsomorphic(Xs4[i], Xs4[j]), false, label = cat("areIsomorphicTest-4-", i, "-", j));
        end if;
    end do;
end do;

# A list of threefolds with the same P-Matrix, but all non-isomorphic to each other
# TODO: Extend this ?
Xs5 := [
    ComplexityOneVariety(PMatrix(2,<-1,-1,-1,-1,2,0;-1,-1,-1,-1,0,2;1,-1,-1,1,0,1;1,1,-1,-1,1,0>), [0,1]),
    ComplexityOneVariety(PMatrix(2,<-1,-1,-1,-1,2,0;-1,-1,-1,-1,0,2;1,-1,-1,1,0,1;1,1,-1,-1,1,0>), [1,4]),
    ComplexityOneVariety(PMatrix(2,<-1,-1,-1,-1,2,0;-1,-1,-1,-1,0,2;1,-1,-1,1,0,1;1,1,-1,-1,1,0>), [-1,4])
];

for i from 1 to nops(Xs5) do
    for j from 1 to nops(Xs5) do
        if i <> j then
            Test(areIsomorphic(Xs5[i], Xs5[j]), false, label = cat("areIsomorphicTest-5-", i, "-", j));
        end if;
    end do;
end do;



