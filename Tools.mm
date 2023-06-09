(*
Some tools.
*)

swapCase := proc(c :: string)
    if c = "EE" then
        "EE";
    elif c = "PE" then
        "EP";
    elif c = "EP" then
        "PE";
    elif c = "PP" then
        "PP";
    else
        error "swapCase only accepts one of the four strings \"EE\", \"PE\", \"EP\" and \"PP\"";
    end if;
end proc;

# Given a list of integers, returns one of the four strings "A", "D", "E" and "X"
platonicityType := proc(ls)
    local sortedLs, x, y, z;
    sortedLs := sort(ls, `>`);
    if not andmap(l -> l = 1, sortedLs[4 .. nops(sortedLs)]) then
        return "X";
    end if;
    x,y,z := sortedLs[1], sortedLs[2], sortedLs[3];
    if z = 1 then
        return "A";
    end if;
    if y = 2 and z = 2 then
        return "D";
    end if;
    if y = 3 and z = 2 then
        return "E";
    else
        return "X";
    end if;
end proc;

(*
Returns the `i`-th canonical basis vector e_i in `n` dimensional space.
Here, the zeroth basis vector is understood to be -(e_1 + ... + e_r)
*)
canonicalBasisVector := proc(n :: integer, i :: integer)
    if n < 0 then
        error "n must be non-negative";
    end if;
    if i < 0 or i > n then
        error "i must range between 0 and n = %1", n;
    end if;
    if i = 0 then
        Vector([-1 $ n]);
    else
        Vector([0 $ i-1, 1, 0 $ n-i]);
    end if;
end proc;

applyPermToList := proc(p :: Perm, ls :: list)
    local i;
    return [seq(ls[(p^(-1))[i]], i = 1 .. nops(ls))];
end proc;

(* Computes all permutations leaving a given list invariant. *)
invariantPermutations := proc(ls :: list, compFun := `=`)
    local indexGroups, a, b, permsList, grp, perms, p, i;
    indexGroups := [ListTools[Categorize]((i,j) -> compFun(ls[i], ls[j]), [seq(1 .. nops(ls))])];
    permsList := [];
    for grp in indexGroups do
        perms := [];
        for a in map(combinat[permute], nops(grp)) do
            p := [seq(1 .. max(grp))];
            for i from 1 to nops(grp) do
                p[grp[i]] := grp[a[i]];
            end do;
            perms := [op(perms), Perm(p)];
        end do;
        permsList := [op(permsList), perms];
    end do; 
    return foldl((p1, p2) -> [seq(seq(a . b, b in p2), a in p1)], [Perm([[]])], op(permsList));
end proc;

sortLexComparison := proc(a :: {numeric, list, Array}, b :: {numeric, list, Array})
    local i;
    if type(a, numeric) and type(b, numeric) then
        return a > b;
    elif (type(a, list) or type(a, Array)) and (type(b, list) or type(b, Array)) then
        if nops(a) > nops(b) then
            return true;
        elif nops(a) < nops(b) then
            return false;
        end if;
        for i from lowerbound(a) to upperbound(a) do
            if sortLexComparison(a[i], b[i]) then
                return true
            elif sortLexComparison(b[i], a[i]) then 
                return false;
            end if;
        end do;
        return false;
    end if;
    error "Both arguments must either be numbers, lists or Arrays";
end proc;

sortLex := proc(ls :: {list, Array})
    if type(ls, 'list'(numeric)) or type(ls, 'Array'(numeric)) then
        sort(ls, `>`);
    else
        sort(map(sortLex, ls), sortLexComparison);
    end if;
end proc;

# For a polytope in two-dimensional space given by a list of vertices,
# check whether the origin is the only interior lattice point.
isAlmostEpsHollow := proc(vs, eps)
    local minX, maxX, minY, maxY, po, x, y;
    minX := min(map(v -> v[1], vs));
    maxX := max(map(v -> v[1], vs));
    minY := min(map(v -> v[2], vs));
    maxY := max(map(v -> v[2], vs));
    po := convhull(op(eps * vs));
    for x from ceil(minX) to floor(maxX) do
        for y from ceil(minY) to floor(maxY) do
            if not (x = 0 and y = 0) and containsrelint(po, [x, y]) then
                return false;
            end if;
        end do;
    end do;
    return true;
end proc;

# Returns the unique intersection points of the two lines in 2D-space that go through
# the points v1 and v2 as well as w1 and w2 respectively (if it exists)
intersectLines2D := proc(v1, v2, w1, w2)
    return LinearSolve(<v2[2] - v1[2], v1[1] - v2[1], v1[1] * v2[2] - v2[1] * v1[2]; w2[2] - w1[2], w1[1] - w2[1], w1[1] * w2[2] - w2[1] * w1[2]>);
end proc;


getSingleDiscrepancy := proc(a, b, v)
    local w;
    w := intersectLines2D(a, b, 0, v);
    (if w[2] = 0 then v[1] / w[1] else v[2] / w[2] end if) - 1; 
end proc;

# Given a 2x2 integer matrix A, computes the normal form of the two-dimensional
# cone spanned by the columns of A. Here, normal form in the sense of
# Proposition 10.1.1 of Cox, Little, Schenck.
# This procedure outputs the parameters d and k, as well as the change of
# basis matrix M such that M . A = <0, d; 1, -k> with d > 0, 0 <= k < d and gcd(d,k) = 1.
cone2DNormalForm := proc(A)

    local x, y, d, l, s, k, sg, M;

    # Get integers 'x' and 'y' such that s * A[1,1] + t * A[2,1]
    igcdex(A[1,1], A[2,1], 'x', 'y');

    d, l := abs(Determinant(A)), x * A[1,2] + y * A[2,2];

    # Performs integer division as in Cox, Little, Schenk (10.1.1)
    s, k := (l + ((-l) mod d)) / d, (-l) mod d;

    sg := sign(Determinant(A));
    M := <-sg * A[2,1], sg * A[1,1]; x + sg * s * A[2,1], y - sg * s * A[1,1]>;
    
    return [d, k, M];

end proc:

# Compute the Hirzebruch-Jung continued fraction expansion of a rational number `x / y`.
hirzebruchJungContinuedFraction := proc(x, y)

    local a, b, res;

    a, b := x, y;
    res := [];

    while b > 0 do
        res := [op(res), (a + ((-a) mod b)) / b];
        a, b := b, (-a) mod b;
    end do;

    return res;

end proc:

# Computes a hilbert basis of a two-dimensional cone whose generators are given by
# the columns of the given matrix `A`. This procedure only outputs the vectors that
# have to be added in order to make a hilbert basis, excluding the columns of A themselves.
# The implementation follows essentially Cox, Little Schenck, Theorem 10.2.3.
hilbertBasisCone2D := proc(A)
    local d, k, M, bs, x, y, a, b, res, i;

    d, k, M := op(cone2DNormalForm(A));
    bs := hirzebruchJungContinuedFraction(d, k);

    x, y := 0, 1;
    a, b := -1, 0;
    res := [];
    for i from 1 to nops(bs) do
        res := [op(res), M^(-1) . <y, -b>];
        x, y := y, bs[i] * y - x;
        a, b := b, bs[i] * b - a;
    end do;

    return res;

end proc:

# Given an integer Matrix `P`, compute the factor group of its image
# The output is given as a list [r, d_1, d_2, ..., d_r], where r encodes
# the rank of the factor group and d_1,...,d_r are the elementary divisors
# encoding the torsion part
imageFactorGroup := proc(P :: Matrix)
    local S, result, i;
    S := SmithForm(P);
    # The rank of the resulting group is the number of zero rows in the smith normal form
    result := [RowDimension(P) - Rank(S)];
    # For the torsion part, we traverse the diagonal of `S` and collect the elementary divisors
    # different from 1.
    for i from 1 to Rank(S) do
        if S[i,i] <> 1 then
            result := [op(result), S[i,i]];
        end if;
    end do;
    return result;
end proc;

# Given an integer Matrix `P`, compute the index of its image in the target domain.
# By convention, this procedure returns `infinity` if `P` is not of full rank.
indexOfImage := proc(P :: Matrix)
    local factorGroup;
    factorGroup := imageFactorGroup(P);
    if factorGroup[1] > 0 then
        return infinity;
    else
        return mul(factorGroup[2 .. nops(factorGroup)]);
    end if;
end proc;

# Computes the column hermite form of a matrix.
# An additional argument specifying the output is also supported
# (same as in LinearAlgebra:-HermiteForm).
dualHermiteForm := proc(A :: Matrix)
    local res;
    res := HermiteForm(Transpose(A), _passed[2.._npassed]);
    if type(res, list) then
        return map(Transpose, res);
    else
        return Transpose(res);
    end if;
end proc;

# Computes the integer kernel of a matrix.
integerKernel := proc(A :: Matrix)
    DeleteColumn(dualHermiteForm(A, output = 'U'), [1 .. Rank(A)]);
end proc;

# Given integer matrices defining sublattices of ZZ^n, compute a basis
# of the intersection of the two sublattices.
# For performance reasons, the output is automatically put into hermite form.
integerIntersectionBasis := proc(A1 :: Matrix, A2 :: Matrix)
    local K, B2;
    K := integerKernel(<A1 | -A2>);
    B2 := DeleteRow(K, [1..ColumnDimension(A1)]);
    return A2 . B2;
end proc;

integerIntersectionBasisList := proc(As :: list(Matrix))
    foldl(integerIntersectionBasis, op(As));
end proc;


