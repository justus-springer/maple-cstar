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
    K := integerKernel(<A1 | -A2>);
    B2 := DeleteRow(K, [1..ColumnDimension(A1)]);
    return A2 . B2;
end proc;

integerIntersectionBasisList := proc(As :: list(Matrix))
    foldl(integerIntersectionBasis, op(As));
end proc;


