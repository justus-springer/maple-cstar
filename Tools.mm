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

sortLexComparison := proc(a :: {numeric, list}, b :: {numeric, list})
    local i;
    if type(a, numeric) and type(b, numeric) then
        return a > b;
    elif type(a, list) and type(b, list) then
        if nops(a) > nops(b) then
            return true;
        elif nops(a) < nops(b) then
            return false;
        end if;
        for i from 1 to nops(a) do
            if sortLexComparison(a[i], b[i]) then
                return true
            elif sortLexComparison(b[i], a[i]) then 
                return false;
            end if;
        end do;
        return false;
    end if;
    error "Both arguments must either be numbers or lists";
end proc;

sortLex := proc(ls :: list)
    if type(ls, list(numeric)) then
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
        return lcm(op(factorGroup[2 .. nops(factorGroup)]));
    end if;
end proc;

