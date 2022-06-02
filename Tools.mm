(*
Some tools.
*)

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

