module PFormat()
    option object;

    export r, ns, n, m, s, dim, classGroupRank;

    export ModuleApply :: static := proc()
        Object(PFormat, _passed);
    end;

    export ModuleCopy :: static := proc(self :: PFormat, proto :: PFormat,
        ns :: list(integer), m :: integer, s :: integer, $)
        local i;
        if nops(ns) < 2 then
            error "r = nops(ns) must be at least 2."
        end if;
        for i in ns do
            if i < 1 then error "each element of ns must be at least 1." end if:
        end do:
        if m < 0 then error "m must be at least 0." end if:
        if s < 0 then error "s must be at least 0." end if:
        self:-r := nops(ns);
        self:-ns := ns;
        self:-n := add(ns);
        self:-m := m;
        self:-s := s;
        self:-dim := s + 1;
        self:-classGroupRank := self:-n + self:-m - (self:-r - 1 + self:-s);
    end;

    (*
    EXPLANATION ABOUT INDEXING.

    There are two different indexing techniques when it comes to labeling columns of a P-Matrix.
    The first one labels them according to the decomposition into leaf cones. Here, a
    column from the L-Block is labeled by a pair (i,j), where 1 <= i <= r and 1 <= j <= ns[i].
    A column from the d-block is labeled by a single number between 1 and m. However, it is
    sometimes useful to instead index all columns by a single number between 1 and n + m instead
    of the two-dimensional labeling. The following two functions translate one labeling to another.
    In the tuple labeling, we use the convention that a tuple (-1, k) refers to the k-th column
    in the d-block, where k is between 1 and m.
    *)

    (*
    Translate a single index 1 <= k <= n into a two-dimensional index (i,j), where 1 <= i <= r
    and 1 <= j <= ns[i]. An index n+1 <= k <= n+m is mapped to (-1, k - n).
    See above for more explanation.
    *)
    export singleToDoubleIndex :: static := proc(self :: PFormat, k_ :: integer)
        local i, k;
        k := k_;
        if k < 1 or k > self:-n + self:-m then
            error "index out of range: k must be between 1 and n + m = %1. Given: %2", self:-n + self:-m, k;
        end if;
        if k > self:-n then
            return -1, k - self:-n;
        end if;
        i := 0;
        while k > 0 do
            i++;
            k -= self:-ns[i];
        end do; 
        return i, k + self:-ns[i];
    end proc:

    (*
    Translate a two-dimensional index (i,j), where 1 <= i <= r and 1 <= j <= ns[i] into a single
    index 1 <= k <= n. A pair (-1, j) is mapped to the single index n + j.
    See above for more explanation.
    *)
    export doubleToSingleIndex :: static := proc(self :: PFormat, i :: integer, j :: integer)
        if i = -1 then
            if j < 1 or j > self:-m then
                error "index out of range: i = -1, hence j must be between 1 and m = %1. Given: %2", self:-m, j;
            end if;
            return self:-n + j;
        end if;
        if i < 1 or i > self:-r then
            error "index out of range: i must either be -1 or between 1 and r = %1. Given: %2.", self:-r, i;
        end if;
        if j < 1 or j > self:-ns[i] then
            error "index out of range: j must be between 1 and ns[%1] = %2. Given: %3.", i, self:-ns[i], j;
        end if;
        return add(self:-ns[1 .. i-1]) + j
    end proc:

    (*
    Checks whether a given `cone` is big with respect to a given PFormat.
    Here, a cone is called big if for every i = 0,...,r-1, we have
    {N, ..., N+ns[i+1]} ∩ cone ≠ {}, where N = ns[1] + ... + ns[i].
    *)
    export isBigCone :: static := proc(self :: PFormat, cone :: set(integer))
        local N, i, k:
        for i from 0 to self:-r - 1 do
            N := add(self:-ns[1..i]);
            if evalb({seq(N + k, k = 1 .. self:-ns[i+1])} intersect cone = { }) then
                return false:
            end if:
        end do:
        return true:
    end proc:

    (*
    Checks whether a given `cone` is big with respect to a given PFormat.
    Here, a cone is called a leaf cone, if there exists an i = 0,...,r-1 such that:
    xcone ⊆ {N+1, ..., N+ns[i+1]} ∪ {n+1,...,n+m},
    where N = ns[1] + ... + ns[i]
    *)
    export isLeafCone :: static := proc(self :: PFormat, cone :: set(integer))
        local N, i, k:
        for i from 0 to self:-r - 1 do
            N := add(self:-ns[1..i]);
            if cone subset ({seq(N + k, k = 1 .. self:-ns[i+1])} union
                {seq(self:-n + k, k = 1 .. self:-m)}) then
                return true:
            end if:
        end do:
        return false:
    end proc:

    (*
    Checks whether a given `cone` is an X-cone with respect to a given P-matrix format.
    Here, a cone is called an X-cone, if it is either a leaf cone or a big cone.
    *)
    export isXCone :: static := proc(self :: PFormat, cone :: set(integer))
        isBigCone(self, cone) or isLeafCone(self, cone):
    end:

    (*
    Given a list of cones, compute the list X-cones which are maximal with respect to a given PFormat.
    *)
    export getMaximalXConesFormat :: static := proc(self :: PFormat, cones :: set(set(integer)))
        local bigCones, nonBigCones, leafCones, leafCone, cone, maxLeafCones, c1, c2, isMaximal, N, i, k;
        # First, we compute the big cones. Note that these are neccessarily maximal and there can
        # be no other maximal big cones.
        bigCones := select(c -> isBigCone(self, c), cones);
        # For the leaf cones, the issue is that there may be maximal leaf cones hiding "inside" the cones.
        # For each of the remaining cones, we compute all the maximal leaf cones it contains.
        leafCones := {};
        nonBigCones := cones minus bigCones;
        for i from 0 to self:-r - 1 do
            N := add(self:-ns[1..i]);
            leafCone := {seq(N + k, k = 1 .. self:-ns[i+1])} union {seq(self:-n + k, k = 1 .. self:-m)};
            leafCones := leafCones union map(c -> c intersect leafCone, nonBigCones);
        end do;
        # Remove the ones that are non-maximal
        maxLeafCones := {};
        for c1 in leafCones do
            isMaximal := true;
            for c2 in bigCones union (leafCones minus {c1}) do
                if c1 subset c2 then
                    isMaximal := false;
                end if;
            end do;
            if isMaximal then
                maxLeafCones := {op(maxLeafCones), c1};
            end if;
        end do;

        return bigCones union maxLeafCones;
    end;

    export `=`::static := proc( l, r, $ )
        if (_npassed <> 2 or not l::PFormat or not r::PFormat) then
           return false;
        end;
        return l:-ns = r:-ns and l:-m = r:-m and l:-s = r:-s;
    end;

    export ModulePrint :: static := proc(self :: PFormat)
        nprintf(cat("PFormat(", self:-ns, ", m = ", self:-m, ", s = ", self:-s, ")"));
    end;

end module: