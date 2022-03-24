ComplexityOne := module()

option package;

export PFormat, PMatrix, TVarOne, FindInDatabase, ImportTVarOneList, ExportTVarOneList;

## TODO: Remove dependency on MDSpackage.
uses LinearAlgebra, MDSpackage, Database[SQLite];

(****************
** CONVENIENCE **
****************)

export applyPermToList := proc(p :: Perm, ls :: list)
    local i;
    return [seq(ls[(p^(-1))[i]], i = 1 .. nops(ls))];
end proc;

(* Computes all permutations leaving a given list invariant. *)
export invariantPermutations := proc(ls :: list, compFun := `=`)
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

(****************
**** PFORMAT ****
****************)

module PFormat()
    option object;

    export r, ns, n, m, s;

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
        if s < 1 then error "s must be at least 1." end if:
        self:-r := nops(ns);
        self:-ns := ns;
        self:-n := add(ns);
        self:-m := m;
        self:-s := s;
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

    (* Helper for `bundleColumnPermutation` *)
    export bundleColumnPermutationApply := proc(self :: PFormat, sigma :: Perm, taus :: list(Perm), rho :: Perm, k :: integer)
        local i, j, newFormat;
        i, j := singleToDoubleIndex(self, k);
        if i = -1 then
            return doubleToSingleIndex(self, -1, rho[j]);
        end if;
        newFormat := PFormat(applyPermToList(sigma, self:-ns), self:-m, self:-s);
        return doubleToSingleIndex(newFormat, sigma[i], taus[sigma[i]][j]);
    end proc;

    (* 
    Bundles the data of an admissible column operation into a single permutation of indices {1 .. n+m}. Here:
    - `sigma` is the permutation of blocks, i.e. a permutation of the numbers {1 .. r},
    - `taus` is a list of permutations, where `taus[i]` is a permutation of the numbers {1 .. ns[i]} and
    - `rho` is the permutation of the last `m` columns, i.e. a permutation of the numbers {1 .. m}
    *)
    export bundleColumnPermutation := proc(self :: PFormat, sigma :: Perm, taus :: list(Perm), rho :: Perm)
        return Perm(map(i -> bundleColumnPermutationApply(self, sigma, taus, rho, i), [seq(1 .. self:-n + self:-m)]));
    end proc;

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

module PMatrix()
    option object;

    # These fields are guaranteed to be filled when a PMatrix is created.
    export format, lss, d, mat, P0;
    export r, ns, n, m, s, dim;

    # These fields are only computed when needed. Use these getters below for these.
    local picardNumber := undefined;
    local Q := undefined;
    local Q0 := undefined;
    local classGroup := undefined;
    local anticanClass := undefined;
    local anitcanCoefficients := undefined;
    local movingCone := undefined;

    export setFormat :: static := proc(self :: PMatrix, f :: PFormat)
        self:-format := f;
        self:-r := f:-r;
        self:-ns := f:-ns;
        self:-n := f:-n;
        self:-m := f:-m;
        self:-s := f:-s;
        self:-dim := f:-s + 1;
    end proc;

    # Check if all columns of P are primitive
    # Throws an error if they are not.
    local assertColumnsPrimitive :: static := proc(self :: PMatrix)
        local i;
        for i from 1 to self:-n + self:-m do
            if igcd(seq(Column(self:-mat, i))) <> 1 then
                error "This is not a P-matrix: The %-1 column is not primitve.", i;
            end if;
        end do;
    end proc;

    # Check if the columns generate QQ^(r+s) as a cone.
    # Throws an error if they do not.
    local assertColumnsGenerateFullCone :: static := proc(self :: PMatrix)
        if poshull(Column(self:-mat, [seq(1..self:-n + self:-m)])) &<> fullcone(self:-r + self:-s - 1) then
            error "This is not a P-matrix. The columns do not generate QQ^(r+s) as a cone.";
        end if;
    end proc;

    export ModuleApply :: static := proc()
        Object(PMatrix, _passed);
    end;

    (*
    This procedure creates a PMatrix from various kinds of data. It supports
    five different input methods:
    (1) format :: PFormat, lss :: list(list(integer)), d :: Matrix.
    (2) lss :: list(list(integer)), d :: Matrix.
    (3) format :: PFormat, P :: Matrix.
    (4) s :: integer, P :: Matrix

    In input method (1) and (2), `lss` is the list of exponent vectors of relations
    of the Cox Ring. These make up the first `r-1` rows of the PMatrix. The lower
    `s` rows have to be provided the matrix `d`. In method (2), the format of the matrix
    is inferred from the given data.

    In methods (3) and (4), the P-Matrix is directly specified. It will be checked
    if this matrix fulfills the properties of a P-Matrix, i.e. if it has the correct shape,
    the columns are primitive and they generate the whole space as a cone. In method (3), the
    expected P-Format is specified and will be checked against P. In method (4), only the dimension
    of the acting torus `s` is specified (the dimension of the overall variety will be `s+1`).
    This is necessary because `s` cannot be uniquely inferred from `P`.

    *)
    export ModuleCopy :: static := proc(self :: PMatrix, proto :: PMatrix)

        local lss, ls, l, rows, rows0, i, j, P, r, s, ns, numZerosBefore;
        if _npassed = 2 then error "not enough arguments." end if;

        if type(_passed[3], 'PFormat') then
            # Input method (1) or (3)

            if _npassed < 4 then
                error "Not enough arguments. Expected input method (1) or (3).";
            end if:

            if type(_passed[4], Matrix) then
                # Input method (3)

                # In this case, we call the function again with input method (4) and check if the
                # inferred format coincides with the given one.
                P := PMatrix(_passed[3]:-s, _passed[4]);
                if P:-ns <> _passed[3]:-ns then
                    error "Expected P-Format does not coincide with the inferred one. Expected: ns = %1. Inferred: ns = %2", _passed[3]:-ns, P:-ns;
                elif P:-m <> _passed[3]:-m then
                    error "Expected P-Format does not coincide with the inferred one. Expected: m = %1. Inferred: m = %2", _passed[3]:-m, P:-m;
                elif P:-s <> _passed[3]:-s then
                    error "Expected P-Format does not coincide with the inferred one. Expected: m = %1. Inferred: m = %2", _passed[3]:-s, P:-s;
                end if;

                setFormat(self, _passed[3]);
                self:-lss := P:-lss;
                self:-mat := P:-mat;
                self:-P0 := P:-P0;
                self:-d := P:-d;

            else
                # Input method (1)
                # PFormat, list(list(integer)), Matrix.
                setFormat(self, _passed[3]);

                if not type(_passed[4], list(list(integer))) then
                    error "Expected 2nd argument to be of type: list(list(integer))";
                end if;

                self:-lss := _passed[4];

                for ls in self:-lss do
                    for l in ls do
                        if l < 1 then
                            error "All entries of `lss` must be greater or equal to 1."
                        end if;
                    end do;
                end do;

                if not type(_passed[5], 'Matrix'(self:-s, self:-n + self:-m, integer)) then
                    error "Expected 3rd argument to be of type: Matrix(%1, %2, integer)", self:-s, self:-n + self:-m;
                end if;
                self:-d := _passed[5];

                # Check if the given ls match the given P-format
                for i from 1 to self:-r do
                    if nops(self:-lss[i]) <> self:-ns[i] then
                        error "length of %-1 vector in lss does not match given P-format. Expected length: %2. Given length: %3.", i, self:-ns[i], nops(self:-lss[i]);
                    end if;
                end do:

                # Construct the P-matrix from the given data
                rows := [seq([seq(-self:-lss[1]), (0 $ add(self:-ns[2..i-1]),
                            seq(self:-lss[i])), (0 $ add(self:-ns[i+1..self:-r]) + self:-m)], i = 2 .. self:-r),
                        seq(:-convert(Row(self:-d, j), list), j = 1 .. self:-s)];
                self:-mat := Matrix(rows);

                # Only the L-block. This gives P0.
                rows0 := [seq([seq(-self:-lss[1]), (0 $ add(self:-ns[2..i-1]),
                            seq(self:-lss[i])), (0 $ add(self:-ns[i..self:-r-1]))], i = 2 .. self:-r)];
                self:-P0 := Matrix(rows0);

                assertColumnsPrimitive(self);
                assertColumnsGenerateFullCone(self);
            end if;

        elif type(_passed[3], list(list(integer))) then
            # Input method (2)
            # lss :: list(list(integer)), d :: Matrix.

            if _npassed < 4 then
                error "Not enough arguments. Expected input: list(list(integer)), Matrix.";
            end if;

            if not type(_passed[4], 'Matrix') then
                error "Expected 2nd argument to be of type: Matrix";
            end if;

            self:-lss := _passed[3];
            self:-ns := map(nops, self:-lss);
            self:-d := _passed[4];
            self:-format := PFormat(self:-ns, ColumnDimension(self:-d) - add(self:-ns), RowDimension(self:-d));
            return PMatrix[ModuleCopy](self, proto, self:-format, self:-lss, self:-d);

        elif type(_passed[3], integer) then
            # Input method (4)
            # s :: integer, P :: Matrix

            if _npassed < 4 then
                error "Not enough arguments. Expected input: integer, Matrix";
            end if;

            s := _passed[3];

            if not type(_passed[4], Matrix) then
                error "Expected 2nd argument to be of type: Matrix";
            end if;

            P := _passed[4];
            self:-mat := P;

            if s > RowDimension(P) - 1 then
                error "The dimension of the acting torus `s` cannot be greater than the number of rows of `P` minus one.";
            end if;

            self:-P0 := SubMatrix(P, [seq(1 .. RowDimension(P) - s)], [seq(1 .. ColumnDimension(P))]);

            r := RowDimension(P) - s + 1;

            for i in P do
                if not type(i, integer) then
                    error "All entries of P must be of type: integer";
                end if;
            end do;

            if RowDimension(P) < r then
                error "P must have at least r = %1 rows", r;
            end if;

            # Get the first vector l1 by looking at the first row of `P`.
            # We move right until the entries are no longer negative.
            ls := [];
            if not (P[1,1] < 0) then
                error "Given matrix is not in P-shape. Expected: P[1,1] < 0. Given: P[1,1] = %1", P[1,1];
            end if;
            i := 1;
            while P[1,i] < 0 do
                ls := [op(ls), -P[1,i]];
                i := i + 1;

                # Needs to be checked before we re-enter the loop
                if i > ColumnDimension(P) then
                    error "Given matrix is not in P-shape. Expected: P[1,%1] > 0. Given: P[1,%1] = %2", i-1, P[1,i-1];
                end if;
            end do;

            lss := [ls];
            ns := [nops(ls), 0 $ r-1]; # Here, we preliminarily fill up with zeros.

            # Check that the first r-1 rows of P all start with -l1
            for i from 1 to r-1 do
                for j from 1 to nops(ls) do
                    if P[i,j] <> P[1,j] then
                        error "Given matrix is not in P-shape. Expected: P[%1,%2] = P[1,%2]. Given: P[%1,%2] = %3 and P[1,%2] = %4", i, j, P[i,j], P[1,j];
                    end if;
                end do;
            end do;

            for i from 1 to r-1 do
                # Note that since we fill up the remaining entries of `ns` with zeros, the following
                # gives the correct result.
                numZerosBefore := add(ns[2..]);
                for j from ns[1] + 1 to ns[1] + numZerosBefore do
                    if P[i,j] <> 0 then
                        error "Given matrix is not in P-shape. Expected: P[%1,%2] = 0. Given: P[%1,%2] = %3", i, j, P[i,j];
                    end if;
                end do;

                if P[i,j] < 1 then
                    error "Given matrix is not in P-shape. Expected: P[%1,%2] > 0. Given: P[%1,%2] = %3", i, j, P[i,j];
                end if;

                ls := [];
                while P[i,j] <> 0 do

                    if P[i,j] < 1 then
                        error "Given matrix is not in P-shape. Expected: P[%1,%2] > 0. Given: P[%1,%2] = %3", i, j, P[i,j];
                    end if;

                    ls := [op(ls), P[i,j]];

                    if i = r-1 and j = ColumnDimension(P) then
                      break; # In this case, we are done.
                    end if;

                    j := j + 1;

                    # Needs to be checked before we re-enter the loop
                    if j > ColumnDimension(P) then
                        error "Given matrix is not in P-shape. Expected: P[%1,%2] = 0. Given: P[%1,%2] = %3", i, j-1, P[i,j-1];
                    end if;

                end do;

                ns[i+1] := nops(ls);
                lss := [op(lss), ls];

            end do;

            self:-lss := lss;

            setFormat(self, PFormat(ns, ColumnDimension(P) - add(ns), RowDimension(P) - r + 1));

            # Check if we really have all-zeros in the upper right block
            for i from 1 to r-1 do
                for j from self:-n + 1 to self:-n + self:-m do
                    if P[i,j] <> 0 then
                        error "Given matrix is not in P-shape. Expected P[%1,%2] = 0. Given: P[%1,%2] = %3", i, j, P[i,j];
                    end if;
                end do;
            end do;

            assertColumnsPrimitive(self);
            assertColumnsGenerateFullCone(self);

            # Finally, we get the d-block
            self:-d := SubMatrix(P, [self:-r .. RowDimension(P)], [1 .. ColumnDimension(P)]);

        end if;

    end;

    export setQ :: static := proc(self :: PMatrix, Q :: Matrix) self:-Q := Q; end proc;

    export setQ0 :: static := proc(self :: PMatrix, Q0 :: Matrix) self:-Q0 := Q0; end proc;

    export setClassGroup :: static := proc(self :: PMatrix, classGroup :: AG) self:-classGroup := classGroup; end proc;

    export setPicardNumber :: static := proc(self :: PMatrix, picardNumber :: integer) self:-picardNumber := picardNumber; end proc;

    export setAnticanCoefficients :: static := proc(self :: PMatrix, anitcanCoefficients) self:-anitcanCoefficients := anitcanCoefficients; end proc;

    export setAnticanClass :: static := proc(self :: PMatrix, anticanClass) self:-anticanClass := anticanClass; end proc;

    export setMovingCone :: static := proc(self :: PMatrix, movingCone :: CONE) self:-movingCone := movingCone; end proc;

    (*
    Gets the degree Matrix (Q) of a P-Matrix. If the degree Matrix has already been computed, it
    is saved and will be used as a cache for any later call of `getQ`. By supplying the option 'forceCompute',
    you can prevent using the cache and recompute the Q-matrix, even if it has already been computed in the past.
    *)
    export getQ :: static := proc(self :: PMatrix)
        local A;
        # Currently, we rely on MDSpackage for these computations.
        # But they should be reimplemented here eventually.
        if type(self:-Q, undefined) or 'forceCompute' in [_passed] then
            A := AGHP2Q(convert(self, matrix));
            setClassGroup(self, AGHdata(A)[2]);
            setQ(self, Matrix(AGHdata(A)[3]));
            setPicardNumber(self, AGdata(self:-classGroup)[3]);
            setQ0(self, DeleteRow(self:-Q, self:-picardNumber + 1 .. RowDimension(self:-Q)));
        end if;
        return self:-Q;
    end;

    export getQ0 :: static := proc(self :: PMatrix)
        if type(self:-Q0, undefined) or 'forceCompute' in [_passed] then
            if not type(self:-Q, undefined) then
                setQ0(self, DeleteRow(self:-Q, self:-picardNumber + 1 .. RowDimension(self:-Q)));
            else
                getQ(self);
            end if;
        end if;
        return self:-Q0;
    end;

    export getClassGroup :: static := proc(self :: PMatrix)
        if type(self:-classGroup, undefined) or 'forceCompute' in [_passed] then
            # Here we do a forceCompute to prevent a situation where Q is defined while classGroup is not.
            getQ(self, 'forceCompute');
        end if;
        return self:-classGroup;
    end;

    export getPicardNumber :: static := proc(self :: PMatrix)
        if type(self:-picardNumber, undefined) then
            setPicardNumber(self, ColumnDimension(self:-mat) - RowDimension(self:-mat));
        end if;
        return self:-picardNumber;
    end;

    export getAnticanCoefficients :: static := proc(self :: PMatrix)
        if type(self:-anitcanCoefficients, undefined) or 'forceCompute' in [_passed] then
            setAnticanCoefficients(self,
                [1 $ self:-n + self:-m] - (self:-r - 2) * [op(self:-lss[1]), 0 $ self:-n + self:-m - self:-ns[1]]);
        end if;
        return self:-anitcanCoefficients;
    end;

    export getAnticanClass :: static := proc(self :: PMatrix)
        local as, i, anticanVec, d;
        if type(self:-anticanClass, undefined) or 'forceCompute' in [_passed] then
            as := getAnticanCoefficients(self);
            anticanVec := add(seq(as[i] * Column(getQ(self), i), i = 1 .. self:-n + self:-m));
            # Some entries in `anticanVec` live in cyclic groups Z/dZ.
            # We normalize these entries, so that each of them is less than `d`.
            for i from getPicardNumber(self) + 1 to RowDimension(getQ(self)) do
                d := AGdata(getClassGroup(self))[4][i - getPicardNumber(self)];
                anticanVec[i] := anticanVec[i] mod d;
            end do;
            setAnticanClass(self, anticanVec);
        end if;
        return self:-anticanClass;
    end;

    (*
    Compute the linear form solving the Gorenstein condition on a given X-cone,
    if there are any.
    *)
    export getLinearFormForXCone :: static := proc(self :: PMatrix, cone :: set(integer))
        local as, u, us, i, j, sol;
        as := getAnticanCoefficients(self);
        us := [seq(u[i], i = 1 .. RowDimension(self:-mat))];
        return solve({seq(DotProduct(us, Column(self:-mat, j)) = as[j], j in cone)});
    end;

    (*
    Computes the gorenstein index of a given X-cone. 
    That is, the smallest positive integer n such that n*K_X is Cartier on the toric orbit 
    defined by the X-cone, where K_X is the anticanonical divisor.
    *)
    export gorensteinIndexForXCone :: static := proc(self :: PMatrix, cone :: set(integer))
        local e;
        return ilcm(seq(denom(rhs(e)), e in getLinearFormForXCone(self, cone)));
    end;

    (*
    Computes the moving cone associated to a P-Matrix.
    This is the intersection of all images of facets of the positive orthant under Q.
    *)
    export getMovingCone :: static := proc(self :: PMatrix)
        local Q0, cones, i, j;
        if type(self:-movingCone, undefined) or 'forceCompute' in [_passed] then
            Q0 := getQ0(self);
            cones := seq(poshull(seq(Column(Q0, i), i in {seq(1 .. ColumnDimension(Q0))} minus {j})), j = 1 .. ColumnDimension(Q0));
            self:-setMovingCone(self, intersection(cones));
        end if;
        return self:-movingCone;
    end;

    (*
    Checks whether there exists a Fano variety coming from this P-Matrix.
    This is equivalent to checking whether the anticanonical class lies in the
    relative interior of the moving cone.
    *)
    export admitsFanoVariety :: static := proc(self :: PMatrix)
        return containsrelint(getMovingCone(self), getAnticanClass(self));
    end;

    (*
    Computes the slope of column in a P-matrix of a C*-surface.
    *)
    export getSlope :: static := proc(self :: PMatrix, i :: integer, j :: integer)
        if self:-s <> 1 then
            error "getSlope is only defined for P-Matrices corresponding to C*-surfaces, i.e. s has to be 1.";
        end if;
        return self:-d[1, doubleToSingleIndex(self:-format, i, j)] / self:-lss[i, j];
    end proc;

    (* ************************
    ** ADMISSIBLE OPERATIONS **
    ***************************)

    (*
    Creates a new P-Matrix by applying an admissible column operation. The admissible operation is given in three parts:
    - `sigma` is a permutation of the blocks, i.e. a permutation of the numbers [1 .. r]
    - `taus` is a list of permutation for the individual blocks. Note that the permutation of
      the blocks is applied *before* the taus, i.e. `taus[i]` is a permutation of [1 .. ns[sigma[i]]].
    - `rho` is a permutation of the last `m` columns, i.e. a permutation of the numbers [1 .. m].
    *)
    export applyAdmissibleColumnOperation := proc(P :: PMatrix, sigma :: Perm, taus :: list(Perm), rho :: Perm)
        local bundledPerm, permutationMatrix, newD, newLss, i, j;
        # First, bundle the column operation into a single permutation of [1 .. n + m]
        # This is done at the level of the format.
        bundledPerm := bundleColumnPermutation(P:-format, sigma, taus, rho);
        # Create the permutation matrix corresponding to `bundledPerm`
        permutationMatrix := Matrix(P:-n + P:-m, P:-n + P:-m, (i,j) -> if j = bundledPerm[i] then 1 else 0 end if);
        # We permute the columns of the `d`-block by multiplying from the right with the permutation matrix.
        newD := P:-d . permutationMatrix;
        # Finally, we create the permuted list of the lss, by first applying `sigma` to the blocks
        # and then applying `tau[i]` to `lss[i]`. This is done simultaneoulsy using `zip`.
        newLss := zip(applyPermToList, taus, applyPermToList(sigma, P:-lss));
        return PMatrix(newLss, newD);
    end proc;

    (*
    Creates a new P-Matrix by applying an admissible column operation of type (1), i.e. 
    permutation of blocks. The parameter `sigma` must be a permutation of the numbers [1 .. r]
    *)
    export applyAdmissibleColumnOperation1 := proc(P :: PMatrix, sigma :: Perm)
        applyAdmissibleColumnOperation(P, sigma, [Perm([]) $ P:-r], Perm([]));
    end proc;

    (*
    Creates a new P-Matrix by applying an admissible column operation of type (2), i.e.
    permutation of column within a given block `i`. The parameter `tau` must be a permutation
    of the numbers [1 .. ns[i]].
    *)
    export applyAdmissibleColumnOperation2 := proc(P :: PMatrix, i :: integer, tau :: Perm)
        applyAdmissibleColumnOperation(P, Perm([]), [Perm([]) $ (i-1), tau, Perm([]) $ P:-r - i], Perm([]));
    end proc;

    (*
    Creates a new P-Matrix by applying an admissible column operation of type (3), i.e.
    permutation of the last m rows. The parameter `rho` must be a permutation of the
    numbers [1 .. m].
    *)
    export applyAdmissibleColumnOperation3 := proc(P :: PMatrix, rho :: Perm)
        applyAdmissibleColumnOperation(P, Perm([]), [Perm([]) $ P:-r], rho);
    end proc;

    (*
    Creates a new P-Matrix by applying an admissible row operation. The row operation is specified by
    an (s x (r-1))-Matrix `A` and an (s x s)-Matrix `B`. The row operation is applied by multiplying `P`
    from the left with an elementary matrix as follows:

    ( E 0 )  *  ( L 0 )  =  (    L      0  )
    ( A B )     ( d d')     ( AL + Bd   Bd')

    *)
    export applyAdmissibleRowOperation := proc(P :: PMatrix, A :: Matrix, B :: Matrix)
        local identityMatrix, zeroMatrix, newP;
        if not type(A, 'Matrix'(P:-s, P:-r - 1, integer)) then
            error "Expected second argument to be of type: Matrix(%1, %2, integer)", P:-s, P:-r - 1;
        end if;
        if not type(B, 'Matrix'(P:-s, P:-s, integer)) then
            error "Expected third argument to be of type: Matrix(%1, %1, integer)", P:-s;
        end if;
        if not abs(Determinant(B)) = 1 then
            error "Third argument must be a unimodular matrix, i.e. its determinant must be 1 or -1";
        end if;
        # An (r-1) x (r-1) identity matrix.
        identityMatrix := Matrix(P:-r - 1, P:-r - 1, shape = diagonal, 1);
        # An (r-1) x s zero matrix.
        zeroMatrix := Matrix(P:-r - 1, P:-s, fill = 0);
        # Multiply P from the left with
        # ( E  0 )
        # ( A  B )
        newP := <<identityMatrix | zeroMatrix>, <A | B>> . P:-mat;
        return PMatrix(P:-s, newP);
    end proc;

    (******************
    ** NORMALIZATION **
    *******************)

    (*
    Removes all redundant blocks from a P-Matrix, i.e. all blocks with a single 1 inside the L-block.
    *)
    export removeRedundantColumns := proc(P :: PMatrix)
        local newLss, redundantIndices, newD;
        newLss := remove(ls -> ls = [1], P:-lss);
        redundantIndices := map(i -> doubleToSingleIndex(P:-format, i, 1), select(i -> P:-lss[i] = [1], [seq(1 .. nops(P:-lss))]));
        newD := DeleteColumn(P:-d, redundantIndices);
        return PMatrix(newLss, newD);
    end proc;

    (*
    Normalizes a P-Matrix by removing redundant columns and applying admissible column permutations.
    The columns of each block in the resulting P-Matrix are ordered descendingly by the ls,
    while the blocks themselves are ordered descendingly by size, and blocks of the same size
    are ordered lexicographically.

    You can obtain the permutation `sigma` of the blocks and the permutations `taus` of the
    individual blocks by specifying the option `'output' = out`, where `out` is a list of the
    names 'normalized', 'sigma' and 'taus'.

    *)
    export normalizePMatrix := proc(P0 :: PMatrix)
        local P1, P2, P3, sortKey, taus_, i, compfun, tau, sigma_, result, str;
        # First, remove any redundant columns
        P1 := removeRedundantColumns(P0);
        # The comparison function to sort the columns by.
        # For surfaces, we sort by slope. In general, we sort by the entries of the lss directly.
        sortKey := proc(P :: PMatrix, i :: integer, j :: integer)
            if P:-s = 1 then
                Column(P:-d, doubleToSingleIndex(P:-format, i, j))[1] / P:-lss[i][j];
            else
                P:-lss[i][j];
            end if;
        end proc;

        # First we sort each individual block.
        taus_ := [];
        for i from 1 to P1:-r do
            tau := Perm(sort([seq(1 .. P1:-ns[i])], (j1, j2) -> sortKey(P1, i, j1) > sortKey(P1, i, j2), 'output' = 'permutation'))^(-1);
            taus_ := [op(taus_), tau];
        end do;
        P2 := applyAdmissibleColumnOperation(P1, Perm([]), taus_, Perm([]));

        # This is the comparison function we will use to sort the blocks themselves.
        # It returns true if the block of index `i1` should occur before the block of index `i2`.
        compfun := proc(P :: PMatrix, i1 :: integer, i2 :: integer)
            if nops(P:-lss[i1]) < nops(P:-lss[i2]) then return false;
            elif nops(P:-lss[i1]) > nops(P:-lss[i2]) then return true;
            else
                for j from 1 to nops(P:-lss[i1]) do
                    if sortKey(P, i1, j) < sortKey(P, i2, j) then return false;
                    elif sortKey(P, i1, j) > sortKey(P, i2, j) then return true;
                    end if;
                end do;
                return true;
            end if;
        end proc;
        sigma_ := Perm(sort([seq(1 .. P2:-r)], (i1, i2) -> compfun(P2, i1, i2), 'output' = 'permutation'))^(-1);
        P3 := applyAdmissibleColumnOperation1(P2, sigma_);

        # Output handling
        if _npassed = 1 then
            return P3;
        end if;
        result := [];
        for i from 2 to _npassed do
            if type(_passed[i], `=`) then
                if lhs(_passed[i]) = 'output' then
                    for str in rhs(_passed[i]) do
                        if str = 'normalized' then
                            result := [op(result), P3];
                        elif str = 'sigma' then
                            result := [op(result), sigma_];
                        elif str = 'taus' then
                            result := [op(result), taus_];
                        end if;
                    end do;
                else
                    error "Unrecognized option: %1", lhs(_passed[i]);
                end if;
            end if;
        end do;
        return result;
    end proc;

    (* *******************
    ** EQUIVALENCE TEST **
    **********************)

    (*
    Checks whether two P-Matrices `P1` and `P2` are equivalent to each other by admissible row operations.
    You can obtain the admissible row operation turning `P1` into `P2` by supplying the parameter `'output' = out`,
    where `out` is a list of the names 'result', 'A', 'B' and 'S'.
    *)
    export areRowEquivalent := proc(P1 :: PMatrix, P2 :: PMatrix)
        local resBool, resA, resB, resS, identityMatrix, zeroMatrix, A_, B_, S_, newP, sol, resultList, str, i, j;
        resBool, resA, resB, resS := false, undefined, undefined, undefined;
        # P-Matrices of different format are never row-equivalent.
        if P1:-format = P2:-format then
            identityMatrix := Matrix(P1:-r - 1, P1:-r - 1, shape = diagonal, 1);
            zeroMatrix := Matrix(P1:-r - 1, P1:-s, fill = 0);
            A_ := Matrix(P1:-s, P1:-r - 1, symbol = 'a');
            B_ := Matrix(P1:-s, P1:-s, symbol = 'b');
            S_ := <<identityMatrix | zeroMatrix>, <A_ | B_>>;
            newP := S_ . P1:-mat;
            sol := isolve({seq(seq(newP[i,j] = P2:-mat[i,j] , j = 1 .. ColumnDimension(P1:-mat)), i = P1:-r .. RowDimension(P1:-mat))});
            if sol = NULL then
                resBool, resA, resB, resS := false, undefined, undefined, undefined;
            else
                if abs(Determinant(subs(sol, B_))) = 1 then
                    resBool := true;
                    resA := subs(sol, A_);
                    resB := subs(sol, B_);
                    resS := subs(sol, S_);
                end if;
            end if;
        end if;
        # Output handling
        if _npassed = 2 then
            return resBool;
        end if;
        resultList := [];
        for i from 2 to _npassed do
            if type(_passed[i], `=`) then
                if lhs(_passed[i]) = 'output' then
                    for str in rhs(_passed[i]) do
                        if str = 'result' then
                            resultList := [op(resultList), resBool];
                        elif str = 'A' then
                            resultList := [op(resultList), resA];
                        elif str = 'B' then
                            resultList := [op(resultList), resB];
                        elif str = 'S' then
                            resultList := [op(resultList), resS];
                        end if;
                    end do;
                else
                    error "Unrecognized option: %1", lhs(_passed[i]);
                end if;
            end if;
        end do;
        return resultList;
    end proc;

    (*
    Checks whether two P-Matrices `P1` and `P2` are eqiuvalent, i.e. can be transformed into each other
    by a series of admissible operations.
    TODO: Output the admissible operations turning `P1` into `P2`.
    *)
    export areEquivalent := proc(P1_ :: PMatrix, P2_ :: PMatrix)
        local P1, P2, Ps, P, i;
        if not evalb('skipNormalization' in [_passed]) then
            # First, normalize P1 and P2
            P1 := normalizePMatrix(P1_);
            P2 := normalizePMatrix(P2_);
        else 
            P1 := P1_;
            P2 := P2_;
        end if;
        # After normalization, the lss must be equal, otherwise `P1` and `P2` can't be equivalent.
        # Note that normalization also removes redundant columns.
        if P1:-lss <> P2:-lss then
            return false;
        end if;
        # Starting from `P1`, we consider all P-Matrices that we can reach via admissible column operations
        # leaving the L-block invariant. We do this in three steps corresponding to the three kinds of admissible column operations.
        # First, we consider all column operations permuting two entire blocks.
        Ps := map(sigma -> applyAdmissibleColumnOperation1(P1, sigma), invariantPermutations(P1:-lss));
        # Second, we consider for each of those P-matrices all further P-matrices we can reach by admissible
        # column permutations within a single block, such that the entire L-block is kept invariant.
        for i from 1 to P1:-r do
            Ps := ListTools[Flatten](map(P -> map(tau -> applyAdmissibleColumnOperation2(P, i, tau), invariantPermutations(P:-lss[i])), Ps));
        end do;
        # Thirdly, for each of those P-matrices, we consider all further P-matrices we reach by admissible
        # column permutations within the last m columns.
        Ps := ListTools[Flatten](map(P -> map(rho -> applyAdmissibleColumnOperation3(P, rho), map(Perm, combinat[permute](P:-m))), Ps));
        # Now, for each of those P-matrices, we check if we can transform them into `P2` by admissible row operations.
        for P in Ps do
            if areRowEquivalent(P, P2) then
                return true;
            end if;
        end do;
        return false;
    end proc;

    export ModulePrint :: static := proc(self :: PMatrix)
        nprintf(cat("PMatrix(", self:-lss, ", m = ", self:-m, ", s = ", self:-s, ")"));
    end;

    export PMatrixInfo :: static := proc(self :: PMatrix)
        local P, Q, picardNumber, classGroup, anticanClass, admitsFano;
        print(P = self:-mat);
        nprintf(cat("Format: ", self:-ns, ", m = ", self:-m, ", s = ", self:-s)); #"
        print(Q = getQ(self));
        print(classGroup = getClassGroup(self));
        print(picardNumber = getPicardNumber(self));
        print(anticanClass = getAnticanClass(self));
        print(admitsFano = admitsFanoVariety(self));
    end;

    export convert :: static := proc(self :: PMatrix, toType, $)
        if toType = ':-Matrix' then
            self:-mat;
        elif toType = ':-matrix' then
            :-convert(self:-mat, matrix);
        else
            error "cannot convert from PMatrix to %1", toType;
        end if;
    end;

end module:

module TVarOne()
    option object;

    export P, Sigma;

    local XCones := undefined;
    local maximalXCones := undefined;
    local gorensteinIndex := undefined;
    local isGorensteinVal := undefined;

    export ModuleApply :: static := proc()
        Object(TVarOne, _passed);
    end;

    (*
    This method creates a T-Variety of complexity one from various kinds of data. It supports three
    different input methods:
    (1) P :: PMatrix, Sigma :: set(set(integer)).
    (2) P :: PMatrix, where the picard number of P is one.
    (3) P :: PMatrix, where P is a P-Matrix of a C* surface.
    *)
    export ModuleCopy :: static := proc(self :: TVarOne, proto :: TVarOne, P :: PMatrix)
        local numColumns, i, j, k, ordered_indices, taus, sigma_plus, sigma_minus, taus_plus, taus_minus;

        self:-P := P;

        if _npassed = 4 then
            # Input method (1)
            if not type(_passed[4], set(set(integer))) then
                error "Expected second argument to be of type: set(set(integer)).";
            end if;
            self:-Sigma := _passed[4];
        else
            if getPicardNumber(P) = 1 then
                # Input method (2)
                # This means, the picard number of X is one.
                numColumns := ColumnDimension(P:-mat);
                self:-Sigma := {seq({seq(1 .. numColumns)} minus {i}, i = 1 .. numColumns)};
            elif P:-s = 1 then
                # Input method (3)
                # This means, X is a C* surface.
                # This implementation is based on Construction 5.4.1.6 of "Cox Rings".
                # It is generalized slightly to support P-Matrices whose columns are
                # not necessarily slope-ordered.
                
                ordered_indices := [seq(sort([seq(1 .. P:-ns[i])], 
                    (j1, j2) -> getSlope(P, i, j1) > getSlope(P, i, j2)), 
                    i = 1 .. P:-r)];
                
                sigma_plus := {seq(doubleToSingleIndex(P:-format, i, ordered_indices[i,1]), i = 1 .. P:-r)};
                sigma_minus := {seq(doubleToSingleIndex(P:-format, i, ordered_indices[i,P:-ns[i]]) , i = 1 .. P:-r)};
                
                taus := {seq(seq({
                    doubleToSingleIndex(P:-format, i, ordered_indices[i,j]), 
                    doubleToSingleIndex(P:-format, i, ordered_indices[i,j+1])}, 
                    j = 1 .. P:-ns[i] - 1), i = 1 .. P:-r)};
                
                if P:-m = 0 then
                  # Case (e-e)
                  self:-Sigma := {sigma_plus} union taus union {sigma_minus};
                elif P:-m = 1 then
                    # Check if the zeros are where they should be
                    for k from 1 to P:-r + P:-s - 2 do
                        if not P:-mat[k, P:-n + P:-m] = 0 then
                            error "Here, s = m = 1, hence this P-matrix should belong to a C*-surface of type (p-e) or (e-p)."
                                  "But the last column is not of the right form for that, see section 5.4.1 of \"Cox Rings\".";
                        end if;
                    end do;
                    if P:-mat[P:-r - 1 + P:-s, P:-n + P:-m] = 1 then
                        # Case (p-e)
                        taus_plus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,1]), P:-n + P:-m} , i = 1 .. P:-r)};
                        self:-Sigma := taus_plus union taus union {sigma_minus};
                    elif P:-mat[P:-r - 1 + P:-s, P:-n + P:-m] = -1 then
                        # Case (e-p)
                        taus_minus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,P:-ns[i]]), P:-n + P:-m} , i = 1 .. P:-r)};
                        self:-Sigma := {sigma_plus} union taus union taus_minus;
                    else
                        error "Here, s = m = 1, hence this P-matrix should belong to a C*-surface of type (p-e) or (e-p)." 
                              "But the last column is not of the right form for that, see section 5.4.1 of \"Cox Rings\".";
                    end if;
                elif P:-m = 2 then
                    # Check if the zeros are where they should be
                    for k from 1 to P:-r + P:-s - 2 do
                        if not (P:-mat[k, P:-n + P:-m - 1] = 0 and P:-mat[k, P:-n + P:-m] = 0) then
                            error "Here, s = 1 and m = 2, hence this P-matrix should belong to a C*-surface of type (p-p)."
                                  "But the last two columns are not of the right form for that, see section 5.4.1 of \"Cox Rings\".";
                        end if;
                    end do;
                    # Check if there is +1 and -1 in the correct places
                    if not (P:-mat[P:-r - 1 + P:-s, P:-n + P:-m - 1] = 1 and P:-mat[P:-r - 1 + P:-s, P:-n + P:-m] = -1) then
                        error "Here, s = 1 and m = 2, hence this P-matrix should belong to a C*-surface of type (p-p)."
                              "But the last two columns are not of the right form for that, see section 5.4.1 of \"Cox Rings\".";
                    end if;
                    # Case (p-p)
                    taus_plus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,1]), P:-n + P:-m - 1} , i = 1 .. P:-r)};
                    taus_minus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,P:-ns[i]]), P:-n + P:-m} , i = 1 .. P:-r)};
                    self:-Sigma := taus_plus union taus union taus_minus;
                else
                    error "Here, s = 1, hence this P-matrix should belong to a C*-surface. But those cannot have m > 2, see section 5.4.1 of \"Cox Rings\".";
                end if;
            else
                error "This PMatrix is neither of Picard number one, nor is it the PMatrix of a surface."
                      "Therefore, you must provide the fan Sigma as input.";
            end if;
        end if;
    end;

    export getXCones :: static := proc(self :: TVarOne)
        if type(self:-XCones, undefined) or 'forceCompute' in [_passed] then
            self:-XCones := select(cone -> isXCone(self:-P:-format, cone), self:-Sigma);
        end if;
        return self:-XCones;
    end;

    export isMaximalXCone :: static := proc(self :: TVarOne, cone :: set(integer))
        local XCones, cone_;
        XCones := getXCones(self);
        for cone_ in XCones do
            if cone subset cone_ and cone <> cone_ then
                return false;
            end if;
        end do;
        return true;
    end;

    export setMaximalXCones :: static := proc(self :: TVarOne, maximalXCones :: set(set(integer))) 
        self:-maximalXCones := maximalXCones;
    end;

    export getMaximalXCones :: static := proc(self :: TVarOne)
        if type(self:-maximalXCones, undefined) or 'forceCompute' in [_passed] then
            setMaximalXCones(self, getMaximalXConesFormat(self:-P:-format, self:-Sigma));
        end if;
        return self:-maximalXCones;
    end;

    export setIsGorensteinVal :: static := proc(self :: TVarOne, isGorensteinVal :: boolean) 
        self:-isGorensteinVal := isGorensteinVal;
    end proc;

    export setGorensteinIndex :: static := proc(self :: TVarOne, gorensteinIndex :: integer) 
        self:-gorensteinIndex := gorensteinIndex;
        setIsGorensteinVal(self, evalb(gorensteinIndex = 1));
    end proc;

    export getGorensteinIndex :: static := proc(self :: TVarOne)
        local cone;
        if type(self:-gorensteinIndex, undefined) or 'forceCompute' in [_passed] then
            setGorensteinIndex(self, ilcm(seq(gorensteinIndexForXCone(self:-P, cone), cone in getMaximalXCones(self))));
        end if;
        return self:-gorensteinIndex;
    end;

    (*
    Checks whether the variety is gorenstein.
    *)
    export isGorenstein :: static := proc(self :: TVarOne)
        local cone;
        if type(self:-isGorensteinVal, undefined) or 'forceCompute' in [_passed] then
            # This method could simply be implemented by computing the gorenstein index and
            # asking if it is equal to one. However, if we are not interested in the gorenstein index
            # and just want to know if it's gorenstein or not, we can just look at the individual
            # gorenstein indices of the X-cones and terminate as soon as we find one X-cone, where it's
            # not equal to one. This has the advantage that if the variety is not Gorenstein, we find that
            # out early, without computing gorenstein indices for *every* X-cone.
            for cone in getMaximalXCones(self) do
                if gorensteinIndexForXCone(self:-P, cone) <> 1 then
                    # Not gorenstein, we are done.
                    setIsGorensteinVal(self, false);
                    break;
                end if;
            end do;
            if type(self:-isGorensteinVal, undefined) then
                setIsGorensteinVal(self, true);
                setGorensteinIndex(self, 1);
            end if;
        end if;
        return self:-isGorensteinVal;
    end;

    (****************************
    *** ADMISSIBLE OPERATIONS ***
    *****************************)

    (*
    Creates a new T-Variety of Complexity One by applying an admissible column operation. 
    See also `PMatrix[applyAdmissibleColumnOperation]`.
    *)
    export applyAdmissibleColumnOperation := proc(X :: TVarOne, sigma :: Perm, taus :: list(Perm), rho :: Perm)
        local newP, bundledPerm, newSigma;
        newP := PMatrix[applyAdmissibleColumnOperation](X:-P, sigma, taus, rho);
        bundledPerm := bundleColumnPermutation(X:-P:-format, sigma, taus, rho);
        newSigma := map(cones -> map(k -> bundledPerm[k], cones), X:-Sigma);
        return TVarOne(newP, newSigma);
    end proc;

    (*
    Creates a new T-Variety of Complexity One by applying an admissible column operation of type (1), i.e. 
    permutation of blocks. See also `PMatrix[applyAdmissibleColumnOperation1]`.
    *)
    export applyAdmissibleColumnOperation1 := proc(X :: TVarOne, sigma :: Perm)
        applyAdmissibleColumnOperation(X, sigma, [Perm([]) $ P:-r], Perm([]));
    end proc;

    (*
    Creates a new T-Variety of Complexity One by applying an admissible column operation of type (2), i.e.
    permutation of column within a given block `i`.See also `PMatrix[applyAdmissibleColumnOperation2]`.
    *)
    export applyAdmissibleColumnOperation2 := proc(X :: TVarOne, i :: integer, tau :: Perm)
        applyAdmissibleColumnOperation(X, Perm([]), [Perm([]) $ (i-1), tau, Perm([]) $ P:-r - i], Perm([]));
    end proc;

    (*
    Creates a new T-Variety of Complexity One by applying an admissible column operation of type (3), i.e.
    permutation of the last m rows. See also `PMatrix[applyAdmissibleColumnOperation3]`.
    *)
    export applyAdmissibleColumnOperation3 := proc(X :: TVarOne, rho :: Perm)
        applyAdmissibleColumnOperation(X, Perm([]), [Perm([]) $ P:-r], rho);
    end proc;

    (*
    Creates a new T-Variety of Complexity One by applying an admissible row operation. 
    See also `PMatrix[applyAdmissibleRowOperation]`.
    *)
    export applyAdmissibleRowOperation := proc(X :: TVarOne, A :: Matrix, B :: Matrix)
        return TVarOne(PMatrix[applyAdmissibleRowOperation](X:-P, A, B), X:-Sigma);
    end proc;

    (*
    Normalizes a T-Variety of Complexity One by removing redundant columns and applying 
    admissible column permutations to the P-Matrix.
    See also `normalizePMatrix`.
    *)
    export normalizeTVarOne := proc(X :: TVarOne)
        local newP, sigma_, taus_, bundledPerm, newSigma;
        newP, sigma_, taus_ := op(normalizePMatrix(X:-P, output = ['normalized', 'sigma', 'taus']));
        bundledPerm := bundleColumnPermutation(X:-P:-format, sigma_, taus_, Perm([]));
        newSigma := map(cones -> map(k -> bundledPerm[k], cones), X:-Sigma);
        return TVarOne(newP, newSigma);
    end proc;

    export ModulePrint :: static := proc(self :: TVarOne)
        nprintf(cat("TVarOne(dim = ", self:-P:-s + 1,
          ", lss = ", self:-P:-lss,
          ", Sigma = ", self:-Sigma));
    end;

    export TVarOneInfo :: static := proc(self :: TVarOne)
        local maximalXCones, gorensteinIndex, gorenstein;
        PMatrixInfo(self:-P);
        print(maximalXCones = getMaximalXCones(self));
        print(gorensteinIndex = getGorensteinIndex(self));
        print(gorenstein = isGorenstein(self));
    end;

end module:

ImportTVarOneList := proc(stmt)
    local columns, INDEX_S, INDEX_P, INDEX_DIMENSION, INDEX_PICARDNUMBER, INDEX_CLASSGROUPTORSION, INDEX_DEGREEMATRIX, INDEX_ANTICANCLASS, INDEX_AMBIENTFAN, INDEX_MAXIMALXCONES, INDEX_GORENSTEININDEX, INDEX_ISGORENSTEIN;
    local i, clm, result, P, X;

    columns := ColumnNames(stmt);
    for i from 0 to nops(columns)-1 do
        clm := columns[i+1];
        if clm = "s" then
            INDEX_S := i;
        elif clm = "P" then
            INDEX_P := i;
        elif clm = "dimension" then
            INDEX_DIMENSION := i;
        elif clm = "picardNumber" then
            INDEX_PICARDNUMBER := i;
        elif clm = "classGroupTorsion" then
            INDEX_CLASSGROUPTORSION := i;
        elif clm = "degreeMatrix" then
            INDEX_DEGREEMATRIX := i;
        elif clm = "anticanClass" then
            INDEX_ANTICANCLASS := i;
        elif clm = "ambientFan" then
            INDEX_AMBIENTFAN := i;
        elif clm = "maximalXCones" then
            INDEX_MAXIMALXCONES := i;
        elif clm = "gorensteinIndex" then
            INDEX_GORENSTEININDEX := i;
        elif clm = "isGorenstein" then
            INDEX_ISGORENSTEIN := i;
        end if;
    end do;

    if not (type(INDEX_S, integer) or type(INDEX_DIMENSION, integer)) then
        error "Not enough data. You must either provide a column \"s\" (dimension of the acting torus) or a column \"dimension\" (the dimension of the variety = s + 1).";
    end if;

    if not type(INDEX_S, integer) then
        error "Not enough data. You must provide a column \"P\" whose entries are the rows of the P-Matrix.";
    end if;

    result := [];
    while Step(stmt) = RESULT_ROW do
        P := PMatrix(Fetch(stmt, INDEX_S), Matrix(parse(Fetch(stmt, INDEX_P))));
        # Read in any already avaliable metadata about P, if available
        if type(INDEX_PICARDNUMBER, integer) then
            setPicardNumber(P, Fetch(stmt, INDEX_PICARDNUMBER));
            if type(INDEX_CLASSGROUPTORSION, integer) then
                setClassGroup(P, createAG(getPicardNumber(P), parse(Fetch(stmt, INDEX_CLASSGROUPTORSION))));
            end if;
        end if;
        if type(INDEX_ANTICANCLASS, integer) then
            setAnticanClass(P, Vector(parse(Fetch(stmt, INDEX_ANTICANCLASS))));
        end if;
        if type(INDEX_DEGREEMATRIX, integer) then
            setQ(P, Matrix(parse(Fetch(stmt, INDEX_DEGREEMATRIX))));
        end if;

        # If a fan is supplied, make a TVarOne.
        if type(INDEX_AMBIENTFAN, integer) then
            X := TVarOne(P, parse(Fetch(stmt, INDEX_AMBIENTFAN)));
            if type(INDEX_MAXIMALXCONES, integer) then
                setMaximalXCones(X, parse(Fetch(stmt, INDEX_MAXIMALXCONES)));
            end if;
            if type(INDEX_GORENSTEININDEX, integer) then
                setGorensteinIndex(X, Fetch(stmt, INDEX_GORENSTEININDEX));
            end if;
            if type(INDEX_ISGORENSTEIN, integer) then
                setIsGorensteinVal(X, evalb(Fetch(stmt, INDEX_ISGORENSTEIN) <> 0));
            end if;
            result := [op(result), X];
        else
            result := [op(result), P];
        end if;
    end do;

    return result;

end proc;

(* 
Attempts to find a given Complexity-1 variety `X` in the table `tableName` of a given SQLite connection `db`.
If it finds the variety, the rowid of the database entry is returned. Otherwise, it returns -1.
*)
FindInDatabase := proc(connection, tableName :: string, X0 :: TVarOne)
    # First, normalize the P-Matrix.
    X := normalizeTVarOne(X0);
    P := X:-P;
    stmt := Prepare(connection, cat("SELECT rowid,s,P FROM ", tableName, " WHERE ",
        "lss = \"", P:-lss, "\" AND ",
        "m = ", P:-m, " AND ",
        "s = ", P:-s, " AND ",
        "picardNumber = ", getPicardNumber(P), " AND ",
        "classGroupTorsion = \"", AGdata(getClassGroup(P))[4], "\" AND ",
        "anticanClass = \"", convert(getAnticanClass(P), list), "\" AND ",
        "gorensteinIndex = ", getGorensteinIndex(X)));
    Ps := ImportTVarOneList(stmt);
    M := FetchAll(stmt):
    rowids := [seq(M[i,1], i = 1 .. RowDimension(M))];
    for i from 1 to nops(rowids) do
        if areEquivalent(P, Ps[i]) then
            return rowids[i];
        end if;
    end do;
    return -1;

end proc;

(*
Given a list of varieties of complexity one `Xs` and a SQLite database connection `db`, this function
inserts those varieties from `Xs` into the database that are not already present.
*)
ExportTVarOneList := proc(connection, tableName :: string, Xs :: list(TVarOne))
    local k, X, P, stmt, i;
    knownCount := 0;

    for k from 1 to nops(Xs) do
        # Normalize the given T-Variety.
        X := normalizeTVarOne(Xs[k]);
        # Only insert `X` if it is not already present.
        if FindInDatabase(connection, tableName, X) = -1 then
            P := X:-P;
            stmt := Prepare(connection, cat("INSERT INTO ", tableName ," VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"));
            Bind(stmt, 1, P:-r, valuetype = "integer");
            Bind(stmt, 2, convert(P:-ns, string), valuetype = "text");
            Bind(stmt, 3, P:-n, valuetype = "integer");
            Bind(stmt, 4, P:-m, valuetype = "integer");
            Bind(stmt, 5, P:-s, valuetype = "integer");
            Bind(stmt, 6, convert(P:-lss, string), valuetype = "text");
            Bind(stmt, 7, convert([seq(convert(Row(P:-mat, i), list), i = 1 .. RowDimension(P:-mat))], string), valuetype = "text");
            Bind(stmt, 8, P:-dim, valuetype = "integer");
            Bind(stmt, 9, getPicardNumber(P), valuetype = "integer");
            Bind(stmt, 10, convert(AGdata(getClassGroup(P))[4], string), valuetype = "text");
            Bind(stmt, 11, convert([seq(convert(Row(getQ(P), i), list), i = 1 .. RowDimension(getQ(P)))], string), valuetype = "text");
            Bind(stmt, 12, convert(convert(getAnticanClass(P), list), string), valuetype = "text");
            Bind(stmt, 13, convert(X:-Sigma, string), valuetype = "text");
            Bind(stmt, 14, convert(getMaximalXCones(X), string), valuetype = "text");
            Bind(stmt, 15, getGorensteinIndex(X), valuetype = "integer");
            Bind(stmt, 16, isGorenstein(X), valuetype = "integer");

            Step(stmt);
            Finalize(stmt);
            if 'logging' in [_passed] then
                print(cat("Variety no. ", k, " inserted."));
            end if;
        else
            knownCount := knownCount + 1;
            if 'logging' in [_passed] then
                print(cat("Variety no. ", k, " already known."));
            end if;
        end if;
    end do;

    print(cat("Supplied ", nops(Xs), " varieties. Already known: ", knownCount, ". Inserted ", nops(Xs) - knownCount, " new varieties into the database."));

end proc;

end module:
