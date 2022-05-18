ComplexityOne := module()

option package;

export PFormat, PMatrix, TVarOne, FindInDatabase, ImportTVarOneList, ImportTVarOne, ExportTVarOneList, performOnDatabase;

## TODO: Remove dependency on MDSpackage.
uses LinearAlgebra, Database[SQLite];

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

(*****************************************************************************
**********************             PFORMAT             ***********************
******************************************************************************)

module PFormat()
    option object;

    export r, ns, n, m, s, dim, picardNumber;

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
        self:-picardNumber := self:-n + self:-m - (self:-r - 1 + self:-s);
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
        return doubleToSingleIndex(newFormat, sigma[i], taus[i][j]);
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

(*****************************************************************************
**********************             PMATRIX             ***********************
******************************************************************************)

module PMatrix()
    option object;

    #########################################################################
    ## These fields are guaranteed to be filled when a PMatrix is created. ##
    #########################################################################

    # Format data. Note that the picardNumber is just n+m - (r+s), hence we consider it
    # part of the format.
    export format, r, ns, n, m, s, dim, picardNumber;

    # The exponents of the relations in the Cox Ring, encoded as a list of lists.
    export lss;
    
    # Matrix data: `mat` is the entire matrix, `d` is only the lower `s` rows, `P0` only 
    # the upper r-1 rows.
    export mat, d, P0;

    ####################################################################################
    ## These fields are only computed when needed. Use these getters below for these. ##
    ####################################################################################

    # The divisor class group, encoded as a list of positive integers, where the first
    # entry is the rank (= picardNumber) and the other entries are the elementary divisors of the torsion part
    local classGroup := undefined;

    # Degree matrix of the Cox Ring
    local Q := undefined;
    
    # Free part of the degree matrix, i.e. the first `picardNumber` rows of `Q`
    local Q0 := undefined;

    # Anticanonical class as a vector in the divisor class group
    local anticanClass := undefined;
    local anitcanCoefficients := undefined;

    # Effective and moving cone in the rational vector space of the divisor class group.
    local effectiveCone := undefined;
    local movingCone := undefined;

    # Says whether there exists a Fano variety having P as its PMatrix.
    # This is equivalent to the anticanonical divisor lying in the moving cone.
    local admitsFanoVal := undefined;

    # Says whether a variety having P as its PMatrix is toric.
    # For irredundant P, this is equivalent to P:-r = 2.
    local isToricVal := undefined;

    # Says whether the PMatrix is irredundant, i.e. has no redundant blocks consisting of
    # a single one. You can use the procedure `removeRedundantBlocks` on a P-Matrix to pass
    # to get an irredundant PMatrix equivalent to the original one.
    local isIrredundantVal := undefined;

    ###########################################################################
    ## These fields are only defined for P-Matrices of surfaces, i.e. s = 1. ##
    ## All of them are computed when the P-Matrix is created.                ##
    ###########################################################################

    # The case of the P-Matrix, i.e. its configuration of elliptic fixed points (E) 
    # vs. parabolic fixed point curves (P). It is exactly one of the fixe strings:
    #                     "EE", "PE", "EP", "PP+", "PP-".
    # The naming follows the literature, e.g. Cox Rings chapter 5.4.1, except that we
    # additionally differentiate between "PP+" and "PP-". The former has [0...1],[0...-1] as
    # last two columns, the latter has [0...-1],[0...1].
    export case := undefined;

    # The slopes of the rays in the fan. This is encoded as a list of lists, following the
    # same indexing as in the `lss`.
    export slopes := undefined;

    # The sum of the biggest resp. smallest slopes in each block.
    export mplus := undefined;
    export mminus := undefined;

    local setFormat :: static := proc(self :: PMatrix, f :: PFormat)
        self:-format := f;
        self:-r := f:-r;
        self:-ns := f:-ns;
        self:-n := f:-n;
        self:-m := f:-m;
        self:-s := f:-s;
        self:-dim := f:-dim;
        self:-picardNumber := f:-picardNumber;
    end proc;

    local setSurfaceData :: static := proc(self :: PMatrix, d :: Matrix, lss :: list(list(integer)))
        local i, j;
        self:-slopes := [seq([seq(d[1,doubleToSingleIndex(self:-format, i, j)] / lss[i,j], j = 1 .. self:-ns[i])], i = 1 .. self:-r)];
        self:-mplus := add([seq(max(self:-slopes[i]), i = 1 .. self:-r)]);
        self:-mminus := add([seq(min(self:-slopes[i]), i = 1 .. self:-r)]);
        if self:-m = 0 then
            self:-case := "EE";
        elif self:-m = 1 then
            if d[1, doubleToSingleIndex(self:-format, -1, 1)] = 1 then
                self:-case := "PE";
            elif d[1, doubleToSingleIndex(self:-format, -1, 1)] = -1 then
                self:-case := "EP";
            end if;
        elif self:-m = 2 then
            if d[1, doubleToSingleIndex(self:-format, -1, 1)] = 1 then
                self:-case := "PP+";
            elif d[1, doubleToSingleIndex(self:-format, -1, 1)] = -1 then
                self:-case := "PP-";
            end if;
        end if;
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

        local lss, ls, l, format, d, rows, rows0, i, j, P, r, s, ns, numZerosBefore;
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

                # If this is a P-Matrix of a surface, set the slopes
                if self:-s = 1 then 
                    setSurfaceData(self, self:-d, self:-lss);
                end if;

                # Construct the P-matrix from the given data
                rows := [seq([seq(-self:-lss[1]), (0 $ add(self:-ns[2..i-1]),
                            seq(self:-lss[i])), (0 $ add(self:-ns[i+1..self:-r]) + self:-m)], i = 2 .. self:-r),
                        seq(:-convert(Row(self:-d, j), list), j = 1 .. self:-s)];
                self:-mat := Matrix(rows);

                # Only the L-block. This gives P0.
                rows0 := [seq([seq(-self:-lss[1]), (0 $ add(self:-ns[2..i-1]),
                            seq(self:-lss[i])), (0 $ add(self:-ns[i..self:-r-1]))], i = 2 .. self:-r)];
                self:-P0 := Matrix(rows0);

                if not 'skipChecks' in [_passed] then
                    assertColumnsPrimitive(self);
                    assertColumnsGenerateFullCone(self);
                end if;
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

            lss := _passed[3];
            d := _passed[4];
            ns := map(nops, lss);
            format := PFormat(ns, ColumnDimension(d) - add(ns), RowDimension(d));
            setFormat(self, format);
            P := PMatrix(format, lss, d);
            self:-lss := P:-lss;
            self:-mat := P:-mat;
            self:-P0 := P:-P0;
            self:-d := P:-d;

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

            if not 'skipChecks' in [_passed] then
                assertColumnsPrimitive(self);
                assertColumnsGenerateFullCone(self);
            end if;

            self:-d := SubMatrix(P, [self:-r .. RowDimension(P)], [1 .. ColumnDimension(P)]);
            
            # If this a P-Matrix of a surface, compute the slopes
            if self:-s = 1 then
                setSurfaceData(self, self:-d, self:-lss);
            end if;
        end if;

    end;

    export setClassGroup :: static := proc(self :: PMatrix, classGroup :: list(integer)) self:-classGroup := classGroup; end proc;

    export setQ :: static := proc(self :: PMatrix, Q :: Matrix) 
        self:-Q := Q; 
        self:-Q0 := DeleteRow(Q, [seq(self:-picardNumber + 1 .. RowDimension(Q))]);
    end proc;

    export setQ0 :: static := proc(self :: PMatrix, Q0 :: Matrix) self:-Q0 := Q0; end proc;

    export setAnticanCoefficients :: static := proc(self :: PMatrix, anitcanCoefficients) self:-anitcanCoefficients := anitcanCoefficients; end proc;

    export setAnticanClass :: static := proc(self :: PMatrix, anticanClass) self:-anticanClass := anticanClass; end proc;

    export setMovingCone :: static := proc(self :: PMatrix, movingCone :: CONE) self:-movingCone := movingCone; end proc;
    
    export setEffectiveCone :: static := proc(self :: PMatrix, effectiveCone :: CONE) self:-effectiveCone := effectiveCone end proc;

    export setAdmitsFanoVal :: static := proc(self :: PMatrix, admitsFanoVal :: boolean) self:-admitsFanoVal := admitsFanoVal; end proc;

    export setIsToricVal :: static := proc(self :: PMatrix, isToricVal :: boolean) self:-isToricVal := isToricVal; end proc;
    
    export setIsIrredundantVal :: static := proc(self :: PMatrix, isIrredundantVal :: boolean) self:-isIrredundantVal := isIrredundantVal; end proc;

    (*
    Compute the smith normal form of the transpose of the matrix.
    From this we can read off the degree matrix.
    *)
    local computeSmithForm :: static := proc(self :: PMatrix)
        local S_, U_, classGroup, i, degreeMatrixTorsion;
        # Compute the Smith normal form.
        S_, U_ := SmithForm(Transpose(self:-mat), output = ['S','U']);
        # The first entry of `classGroup` is the picard number, i.e. the rank of the free part
        classGroup := [self:-picardNumber];
        # For the torsion part, we traverse the diagonal of `S` and add all entries that are not equal to one
        # Note that they are already given in ascending order by `SmithForm`
        for i from 1 to ColumnDimension(S_) do
            if S_[i,i] <> 1 then
                classGroup := [op(classGroup), S_[i,i]];
            end if;
        end do;
        setClassGroup(self, classGroup);
        # Now read off the degree matrix from `U`.
        setQ0(self, IntegerRelations[LLL](DeleteRow(U_, [seq(1 .. ColumnDimension(S_))])));
        degreeMatrixTorsion := Matrix(nops(classGroup) - 1, self:-n + self:-m,
            [seq(map(x -> x mod classGroup[i+1], :-convert(Row(U_, ColumnDimension(S_) - (nops(classGroup) - 1) + 1), list)), i = 1 .. nops(classGroup) - 1)]);
        setQ(self, Matrix(self:-picardNumber + nops(classGroup) - 1, self:-n + self:-m, [[self:-Q0],[degreeMatrixTorsion]]));
    end;

    export getClassGroup :: static := proc(self :: PMatrix)
        if type(self:-classGroup, undefined) or 'forceCompute' in [_passed] then
            computeSmithForm(self);
        end if;
        return self:-classGroup;
    end;

    export getQ :: static := proc(self :: PMatrix)
        if type(self:-Q, undefined) or 'forceCompute' in [_passed] then
            computeSmithForm(self);
        end if;
        return self:-Q;
    end;

    export getQ0 :: static := proc(self :: PMatrix)
        if type(self:-Q0, undefined) or 'forceCompute' in [_passed] then
            computeSmithForm(self);
        end if;
        return self:-Q0;
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
            for i from 1 to nops(getClassGroup(self)) - 1 do
                anticanVec[self:-picardNumber + i] := anticanVec[self:-picardNumber + i] mod getClassGroup(self)[i+1];
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

    export getEffectiveCone :: static := proc(self :: PMatrix)
        if type(self:-effectiveCone, undefined) or 'forceCompute' in [_passed] then
            setEffectiveCone(self, poshull(Column(getQ0(self), [seq(1 .. self:-n + self:-m)])));
        end if;
        return self:-effectiveCone;
    end proc;

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
    export admitsFano :: static := proc(self :: PMatrix)
        if type(self:-admitsFanoVal, undefined) or 'forceCompute' in [_passed] then
            setAdmitsFanoVal(self, containsrelint(getMovingCone(self), getAnticanClass(self)));
        end if;
        return self:-admitsFanoVal;
    end;

    export isToric :: static := proc(P :: PMatrix)
        if type(P:-isToricVal, undefined) or 'forceCompute' in [_passed] then
            setIsToricVal(P, evalb(removeRedundantBlocks(P):-r = 2));
        end if;
        return P:-isToricVal;
    end proc;

    export isIrredundant :: static := proc(P :: PMatrix)
        local redundantIndices;
        if type(P:-isIrredundantVal, undefined) or 'forceCompute' in [_passed] then
            redundantIndices := select(i -> P:-lss[i] = [1], [seq(1 .. nops(P:-lss))]);
            setIsIrredundantVal(P, evalb(P:-r = 2 or redundantIndices = []));
        end if;
        return P:-isIrredundantVal;        
    end proc;

    (*
    Computes the local class group of the variety associated to a P-Matrix.
    The first argument is the P-Matrix. The second argument is a cone corresponding
    to a toric orbit, encoded as a list of integers.
    The output is encoded as a list of integers, where the first entry is the rank of
    the free part and all remaining entries are the elementary divisors of the torsion part.
    *)
    export getLocalClassGroup :: static := proc(self :: PMatrix, cone :: set(integer))
        local Psub, S_, U_, localRank, localClassGroup, i;
        # We need to compute a smith form of a submatrix of the P-Matrix, where
        # we removed all columns that do not occur on `cone`.
        Psub := DeleteColumn(self:-mat, [op({seq(1 .. self:-n + self:-m)} minus cone)]);
        S_, U_ := SmithForm(Transpose(Psub), output = ['S', 'U']);
        # The rank of the local class group
        localRank := max(0, RowDimension(S_) - ColumnDimension(S_));
        localClassGroup := [localRank];
        # For the torsion part, we traverse the diagonal of `S` and add all entries that are not equal to one
        # Note that they are already given in ascending order by `SmithForm`
        for i from 1 to min(ColumnDimension(S_), RowDimension(S_)) do
            if S_[i,i] <> 1 then
                localClassGroup := [op(localClassGroup), S_[i,i]];
            end if;
        end do;
        return localClassGroup;
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
    export applyAdmissibleColumnOperation :: static := proc(P :: PMatrix, sigma :: Perm, taus :: list(Perm), rho :: Perm)
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
    export applyAdmissibleColumnOperation1 :: static := proc(P :: PMatrix, sigma :: Perm)
        applyAdmissibleColumnOperation(P, sigma, [Perm([]) $ P:-r], Perm([]));
    end proc;

    (*
    Creates a new P-Matrix by applying an admissible column operation of type (2), i.e.
    permutation of column within a given block `i`. The parameter `tau` must be a permutation
    of the numbers [1 .. ns[i]].
    *)
    export applyAdmissibleColumnOperation2 :: static := proc(P :: PMatrix, i :: integer, tau :: Perm)
        applyAdmissibleColumnOperation(P, Perm([]), [Perm([]) $ (i-1), tau, Perm([]) $ P:-r - i], Perm([]));
    end proc;

    (*
    Creates a new P-Matrix by applying an admissible column operation of type (3), i.e.
    permutation of the last m rows. The parameter `rho` must be a permutation of the
    numbers [1 .. m].
    *)
    export applyAdmissibleColumnOperation3 :: static := proc(P :: PMatrix, rho :: Perm)
        applyAdmissibleColumnOperation(P, Perm([]), [Perm([]) $ P:-r], rho);
    end proc;

    (*
    Creates a new P-Matrix by applying an admissible row operation. The row operation is specified by
    an (s x (r-1))-Matrix `A` and an (s x s)-Matrix `B`. The row operation is applied by multiplying `P`
    from the left with an elementary matrix as follows:

    ( E 0 )  *  ( L 0 )  =  (    L      0  )
    ( A B )     ( d d')     ( AL + Bd   Bd')

    *)
    export applyAdmissibleRowOperation :: static := proc(P :: PMatrix, C :: Matrix, T :: Matrix)
        local identityMatrix, zeroMatrix, newP;
        if not type(C, 'Matrix'(P:-s, P:-r - 1, integer)) then
            error "Expected second argument to be of type: Matrix(%1, %2, integer)", P:-s, P:-r - 1;
        end if;
        if not type(T, 'Matrix'(P:-s, P:-s, integer)) then
            error "Expected third argument to be of type: Matrix(%1, %1, integer)", P:-s;
        end if;
        if not abs(Determinant(T)) = 1 then
            error "Third argument must be a unimodular matrix, i.e. its determinant must be 1 or -1";
        end if;
        # An (r-1) x (r-1) identity matrix.
        identityMatrix := Matrix(P:-r - 1, P:-r - 1, shape = diagonal, 1);
        # An (r-1) x s zero matrix.
        zeroMatrix := Matrix(P:-r - 1, P:-s, fill = 0);
        # Multiply P from the left with
        # ( E  0 )
        # ( A  B )
        newP := <<identityMatrix | zeroMatrix>, <C | T>> . P:-mat;
        return PMatrix(P:-s, newP);
    end proc;

    (******************
    ** NORMALIZATION **
    *******************)

    (*
    Removes a single redundant column from a P-Matrix.
    *)
    export removeSingleRedundantBlock :: static := proc(P_ :: PMatrix, i0 :: integer)
        local P, A, B, newLss, newD, newFormat, i;
        P := P_;
        # First, we have to apply admissible row operations to achieve all zeros in the d-block under
        # the redundant columns. We construct the A-Matrix necessary for this.
        # This step is necessary to ensure the columns of the resulting P-Matrix still generate the whole
        # space as a cone.
        if i0 = 1 then
            A := Matrix(P:-s, P:-r - 1, (k,l) -> if l = 1 then P:-d[k, add(P:-ns[1 .. i0-1]) + 1] else 0 end if);
        else
            A := Matrix(P:-s, P:-r - 1, (k,l) -> if l = i0 - 1 then -P:-d[k, add(P:-ns[1 .. i0-1]) + 1] else 0 end if);
        end if;
        # Just the identity
        B := Matrix(P:-s, P:-s, shape = diagonal, 1);
        P := applyAdmissibleRowOperation(P, A, B);
        # Now construct the new P-Matrix data
        newLss := [seq(P:-lss[i], i in {seq(1 .. P:-r)} minus {i0})];
        newFormat := PFormat([seq(P:-ns[i], i in {seq(1 .. P:-r)} minus {i0})], P:-m, P:-s);
        newD := DeleteColumn(P:-d, add(P:-ns[1 .. i0-1]) + 1);
        return PMatrix(newFormat, newLss, newD);
    end proc;

    (*
    Removes all redundant blocks from a P-Matrix, i.e. all blocks with a single 1 inside the L-block.
    Note that this can change the lower `s` rows of a P-Matrix, as admissible row operations may be 
    necessary to achieve all zeros in the columns under the redundant blocks.
    *)
    export removeRedundantBlocks :: static := proc(P :: PMatrix)
        local redundantIndices;
        redundantIndices := select(i -> P:-lss[i] = [1], [seq(1 .. nops(P:-lss))]);
        # If there is a redundant block and we still have more than two blocks, remove it.
        if nops(redundantIndices) > 0 and P:-r > 2 then
            # Recursive call
            return removeRedundantBlocks(removeSingleRedundantBlock(P, redundantIndices[1]));
        else
            return P;
        end if;
    end proc;

    (*
    Applies admissible operations to sort the columns of a P-Matrix by the entries in the L-block.
    The columns of block are sorted descendingly. The blocks themselves are sorted first descendingly
    according to size, and within the same size lexicographically by the lss.
    *)
    export sortColumnsByLss :: static := proc(P0 :: PMatrix)
        local P1, P2, P3, sortKey, taus_, i, compfun, tau, sigma_, result, str;

        # First, remove any redundant columns
        P1 := removeRedundantBlocks(P0);
        # Sort each individual block
        taus_ := [];
        for i from 1 to P1:-r do
            tau := Perm(sort([seq(1 .. P1:-ns[i])], (j1, j2) -> P1:-lss[i][j1] > P1:-lss[i][j2], 'output' = 'permutation'))^(-1);
            taus_ := [op(taus_), tau];
        end do;
        P2 := applyAdmissibleColumnOperation(P1, Perm([]), taus_, Perm([]));

        # This is the comparison function we will use to sort the blocks themselves.
        # It returns true if the block of index `i1` should occur before the block of index `i2`.
        compfun := proc(P :: PMatrix, i1 :: integer, i2 :: integer)
            local j;
            if nops(P:-lss[i1]) < nops(P:-lss[i2]) then return false;
            elif nops(P:-lss[i1]) > nops(P:-lss[i2]) then return true;
            else
                for j from 1 to nops(P:-lss[i1]) do
                    if P:-lss[i1][j] < P:-lss[i2][j] then return false;
                    elif P:-lss[i1][j] > P:-lss[i2][j] then return true;
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

    (*
    Sort the columns of a P-Matrix of a K*-surface descendingly by the values of the alphas (adjusted slopes).
    This is experimental and only works for surfaces at the moment.
    *)
    export sortColumnsByAdjustedSlopes :: static := proc(P0 :: PMatrix)
        local P1, P2, P3, sortKey, taus_, i, compfun, tau, sigma_, result, str;

        # First, remove any redundant columns
        P1 := removeRedundantBlocks(P0);

        # The sortkey to sort the columns by
        sortKey := proc(P :: PMatrix, i :: integer, j :: integer)
            local i_;
            P:-d[1,doubleToSingleIndex(P:-format, i, j)] / P:-lss[i][j]
                + add([seq(floor(P:-d[1,doubleToSingleIndex(P:-format, i_, 1)] / P:-lss[i_][1]), i_ in ({seq(1..P:-r)} minus {i}))]);
        end proc;

        # Sort each individual block
        taus_ := [];
        for i from 1 to P1:-r do
            tau := Perm(sort([seq(1 .. P1:-ns[i])], (j1, j2) -> sortKey(P1, i, j1) > sortKey(P1, i, j2), 'output' = 'permutation'))^(-1);
            taus_ := [op(taus_), tau];
        end do;
        P2 := applyAdmissibleColumnOperation(P1, Perm([]), taus_, Perm([]));

        # This is the comparison function we will use to sort the blocks themselves.
        # It returns true if the block of index `i1` should occur before the block of index `i2`.
        compfun := proc(P :: PMatrix, i1 :: integer, i2 :: integer)
            local j;
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
    where `out` is a list of the names 'result', 'C', 'T' and 'S'.
    *)
    export areRowEquivalent :: static := proc(P1 :: PMatrix, P2 :: PMatrix)
        local resBool, resC, resT, resS, identityMatrix, zeroMatrix, C_, T_, S_, newP, sol, resultList, str, i, j;
        resBool, resC, resT, resS := false, undefined, undefined, undefined;
        
        # P-Matrices of different tower structure (lss) are never row-equivalent.
        if P1:-lss = P2:-lss then
            identityMatrix := Matrix(P1:-r - 1, P1:-r - 1, shape = diagonal, 1);
            zeroMatrix := Matrix(P1:-r - 1, P1:-s, fill = 0);
            C_ := Matrix(P1:-s, P1:-r - 1, symbol = 'a');
            T_ := Matrix(P1:-s, P1:-s, symbol = 'b');
            S_ := <<identityMatrix | zeroMatrix>, <C_ | T_>>;
            newP := S_ . P1:-mat;
            sol := isolve({seq(seq(newP[i,j] = P2:-mat[i,j] , j = 1 .. ColumnDimension(P1:-mat)), i = P1:-r .. RowDimension(P1:-mat))});
            if sol = NULL then
                resBool, resC, resT, resS := false, undefined, undefined, undefined;
            else
                if abs(Determinant(subs(sol, T_))) = 1 then
                    resBool := true;
                    resC := subs(sol, C_);
                    resT := subs(sol, T_);
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
                        elif str = 'C' then
                            resultList := [op(resultList), resC];
                        elif str = 'T' then
                            resultList := [op(resultList), resT];
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

    Currently, we are working with two different algorithms for equivalence checking.
    The first one only works for surfaces, but is a lot more efficient (sort the columns by the adjusted slopes)
    The second one works in general, but has exponential running time in the number of blocks (roughly speaking).
    Eventually, the first one should be generalized to get an efficient algorithm in general, but this will require some
    theoretical work first.

    *)
    export areEquivalent :: static := proc(P1_ :: PMatrix, P2_ :: PMatrix)
        local Ps, P, i, P1, P2, P11, P12, newP1, sol, j;

        # First, remove any redundant blocks.
        P1 := removeRedundantBlocks(P1_);
        P2 := removeRedundantBlocks(P2_);

        # After removing redundant blocks, the number of blocks must coincide.
        if P1:-r <> P2:-r then
            return false;
        end if;

        # Varieties of different dimension can't be isomorphic.
        if P1:-s <> P2:-s then
            return false;
        end if; 

        if P1:-r = 2 then
            # TORIC CASE
            # In this case, the question is purely toric. We need to find an invertible integer 
            # matrix `S` sending P1 to P2.
            newP1 := Matrix(RowDimension(P1:-mat), RowDimension(P1:-mat), symbol = 'x') . P1:-mat;
            sol := isolve({seq(seq(newP1[i,j] = P2:-mat[i,j] , j = 1 .. ColumnDimension(P1:-mat)), i = 1 .. RowDimension(P1:-mat))});
            return evalb(sol <> NULL);
        elif P1:-s = 1 then
            # SURFACE CASE
            # We sort the columns by adjusted slopes
            if not 'skipSorting' in [_passed] then
                P1 := sortColumnsByAdjustedSlopes(P1);
                P2 := sortColumnsByAdjustedSlopes(P2);
            end if;
            P11 := P1;
            P12 := sortColumnsByAdjustedSlopes(applyAdmissibleRowOperation(P1, Matrix([[0 $ P1:-r - 1]]), Matrix([[-1]])));
            # After sorting by adjusted slopes, P1 and P2 are equivalent if and only if they are row-equivalent
            # i.e. no need to worry about column permutations any more.
            # (With the caveat that we don't fix the direction of the P-Matrix, which is why we work both with P1 and its negative at the same time)
            return areRowEquivalent(P11, P2) or areRowEquivalent(P12, P2);
        else
            # GENERAL CASE
            # First, sort the columns by the lss.
            if 'skipSorting' in [_passed] then
                P1 := P1_;
                P2 := P2_;
            else 
                P1 := sortColumnsByLss(P1_);
                P2 := sortColumnsByLss(P2_);
            end if;
            # After sorting, the lss must be equal, otherwise `P1` and `P2` can't be equivalent.
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
        end if;
    end proc;

    export ModulePrint :: static := proc(self :: PMatrix)
        nprintf(cat("PMatrix(", self:-lss, ", m = ", self:-m, ", s = ", self:-s, ")"));
    end;

    export PMatrixInfo :: static := proc(self :: PMatrix)
        local P, Q, n, m, picardNumber, classGroup, anticanClass, admitsFano, i;
        print(P = self:-mat);
        print([seq(cat(n,i), i = 0 .. self:-r - 1), m] = [seq(self:-ns[i], i = 1 .. self:-r), self:-m]);
        print(Q = getQ(self));
        print(classGroup = getClassGroup(self));
        print(picardNumber = self:-picardNumber);
        print(anticanClass = getAnticanClass(self));
        print(effectiveConeRays = rays(getEffectiveCone(self)));
        print(movingConeRays = rays(getMovingCone(self)));
        print(admitsFano = PMatrix[admitsFano](self));
    end;

    export convert :: static := proc(self :: PMatrix, toType, $)
        if toType = ':-Matrix' then
            self:-mat;
        elif toType = 'linalg:-matrix' or toType = 'matrix' then
            :-convert(self:-mat, matrix);
        else
            error "cannot convert from PMatrix to %1", toType;
        end if;
    end;

end module:


(*****************************************************************************
**********************             TVARONE             ***********************
******************************************************************************)

module TVarOne()
    option object;

    # These fields are guaranteed to be filled when a TVarOne is created.
    export P;
    export Sigma := undefined;
    export A := undefined;

    # Variables and relations in the Cox Ring
    export variables, monomials, relations;

    # These fields are only computed when needed. Use the getters below for them.
    local maximalXCones := undefined;
    local gorensteinIndex := undefined;
    local isGorensteinVal := undefined;
    local ampleCone := undefined;
    local isFanoVal := undefined;

    # The following fields are only defined for K*-surfaces, i.e. when s = 1.
    local intersectionTable := undefined;
    local anticanonicalSelfIntersection := undefined;

    export ModuleApply :: static := proc()
        Object(TVarOne, _passed);
    end;

    local setRelations :: static := proc(self :: TVarOne, P :: PMatrix, A :: Matrix)
        local f, lss, i, j, alpha;
        f := P:-format;
        lss := P:-lss;
        self:-variables := [seq(seq(T[i,j], j = 1 .. f:-ns[i]), i = 1 .. f:-r), seq(S[i], i = 1 .. f:-m)];
        self:-monomials := [seq(mul([seq(T[i,j] ^ lss[i][j], j = 1 .. f:-ns[i])]), i = 1 .. f:-r)];

        alpha := (i,j) -> Determinant(Matrix([Column(A, [i,j])]));
        # TODO: Add support for an optional parameter A during creation of the P-Matrix to modify the coefficients in the relations.
        self:-relations := [seq(alpha(i+1,i+2) * self:-monomials[i] + alpha(i+2, i) * self:-monomials[i+1] + alpha(i, i+1) * self:-monomials[i+2], i = 1 .. f:-r - 2)];
    end proc;

    (*
    This method creates a T-Variety of complexity one from various kinds of data. 
    The standard input method looks like this:

    P :: PMatrix, Sigma :: {set(set(integer)), Vector, list}, A :: Matrix.

    In some cases (see below), `Sigma` is optional. `A` is always optional.

    The parameter `Sigma` encodes the fan of the ambient toric variety. It can either be given directly as a set of
    cones, where each cone is encoded by the set of column indices of P it contains, or it can be given by 
    a weight in the rational vector space associated to the class group of P lying in the moving cone of P. 
    In the latter case, we take the fan to be the gale dual of the bunch of orbit cones generated by that weight.

    The parameter `Sigma` is optional in the following three cases:

    (1) s = 1. In this case, we have a C*-surface and there is only one possible minimal ambient fan, which we
        can explicitly construct from P (see Cox Rings 5.4.1.6).
    (2) picardNumber = 1. In this case, the number of rays of the fan is one more than the ambient dimension, hence
        there is only one complete fan having P as generator matrix.
    (3) P admits a Fano variety. In this case, the anticanonical class is contained in the moving cone of P, so we can 
        use it to define a bunch of cones and hence an ambient fan. The resulting variety will be the unique Fano variety
        having P as its PMatrix. Note however, that there are other possible non-Fano varieties with P as PMatrix.

    The parameter `A` is the coefficient matrix for the trinomial euqations defining the variety.
    It is an optional parameter. If it is not proved, the standard coefficient matrix is used, hence we
    have one's everywhere in the defining trinomials.

    *)
    export ModuleCopy :: static := proc(self :: TVarOne, proto :: TVarOne, P :: PMatrix)
        local numColumns, i, j, k, ordered_indices, taus, sigma_plus, sigma_minus, taus_plus, taus_minus, w, candidates, minimalBunchCones, cone;

        self:-P := P;

        # Check the arguments, set Sigma and A if present.
        if _npassed > 3 then
            if type(_passed[4], Matrix) then
                # No Sigma provided, coefficient Matrix is fourth argument.
                if RowDimension(_passed[4]) <> 2 or ColumnDimension(_passed[4]) <> P:-r then
                    error "The coefficient matrix must be a (2 x r)-Matrix. Here, r = %1", P:-r;
                end if;
                self:-A := _passed[4];
            elif type(_passed[4], set(set(integer))) then
                # Fan Sigma provided
                self:-Sigma := _passed[4];
            elif type(_passed[4], {Vector, list(integer)}) then
                # Weight w provided. Compute the associated fan.
                w := _passed[4];

                if not containsrelint(getMovingCone(P), w) then
                    error "The given weight w = %1 does not lie in the moving cone of P.", convert(w, list);
                end if;

                # TODO: Make this computation more efficient by incrementally searching only the minimal cones containing w.
                candidates := [op(combinat[powerset]({seq(1 .. P:-n + P:-m)}) minus {{}})];
                minimalBunchCones := {};
                for i from 1 to nops(candidates) do
                    cone := candidates[i];
                    # If we already have a cone that's contained in this one, skip it.
                    if select(c -> c subset cone, minimalBunchCones) = {} then
                        if containsrelint(poshull(Column(getQ0(P), [op(cone)])), w) then
                            minimalBunchCones := {op(minimalBunchCones), cone};
                        end if;
                    end if;
                end do;
                # We dualize to get the maximal cones of the associated fan.
                self:-Sigma := map(c -> {seq(1 .. P:-n + P:-m)} minus c, minimalBunchCones);
            else
                error "Expected second argument to be of type set(set(integer), Vector, list or Matrix";
            end if;

            if _npassed > 4 then
                if type(_passed[5], Matrix) then
                    if RowDimension(_passed[5]) <> 2 or ColumnDimension(_passed[5]) <> P:-r then
                        error "The coefficient matrix must be a (2 x r)-Matrix. Here, r = %1", P:-r;
                    end if;
                    self:-A := _passed[5];
                else
                    error "Expected third argument to be of type Matrix";
                end if;
            end if;
        end if;

        # If no coefficient matrix has been provided, use the standard one.
        if type(self:-A, undefined) then
            self:-A := Matrix(2, P:-r, [[1, 0, -1 $ P:-r - 2], [0, 1, seq(-i, i = 1 .. P:-r - 2)]]);
        end if;

        # Construct the trinomial relations in the Cox Ring
        setRelations(self, P, self:-A);

        # If no Sigma has been provided, check if we are in one of the allowed cases (1)-(3)
        # and compute it.
        if type(self:-Sigma, undefined) then
            if P:-s = 1 then
                # Here, X is a C*-surface.
                # In this case, the moving cone equals the semiample cone, hence there is only one
                # possible fan for the P-Matrix. We can write it down explicitly, following Construction 5.4.1.6 of "Cox Rings".
                # Note that we do not require the P-Matrix to be slope-ordered.
                
                ordered_indices := [seq(sort([seq(1 .. P:-ns[i])], 
                    (j1, j2) -> P:-slopes[i, j1] > P:-slopes[i, j2]), 
                    i = 1 .. P:-r)];
                
                sigma_plus := {seq(doubleToSingleIndex(P:-format, i, ordered_indices[i,1]), i = 1 .. P:-r)};
                sigma_minus := {seq(doubleToSingleIndex(P:-format, i, ordered_indices[i,P:-ns[i]]) , i = 1 .. P:-r)};
                
                taus := {seq(seq({
                    doubleToSingleIndex(P:-format, i, ordered_indices[i,j]), 
                    doubleToSingleIndex(P:-format, i, ordered_indices[i,j+1])}, 
                    j = 1 .. P:-ns[i] - 1), i = 1 .. P:-r)};
                
                if P:-case = "EE" then
                    self:-Sigma := {sigma_plus} union taus union {sigma_minus};
                elif P:-case = "PE" then
                    taus_plus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,1]), P:-n + P:-m} , i = 1 .. P:-r)};
                    self:-Sigma := taus_plus union taus union {sigma_minus};
                elif P:-case = "EP" then
                    taus_minus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,P:-ns[i]]), P:-n + P:-m} , i = 1 .. P:-r)};
                    self:-Sigma := {sigma_plus} union taus union taus_minus;
                elif P:-case = "PP+" then
                    taus_plus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,1]), P:-n + P:-m - 1} , i = 1 .. P:-r)};
                    taus_minus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,P:-ns[i]]), P:-n + P:-m} , i = 1 .. P:-r)};
                    self:-Sigma := taus_plus union taus union taus_minus;
                elif P:-case = "PP-" then
                    taus_plus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,1]), P:-n + P:-m} , i = 1 .. P:-r)};
                    taus_minus := {seq({doubleToSingleIndex(P:-format, i, ordered_indices[i,P:-ns[i]]), P:-n + P:-m - 1} , i = 1 .. P:-r)};
                    self:-Sigma := taus_plus union taus union taus_minus;
                else
                    error "Invalid case. P:-case should be one of the five strings: \"EE\", \"PE\", \"EP\", \"PP+\", \"PP-\"";
                end if;

            elif P:-picardNumber = 1 then
                # In this case, the picard number is one, hence there is only one complete fan in the lattice 
                # containing the columns of the P-Matrix as rays.
                numColumns := ColumnDimension(P:-mat);
                self:-Sigma := {seq({seq(1 .. numColumns)} minus {i}, i = 1 .. numColumns)};
            elif admitsFano(P) then
                # In this case, we call the procedure again with the anticanonical class as the weight.
                return TVarOne[ModuleCopy](self, proto, P, getAnticanClass(P));
            else
                error "This PMatrix is neither of Picard number one, nor is it a surface, nor does it admit a Fano variety."
                      "Therefore, you must provide a fan Sigma or a weight w as input.";
            end if;
        end if;
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
    end proc;

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
    end proc;

    export setAmpleCone :: static := proc(self :: TVarOne, ampleCone :: CONE) self:-ampleCone := ampleCone; end proc;

    export getAmpleCone :: static := proc(self :: TVarOne)
        local cone;
        if type(self:-ampleCone, undefined) or 'forceCompute' in [_passed] then
            setAmpleCone(self, intersection(seq(poshull(Column(getQ0(self:-P), [op({seq(1 .. ColumnDimension(getQ0(self:-P)))} minus cone)])), cone in getMaximalXCones(self))));
        end if;
        return self:-ampleCone;
    end proc;

    export setIsFanoVal :: static := proc(self :: TVarOne, isFanoVal :: boolean) self:-isFanoVal := isFanoVal; end proc;

    export isFano :: static := proc(self :: TVarOne)
        if type(self:-isFanoVal, undefined) or 'forceCompute' in [_passed] then
            setIsFanoVal(self, containsrelint(getAmpleCone(self), getAnticanClass(self:-P)));
        end if;
        return self:-isFanoVal
    end proc;

    export setintersectionTable :: static := proc(self :: TVarOne, intersectionTable :: Array) self:-intersectionTable := intersectionTable; end proc;

    export getIntersectionTable :: static := proc(X :: TVarOne)   
        local P, res, ordered_indices, i, mcal, newM, j, j_, j1, j2, k1, k2, k, i1, i2, kplus, kminus;

        if type(X:-intersectionTable, undefined) or 'forceCompute' in [_passed] then
            P := X:-P;
            if X:-P:-s > 1 then
                error "Intersection numbers are only defined for K*-surfaces, i.e. s = 1.";
            end if;

            # We use the formulas of Chapter 6 in the paper: "Log del pezzo surfaces with torus action" by Hättig, Hausen, Hummel
            # Note however that we do not require P to be slope-ordered, hence we have to generalize the formulas
            # a bit by translating the indices through the `ordered_indices` array.

            res := Array(1 .. ColumnDimension(P:-mat), 1 .. ColumnDimension(P:-mat), fill = 0);
            ordered_indices := [seq(sort([seq(1 .. P:-ns[i])], 
                        (j1, j2) -> P:-slopes[i, j1] > P:-slopes[i, j2]), 
                        i = 1 .. P:-r)];
            
            # First, compute calligraphic m. This is encoded as a list of arrays, since we want to have
            # the indexing go from 0 to n_i.
            mcal := [];
            for i from 1 to P:-r do
                newM := Array(0 .. P:-ns[i], fill = 0);
                if P:-case = "EE" or P:-case = "EP" then
                    # There is an elliptic fixed point x^+
                    newM[0] := - 1 / P:-mplus;
                end if;
                if P:-case = "EE" or P:-case = "PE" then
                    # There is an elliptic fixed point x^-
                    newM[P:-ns[i]] := 1 / P:-mminus;
                end if;

                for j from 1 to P:-ns[i] - 1 do
                    newM[j] := 1 / (P:-slopes[i,ordered_indices[i,j]] - P:-slopes[i,ordered_indices[i,j+1]]);
                end do;

                mcal := [op(mcal), newM];
            end do;

            # Compute intersection numbers of two adjacent rays in the leaves.
            # This is independent of the case of P
            for i from 1 to P:-r do
                for j_ from 1 to P:-ns[i] - 1 do
                    j1 := ordered_indices[i,j_];
                    j2 := ordered_indices[i,j_+1];
                    k1 := doubleToSingleIndex(P:-format, i, j1);
                    k2 := doubleToSingleIndex(P:-format, i, j2);
                    res[k1,k2] := 1 / (P:-lss[i,j1] * P:-lss[i,j2]) * mcal[i][j_];
                    res[k2,k1] := res[k1,k2];
                end do;
            end do:

            # Compute the self intersection numbers of rays in the leaves
            # This is independent of the case of P
            for i from 1 to P:-r do
                for j_ from 1 to P:-ns[i] do
                    j := ordered_indices[i,j_];
                    k := doubleToSingleIndex(P:-format, i, j);
                    res[k,k] := - 1 / P:-lss[i][j]^2 * (mcal[i][j_ - 1] + mcal[i][j_]);
                end do;
            end do:

            if P:-case = "EE" or P:-case = "EP" then
                ##########################################  
                ## There is an elliptic fixed point x^+ ##
                ##########################################

                # Intersection numbers of the highest rays in each block with each other.
                for i1 from 1 to P:-r do
                    for i2 in {seq(1 .. P:-r)} minus {i1} do
                        j1 := ordered_indices[i1,1];
                        j2 := ordered_indices[i2,1];
                        k1 := doubleToSingleIndex(P:-format, i1, j1);
                        k2 := doubleToSingleIndex(P:-format, i2, j2);
                        if P:-ns[i1] = 1 and P:-ns[i2] = 1 then
                            res[k1,k2] := - 1 / (P:-lss[i1,j1] * P:-lss[i2,j2]) * (mcal[i1][0] + mcal[i1][1]);
                        else 
                            res[k1,k2] := - 1 / (P:-lss[i1,j1] * P:-lss[i2,j2]) * mcal[i1][0];
                        end if;
                    end do;
                end do;
            else
                ################################################
                ## There is a parabolic fixed point curce D^+ ##
                ################################################

                # kplus is the index of the divisor D^+. This depends on the case.
                if P:-case = "PE" or P:-case = "PP+" then
                    kplus := doubleToSingleIndex(P:-format, -1, 1);
                elif P:-case = "PP-" then
                    kplus := doubleToSingleIndex(P:-format, -1, 2);
                end if;

                # Intersection numbers of the highest rays of each block with D^+
                for i from 1 to P:-r do
                    j := ordered_indices[i,1];
                    k := doubleToSingleIndex(P:-format, i, j);
                    res[k,kplus] := 1 / P:-lss[i][j];
                    res[kplus,k] := res[k,kplus];
                end do;

                # Self intersection of D^+
                res[kplus, kplus] := - P:-mplus;

            end if;

            if P:-case = "EE" or P:-case = "PE" then
                ##########################################  
                ## There is an elliptic fixed point x^- ##
                ##########################################

                # Intersection numbers of the lowest rays in each block with each other.
                for i1 from 1 to P:-r do
                    for i2 in {seq(1 .. P:-r)} minus {i1} do
                        j1 := ordered_indices[i1,P:-ns[i1]];
                        j2 := ordered_indices[i2,P:-ns[i2]];
                        k1 := doubleToSingleIndex(P:-format, i1, j1);
                        k2 := doubleToSingleIndex(P:-format, i2, j2);
                        if P:-ns[i1] = 1 and P:-ns[i2] = 1 then
                            res[k1,k2] := - 1 / (P:-lss[i1,j1] * P:-lss[i2,j2]) * (mcal[i1][0] + mcal[i1][1]);
                        else 
                            res[k1,k2] := - 1 / (P:-lss[i1,j1] * P:-lss[i2,j2]) * mcal[i1][P:-ns[i1]];
                        end if;
                    end do;
                end do;
            else 
                ################################################
                ## There is a parabolic fixed point curce D^- ##
                ################################################

                # kminus is the index of the divisor D^+. This depends on the case.
                if P:-case = "EP" or P:-case = "PP-" then
                    kminus := doubleToSingleIndex(P:-format, -1, 1);
                elif P:-case = "PP+" then
                    kminus := doubleToSingleIndex(P:-format, -1, 2);
                end if;

                # Intersection numbers of the highest rays of each block with D^-
                for i from 1 to P:-r do
                    j := ordered_indices[i,P:-ns[i]];
                    k := doubleToSingleIndex(P:-format, i, j);
                    res[k,kminus] := 1 / P:-lss[i][j];
                    res[kminus,k] := res[k,kminus];
                end do;

                # Self intersection of D^-
                res[kminus, kminus] := P:-mminus;

            end if;

            setintersectionTable(X, res);
        
        end if;

        return X:-intersectionTable;

    end proc;

    (*
    Compute the intersection number of any two divisors, given as linear combinations
    of the primitive invariant divisors D_X^{ij} and D_X^{\pm}.
    The input data are lists of integers encoding the coefficients this linear combination.
    *)
    export intersectionNumber :: static := proc(X :: TVarOne, D1 :: list(integer), D2 :: list(integer))
        local k1, k2;
        
        if nops(D1) <> X:-P:-n + X:-P:-m or nops(D2) <> X:-P:-n + X:-P:-m then
            error "The list of integers encoding the divisor must have length n + m = %1", X:-P:-n + X:-P:-m;
        end if;

        return add([seq(seq(D1[k1] * D2[k2] * getIntersectionTable(X)[k1,k2], k2 = 1 .. nops(D2)), k1 = 1 .. nops(D1))]);

    end proc;

    export setAnticanonicalSelfIntersection :: static := proc(self :: TVarOne, anticanonicalSelfIntersection) self:-anticanonicalSelfIntersection := anticanonicalSelfIntersection; end proc;

    export getAnticanonicalSelfIntersection :: static := proc(X :: TVarOne)
        if type(X:-anticanonicalSelfIntersection, undefined) or 'forceCompute' in [_passed] then
            setAnticanonicalSelfIntersection(X, intersectionNumber(X, getAnticanCoefficients(X:-P), getAnticanCoefficients(X:-P)));
        end if;
        return X:-anticanonicalSelfIntersection;
    end proc;


    (****************************
    *** ADMISSIBLE OPERATIONS ***
    *****************************)

    export removeSingleRedundantBlock :: static := proc(X :: TVarOne, i0 :: integer)
        local newP, oldToNewIndex, newSigma;
        newP := PMatrix[removeSingleRedundantBlock](X:-P, i0);
        oldToNewIndex := k -> if k < add(X:-P:-ns[1 .. i0 - 1]) + 1 then k else k - 1 end if;
        newSigma := map(cones -> map(oldToNewIndex, cones), X:-Sigma);
        return TVarOne(newP, newSigma);
    end proc;

    export removeRedundantBlocks :: static := proc(X :: TVarOne)
        local P, redundantIndices;
        P := X:-P;
        redundantIndices := select(i -> P:-lss[i] = [1], [seq(1 .. nops(P:-lss))]);
        # If there is a redundant block and we still have more than two blocks, remove it.
        if nops(redundantIndices) > 0 and P:-r > 2 then
            # Recursive call
            return removeRedundantBlocks(removeSingleRedundantBlock(X, redundantIndices[1]));
        else
            return X;
        end if;
    end proc;

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
    export applyAdmissibleRowOperation := proc(X :: TVarOne, C :: Matrix, T :: Matrix)
        return TVarOne(PMatrix[applyAdmissibleRowOperation](X:-P, C, T), X:-Sigma);
    end proc;

    (*
    See also `sortColumnsByLss` for PMatrix.
    *)
    export sortColumnsByLss := proc(X0 :: TVarOne)
        local X, newP, sigma_, taus_, bundledPerm, newSigma, newX;
        # First remove redundant columns. Note that this also changes the fan.
        X := removeRedundantBlocks(X0);
        # Now, sort the columns using the method from `PMatrix`.
        newP, sigma_, taus_ := op(PMatrix[sortColumnsByLss](X:-P, output = ['normalized', 'sigma', 'taus']));
        # Adjust the fan accordingly, by applying the permutation given from the sorting.
        bundledPerm := bundleColumnPermutation(X:-P:-format, sigma_, taus_, Perm([]));
        newSigma := map(cones -> map(k -> bundledPerm[k], cones), X:-Sigma);
        newX := TVarOne(newP, newSigma);
        return newX;
    end proc;

    export ModulePrint :: static := proc(self :: TVarOne)
        nprintf(cat("TVarOne(dim = ", self:-P:-s + 1,
          ", lss = ", self:-P:-lss,
          ", Sigma = ", self:-Sigma));
    end;

    export TVarOneInfo :: static := proc(self :: TVarOne)
        local P, i, relations, maximalXCones, Q, classGroup, picardNumber, anticanClass, effectiveConeRays, movingConeRays, ampleConeRays, isFano, gorensteinIndex, intersectionTable, anticanonicalSelfIntersection;
        print(P = self:-P:-mat);
        print([seq(cat(n,i), i = 0 .. self:-P:-r - 1), m] = [seq(self:-P:-ns[i], i = 1 .. self:-P:-r), self:-P:-m]);
        print(relations = self:-relations);
        print(maximalXCones = getMaximalXCones(self));
        print(Q = getQ(self:-P));
        print(classGroup = getClassGroup(self:-P));
        print(picardNumber = self:-P:-picardNumber);
        print(anticanClass = getAnticanClass(self:-P));
        print(effectiveConeRays = rays(getEffectiveCone(self:-P)));
        print(movingConeRays = rays(getMovingCone(self:-P)));
        print(ampleConeRays = rays(getAmpleCone(self)));
        print(isFano = TVarOne[isFano](self));
        print(gorensteinIndex = getGorensteinIndex(self));
        if self:-P:-s = 1 then
            print(intersectionTable = getIntersectionTable(self));
            print(anticanonicalSelfIntersection = getAnticanonicalSelfIntersection(self));
        end if;
    end;

end module:

ImportTVarOneList := proc(stmt)
    local columns, INDEX_S, INDEX_P, INDEX_DIMENSION, INDEX_PICARDNUMBER, INDEX_CLASSGROUP, INDEX_DEGREEMATRIX, INDEX_ANTICANCLASS, INDEX_AMBIENTFAN, INDEX_MAXIMALXCONES, INDEX_GORENSTEININDEX, INDEX_ISGORENSTEIN;
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
        elif clm = "classGroup" then
            INDEX_CLASSGROUP := i;
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
        P := PMatrix(Fetch(stmt, INDEX_S), Matrix(parse(Fetch(stmt, INDEX_P))), 'skipChecks');
        # Read in any already avaliable metadata about P, if available
        if type(INDEX_PICARDNUMBER, integer) then
            setPicardNumber(P, Fetch(stmt, INDEX_PICARDNUMBER));
            if type(INDEX_CLASSGROUP, integer) then
                setClassGroup(P, parse(Fetch(stmt, INDEX_CLASSGROUP)));
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
It returns the list of rowids in the database where the P-Matrix is equivalent to the given one. 
If the P-Matrix does not occur, it returns the empty list.
*)
FindInDatabase := proc(connection, tableName :: string, X :: TVarOne)
    local P, stmtString, stmt, Ps, M, rowids, i, resids;
    P := X:-P;
    stmtString := cat("SELECT rowid,s,P FROM ", tableName, " WHERE ",
        "m = ", P:-m, " AND ",
        "s = ", P:-s, " AND ",
        "classGroup = \"", getClassGroup(X:-P), "\" AND ",
        "gorensteinIndex = ", getGorensteinIndex(X));
    # Use the 'lss' as search criterion, only if X is non-toric.
    if not isToric(P) then
        stmtString := cat(stmtString, " AND orderedLss = \"", sortColumnsByLss(P):-lss, "\"");
    end if;
    
    stmt := Prepare(connection, stmtString);
    Ps := ImportTVarOneList(stmt);
    M := FetchAll(stmt):
    rowids := [seq(M[i,1], i = 1 .. RowDimension(M))];
    resids := [];
    for i from 1 to nops(rowids) do
        if areEquivalent(P, Ps[i]) then
            resids := [op(resids), rowids[i]];
        end if;
    end do;
    return resids;
end proc;

(*
Imports a single TVarOne from a database with a given rowid.
To import multiple varieties using an arbitrary SQLite query, use `ImportTVarOneList`
*)
ImportTVarOne := proc(db, tableName :: string, rowid :: integer)
    ImportTVarOneList(Prepare(db, cat("SELECT * FROM ", tableName, " WHERE rowid = ", rowid)))[1];
end proc;

(*
Given a list of varieties of complexity one `Xs` and a SQLite database connection `db`, this function
inserts those varieties from `Xs` into the database that are not already present.
*)
ExportTVarOneList := proc(connection, tableName :: string, Xs :: list(TVarOne))
    local k, X, P, stmt, i, knownCount;
    knownCount := 0;

    for k from 1 to nops(Xs) do
        X := Xs[k];
        # Normalize the given T-Variety.
        # DISABLED FOR NOW
        # X := normalizeTVarOne(Xs[k]);
        
        # Only insert `X` if it is not already present.
        # Can be disabled by the 'noChecks' option
        if 'noChecks' in [_passed] or FindInDatabase(connection, tableName, X) = [] then
            P := X:-P;
            stmt := Prepare(connection, cat("INSERT INTO ", tableName ," VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"));
            Bind(stmt, 1, P:-r, valuetype = "integer");
            Bind(stmt, 2, convert(P:-ns, string), valuetype = "text");
            Bind(stmt, 3, P:-n, valuetype = "integer");
            Bind(stmt, 4, P:-m, valuetype = "integer");
            Bind(stmt, 5, P:-s, valuetype = "integer");
            Bind(stmt, 6, convert(P:-lss, string), valuetype = "text");
            Bind(stmt, 7, convert([seq(convert(Row(P:-mat, i), list), i = 1 .. RowDimension(P:-mat))], string), valuetype = "text");
            Bind(stmt, 8, P:-dim, valuetype = "integer");
            Bind(stmt, 9, P:-picardNumber, valuetype = "integer");
            Bind(stmt, 10, convert(getClassGroup(P), string), valuetype = "text");
            Bind(stmt, 11, convert([seq(convert(Row(getQ(P), i), list), i = 1 .. RowDimension(getQ(P)))], string), valuetype = "text");
            Bind(stmt, 12, convert(convert(getAnticanClass(P), list), string), valuetype = "text");
            Bind(stmt, 13, convert(X:-Sigma, string), valuetype = "text");
            Bind(stmt, 14, convert(getMaximalXCones(X), string), valuetype = "text");
            Bind(stmt, 15, getGorensteinIndex(X), valuetype = "integer");
            Bind(stmt, 16, isGorenstein(X), valuetype = "integer");
            Bind(stmt, 17, convert(sortColumnsByLss(P):-lss, string), valuetype = "text");
            Bind(stmt, 18, convert(rays(getEffectiveCone(P)), string), valuetype = "text");
            Bind(stmt, 19, convert(rays(getMovingCone(P)), string), valuetype = "text");
            Bind(stmt, 20, convert(rays(getAmpleCone(X)), string), valuetype = "text");
            Bind(stmt, 21, isFano(X), valuetype = "integer");
            Bind(stmt, 22, convert(X:-variables, string), valuetype = "text");
            Bind(stmt, 23, convert(X:-monomials, string), valuetype = "text");
            Bind(stmt, 24, convert(X:-relations, string), valuetype = "text");
            Bind(stmt, 25, convert(convert(map(x -> [numer(x), denom(x)], getIntersectionTable(X)), list, nested), string), valuetype = "text");
            Bind(stmt, 26, convert([numer(getAnticanonicalSelfIntersection(X)), denom(getAnticanonicalSelfIntersection(X))], string), valuetype = "text");
            Bind(stmt, 27, convert(getAnticanonicalSelfIntersection(X), hfloat), valuetype = "float");
            Bind(stmt, 28, isToric(P), valuetype = "integer");
            Bind(stmt, 29, isIrredundant(P), valuetype = "integer");


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

    if 'logging' in [_passed] then
        if 'noChecks' in [_passed] then
            print(cat("Written ", nops(Xs), "varieties to the database"));
        else
            print(cat("Supplied ", nops(Xs), " varieties. Already known: ", knownCount, ". Inserted ", nops(Xs) - knownCount, " new varieties into the database."));
        end if;
    end if;
    
end proc;

(*
Helper function to perform an operation on every variety in the database, without loading them into memory all at once
(which can be unpractical, when the database is large).
The parameter `f` is a function taking a `TVarOne` as first argument and an integer as second argument, which is the
position the variety occurs in the database.
*)
performOnDatabase := proc(db, tableName :: string, f, step := 1000)
    local numberOfEntries, numberOfSteps, i, offset, Xs, j;
    # First, count the number of database entries
    numberOfEntries := FetchAll(Prepare(db, cat("SELECT COUNT(*) FROM ", tableName)))[1,1];
    numberOfSteps := ceil(numberOfEntries / step);
    for i from 1 to numberOfSteps do 
        offset := (i - 1) * step;
        Xs := ImportTVarOneList(Prepare(db, cat("SELECT * FROM ", tableName, " LIMIT ", step, " OFFSET ", offset)));
        for j from 1 to nops(Xs) do
            f(Xs[j], offset + j);
        end do;
    end do;
end proc;

end module:
