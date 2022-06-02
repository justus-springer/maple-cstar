
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

    # The maximum resp. minimum slopes in each block
    export maximumSlopes := undefined;
    export minimumSlopes := undefined;

    # The sum of the maximum resp. minimum slopes
    export mplus := undefined;
    export mminus := undefined;

    export mplusFloor := undefined;
    export mminusCeil := undefined;

    export betasPlus := undefined;
    export betasMinus := undefined;

    # The orientation of the P-Matrix. It is either +1, -1 or 0.
    export orientation := undefined;

    export magicInvariant := undefined;

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
        self:-maximumSlopes := map(max, self:-slopes);
        self:-minimumSlopes := map(min, self:-slopes);
        self:-mplus := add(self:-maximumSlopes);
        self:-mminus := add(self:-minimumSlopes);
        self:-mplusFloor := add(map(floor, self:-maximumSlopes));
        self:-mminusCeil := add(map(ceil, self:-minimumSlopes));
        self:-betasPlus := [seq([seq(self:-slopes[i,j] - floor(self:-maximumSlopes[i]), j = 1 .. self:-ns[i])], i = 1 .. self:-r)];
        self:-betasMinus := [seq([seq(self:-slopes[i,j] - ceil(self:-minimumSlopes[i]), j = 1 .. self:-ns[i])], i = 1 .. self:-r)];
        self:-magicInvariant := [sortLex([self:-mplus, -self:-mminus]), sortLex([self:-betasPlus, -self:-betasMinus])];
        
        if self:-m = 0 then
            self:-case := "EE";
            if self:-mplus > - self:-mminus then
                self:-orientation := 1;
            elif self:-mplus < -self:-mminus then
                self:-orientation := -1;
            else    
                if sortLexComparison(sortLex(self:-betasPlus), sortLex(-self:-betasMinus)) then
                    self:-orientation := 1;
                elif sortLexComparison(sortLex(-self:-betasMinus), sortLex(self:-betasPlus)) then
                    self:-orientation := -1;
                else
                    self:-orientation := 0;
                end if;
            end if;
        elif self:-m = 1 then
            if d[1, doubleToSingleIndex(self:-format, -1, 1)] = 1 then
                self:-case := "PE";
                self:-orientation := 1;
            elif d[1, doubleToSingleIndex(self:-format, -1, 1)] = -1 then
                self:-case := "EP";
                self:-orientation := -1;
            end if;
        elif self:-m = 2 then
            if d[1, doubleToSingleIndex(self:-format, -1, 1)] = 1 then
                self:-case := "PP+";
                self:-orientation := 1;
            elif d[1, doubleToSingleIndex(self:-format, -1, 1)] = -1 then
                self:-case := "PP-";
                self:-orientation := -1;
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
    (5) mplusFloor :: integer, betasPlus :: list(list(fraction)), case :: string

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

    TODO: Explain method (5).

    *)
    export ModuleCopy :: static := proc(self :: PMatrix, proto :: PMatrix)

        local lss, ls, l, format, d, rows, rows0, i, j, P, r, s, ns, numZerosBefore, mplusFloor, betasPlus, case;
        
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

                # If this is a P-Matrix of a surface, set the slopes
                if self:-s = 1 then 
                    setSurfaceData(self, self:-d, self:-lss);
                end if;
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

        elif type(_passed[3], integer) and type(_passed[4], Matrix) then
            # Input method (4)
            # s :: integer, P :: Matrix

            s := _passed[3];
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
        elif type(_passed[3], integer) and type(_passed[4], list) and type(_passed[5], string) then
            # Input method (5)
            # mplusFloor :: integer, betasPlus :: list(list(fraction)), case :: string

            mplusFloor := _passed[3];
            betasPlus := _passed[4];
            case := _passed[5];

            # TODO: Proper input checking

            ns := map(nops, betasPlus);

            self:-lss := map(betas -> map(beta -> denom(beta), betas), betasPlus);
            self:-d := Matrix(1, add(ns),
                [[seq(self:-lss[1,j] * (betasPlus[1,j] + mplusFloor), j = 1 .. ns[1]),
                  seq(seq(self:-lss[i,j] * betasPlus[i,j], j = 1 .. ns[i]), i = 2 .. nops(ns))]]);

            if case = "PE" then
                self:-d := <self:-d | 1>;
                setFormat(self, PFormat(ns, 1, 1));
            elif case = "EP" then
                self:-d := <self:-d | -1>;
                setFormat(self, PFormat(ns, 1, 1));
            elif case = "PP+" then
                self:-d := <self:-d | 1 | -1>;
                setFormat(self, PFormat(ns, 2, 1));
            elif case = "PP-" then
                self:-d := <self:-d | -1 | 1>;
                setFormat(self, PFormat(ns, 2, 1));
            elif case = "EE" then
                setFormat(self, PFormat(ns, 0, 1));
            else
                error "The case must be one of the five strings: \"EE\", \"EP\", \"PE\", \"PP+\" and \"PP-\".";
            end if; 

            P := PMatrix(self:-format, self:-lss, self:-d);
            self:-mat := P:-mat;
            self:-P0 := P:-P0;
            setSurfaceData(self, self:-d, self:-lss);

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

    export applyAdmissibleOperation :: static := proc(P :: PMatrix, a :: AdmissibleOperation)
        if P:-format <> a:-formatFrom then
            error "Format of the P-Matrix does not coincide with format of the admissible operation.";
        end if;
        return PMatrix(a:-formatTo, a:-S . P:-mat . a:-bundledPermutationMatrix^(-1));
    end proc;

    (******************
    ** NORMALIZATION **
    *******************)

    (*
    Removes a single redundant column from a P-Matrix.
    *)
    export removeSingleRedundantBlock :: static := proc(P_ :: PMatrix, i0 :: integer)
        local P, C, T, a, newLss, newD, newFormat, i;
        P := P_;
        # First, we have to apply admissible row operations to achieve all zeros in the d-block under
        # the redundant columns. We construct the A-Matrix necessary for this.
        # This step is necessary to ensure the columns of the resulting P-Matrix still generate the whole
        # space as a cone.
        if i0 = 1 then
            C := Matrix(P:-s, P:-r - 1, (k,l) -> if l = 1 then P:-d[k, add(P:-ns[1 .. i0-1]) + 1] else 0 end if);
        else
            C := Matrix(P:-s, P:-r - 1, (k,l) -> if l = i0 - 1 then -P:-d[k, add(P:-ns[1 .. i0-1]) + 1] else 0 end if);
        end if;
        T := Matrix(P:-s, P:-s, shape = diagonal, 1);
        a := AdmissibleOperation[FromRowOperation](P:-format, C, T);
        P := applyAdmissibleOperation(P,a);
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
    Computes the admissible operation necessary to sort the columns of a P-Matrix by the entries in the L-block.
    The columns of block are sorted descendingly. The blocks themselves are sorted first descendingly
    according to size, and within the same size lexicographically by the lss.
    *)
    export sortColumnsByLssOperation :: static := proc(P :: PMatrix)
        local taus, newLss, i, compfun, tau, sigma;

        # Sort each individual block
        taus := [];
        for i from 1 to P:-r do
            tau := Perm(sort([seq(1 .. P:-ns[i])], (j1, j2) -> P:-lss[i][j1] > P:-lss[i][j2], 'output' = 'permutation'))^(-1);
            taus := [op(taus), tau];
        end do;
        newLss := map(i -> applyPermToList(taus[i], P:-lss[i]), [seq(1 .. P:-r)]);

        # This is the comparison function we will use to sort the blocks themselves.
        # It returns true if the block of index `i1` should occur before the block of index `i2`.
        compfun := proc(i1 :: integer, i2 :: integer)
            local j;
            if nops(newLss[i1]) < nops(newLss[i2]) then return false;
            elif nops(newLss[i1]) > nops(newLss[i2]) then return true;
            else
                for j from 1 to nops(newLss[i1]) do
                    if newLss[i1][j] < newLss[i2][j] then return false;
                    elif newLss[i1][j] > newLss[i2][j] then return true;
                    end if;
                end do;
                return true;
            end if;
        end proc;
        sigma := Perm(sort([seq(1 .. P:-r)], (i1, i2) -> compfun(i1, i2), 'output' = 'permutation'))^(-1);
        
        return AdmissibleOperation[FromPermutations](P:-format, sigma, taus, Perm([]));

    end proc;

    export sortColumnsByLss :: static := proc(P :: PMatrix)
        applyAdmissibleOperation(P, sortColumnsByLssOperation(P));
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
        local resBool, resOperation, identityMatrix, zeroMatrix, C, T, S, newP, sol, resultList, str, i, j;

        identityMatrix := Matrix(P1:-r - 1, P1:-r - 1, shape = diagonal, 1);
        zeroMatrix := Matrix(P1:-r - 1, P1:-s, fill = 0);
        C := Matrix(P1:-s, P1:-r - 1, symbol = 'x');
        T := Matrix(P1:-s, P1:-s, symbol = 'y');
        S := <<identityMatrix | zeroMatrix>, <C | T>>;
        newP := S . P1:-mat;
        sol := isolve({seq(seq(newP[i,j] = P2:-mat[i,j] , j = 1 .. ColumnDimension(P1:-mat)), i = P1:-r .. RowDimension(P1:-mat))});
        
        if sol <> NULL and abs(Determinant(subs(sol, T))) = 1 then
            return AdmissibleOperation[FromRowOperation](P1:-format, subs(sol, C), subs(sol, T));             
        end if;
        
        return false;

    end proc;

    (*
    Computes all admissible column operations leaving the L-block of a P-Matrix invariant.
    This is an important step for checking if two P-Matrices are equivalent by admissible operations.
    *)
    export invariantAdmissibleOperations :: static := proc(P0 :: PMatrix)
    
        local a0, P, admOps, i, taus, rhos;

        a0 := sortColumnsByLssOperation(P0);
        P := sortColumnsByLss(P0);

        admOps := map(sigma -> AdmissibleOperation[FromSigma](P:-format, sigma), invariantPermutations(P:-lss));

        for i from 1 to P:-r do
            taus := map(tau -> AdmissibleOperation[FromSingleTau](P:-format, i, tau), invariantPermutations(P:-lss[i]));
            admOps := map(a -> op(map(tau -> compose(a, tau), taus)), admOps);
        end do;

        rhos := map(rho -> AdmissibleOperation[FromRho](P:-format, Perm(rho)), combinat[permute](P:-m));
        admOps := map(a -> op(map(rho -> compose(a, rho), rhos)), admOps);
        
        admOps := map(a -> compose(compose(a0, a), inverse(a0)), admOps);

        return admOps;

    end proc;

    (*
    Computes all admissible operations turning `P1_` into `P2_`. If `P1_` and `P2_` are not equivalent
    by admissible operations, this returns the empty list.
    *)
    export areEquivalentOperations :: static := proc(P1_ :: PMatrix, P2_ :: PMatrix)
        local Ps, P, i, P1, P2, P11, P12, a0, admOps, resultOps, a, rowOp;

        # First, remove any redundant blocks.
        P1 := removeRedundantBlocks(P1_);
        P2 := removeRedundantBlocks(P2_);

        # After removing redundant blocks, the number of blocks must coincide.
        if P1:-r <> P2:-r then
            return [];
        end if;

        # Varieties of different dimension can't be isomorphic.
        if P1:-s <> P2:-s then
            return [];
        end if; 
        
        # Composing the sorting opreation of P1 with the inverse sorting operation of P2, we
        # should obtain an admissible operation sending the L-block of P1 to the L-block of P2.
        a0 := compose(sortColumnsByLssOperation(P1), inverse(sortColumnsByLssOperation(P2)));
        
        # If the L-block of P1 does not coincide with the L-block of P2 after sorting, then
        # the P-Matrices can't be equivalent.
        if applyAdmissibleOperation(P1, a0):-lss <> P2:-lss then
            return [];
        end if;

        # By composing `a0` with all invariant operations of P1, we obtain *all* admissible operations
        # sending the L-block of P1 to the L-block of P2.
        admOps := map(a -> compose(a, a0), invariantAdmissibleOperations(P1));

        # For each of these operations, we check if there is an admissible row operation turning the matrix into P2.
        resultOps := [];
        for a in admOps do
            rowOp := areRowEquivalent(applyAdmissibleOperation(P1, a), P2);
            if type(rowOp, AdmissibleOperation) then
                resultOps := [op(resultOps), compose(a, rowOp)];
            end if;
        end do;

        return resultOps;

    end proc;

    export areEquivalent :: static := proc(P1 :: PMatrix, P2 :: PMatrix)
        if P1:-s = 1 and P2:-s = 1 then
            # If we are in the surface case, the magic invariant is all we need
            return evalb(P1:-magicInvariant = P2:-magicInvariant);
        end if;
        # In general, we need to to the hard work of determining all admissible operations
        # turning P1 into P2.
        return areEquivalentOperations(P1, P2) <> [];
    end proc;

    (*
    Computes the normal form for a P-Matrix of a surface.
    *)
    export normalForm :: static := proc(P :: PMatrix)
        if P:-orientation = -1 then
            PMatrix(- P:-mminusCeil, sortLex(- P:-betasMinus), P:-case);
        else
            PMatrix(P:-mplusFloor, sortLex(P:-betasPlus), P:-case);
        end if;
    end proc;

    (*
    Converts a cone given as a list of column indices into a CONE from the convex package.
    *)
    export intSetConeToConvexCone :: static := proc(P :: PMatrix, cone :: set(integer))
        return poshull(op(map(i -> Column(P:-mat, i), cone)));
    end proc;

    (*
    Converts a CONE from the convex package into a list of column indices for a P-Matrix.
    *)
    export convexConeToIntSetCone :: static := proc(P :: PMatrix, cone :: CONE)
        local cols;
        cols := map(v -> :-convert(v, list), [Column(P:-mat, [seq(1 .. P:-n + P:-m)])]);
        return map(ray -> ListTools[Search](ray, cols), {op(rays(cone))});
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