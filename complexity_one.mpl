ComplexityOne := module()

option package;

export PFormat, PMatrix, TVarOne, parseCSV;

## TODO: Remove dependency on MDSpackage.
uses LinearAlgebra, MDSpackage;

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
    Checks whether a given cone `cone` is big with respect to a given PFormat.
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
    Checks whether a given `cone` is an X-cone with respect to a given P-matrix format `nL, m`.
    Here, a cone is called an X-cone, if it is either a leaf cone or a big cone.
    *)
    export isXCone :: static := proc(self :: PFormat, cone :: set(integer))
        isBigCone(self, cone) or isLeafCone(self, cone):
    end:

    export ModulePrint :: static := proc(self :: PFormat)
        nprintf(cat("PFormat(", self:-ns, ", m = ", self:-m, ", s = ", self:-s, ")"));
    end;

end module:

module PMatrix()
    option object;

    # These fields are guaranteed to be filled when a PMatrix is created.
    export format, lss, d, mat;
    export r, ns, n, m, s;

    # These fields are only computed when needed. Use these getters below for these.
    local picardNumber := undefined;
    local Q := undefined;
    local Q0 := undefined;
    local classGroup := undefined;
    local anticanClass := undefined;
    local anitcanCoefficients := undefined;

    export setFormat :: static := proc(self :: PMatrix, f :: PFormat)
        self:-format := f;
        self:-r := f:-r;
        self:-ns := f:-ns;
        self:-n := f:-n;
        self:-m := f:-m;
        self:-s := f:-s;
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
    (4) r :: integer, P :: Matrix
    (5) P :: Matrix (NOT YET IMPLEMMENTED)

    In input method (1) and (2), `lss` is the list of exponent vectors of relations
    of the Cox Ring. These make up the first `r-1` rows of the PMatrix. The lower
    `s` rows have to be provided the matrix `d`. In method (2), the format of the matrix
    is inferred from the given data.

    In methods (3), (4) and (5), the P-Matrix is directly specified. It will be checked
    if this matrix fulfills the properties of a P-Matrix, i.e. if it has the correct shape,
    the columns are primitive and they generate the whole space as a cone. In method (3), the
    expected P-Format is specified and will be checked against P. In method (4), only the number of
    exponent vectors `r` is specified, all other data is inferred from `P`. In method (5), only the
    matrix is provided. Note that this comes with some ambiguity: If we don't know how many
    exponent vectors there are, hence we don't know where the upper block of `P` ends and the d-block 
    starts. In this case, we can only make an educated guess based on the entries of `P`.

    TODO: Explain more about input method (5).

    *)
    export ModuleCopy :: static := proc(self :: PMatrix, proto :: PMatrix)

        local lss, ls, l, rows, i, j, P, r, ns, numZerosBefore;
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
                P := PMatrix(_passed[3]:-r, _passed[4]);
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
                            seq(self:-lss[i])), (0 $ add(self:-ns[i..self:-r-1]) + self:-m)], i = 2 .. self:-r),
                        seq(:-convert(Row(self:-d, j), list), j = 1 .. self:-s)];
                self:-mat := Matrix(rows);

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
            # r :: integer, P :: Matrix
        
            if _npassed < 4 then
                error "Not enough arguments. Expected input: integer, Matrix";
            end if;

            r := _passed[3];
            if r < 2 then
                error "r must be at least 2."
            end if;

            if not type(_passed[4], Matrix) then
                error "Expected 2nd argument to be of type: Matrix";
            end if;

            P := _passed[4];
            self:-mat := P;

            for i in P do
                if not type(i, integer) then
                    error "All entries of P must be of type: integer";
                end if;
            end do;

            if RowDimension(P) < r then
                error "P must have at least r = %1 rows", r;
            end if;

            # Get the first vector l1 by looking at the second row of P
            ls := [];
            i := 1;
            while P[2,i] <> 0 do
                if P[2,i] > -1 then
                    error " Given matrix is not in P-shape. Expected: P[2,%1] < 0. Given: P[2,%1] = %2", i, P[2,i];
                end if;

                ls := [op(ls), -P[2,i]];
                i := i + 1;

                # Needs to be checked before we re-enter the loop
                if i > ColumnDimension(P) then
                    error "Given matrix is not in P-shape. Expected: P[2,%1] = 0. Given: P[2,%1] = %2", i-1, P[2,i-1];
                end if;
            end do;

            lss := [ls];
            ns := [nops(ls), 0 $ r-1]; # Here, we preliminarily fill up with zeros.

            # Check that the first r-1 rows of P all start with -l1
            for i from 1 to r-1 do
                for j from 1 to nops(ls) do
                    if P[i,j] <> P[2,j] then
                        error "Given matrix is not in P-shape. Expected: P[%1,%2] = P[2,%2]. Given: P[%1,%2] = %3 and P[2,%2] = %4", i, j, P[i,j], P[2,j];
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

    export getQ :: static := proc(self :: PMatrix)
        local A;
        # Currently, we rely on MDSpackage for these computations.
        # But they should be reimplemented here eventually.
        if not type(self:-Q, undefined) then
            self:-Q
        else
            A := AGHP2Q(convert(self, matrix));
            self:-classGroup := AGHdata(A)[2];
            self:-Q := Matrix(AGHdata(A)[3]);
            self:-picardNumber := AGdata(self:-classGroup)[3];
            self:-Q0 := DeleteRow(self:-Q, self:-picardNumber + 1 .. RowDimension(self:-Q));
            self:-Q;
        end if;
    end;

    export getQ0 :: static := proc(self :: PMatrix)
        if not type(self:-Q0, undefined) then
            self:-Q0;
        else
            getQ(self);
            getQ0(self);
        end if;
    end;

    export getClassGroup :: static := proc(self :: PMatrix)
        if not type(self:-classGroup, undefined) then
            self:-classGroup;
        else
            getQ(self);
            getClassGroup(self);
        end if;
    end;

    export getPicardNumber :: static := proc(self :: PMatrix)
        if not type(self:-picardNumber, undefined) then
            self:-picardNumber
        else
            getQ(self);
            getPicardNumber(self);
        end if;
    end;

    export getAnticanCoefficients :: static := proc(self :: PMatrix)
        if not type(self:-anitcanCoefficients, undefined) then
            self:-anitcanCoefficients;
        else
            self:-anitcanCoefficients :=
                [1 $ self:-n + self:-m] - (self:-r - 2) * [op(self:-lss[1]), 0 $ self:-n + self:-m - self:-ns[1]]
        end if;
    end;

    export getAnticanClass :: static := proc(self :: PMatrix)
        local as, i, anticanVec, d;
        if not type(self:-anticanClass, undefined) then
            self:-anticanClass;
        else
            as := getAnticanCoefficients(self);
            anticanVec := add(seq(as[i] * Column(getQ(self), i), i = 1 .. self:-n + self:-m));
            # Some entries in `anticanVec` live in cyclic groups Z/dZ.
            # We normalize these entries, so that each of them is less than `d`.
            for i from getPicardNumber(self) + 1 to RowDimension(getQ(self)) do
                d := AGdata(getClassGroup(self))[4][i - getPicardNumber(self)];
                anticanVec[i] := anticanVec[i] mod d;
            end do;
            return anticanVec;
        end if;
    end;

    (*
    Checks whether the variety of a given PMatrix satisfies the Gorenstein condition
    with respect to a given X-cone. The Gorenstein condition is satisfied if
    the system of equations {<u,v_i> = a_i, i in cone} has an integer solution
    vector u. Here, v_i is the ith column of P and the a_i are the coefficents of the
    anticanonical class.
    *)
    export isGorensteinForXCone :: static := proc(self :: PMatrix, cone :: set(integer))
        local as, u, us, i, j, sol;
        as := getAnticanCoefficients(self);
        us := [seq(u[i], i = 1 .. RowDimension(self:-mat))];
        sol := isolve({seq(DotProduct(us, Column(self:-mat, j)) = as[j], j in cone)});
        return evalb(sol <> NULL);
    end;

    export ModulePrint :: static := proc(self :: PMatrix)
        nprintf(cat("PMatrix(", self:-lss, ", m = ", self:-m, ", s = ", self:-s, ")"));
    end;

    export PMatrixInfo :: static := proc(self :: PMatrix)
        local P, Q, picardNumber, classGroup, anticanClass;
        print(P = self:-mat);
        nprintf(cat("Format: ", self:-ns, ", m = ", self:-m, ", s = ", self:-s)); #"
        print(Q = getQ(self));
        print(classGroup = getClassGroup(self));
        print(picardNumber = getPicardNumber(self));
        print(anticanClass = getAnticanClass(self));
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

    export ModuleApply :: static := proc()
        Object(TVarOne, _passed);
    end;

    export ModuleCopy :: static := proc(self :: TVarOne, proto :: TVarOne,
        P :: PMatrix, Sigma :: list(set(integer)))
        self:-P := P;
        self:-Sigma := Sigma;
    end;

    export getXCones :: static := proc(self :: TVarOne)
        if not type(self:-XCones, undefined) then
            self:-XCones;
        else
            select(cone -> isXCone(self:-P:-format, cone), self:-Sigma);
        end if;
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

    export getMaximalXCones :: static := proc(self :: TVarOne)
        select(cone -> isMaximalXCone(self, cone), getXCones(self));
    end;

    export isGorenstein :: static := proc(self :: TVarOne)
      andmap(cone -> isGorensteinForXCone(self:-P, cone), getMaximalXCones(self));
    end;

    export ModulePrint :: static := proc(self :: TVarOne)
        nprintf(cat("TVarOne(dim = ", self:-P:-s + 1,
          ", picardNum = ", getPicardNumber(self:-P),
          ", w = ", convert(getAnticanClass(self:-P), list), ")")); # "
    end;

    export TVarOneInfo :: static := proc(self :: TVarOne)
        local maximalXCones, isGorenstein;
        PMatrixInfo(self:-P);
        print(maximalXCones = getMaximalXCones(self));
        print(isGorenstein = :-isGorenstein(self));
    end;

end module:


parseCSV := proc(fn :: string)
    local CSVData, index_row, i, INDEX_P, INDEX_FORMAT, INDEX_Q, INDEX_ANTICAN, row, Ps, P;
    local evaluated_Q, expected_Q, evaluated_antican, expected_antican;

    CSVData := ImportMatrix(fn, delimiter = ";");
    index_row := CSVData[1];

    # First, find out index numbers of the fields by looking at the first row.
    for i from 1 to Dimension(index_row) do
        if index_row[i] = "P-Format" or index_row[i] = "Format" then
            INDEX_FORMAT := i;
        elif index_row[i] = "Generator matrix" or index_row[i] = "P" then
            INDEX_P := i;
        elif index_row[i] = "Degree matrix" or index_row[i] = "Q" then
            INDEX_Q := i;
        elif index_row[i] = "Anticanonical class" or index_row[i] = "Antican class" then
            INDEX_ANTICAN := i;
        end if;
    end do;

    if not type(INDEX_FORMAT, integer) then
        error "Not enough data. You must provide a column for the format of the P-Matrix. Please call it \"P-Format\" or \"Format\"";
    end if;

    (*
    if not type(INDEX_P, integer) and not (type(INDEX_LSS, integer) or type(INDEX_LSS, integer)) then
        error "Not enough data. You must either provide a column with the P-Matrix (called \"P\" or \"Generator Matrix\"), or a column for the exponent vectors (called \"Exponent vectors\") and a column with the d-Block (called \"d\" or \"d-Block\")";
    end if;
    *)

    if not type(INDEX_FORMAT, integer) then
        error "Not enough data. You must provide a column for the P-Matrix. Please call it \"P\" or \"Generator Matrix\"";
    end if;


    Ps := [];
    for i from 2 to RowDimension(CSVData) do
        row := CSVData[i];
        try
            format := PFormat(parse(row[INDEX_FORMAT]));
            P := PMatrix(format, Matrix(parse(row[INDEX_P])));
        catch:
            error "Error while reading line %1: %2", i, StringTools[FormatMessage](lastexception[2 .. -1]);
        end try;
        Ps := [op(Ps), P];

        if not 'nochecks' in { _passed } then
            # Test if the degree matrix is as expected
            expected_Q := parse(row[INDEX_Q]):
            evaluated_Q := map(x -> convert(x, list), [Row(getQ(P), [seq(1 .. RowDimension(getQ(P)))])]):
            if expected_Q <> evaluated_Q then
                error "In %-1 row: Degree matrix check failed. Evaluated: %2. Given: %3.", i, evaluated_Q, expected_Q;
            end if;

            # Test if the anticanonical class is as expected
            expected_antican := parse(row[INDEX_ANTICAN]):
            evaluated_antican := convert(getAnticanClass(P), list);
            if expected_antican <> evaluated_antican then
                error "In %-1 row: Antican check failed. Evaluated: %2. Given: %3.", i, evaluated_antican, expected_antican;
            end if;
        end if;

    end do;

    return Ps;

end proc;

end module:
