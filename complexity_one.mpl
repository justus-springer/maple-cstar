ComplexityOne := module()

option package;

export PFormat, PMatrix, TVarOne;

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

    export ModuleApply :: static := proc()
        Object(PMatrix, _passed);
    end;

    (*
    This method is for creates a PMatrix from various kinds of data. It supports
    four different input methods:
    (1) format :: PFormat, lss :: list(list(integer)), d :: Matrix.
    (2) lss :: list(list(integer)), d :: Matrix.
    (3) r :: integer, P :: Matrix (NOT YET IMPLEMMENTED)
    (4) P :: Matrix (NOT YET IMPLEMMENTED)

    In input method (1) and (2), `lss` is the list of exponent vectors of relations
    of the Cox Ring. These make up the first `r-1` rows of the PMatrix. The lower
    `s` rows are given in the matrix `d`. In method (2), the format of the matrix
    is inferred from the given data.

    In methods (3) and (4), the P-Matrix is directly specified and all other data
    is inferred.

    TODO: Explain why method (4) is bad.

    *)
    export ModuleCopy :: static := proc(self :: PMatrix, proto :: PMatrix)

        local ls, l, rows, i, j, P;
        if _npassed = 2 then error "not enough arguments." end if;

        if type(_passed[3], 'PFormat') then
            # Input method (1)
            # format :: PFormat, lss :: list(list(integer)), d :: Matrix.

            if _npassed < 5 then
                error "Not enough arguments. Expected input: PFormat, list(list(integer)), Matrix.";
            end if:

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
            self:-Q := Matrix(AGHdata(A)[3]);
            self:-picardNumber := AGdata(AGHdata(A)[2])[3];
            self:-Q0 := DeleteRow(self:-Q, self:-picardNumber + 1 .. RowDimension(self:-Q));
            self:-Q;
        end if;
    end;

    export getQ0 :: static := proc(self :: PMatrix)
        if not type(self:-Q0, undefined) then
            self:-Q0
        else
            getQ(self);
            getQ0(self);
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
        local Q0, as, i;
        if not type(self:-anticanClass, undefined) then
            self:-anticanClass;
        else
            Q0 := getQ0(self);
            as := getAnticanCoefficients(self);
            self:-anticanClass := add(seq(as[i] * Column(Q0, i), i = 1 .. self:-n + self:-m));
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
        as := getAnticanCoefficients(P);
        us := [seq(u[i], i = 1 .. RowDimension(P:-mat))];
        sol := isolve({seq(DotProduct(us, Column(P:-mat, j)) = as[j], j in cone)});
        return evalb(sol <> NULL);
    end;

    export ModulePrint :: static := proc(self :: PMatrix)
        nprintf(cat("PMatrix(", self:-lss, ", m = ", self:-m, ", s = ", self:-s, ")"));
    end;

    export PMatrixInfo :: static := proc(self :: PMatrix)
        local P, Q, picardNumber, anticanClass;
        print(P = self:-mat);
        nprintf(cat("Format: ", self:-ns, ", m = ", self:-m, ", s = ", self:-s)); #"
        print(Q = getQ(self));
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

end module:
