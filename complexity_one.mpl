ComplexityOne := module()

option package;

export PFormat, PMatrix, isBigCone, isLeafCone, isXCone, isGorensteinForXCone, TVarOne;

## TODO: Remove dependency on MDSpackage.
uses LinearAlgebra, MDSpackage;

module PFormat()
    option object;

    export r, ns, n, m, s;

    export ModuleApply :: static := proc()
        Object(PFormat, _passed);
    end;

    export ModuleCopy :: static := proc(self :: PFormat, proto :: PFormat,
        ns_ :: list(integer), m_ :: integer, s_ :: integer, $)
        local i;
        for i in ns_ do
            if i < 1 then error "each element of ns must be at least 1." end if:
        end do:
        if m_ < 0 then error "m must be at least 0." end if:
        if s_ < 1 then error "s must be at least 1." end if:
        self:-r := nops(ns_);
        self:-ns := ns_;
        self:-n := add(ns_);
        self:-m := m_;
        self:-s := s_;
    end;

    export ModulePrint :: static := proc(self :: PFormat)
        nprintf(cat("(ns = ", self:-ns, ", m = ", self:-m, ", s = ", self:-s, ")")); # "
    end;

end module:

(*
Checks whether a given cone `cone` is big with respect to a given PFormat.
Here, a cone is called big if for every i = 0,...,r-1, we have {N, ..., N+ns[i+1]} ∩ cone ≠ {},
where N = ns[1] + ... + ns[i]
*)
isBigCone := proc(format :: PFormat, cone :: set(integer))
    local r, n, ns, N, i, k:
    r := format:-r;
    n := format:-n;
    ns := format:-ns;
    for i from 0 to r - 1 do
        N := add(ns[1..i]);
        if evalb({seq(N + k, k = 1 .. ns[i+1])} intersect cone = { }) then
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
isLeafCone := proc(format :: PFormat, cone :: set(integer))
    local r, n, ns, m, N, i, k:
    r := format:-r;
    n := format:-n;
    ns := format:-ns;
    m := format:-m;
    for i from 0 to r - 1 do
        N := add(ns[1..i]);
        if cone subset ({seq(N + k, k = 1 .. ns[i+1])} union {seq(n + k, k = 1 .. m)}) then
            return true:
        end if:
    end do:
    return false:
end proc:

(*
Checks whether a given `cone` is an X-cone with respect to a given P-matrix format `nL, m`.
Here, a cone is called an X-cone, if it is either a leaf cone or a big cone.
*)
isXCone := proc(format :: PFormat, cone :: set(integer))
    isBigCone(format, cone) or isLeafCone(format, cone):
end:

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

    export ModuleCopy :: static := proc(self :: PMatrix, proto :: PMatrix)

        local rows, i, j, P;
        if _npassed = 2 then error "not enough arguments." end if;

        # Input method (1)
        # format, ls, d
        if type(_passed[3], 'PFormat') then
            if _npassed < 5 then
                error "Not enough arguments. Expected input: PFormat, list(list(integer)), Matrix, Matrix.";
            end if:

            if not type(_passed[3], 'PFormat') then
                error "Expected 2nd argument to be of type: PFormat.";
            end if;
            setFormat(self, _passed[3]);

            if not type(_passed[4], list(list(integer))) then
                error "Expected 3rd argument to be of type: list(list(integer))";
            end if;
            self:-lss := _passed[4];

            if not type(_passed[5], 'Matrix'(self:-s, self:-n + self:-m, integer)) then
                error "Expected 4th argument to be of type: Matrix(%1, %2, integer)", self:-s, self:-n + self:-m;
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

        else
            error "Wrong arguments.";
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

    export ModulePrint :: static := proc(self :: PMatrix)
        ## TODO: Make this work as expected.
        print(self:-mat);
        print(self:-format);
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

(*
isGorensteinForXCone := proc(P :: PMatrix, cone :: set(integer))
    local as, u, i, j;
    as := getAnticanCoefficients(P);
    us := [seq(u[i], i = 1 .. nops(cone))];
    isolve({seq(DotProduct(us, Column(P:-mat, j)) = as[j], j in cone)});
end;
*)

module TVarOne()
    option object;

    export P, Sigma;

    local SigmaX := undefined;
    local SigmaXMax := undefined;

    export ModuleApply :: static := proc()
        Object(TVarOne, _passed);
    end;

    export ModuleCopy :: static := proc(self :: TVarOne, proto :: TVarOne,
        P :: PMatrix, Sigma :: list(set(integer)))
        self:-P := P;
        self:-Sigma := Sigma;
    end;

    export getSigmaX :: static := proc(self :: TVarOne)
        if not type(self:-SigmaX, undefined) then
            self:-SigmaX;
        else
            select(cone -> isXCone(self:-P:-format, cone), self:-Sigma);
        end if;
    end;

    export isMaximalXCone :: static := proc(self :: TVarOne, cone :: set(integer))
        local SigmaX, cone_;
        SigmaX := getSigmaX(self);
        for cone_ in SigmaX do
            if cone subset cone_ and cone <> cone_ then
                return false;
            end if;
        end do;
        return true;
    end;

    export getSigmaXMax :: static := proc(self :: TVarOne)
        select(cone -> isMaximalXCone(self, cone), getSigmaX(self));
    end;

    export ModulePrint :: static := proc(self :: TVarOne)
        ## TODO: Make this work as expected.
        print(self:-P);
        print(self:-Sigma);
    end;

end module:

end module:
