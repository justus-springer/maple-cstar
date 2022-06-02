
module AdmissibleOperation()
    option object;

    export formatFrom, sigma, taus, rho, C, T, U, ds;
    export formatTo, bundledPerm, bundledPermutationMatrix, sigmaPermutationMatrix, B, S, D;

    export ModuleApply :: static := proc()
        Object(AdmissibleOperation, _passed);
    end proc;

    export ModuleCopy :: static := proc(self :: AdmissibleOperation, proto :: AdmissibleOperation,
        formatFrom :: PFormat, sigma :: Perm, taus :: list(Perm), rho :: Perm, C :: Matrix, T :: Matrix, U :: Matrix, ds :: list(complex), $)

        local bundledPermApply, zeroMatrix, i;

        # Copy given data
        self:-formatFrom := formatFrom;
        self:-sigma := sigma;
        self:-taus := taus;
        self:-rho := rho;
        self:-C := C;
        self:-T := T;
        self:-U := U;
        self:-ds := ds;

        # Input checking

        if nops(taus) <> formatFrom:-r then
            error "Expected `taus` to be a list of length r = %1", formatFrom:-r;
        end if;

        if not type(C, 'Matrix'(formatFrom:-s, formatFrom:-r - 1, integer)) then
            error "Expected C to be of type Matrix(%1,%2,integer)", formatFrom:-s, formatFrom:-r - 1;
        end if;

        if not type(T, 'Matrix'(formatFrom:-s, formatFrom:-s, integer)) then
            error "Expected T to be of type Matrix(%1,%1,integer)", formatFrom:-s;
        end if;

        if not abs(Determinant(T)) = 1 then
            error "T must be a unimodular matrix, i.e. its determinant must be 1 or -1";
        end if;

        if not type(U, 'Matrix'(2, 2, complex)) then
            error "U must be a 2x2 matrix.";
        end if;

        if Determinant(U) = 0 then
            error "U must be invertible.";
        end if;

        if nops(ds) <> formatFrom:-r then
            error "Expected `ds` to be a list of length r = %1 with nonzero complex entries.", formatFrom:-r;
        end if;

        for i from 1 to formatFrom:-r do
            if ds[i] = 0 then
                error "entries of `ds` cannot be zero.";
            end if;
        end do;

        # Construct missing data
        self:-formatTo := PFormat(applyPermToList(sigma, formatFrom:-ns), formatFrom:-m, formatFrom:-s);
        bundledPermApply := proc(k :: integer)
            local i, j;
            i, j := singleToDoubleIndex(formatFrom, k);
            if i = -1 then
                return doubleToSingleIndex(self:-formatTo, -1, rho[j]);
            end if;
            return doubleToSingleIndex(self:-formatTo, sigma[i], taus[i][j]);
        end proc;
        self:-bundledPerm := Perm(map(bundledPermApply, [seq(1 .. formatFrom:-n + formatFrom:-m)]));

        self:-bundledPermutationMatrix := Matrix(formatFrom:-n + formatFrom:-m, formatFrom:-n + formatFrom:-m, 
            (i,j) -> if i = self:-bundledPerm[j] then 1 else 0 end if);
    
        self:-B := Matrix(formatFrom:-r - 1, formatFrom:-r - 1,
            (i,j) -> if sigma[j+1] = 1 then -1 
                     elif sigma[j+1] = i+1 then 1
                     else 0 end if);
        
        zeroMatrix := Matrix(formatFrom:-r - 1, formatFrom:-s, fill = 0);
        self:-S := <<self:-B | zeroMatrix>, <C | T>>;

        self:-sigmaPermutationMatrix := Matrix(formatFrom:-r, formatFrom:-r, 
            (i,j) -> if i = sigma[j] then 1 else 0 end if);

        self:-D := Matrix(formatFrom:-r, formatFrom:-r, Vector(ds), shape = diagonal);
        
    end proc;

    ##############################
    ## Alternative constructors ##
    ##############################

    export Identity :: static := proc(format :: PFormat)
        AdmissibleOperation(format, Perm([]), [Perm([]) $ format:-r], Perm([]), Matrix(format:-s, format:-r - 1, fill = 0), Matrix(format:-s, format:-s, shape = diagonal, 1), Matrix([[1,0],[0,1]]), [1 $ format:-r]);
    end proc;

    export POperation :: static := proc(formatFrom :: PFormat, sigma :: Perm, taus :: list(Perm), rho :: Perm, C :: Matrix, T :: Matrix)
        AdmissibleOperation(formatFrom, sigma, taus, rho, C, T, Matrix([[1,0],[0,1]]), [1 $ formatFrom:-r]);
    end proc;

    export ColumnOperation :: static := proc(formatFrom :: PFormat, sigma :: Perm, taus :: list(Perm), rho :: Perm)
        AdmissibleOperation[POperation](formatFrom, sigma, taus, rho, Matrix(formatFrom:-s, formatFrom:-r - 1, fill = 0), Matrix(formatFrom:-s, formatFrom:-s, shape = diagonal, 1));
    end proc;

    export BlockPermutation :: static := proc(formatFrom, sigma :: Perm)
        AdmissibleOperation[ColumnOperation](formatFrom, sigma, [Perm([]) $ formatFrom:-r], Perm([]));
    end proc;

    export SingleBlockSwap :: static := proc(formatFrom :: PFormat, i1 :: integer, i2 :: integer)
        AdmissibleOperation[BlockPermutation](formatFrom, Perm([[i1, i2]]));
    end proc;

    export SingleBlockPermutation :: static := proc(formatFrom, i :: integer, tau :: Perm)
        AdmissibleOperation[ColumnOperation](formatFrom, Perm([]), [Perm([]) $ i - 1, tau, Perm([]) $ P:-r - i], Perm([]));
    end proc;

    export LastColumnsPermutation :: static := proc(formatFrom, rho :: Perm)
        AdmissibleOperation[ColumnOperation](formatFrom, Perm([]), [Perm([]) $ formatFrom:-r], rho);
    end proc;

    export SingleColumnSwap :: static := proc(formatFrom)

    end proc;

    export RowOperation :: static := proc(formatFrom :: PFormat, C :: Matrix, T :: Matrix)
        AdmissibleOperation[POperation](formatFrom, Perm([]), [Perm([]) $ formatFrom:-r], Perm([]), C, T);
    end proc;

    export AOperation :: static := proc(formatFrom, U :: Matrix, ds :: list(complex))
        AdmissibleOperation(formatFrom, Perm([]), [Perm([]) $ formatFrom:-r], Perm([]), Matrix(formatFrom:-s, formatFrom:-r - 1, fill = 0), Matrix(formatFrom:-s, formatFrom:-s, shape = diagonal, 1), U, ds);
    end proc;

    export UOperation :: static := proc(formatFrom, U :: Matrix)
        AdmissibleOperation[AOperation](formatFrom, U, [1 $ formatFrom:-r]);
    end proc;

    export ScaleOperation :: static := proc(formatFrom, ds :: list(complex))
        AdmissibleOperation[AOperation](formatFrom, Matrix([[1,0],[0,1]]), ds);
    end proc;

    export SingleScaleOperation :: static := proc(formatFrom, i :: integer, d :: complex)
        AdmissibleOperation[ScaleOperation](formatFrom, [1 $ i - 1, d, 1 $ formatFrom:-r - i]);
    end proc;

    (*
    Composition of two admissible operations. Here, `a1` is applied before `a2`.
    *)
    export compose :: static := proc(a1 :: AdmissibleOperation, a2 :: AdmissibleOperation)
        local sigma, taus, rho, C, T, U, ds, i;
        if a1:-formatTo <> a2:-formatFrom then
            error "Can't compose these admissible operations, their formats do not match up.";
        end if;
        
        sigma := Perm:-perm_mul(a1:-sigma, a2:-sigma);
        taus := [seq(Perm:-perm_mul(a1:-taus[i], a2:-taus[a1:-sigma[i]]) , i = 1 .. a2:-formatFrom:-r)];
        rho := Perm:-perm_mul(a1:-rho, a2:-rho);
        C := a2:-C . a1:-B + a2:-T . a1:-C;
        T := a2:-T . a1:-T;
        U := a2:-U . a1:-U;
        ds := zip((d1, d2) -> d1 * d2, a1:-ds, a2:-ds);
        AdmissibleOperation(a1:-formatFrom, sigma, taus, rho, C, T, U, ds);
    end proc;

    export inverse :: static := proc(a :: AdmissibleOperation)
        local sigma, taus, rho, C, T, U, ds;

        sigma := a:-sigma^(-1);
        taus := map(tau -> tau^(-1), applyPermToList(a:-sigma, a:-taus));
        rho := a:-rho^(-1);
        T := a:-T^(-1);
        C := - a:-T^(-1) . a:-C . a:-B^(-1);
        U := a:-U^(-1);
        ds := map(d -> 1 / d, a:-ds);

        AdmissibleOperation(a:-formatTo, sigma, taus, rho, C, T, U, ds);

    end proc;

    export ModulePrint :: static := proc(a :: AdmissibleOperation)
        nprintf(cat("AdmissibleOperation(sigma = ", a:-sigma, ", taus = ", a:-taus, ", rho = ", a:-rho, ", C = ", a:-C, ", T = ", a:-T, ")"));
    end;

end module;