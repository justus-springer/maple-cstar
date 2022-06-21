
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

        if nops(taus) <> formatFrom:-r + 1 then
            error "Expected `taus` to be a list of length r+1 = %1", formatFrom:-r+1;
        end if;

        if not type(C, 'Matrix'(formatFrom:-s, formatFrom:-r, integer)) then
            error "Expected C to be of type Matrix(%1,%2,integer)", formatFrom:-s, formatFrom:-r;
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

        if nops(ds) <> formatFrom:-r + 1 then
            error "Expected `ds` to be a list of length r+1 = %1 with nonzero complex entries.", formatFrom:-r+1;
        end if;

        for i from 0 to formatFrom:-r do
            if ds[i+1] = 0 then
                error "entries of `ds` cannot be zero.";
            end if;
        end do;

        # Construct missing data
        self:-formatTo := PFormat(applyPermToList(sigma, convert(formatFrom:-ns, list)), formatFrom:-m, formatFrom:-s);
        bundledPermApply := proc(k :: integer)
            local i, j;
            i, j := singleToDoubleIndex(formatFrom, k);
            if i = -1 then
                return doubleToSingleIndex(self:-formatTo, -1, rho[j]);
            end if;
            # Some annoying shifting necessary because permutations start counting at one, while we count blocks from zero
            return doubleToSingleIndex(self:-formatTo, sigma[i+1]-1, taus[i+1][j]);
        end proc;
        self:-bundledPerm := Perm(map(bundledPermApply, [seq(1 .. formatFrom:-numCols)]));

        self:-bundledPermutationMatrix := Matrix(formatFrom:-numCols, formatFrom:-numCols, 
            (i,j) -> if i = self:-bundledPerm[j] then 1 else 0 end if);
    
        self:-B := Matrix(formatFrom:-r, formatFrom:-r,
            (i,j) -> if sigma[j+1] = 1 then -1 
                     elif sigma[j+1] = i+1 then 1
                     else 0 end if);
        
        zeroMatrix := Matrix(formatFrom:-r, formatFrom:-s, fill = 0);
        self:-S := <<self:-B | zeroMatrix>, <C | T>>;

        self:-sigmaPermutationMatrix := Matrix(formatFrom:-r + 1, formatFrom:-r + 1, 
            (i,j) -> if i = sigma[j] then 1 else 0 end if);

        self:-D := Matrix(formatFrom:-r + 1, formatFrom:-r + 1, Vector(ds), shape = diagonal);
        
    end proc;

    ##############################
    ## Alternative constructors ##
    ##############################

    export Identity :: static := proc(format :: PFormat)
        AdmissibleOperation(format, Perm([]), [Perm([]) $ format:-r + 1], Perm([]), Matrix(format:-s, format:-r, fill = 0), Matrix(format:-s, format:-s, shape = diagonal, 1), Matrix([[1,0],[0,1]]), [1 $ format:-r + 1]);
    end proc;

    export OnP :: static := proc(formatFrom :: PFormat, sigma :: Perm, taus :: list(Perm), rho :: Perm, C :: Matrix, T :: Matrix)
        AdmissibleOperation(formatFrom, sigma, taus, rho, C, T, Matrix([[1,0],[0,1]]), [1 $ formatFrom:-r + 1]);
    end proc;

    export FromPermutations :: static := proc(formatFrom :: PFormat, sigma :: Perm, taus :: list(Perm), rho :: Perm)
        AdmissibleOperation[OnP](formatFrom, sigma, taus, rho, Matrix(formatFrom:-s, formatFrom:-r, fill = 0), Matrix(formatFrom:-s, formatFrom:-s, shape = diagonal, 1));
    end proc;

    export FromSigma :: static := proc(formatFrom, sigma :: Perm)
        AdmissibleOperation[FromPermutations](formatFrom, sigma, [Perm([]) $ formatFrom:-r + 1], Perm([]));
    end proc;

    export FromSingleBlockSwap :: static := proc(formatFrom :: PFormat, i1 :: integer, i2 :: integer)
        AdmissibleOperation[FromSigma](formatFrom, Perm([[i1+1, i2+1]]));
    end proc;

    export FromTaus :: static := proc(formatFrom, taus :: list(Perm))
        AdmissibleOperation[FromPermutations](formatFrom, Perm([]), taus, Perm([]));
    end proc;

    export FromSingleTau :: static := proc(formatFrom, i :: integer, tau :: Perm)
        AdmissibleOperation[FromTaus](formatFrom, [Perm([]) $ i, tau, Perm([]) $ formatFrom:-r + 1 - i]);
    end proc;

    export FromRho :: static := proc(formatFrom, rho :: Perm)
        AdmissibleOperation[FromPermutations](formatFrom, Perm([]), [Perm([]) $ formatFrom:-r + 1], rho);
    end proc;

    export FromSingleColumnSwap :: static := proc(formatFrom, i :: integer, j1 :: integer, j2 :: integer)
        AdmissibleOperation[FromSingleTau](formatFrom, i, Perm([[j1, j2]]));
    end proc;

    export FromRowOperation :: static := proc(formatFrom :: PFormat, C :: Matrix, T :: Matrix)
        AdmissibleOperation[OnP](formatFrom, Perm([]), [Perm([]) $ formatFrom:-r + 1], Perm([]), C, T);
    end proc;

    export OnA :: static := proc(formatFrom, U :: Matrix, ds :: list(complex))
        AdmissibleOperation(formatFrom, Perm([]), [Perm([]) $ formatFrom:-r + 1], Perm([]), Matrix(formatFrom:-s, formatFrom:-r, fill = 0), Matrix(formatFrom:-s, formatFrom:-s, shape = diagonal, 1), U, ds);
    end proc;

    export FromU :: static := proc(formatFrom, U :: Matrix)
        AdmissibleOperation[OnA](formatFrom, U, [1 $ formatFrom:-r + 1]);
    end proc;

    export FromScaling :: static := proc(formatFrom, ds :: list(complex))
        AdmissibleOperation[OnA](formatFrom, Matrix([[1,0],[0,1]]), ds);
    end proc;

    export FromSingleScaling :: static := proc(formatFrom, i :: integer, d :: complex)
        AdmissibleOperation[FromScaling](formatFrom, [1 $ i, d, 1 $ formatFrom:-r + 1 - i]);
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
        taus := [seq(Perm:-perm_mul(a1:-taus[i], a2:-taus[a1:-sigma[i]]) , i = 1 .. a2:-formatFrom:-r + 1)];
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

    export `=`::static := proc( l, r, $ )
        if (_npassed <> 2 or not type(l, AdmissibleOperation) or not type(r, AdmissibleOperation)) then
           return false;
        end;
        return 
            l:-formatFrom = r:-formatFrom and
            l:-sigma = r:-sigma and 
            l:-taus = r:-taus and 
            l:-rho = r:-rho and
            Equal(l:-C, r:-C) and
            Equal(l:-T, r:-T) and
            Equal(l:-U, r:-U) and
            l:-ds = r:-ds;
    end;

    export ModulePrint :: static := proc(a :: AdmissibleOperation)
        local str, isFirstEntry;
        str := "AdmissibleOperation(";

        if a:-sigma <> Perm([]) then
            str := cat(str, "sigma = ", convert(a:-sigma, string), ", ");
        end if;
        if a:-taus <> [Perm([]) $ a:-formatFrom:-r + 1] then
            tausString := cat("[", convert(a:-taus[1], string));
            for i from 2 to a:-formatFrom:-r + 1 do
                tausString := cat(tausString, ", ", convert(a:-taus[i], string));
            end do;
            tausString := cat(tausString, "]");
            str := cat(str, "taus = ", tausString, ", ");
        end if;
        if a:-rho <> Perm([]) then
            str := cat(str, "rho = ", convert(a:-rho, string), ", ");
        end if;
        if not Equal(a:-C, Matrix(a:-formatFrom:-s, a:-formatFrom:-r, fill = 0)) then
            str := cat(str, "C = ", convert(convert(a:-C, list, nested), string), ", ");
        end if;
        if not Equal(a:-T, Matrix(a:-formatFrom:-s, a:-formatFrom:-s, shape = diagonal, 1)) then
            str := cat(str, "T = ", convert(convert(a:-T, list, nested), string), ", ");
        end if;
        if not Equal(a:-U, Matrix([[1,0],[0,1]])) then
            str := cat(str, "U = ", convert(convert(a:-U, list, nested), string), ", ");
        end if;
        if a:-ds <> [1 $ a:-formatFrom:-r + 1] then
            str := cat(str, "ds = ", convert(a:-ds, string), ", ");
        end if;

        if a = AdmissibleOperation[Identity](a:-formatFrom) then
            str := cat(str, "1");
        else
            # Cut the last comma
            str := str[1 .. length(str) - 2];
        end if;
        str := cat(str, ")");

        nprintf(str);
    end;



end module;