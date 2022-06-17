
module ComplexityOneVariety()
    option object;

    # These fields are guaranteed to be filled when a ComplexityOneVariety is created.
    export P;
    export Sigma := undefined;
    export A := undefined;

    # Variables and relations in the Cox Ring
    export variables := undefined;
    export monomials := undefined;
    export relations := undefined;

    ###################################################################################
    ### These fields are only computed when needed. Use the getters below for them. ###
    ###################################################################################

    # The maximal X-cones of the fan Sigma, i.e. the maximal cones of the fan of the minimal
    # ambient toric variety
    local maximalXCones := undefined;

    # The fan of the minimal ambient toric variety, stored as a `FAN` from the convex package
    local minimalAmbientFan := undefined;

    # Says whether the variety is Q-factorial, i.e. if every weil divisor admits a Cartier multiple.
    # This is equivalent to whether the minimal ambient toric variety is Q-factorial.
    local isQfactorialVal := undefined;

    # A list of local picard indices, one for each X-cone. Might contain infinite values
    # if the variety is not Q-factorial.
    local localPicardIndices := undefined;

    # The picard index of the variety, i.e. the index of the picard group in the divisor class group.
    # It is given as the least common multiple of all local picard indices.
    # It is infinity if and only if the variety is not Q-factorial
    local picardIndex := undefined;

    # Says whether the variety is factorial, i.e. if every weil divisor is Cartier
    # This is equivalent to whether the minimal ambient toric variety is smooth.
    # It is also euivalent to the picard index being one.
    local isFactorialVal := undefined;

    # Says whether the variety is Q-gorenstein, i.e. if the anticanonical class admits a Cartier multiple.
    local isQgorensteinVal := undefined;

    # A list of local gorenstein indices, one for each X-cone. Might contain infinite values
    # if the variety is not Q-gorenstein.
    local localGorensteinIndices := undefined;

    # The gorenstein index of the variety.
    # It is given as the least common multiple of all local gorenstein indices.
    # It is infinity if and only if the variety is not Q-gorenstein.
    local gorensteinIndex := undefined;

    # Says whether the variety is gorenstein, i.e. if the anticanonical class is Cartier.
    # This is equivalent to the gorenstein index being one.
    local isGorensteinVal := undefined;

    local ampleCone := undefined;

    local isFanoVal := undefined;

    ###############################################################################
    ### The following fields are only defined for K*-surfaces, i.e. when s = 1. ###
    ###############################################################################
    local intersectionTable := undefined;
    local anticanonicalSelfIntersection := undefined;

    # This field is only filled if the variety is created as some resolution of singularities from another variety
    export exceptionalDivisorsIndices := undefined;

    export ModuleApply :: static := proc()
        Object(ComplexityOneVariety, _passed);
    end;

    export setCoefficientMatrix :: static := proc(self :: ComplexityOneVariety, A :: Matrix)
        local P, f, lss, i, j, alpha;
        P := self:-P;
        f := P:-format;
        lss := P:-lss;
        self:-variables := [seq(seq(T[i,j], j = 1 .. f:-ns[i]), i = 1 .. f:-r), seq(S[i], i = 1 .. f:-m)];
        self:-monomials := [seq(mul([seq(T[i,j] ^ lss[i][j], j = 1 .. f:-ns[i])]), i = 1 .. f:-r)];

        alpha := (i,j) -> Determinant(Matrix([Column(A, [i,j])]));
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
    (2) classGroupRank = 1. In this case, the number of rays of the fan is one more than the ambient dimension, hence
        there is only one complete fan having P as generator matrix.
    (3) P admits a Fano variety. In this case, the anticanonical class is contained in the moving cone of P, so we can 
        use it to define a bunch of cones and hence an ambient fan. The resulting variety will be the unique Fano variety
        having P as its PMatrix. Note however, that there are other possible non-Fano varieties with P as PMatrix.

    The parameter `A` is the coefficient matrix for the trinomial euqations defining the variety.
    It is an optional parameter. If it is not provided, it will be left undefined. It can always be
    provided later with the `setCoefficientMatrix` procedure.

    *)
    export ModuleCopy :: static := proc(self :: ComplexityOneVariety, proto :: ComplexityOneVariety, P :: PMatrix)
        local numColumns, i, j, k, ordered_indices, taus, sigma_plus, sigma_minus, taus_plus, taus_minus, w, candidates, minimalBunchCones, cone, A;

        self:-P := P;

        # Check the arguments, set Sigma and A if present.
        if _npassed > 3 then
            if type(_passed[4], Matrix) then
                # No Sigma provided, coefficient Matrix is fourth argument.
                A := _passed[4];
                if RowDimension(A) <> 2 or ColumnDimension(A) <> P:-r then
                    error "The coefficient matrix must be a (2 x r)-Matrix. Here, r = %1", P:-r;
                end if;
                for i from 1 to P:-r do
                    for j from i+1 to P:-r do
                        if Determinant(Matrix([Column(A, [i,j])])) = 0 then
                            error "The columns of the coefficient matrix must be linearly independent. Here the %1-th and the %2-th columns are linearly dependent.", i, j;
                        end if;
                    end do;
                end do;
                self:-A := A;
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
                        if containsrelint(poshull(Column(getDegreeMatrixFree(P), [op(cone)])), w) then
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
                    A := _passed[5];
                    if RowDimension(A) <> 2 or ColumnDimension(A) <> P:-r then
                        error "The coefficient matrix must be a (2 x r)-Matrix. Here, r = %1", P:-r;
                    end if;
                    for i from 1 to P:-r do
                        for j from i+1 to P:-r do
                            if Determinant(Matrix([Column(A, [i,j])])) = 0 then
                                error "The columns of the coefficient matrix must be linearly independent. Here the %1-th and the %2-th columns are linearly dependent.", i, j;
                            end if;
                        end do;
                    end do;
                    self:-A := A;
                else
                    error "Expected third argument to be of type Matrix";
                end if;
            end if;
        end if;

        if not type(self:-A, undefined) then
            setCoefficientMatrix(self, self:-A);
        end if;

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

            elif P:-classGroupRank = 1 then
                # In this case, the picard number is one, hence there is only one complete fan in the lattice 
                # containing the columns of the P-Matrix as rays.
                numColumns := ColumnDimension(P:-mat);
                self:-Sigma := {seq({seq(1 .. numColumns)} minus {i}, i = 1 .. numColumns)};
            elif admitsFano(P) then
                # In this case, we call the procedure again with the anticanonical class as the weight.
                return ComplexityOneVariety[ModuleCopy](self, proto, P, getAnticanonicalClass(P));
            else
                error "This PMatrix is neither of Picard number one, nor is it a surface, nor does it admit a Fano variety."
                      "Therefore, you must provide a fan Sigma or a weight w as input.";
            end if;
        end if;
    end;

    export setMaximalXCones :: static := proc(self :: ComplexityOneVariety, maximalXCones :: set(set(integer))) 
        self:-maximalXCones := maximalXCones;
    end;

    export getMaximalXCones :: static := proc(self :: ComplexityOneVariety)
        if type(self:-maximalXCones, undefined) or 'forceCompute' in [_passed] then
            setMaximalXCones(self, getMaximalXConesFormat(self:-P:-format, self:-Sigma));
        end if;
        return self:-maximalXCones;
    end;

    export setMinimalAmbientFan :: static := proc(self :: ComplexityOneVariety, minimalAmbientFan :: FAN) 
        self:-minimalAmbientFan := minimalAmbientFan;
    end;

    export getMinimalAmbientFan :: static := proc(self :: ComplexityOneVariety)
        if type(self:-minimalAmbientFan, undefined) or 'forceCompute' in [_passed] then
            setMinimalAmbientFan(self, fan(op(map(c -> intSetConeToConvexCone(self:-P, c), getMaximalXCones(self)))));
        end if;
        return self:-minimalAmbientFan;
    end;

    export setIsQfactorialVal :: static := proc(self :: ComplexityOneVariety, isQfactorialVal :: boolean)
        self:-isQfactorialVal := isQfactorialVal
    end proc;

    export isQfactorial :: static := proc(self :: ComplexityOneVariety)
        if type(self:-isQfactorialVal, undefined) or 'forceCompute' in [_passed] then
            setIsQfactorialVal(self, issimplicial(getMinimalAmbientFan(self)));
        end if;
        return self:-isQfactorialVal;
    end proc;

    export setIsFactorialVal :: static := proc(self :: ComplexityOneVariety, isFactorialVal :: boolean)
        self:-isFactorialVal := isFactorialVal
    end proc;

    export isFactorial :: static := proc(self :: ComplexityOneVariety)
        if type(self:-isFactorialVal, undefined) or 'forceCompute' in [_passed] then
            setIsFactorialVal(self, isregular(getMinimalAmbientFan(self)));
        end if;
        return self:-isFactorialVal;
    end proc;

    export setLocalPicardIndices :: static := proc(self :: ComplexityOneVariety, localPicardIndices)
        self:-localPicardIndices := localPicardIndices;
    end proc;

    export getLocalPicardIndices :: static := proc(self :: ComplexityOneVariety)
        if type(self:-localPicardIndices, undefined) or 'forceCompute' in [_passed] then
            setLocalPicardIndices(self, map(c -> localPicardIndex(self:-P, c), convert(getMaximalXCones(self), list)));
        end if;
        return self:-localPicardIndices;
    end proc;

    export setPicardIndex :: static := proc(self :: ComplexityOneVariety, picardIndex)
        self:-picardIndex := picardIndex;
    end proc;

    export getPicardIndex :: static := proc(self :: ComplexityOneVariety)
        if type(self:-picardIndex, undefined) or 'forceCompute' in [_passed] then
            setPicardIndex(self, lcm(op(getLocalPicardIndices(self))));
        end if;
        return self:-picardIndex;
    end proc;

    export getLocalCartierIndices :: static := proc(self :: ComplexityOneVariety, D :: list(integer))
        map(c -> localCartierIndex(self:-P, c, D), convert(getMaximalXCones(self), list));
    end proc;

    export getCartierIndex :: static := proc(self :: ComplexityOneVariety, D :: list(integer))
        lcm(op(getLocalCartierIndices(self, D)));
    end proc;

    export setLocalGorensteinIndices :: static := proc(self :: ComplexityOneVariety, localGorensteinIndices :: list)
        self:-localGorensteinIndices := localGorensteinIndices;
    end proc;

    export getLocalGorensteinIndices :: static := proc(self :: ComplexityOneVariety)
        if type(self:-localGorensteinIndices, undefined) or 'forceCompute' in [_passed] then
            setLocalGorensteinIndices(self, map(c -> localGorensteinIndex(self:-P, c), convert(getMaximalXCones(self), list)));
        end if;
        return self:-localGorensteinIndices;
    end proc;

    export setGorensteinIndex :: static := proc(self :: ComplexityOneVariety, gorensteinIndex) 
        self:-gorensteinIndex := gorensteinIndex;
    end proc;

    export getGorensteinIndex :: static := proc(self :: ComplexityOneVariety)
        local cone;
        if type(self:-gorensteinIndex, undefined) or 'forceCompute' in [_passed] then
            setGorensteinIndex(self, lcm(op(getLocalGorensteinIndices(self))));
        end if;
        return self:-gorensteinIndex;
    end proc;

    export setIsGorensteinVal :: static := proc(self :: ComplexityOneVariety, isGorensteinVal :: boolean) 
        self:-isGorensteinVal := isGorensteinVal;
    end proc;

    export isGorenstein :: static := proc(self :: ComplexityOneVariety)
        local cone;
        if type(self:-isGorensteinVal, undefined) or 'forceCompute' in [_passed] then
            setIsGorensteinVal(self, evalb(getGorensteinIndex(self) = 1));
        end if;
        return self:-isGorensteinVal;
    end proc;

    export setIsQgorensteinVal :: static := proc(self :: ComplexityOneVariety, isQgorensteinVal :: boolean)
        self:-isQgorensteinVal := isQgorensteinVal;
    end proc;

    export isQgorenstein :: static := proc(self :: ComplexityOneVariety)
        if type(self:-isQgorensteinVal, undefined) or 'forceCompute' in [_passed] then
            setIsQgorensteinVal(self, evalb(getGorensteinIndex(self) < infinity));
        end if;
        return self:-isQgorensteinVal;
    end proc;

    export setAmpleCone :: static := proc(self :: ComplexityOneVariety, ampleCone :: CONE) self:-ampleCone := ampleCone; end proc;

    export getAmpleCone :: static := proc(self :: ComplexityOneVariety)
        local cone;
        if type(self:-ampleCone, undefined) or 'forceCompute' in [_passed] then
            setAmpleCone(self, intersection(seq(poshull(Column(getDegreeMatrixFree(self:-P), [op({seq(1 .. ColumnDimension(getDegreeMatrixFree(self:-P)))} minus cone)])), cone in getMaximalXCones(self))));
        end if;
        return self:-ampleCone;
    end proc;

    export setIsFanoVal :: static := proc(self :: ComplexityOneVariety, isFanoVal :: boolean) self:-isFanoVal := isFanoVal; end proc;

    export isFano :: static := proc(self :: ComplexityOneVariety)
        if type(self:-isFanoVal, undefined) or 'forceCompute' in [_passed] then
            setIsFanoVal(self, containsrelint(getAmpleCone(self), getAnticanonicalClass(self:-P)));
        end if;
        return self:-isFanoVal
    end proc;

    export setintersectionTable :: static := proc(self :: ComplexityOneVariety, intersectionTable :: Array) self:-intersectionTable := intersectionTable; end proc;

    export getIntersectionTable :: static := proc(X :: ComplexityOneVariety)   
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
    export intersectionNumber :: static := proc(X :: ComplexityOneVariety, D1 :: list(integer), D2 :: list(integer))
        local k1, k2;
        
        if nops(D1) <> X:-P:-n + X:-P:-m or nops(D2) <> X:-P:-n + X:-P:-m then
            error "The list of integers encoding the divisor must have length n + m = %1", X:-P:-n + X:-P:-m;
        end if;

        return add([seq(seq(D1[k1] * D2[k2] * getIntersectionTable(X)[k1,k2], k2 = 1 .. nops(D2)), k1 = 1 .. nops(D1))]);

    end proc;

    export setAnticanonicalSelfIntersection :: static := proc(self :: ComplexityOneVariety, anticanonicalSelfIntersection) self:-anticanonicalSelfIntersection := anticanonicalSelfIntersection; end proc;

    export getAnticanonicalSelfIntersection :: static := proc(X :: ComplexityOneVariety)
        if type(X:-anticanonicalSelfIntersection, undefined) or 'forceCompute' in [_passed] then
            setAnticanonicalSelfIntersection(X, intersectionNumber(X, getAnticanCoefficients(X:-P), getAnticanCoefficients(X:-P)));
        end if;
        return X:-anticanonicalSelfIntersection;
    end proc;


    (****************************
    *** ADMISSIBLE OPERATIONS ***
    *****************************)

    export removeSingleRedundantBlock :: static := proc(X :: ComplexityOneVariety, i0 :: integer)
        local newP, oldToNewIndex, newSigma;
        newP := PMatrix[removeSingleRedundantBlock](X:-P, i0);
        oldToNewIndex := k -> if k < add(X:-P:-ns[1 .. i0 - 1]) + 1 then k else k - 1 end if;
        newSigma := map(cones -> map(oldToNewIndex, cones), X:-Sigma);
        return ComplexityOneVariety(newP, newSigma);
    end proc;

    export removeErasableBlocks :: static := proc(X :: ComplexityOneVariety)
        local P, redundantIndices;
        P := X:-P;
        redundantIndices := select(i -> P:-lss[i] = [1], [seq(1 .. nops(P:-lss))]);
        # If there is a redundant block and we still have more than two blocks, remove it.
        if nops(redundantIndices) > 0 and P:-r > 2 then
            # Recursive call
            return removeErasableBlocks(removeSingleRedundantBlock(X, redundantIndices[1]));
        else
            return X;
        end if;
    end proc;

    export applyAdmissibleOperation :: static := proc(X :: ComplexityOneVariety, a :: AdmissibleOperation)
        local newP, newSigma, newA;
        newP := PMatrix[applyAdmissibleOperation](X:-P, a);
        newSigma := map(cones -> map(k -> a:-bundledPerm[k], cones), X:-Sigma);
        if type(X:-A, undefined) then
            ComplexityOneVariety(newP, newSigma);
        else
            newA := a:-U . X:-A . a:-D . a:-sigmaPermutationMatrix^(-1);
            ComplexityOneVariety(newP, newSigma, newA);
        end if;

    end proc;

    export sortColumnsByLss :: static := proc(X :: ComplexityOneVariety)
       applyAdmissibleOperation(X, PMatrix[sortColumnsByLssOperation](X:-P));
    end proc;

    export standardizeCoefficientMatrixOperation :: static := proc(X :: ComplexityOneVariety)
        local U1, newA, ds, U2, i;
        U1 := Matrix([Column(X:-A,[1,2])])^(-1);
        newA := U1 . X:-A;
        ds := [1, newA[2,3] / newA[1,3], seq(-1 / newA[1,i], i = 3 .. X:-P:-r)];
        U2 := <1,0;0,newA[1,3]/newA[2,3]>;
        AdmissibleOperation[OnA](X:-P:-format, U2 . U1, ds);
    end proc;

    export standardizeCoefficientMatrix :: static := proc(X :: ComplexityOneVariety)
        applyAdmissibleOperation(X, standardizeCoefficientMatrixOperation(X));
    end proc;
    
    export areCoefficientMatricesEquivalentOperation :: static := proc(X1 :: ComplexityOneVariety, X2 :: ComplexityOneVariety)
        local a1, a2;
        a1 := standardizeCoefficientMatrixOperation(X1);
        a2 := standardizeCoefficientMatrixOperation(X2);
        if Equal(applyAdmissibleOperation(X1, a1):-A, applyAdmissibleOperation(X2, a2):-A) then
            return compose(a1, inverse(a2));
        end if;
        return NULL;
    end proc;

    export areCoefficientMatricesEquivalent :: static := proc(X1 :: ComplexityOneVariety, X2 :: ComplexityOneVariety)
        type(areCoefficientMatricesEquivalentOperation(X1, X2), AdmissibleOperation);
    end proc;

    export areIsomorphicOperation :: static := proc(X1_ :: ComplexityOneVariety, X2_ :: ComplexityOneVariety)
        local X1, X2, a, newX1, coefficientOp;
        X1 := removeErasableBlocks(X1_);
        X2 := removeErasableBlocks(X2_);

        for a in PMatrix[areEquivalentOperations](X1:-P, X2:-P) do
            newX1 := applyAdmissibleOperation(X1, a);
            if getMaximalXCones(newX1) = getMaximalXCones(X2) then
                # If there are no coefficient matrices, we are done.
                if type(X1:-A, undefined) or type(X2:-A, undefined) then
                    return a;
                end if;
                # Else we check if the coefficient matrices can be turned into each other by an
                # admissible operation on the coefficient matrix.
                coefficientOp := areCoefficientMatricesEquivalentOperation(newX1, X2);
                if type(coefficientOp, AdmissibleOperation) then
                    return compose(a, coefficientOp);
                end if;
            end if;
        end do;

        return NULL;
    end proc;

    export areIsomorphic :: static := proc(X1 :: ComplexityOneVariety, X2 :: ComplexityOneVariety)
        # If we are in the surface case and there is no A-Matrix to worry about, we can use the faster test.
        if X1:-P:-s = 1 and X2:-P:-s = 1 and (type(X1:-A, undefined) or type(X2:-A, undefined)) then
            return areEquivalent(X1:-P, X2:-P);
        end if;
        type(areIsomorphicOperation(X1,X2), AdmissibleOperation);
    end proc;

    (*
    Construct the tropical resolution of a ComplexityOneVariety.
    *)
    export tropicalResolution :: static := proc(X :: ComplexityOneVariety)
        local P, tropicalSheetsGenerators, tropicalCones, newCones, Pcols, newPMatrixColumns, newP, newSigma, k, i;

        P := X:-P;

        # The tropical cones are the cones rho_i x QQ^s, where rho_i is generated by the canonical basis vector e_i
        # and e_0 = -(e_1 + ... + e_r). 
        tropicalSheetsGenerators := [seq([0 $ (P:-r - 1 + k - 1), 1, 0 $ (P:-s - k)], k = 1 .. P:-s), seq([0 $ (P:-r -1 + k - 1), -1, 0 $ (P:-s - k)], k = 1 .. P:-s)];
        tropicalCones := [poshull([-1 $ (P:-r - 1), 0 $ P:-s], op(tropicalSheetsGenerators)), seq(poshull([0 $ (i-1), 1, 0 $ (P:-r - 1 - i + P:-s)], op(tropicalSheetsGenerators)), i = 1 .. P:-r - 1)];
        
        # To get the new cones of our fan, we intersect each X-Cone with each of the tropical cones.
        # This might yield new rays that we append to our P-Matrix later
        newCones := map(c1 -> op(map(c2 -> intersection(c1,c2), tropicalCones)), map(c -> intSetConeToConvexCone(X:-P, c), getMaximalXCones(X)));
        
        # We construct newP, which arises from P by appending all the columns that arised during the intersection with the tropical cones.
        Pcols := map(v -> :-convert(v, list), [Column(P:-mat, [seq(1 .. P:-n + P:-m)])]);
        newPMatrixColumns := select(ray -> not (ray in Pcols), map(c -> op(rays(c)), newCones));
        newP := PMatrix(P:-s, <P:-mat | Transpose(Matrix([op(newPMatrixColumns)]))>);

        # Convert back from CONEs to lists of integers to get out new fan
        newSigma := map(c -> convexConeToIntSetCone(newP, c), newCones);

        return ComplexityOneVariety(newP, newSigma);

    end proc;

    export canonicalResolution :: static := proc(X0 :: ComplexityOneVariety)
        local X, P, newFan, Pcols, allRays, newRayBlocks, newPColumns, newP, newSigma, newX, i, exceptDivIndices;

        X := tropicalResolution(X0);
        P := X:-P;
        newFan := regularsubdiv(getMinimalAmbientFan(X));

        Pcols := map(v -> :-convert(v, list), [Column(P:-mat, [seq(1 .. P:-n + P:-m)])]);
        allRays := {op(map(c -> op(rays(c)), maximal(newFan)))};
        newRayBlocks := [select(ray -> not (ray in Pcols) and andmap(k -> ray[k] = -1, [seq(1 .. P:-r - 1)]), allRays)];
        for i from 2 to P:-r do
            newRayBlocks := [op(newRayBlocks), select(ray -> not (ray in Pcols) and ray[i-1] = 1 and andmap(k -> ray[k] = 0, {seq(1 .. P:-r - 1)} minus {i-1}), allRays)];
        end do;

        newPColumns := [];
        exceptDivIndices := [];
        for i from 1 to P:-r do
            newPColumns := [op(newPColumns), op(map(v -> :-convert(v, list), [Column(P:-mat, [seq(add(P:-ns[1 .. i-1]) + 1 .. add(P:-ns[1 .. i]))])]))];
            exceptDivIndices := [op(exceptDivIndices), seq(nops(newPColumns) + 1 .. nops(newPColumns) + nops(newRayBlocks[i]))];
            newPColumns := [op(newPColumns), op(newRayBlocks[i])];
        end do;
        newPColumns := [op(newPColumns), op(map(v -> :-convert(v, list), [Column(P:-mat, [seq(P:-n + 1 .. P:-n + P:-m)])]))];

        newP := PMatrix(P:-s, Transpose(Matrix(newPColumns)));
        newSigma := {op(map(c -> convexConeToIntSetCone(newP, c), maximal(newFan)))};
        newX := ComplexityOneVariety(newP, newSigma);
        newX:-exceptionalDivisorsIndices := exceptDivIndices;
        return newX;

    end proc;

    export minimalResolution :: static := proc(X0 :: ComplexityOneVariety)
        local X, P, contracticleIndices;
        X := canonicalResolution(X0);
        contracticleIndices := select(i -> getIntersectionTable(X)[i,i] = -1, X:-exceptionalDivisorsIndices);
        while contracticleIndices <> [] do
            X := ComplexityOneVariety(PMatrix(1, DeleteColumn(X:-P:-mat, contracticleIndices)));
            contracticleIndices := select(i -> getIntersectionTable(X)[i,i] = -1, X:-exceptionalDivisorsIndices);
        end do;
        return X;
    end proc;

    export ModulePrint :: static := proc(self :: ComplexityOneVariety)
        nprintf(cat("ComplexityOneVariety(dim = ", self:-P:-s + 1,
          ", lss = ", self:-P:-lss,
          ", Sigma = ", self:-Sigma));
    end;

    export ComplexityOneVarietyInfo :: static := proc(self :: ComplexityOneVariety)
        local P, i, relations, maximalXCones, Q, classGroup, classGroupRank, anticanonicalClass, effectiveConeRays, movingConeRays, ampleConeRays, isFano, gorensteinIndex, intersectionTable, anticanonicalSelfIntersection;
        print(P = self:-P:-mat);
        print([seq(cat(n,i), i = 0 .. self:-P:-r - 1), m] = [seq(self:-P:-ns[i], i = 1 .. self:-P:-r), self:-P:-m]);
        print(relations = self:-relations);
        print(maximalXCones = getMaximalXCones(self));
        print(Q = getDegreeMatrix(self:-P));
        print(classGroup = getClassGroup(self:-P));
        print(classGroupRank = self:-P:-classGroupRank);
        print(anticanonicalClass = getAnticanonicalClass(self:-P));
        print(effectiveConeRays = rays(getEffectiveCone(self:-P)));
        print(movingConeRays = rays(getMovingCone(self:-P)));
        print(ampleConeRays = rays(getAmpleCone(self)));
        print(isFano = ComplexityOneVariety[isFano](self));
        print(gorensteinIndex = getGorensteinIndex(self));
        if self:-P:-s = 1 then
            print(intersectionTable = getIntersectionTable(self));
            print(anticanonicalSelfIntersection = getAnticanonicalSelfIntersection(self));
        end if;
    end;

end module: