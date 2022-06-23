
ImportComplexityOneVarietyList := proc(stmt)
    local columns, INDEX_S, INDEX_P, INDEX_DIMENSION, INDEX_CLASSGROUP, INDEX_DEGREEMATRIX, INDEX_ANTICANCLASS, INDEX_AMBIENTFAN, INDEX_MAXIMALXCONES, INDEX_GORENSTEININDEX, INDEX_ISGORENSTEIN;
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
        if type(INDEX_CLASSGROUP, integer) then
            setClassGroup(P, parse(Fetch(stmt, INDEX_CLASSGROUP)));
        end if;
        if type(INDEX_ANTICANCLASS, integer) then
            setAnticanonicalClass(P, Vector(parse(Fetch(stmt, INDEX_ANTICANCLASS))));
        end if;
        if type(INDEX_DEGREEMATRIX, integer) then
            setDegreeMatrix(P, Matrix(parse(Fetch(stmt, INDEX_DEGREEMATRIX))));
        end if;

        # If a fan is supplied, make a ComplexityOneVariety.
        if type(INDEX_AMBIENTFAN, integer) then
            X := ComplexityOneVariety(P, parse(Fetch(stmt, INDEX_AMBIENTFAN)));
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
FindInDatabase := proc(connection, tableName :: string, X :: ComplexityOneVariety)
    local P, stmtString, stmt, Xs, M, rowids, i, resids;
    P := X:-P;
    stmtString := cat("SELECT rowid,* FROM ", tableName, " WHERE ",
        "m = ", P:-m, " AND ",
        "s = ", P:-s, " AND ",
        "classGroup = \"", getClassGroup(X:-P), "\" AND ",
        "gorensteinIndex = ", getGorensteinIndex(X));
    # Use the 'lss' as search criterion, only if X is non-toric.
    if not isToric(P) then
        stmtString := cat(stmtString, " AND orderedLss = \"", sortColumnsByLss(P):-lss, "\"");
    end if;
    
    stmt := Prepare(connection, stmtString);
    Xs := ImportComplexityOneVarietyList(stmt);
    M := FetchAll(stmt):
    rowids := [seq(M[i,1], i = 1 .. RowDimension(M))];
    resids := [];
    for i from 1 to nops(rowids) do
        if areIsomorphic(X, Xs[i]) then
            resids := [op(resids), rowids[i]];
        end if;
    end do;
    return resids;
end proc;

(*
Imports a single ComplexityOneVariety from a database with a given rowid.
To import multiple varieties using an arbitrary SQLite query, use `ImportComplexityOneVarietyList`
*)
ImportComplexityOneVariety := proc(db, tableName :: string, rowid :: integer)
    ImportComplexityOneVarietyList(Prepare(db, cat("SELECT * FROM ", tableName, " WHERE rowid = ", rowid)))[1];
end proc;

(*
Given a list of varieties of complexity one `Xs` and a SQLite database connection `db`, this function
inserts those varieties from `Xs` into the database that are not already present.
*)
ExportComplexityOneVarietyList := proc(connection, tableName :: string, Xs :: list(ComplexityOneVariety))
    local k, X, P, stmt, i, knownCount;
    knownCount := 0;

    for k from 1 to nops(Xs) do
        X := Xs[k];
        
        # Only insert `X` if it is not already present.
        # Can be disabled by the 'noChecks' option
        if 'noChecks' in [_passed] or FindInDatabase(connection, tableName, X) = [] then
            P := X:-P;
            stmt := Prepare(connection, cat("INSERT INTO ", tableName ," VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"));
            Bind(stmt, 1, P:-r, valuetype = "integer");
            Bind(stmt, 2, convert(P:-ns, string), valuetype = "text");
            Bind(stmt, 3, P:-n, valuetype = "integer");
            Bind(stmt, 4, P:-m, valuetype = "integer");
            Bind(stmt, 5, P:-s, valuetype = "integer");
            Bind(stmt, 6, convert(P:-lss, string), valuetype = "text");
            Bind(stmt, 7, convert([seq(convert(Row(P:-mat, i), list), i = 1 .. RowDimension(P:-mat))], string), valuetype = "text");
            Bind(stmt, 8, P:-dim, valuetype = "integer");
            Bind(stmt, 9, P:-classGroupRank, valuetype = "integer");
            Bind(stmt, 10, convert(getClassGroup(P), string), valuetype = "text");
            Bind(stmt, 11, convert([seq(convert(Row(getDegreeMatrix(P), i), list), i = 1 .. RowDimension(getDegreeMatrix(P)))], string), valuetype = "text");
            Bind(stmt, 12, convert(convert(getAnticanonicalClass(P), list), string), valuetype = "text");
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
            Bind(stmt, 25, convert(convert(map(x -> [numer(x), denom(x)], getintersectionMatrix(X)), list, nested), string), valuetype = "text");
            Bind(stmt, 26, convert([numer(getAnticanonicalSelfIntersection(X)), denom(getAnticanonicalSelfIntersection(X))], string), valuetype = "text");
            Bind(stmt, 27, convert(getAnticanonicalSelfIntersection(X), hfloat), valuetype = "float");
            Bind(stmt, 28, isToric(P), valuetype = "integer");
            Bind(stmt, 29, isIrredundant(P), valuetype = "integer");
            Bind(stmt, 30, isQfactorial(X), valuetype = "integer");
            Bind(stmt, 31, isFactorial(X), valuetype = "integer");
            Bind(stmt, 32, convert(getLocalPicardIndices(X), string), valuetype = "text");
            Bind(stmt, 33, getPicardIndex(X), valuetype = "integer");
            Bind(stmt, 34, isQgorenstein(X), valuetype = "integer");
            Bind(stmt, 35, convert(getLocalGorensteinIndices(X), string), valuetype = "text");
            Bind(stmt, 36, convert(getLocalGorensteinQuotients(X), string), valuetype = "text");
            Bind(stmt, 37, getGorensteinQuotient(X), valuetype = "integer");
            Bind(stmt, 38, X:-P:-case, valuetype = "text");
            Bind(stmt, 39, X:-P:-orientation, valuetype = "integer");

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

UpdateInDatabase := proc(connection, tableName :: string, rowid :: integer, X :: ComplexityOneVariety)
    local P, stmt;
    P := X:-P;    
    stmt := cat("UPDATE ", tableName, " SET ",
        "r = ", P:-r, ", ",
        "ns = \"", P:-ns, "\", ",
        "n = ", P:-n, ", ",
        "m = ", P:-m, ", ",
        "s = ", P:-s, ", ",
        "lss = \"", P:-lss, "\", ",
        "P = \"", convert(P:-mat, list, nested), "\", ",
        "dimension = ", P:-dim, ", ",
        "classGroupRank = ", P:-classGroupRank, ", ",
        "classGroup = \"", getClassGroup(P), "\", ",
        "degreeMatrix = \"", convert(getDegreeMatrix(P), list, nested), "\", ",
        "anticanClass = \"", convert(getAnticanonicalClass(P), list), "\", ",
        "ambientFan = \"", X:-Sigma, "\", ",
        "maximalXCones = \"", getMaximalXCones(X), "\", ",
        "gorensteinIndex = ", getGorensteinIndex(X), ", ",
        "isGorenstein = ", if isGorenstein(X) then 1 else 0 end if, ", ",
        "orderedLss = \"", sortColumnsByLss(P):-lss, "\", ",
        "effectiveCone = \"", rays(getEffectiveCone(P)), "\", ",
        "movingCone = \"", rays(getMovingCone(P)), "\", ",
        "ampleCone = \"", rays(getAmpleCone(X)), "\", ",
        "isFano = ", if isFano(X) then 1 else 0 end if, ", ",
        "variables = \"", X:-variables, "\", ",
        "monomials = \"", X:-monomials, "\", ",
        "relations = \"", X:-relations, "\", ",
        "intersectionTable = \"", convert(map(x -> [numer(x), denom(x)], getintersectionMatrix(X)), list, nested), "\", ",
        "anticanonicalSelfIntersectionFraction = \"", [numer(getAnticanonicalSelfIntersection(X)), denom(getAnticanonicalSelfIntersection(X))], "\", ",
        "anticanonicalSelfIntersectionFloat = ", convert(getAnticanonicalSelfIntersection(X), hfloat), ", ",
        "isToric = ", if isToric(P) then 1 else 0 end if, ", ",
        "isIrredundant = ", if isIrredundant(P) then 1 else 0 end if, ", ",
        "isQfactorial = ", if isQfactorial(X) then 1 else 0 end if, ", ",
        "isFactorial = ", if isFactorial(X) then 1 else 0 end if, ", ",
        "localPicardIndices = \"", getLocalPicardIndices(X), "\", ",
        "picardIndex = ", getPicardIndex(X), ", ",
        "isQgorenstein = ", if isQgorenstein(X) then 1 else 0 end if, ", ",
        "localGorensteinIndices = \"", getLocalGorensteinIndices(X), "\", ",
        "localGorensteinQuotients = \"", getLocalGorensteinQuotients(X), "\", ",
        "gorensteinQuotient = ", getGorensteinQuotient(X), ", ",
        "case_ = \"", P:-case, "\", ",
        "orientation = ", P:-orientation, " ",
        "WHERE rowid = ", rowid
        );
    Execute(db, stmt);
end proc;

(*
Helper function to perform an operation on every variety in the database, without loading them into memory all at once
(which can be unpractical, when the database is large).
The parameter `f` is a function taking a `ComplexityOneVariety` as first argument and an integer as second argument, which is the
position the variety occurs in the database.
*)
performOnDatabase := proc(db, tableName :: string, f, step := 1000)
    local numberOfEntries, numberOfSteps, i, offset, Xs, j;
    # First, count the number of database entries
    numberOfEntries := FetchAll(Prepare(db, cat("SELECT COUNT(*) FROM ", tableName)))[1,1];
    numberOfSteps := ceil(numberOfEntries / step);
    for i from 1 to numberOfSteps do 
        offset := (i - 1) * step;
        Xs := ImportComplexityOneVarietyList(Prepare(db, cat("SELECT * FROM ", tableName, " LIMIT ", step, " OFFSET ", offset)));
        for j from 1 to nops(Xs) do
            f(Xs[j], offset + j);
        end do;
    end do;
end proc;