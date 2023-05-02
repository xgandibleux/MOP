# ==============================================================================
# ==============================================================================

println("""\nA study of algorithms for computing a lower bound set -------\n""")

global verbose    = true
global graphic    = true
global exact      = true
global showX      = true
global experiment = false

print("  verbose.....: "); verbose    ? println("yes") : println("no") 
print("  graphics....: "); graphic    ? println("yes") : println("no") 
print("  exact.......: "); exact      ? println("yes") : println("no") 
print("  showX.......: "); showX      ? println("yes") : println("no") 
print("  experiment..: "); experiment ? println("yes") : println("no") 


println("\n\n-) Activate the required packages\n")
using JuMP, GLPK, MathOptInterface
import MultiObjectiveAlgorithms as MOA
using Printf, PyPlot
println("  Done\n")


# ==============================================================================
# ==============================================================================
# collect the un-hidden filenames available in a given directory

function getfname(target)
    # target : string := chemin + nom du repertoire ou se trouve les instances

    # positionne le currentdirectory dans le repertoire cible
    #cd(joinpath(homedir(),target))
    savePwd = pwd()
    newPwd = pwd() * "/" * target
    cd(newPwd)
    # retourne le repertoire courant
    println("pwd = ", pwd())

    # recupere tous les fichiers se trouvant dans le repertoire data
    allfiles = readdir()

    # vecteur booleen qui marque les noms de fichiers valides
    flag = trues(size(allfiles))

    k=1
    for f in allfiles
        # traite chaque fichier du repertoire
        if f[1] != '.'
            # pas un fichier cache => conserver
            println("fname = ", f)
        else
            # fichier cache => supprimer
            flag[k] = false
        end
        k = k+1
    end

    cd(savePwd)

    # extrait les noms valides et retourne le vecteur correspondant
    finstances = allfiles[flag]
    return finstances
end


# ==============================================================================
# Parse an instance of 2SPA problem (bi-objective partionning problem) from a file

function parse2SPA(fname::String)

    f = open(fname)    
    nbctr, nbvar = parse.(Int, split(readline(f))) # nombre de contraintes , nombre de variables
    A = zeros(Int, nbctr, nbvar)                   # matrice des contraintes
    C = zeros(Int, 2,nbvar)                        # matrice des couts
    nb = zeros(Int, nbvar)
    for j in 1:nbvar
        flag = 1
        for valeur in split(readline(f))
            if flag == 1
                C[1,j] = parse(Int, valeur)
                flag +=1
            elseif flag == 2
                C[2,j] = parse(Int, valeur)
                flag +=1
            elseif flag == 3
                nb[j] = parse(Int, valeur)
                flag +=1
            else
                i = parse(Int, valeur)
                A[i,j] = 1
            end
        end
    end
    close(f)
    return C, A
end


# ==============================================================================
# Create (set) a 2SPA model (bi-objective partionning problem) with JuMP/MOA 
# from a instance file,  with all variables :Bin or :Con

function set2SPA(C::Array{Int,2}, A::Array{Int,2}, varType::Symbol)

    nbctr, nbvar = size(A)
    m2SPA = Model()
    if varType == :Bin
        @variable(m2SPA, x[1:nbvar], Bin)
    elseif varType == :Con
        @variable(m2SPA, 0.0 <= x[1:nbvar] <= 1.0 )
    else
        @assert false "error: unknow type to set for the variables"        
    end
    @expression(m2SPA, obj1, sum((C[1,i])*x[i] for i in 1:nbvar))
    @expression(m2SPA, obj2, sum((C[2,i])*x[i] for i in 1:nbvar))
    @objective(m2SPA, Min, [obj1, obj2])
    @constraint(m2SPA, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1) 

    return m2SPA
end


# ==============================================================================
# Create (load) a 2SPA model (bi-objective partionning problem) with JuMP/MOA 
# from a instance file,  with all variables set to :Bin

function load2SPA(fname::String)

    C, A = parse2SPA(fname)
    m2SPA = set2SPA(C, A, :Bin)

    return m2SPA
end


# ==============================================================================
# compute S_N for a 2SPA (JuMP/MOA model) for methods and options selected

function solve2SPA(  
    m2SPA::Model,          # a 2SPA model
    solverMIP::Symbol,     # MIP solver to use (:GLPK or :Gurobi)
    methodMOA::Symbol,     # MOA method to use (:EpsilonConstraint or :Dichotomy or :Lexicographic)
    varType::Symbol;       # nature of variables (:Bin or :Con)
    nbPoints::Int64=0      # number of points to compute (optional parameter with default value)
    )

    println("  $solverMIP  $methodMOA  $varType  $nbPoints")
    # ---- precondition
    if methodMOA==:EpsilonConstraint && varType==:Con && nbPoints==0
        @assert false "warning: no reasonable to probe all YN with ϵ-constraint and 0≤x≤1"        
    end        

    # ---- Get the dimensions of the 2SPA instance to solve
    nbvar = num_variables(m2SPA)
    nbctr = num_constraints(m2SPA, AffExpr, MathOptInterface.EqualTo{Float64})
    start = time()

    # ---- Relax the integrality constraints on variables if varType == :Con
    if varType == :Con
        undo_relax = relax_integrality(m2SPA)
    end

    # ---- Setting the solver
    if solverMIP == :GLPK
        set_optimizer(m2SPA, () -> MOA.Optimizer(GLPK.Optimizer))
    elseif solverMIP == :Gurobi
        set_optimizer(m2SPA, () -> MOA.Optimizer(Gurobi.Optimizer))
    else
        @assert false "error: unavailable MIP solver requested"
    end   

    set_silent(m2SPA)

    if methodMOA == :EpsilonConstraint
        set_optimizer_attribute(m2SPA, MOA.Algorithm(), MOA.EpsilonConstraint())
    elseif methodMOA == :Dichotomy
        set_optimizer_attribute(m2SPA, MOA.Algorithm(), MOA.Dichotomy())
    elseif methodMOA == :Lexicographic
        set_optimizer_attribute(m2SPA, MOA.Algorithm(), MOA.Lexicographic())
    else
        @assert false "error: unavailable MOA method requested"
    end

    if nbPoints > 0
        set_optimizer_attribute(m2SPA, MOA.SolutionLimit(), nbPoints)
    end

    optimize!(m2SPA)

    # ---- Querying the results
    #@show solution_summary(m2SPA)
    cardSN = result_count(m2SPA)
    verbose ? println("  cardSN = $cardSN") : nothing

    SN = Array{Number}(undef,2,cardSN) # matrix of performances
    SE = Array{Array{Number}}(undef,1,cardSN) # vector of decisions
    sumNbFrac = 0.0
    for i in 1:cardSN

        if varType == :Bin
            SN[1,i] = convert(Int64,value(m2SPA[:obj1]; result = i))
            SN[2,i] = convert(Int64,value(m2SPA[:obj2]; result = i))
            verbose ? @printf("  %3d: z=[%6d,%6d] | ", i, SN[1,i], SN[2,i]) : nothing
            SE[i] = [convert(Int64,value(m2SPA[:x][j]; result = i)) for j in 1:nbvar]
            verbose && showX ? print("x= ", SE[i], "  ") : nothing            
        else
            SN[1,i] = value(m2SPA[:obj1]; result = i)
            SN[2,i] = value(m2SPA[:obj2]; result = i)
            verbose ? @printf("  %3d: z=[%9.2f,%9.2f] | ", i, SN[1,i], SN[2,i]) : nothing
            SE[i] = [value(m2SPA[:x][j]; result = i) for j in 1:nbvar]
            verbose && showX ? print("x= ", SE[i], "  ") : nothing
        end

        nbUns,nbFrac = examineVectorVariables( value.(m2SPA[:x]; result = i) )
        sumNbFrac += nbFrac
        verbose ? print("nbOne=$nbUns  nbFrac=$nbFrac  \n") : nothing

    end

    verbose ? println("  nbVar = $nbvar  moyNbFrac = ",sumNbFrac/cardSN, " ⇒ ", round(100*sumNbFrac/cardSN/nbvar, digits=2),"%") : nothing

    # ---- Restore the integrality constraints on variables if varType == :Con    
    if varType == :Con
        undo_relax()
    end

    elapsedTime = time()-start
    println("  Elapsed time: $(round(elapsedTime,digits=3))s \n\n ")

    return SN, SE
end


# ==============================================================================
# Compute the performance on 2 objectives of a solution x

function evaluerSolution(
    x::Vector{Float64}, 
    C::Array{Int,2}
    )

    z1 = 0.0; z2 = 0.0
    for i in 1:length(x)
        z1 += x[i] * C[1,i]
        z2 += x[i] * C[2,i]
    end
    return round(z1, digits=2), round(z2, digits=2)
end

# ==============================================================================
# Compute the performance on 2 objectives of a solution x

function examineVectorVariables(x)

    nbvar = length(x)
    nbUns = 0
    nbFrac = 0
    for j in 1:nbvar
        v = value(x[j])
        if round(v,digits=6) == 1.0
            nbUns += 1
        end
        if round(v,digits=6) != 0.0 && round(v,digits=6) != 1.0
            nbFrac += 1
        end
    end
    return nbUns,nbFrac
end

# ==============================================================================
# Calcul des generateurs avec une ϵ-contrainte (ou z1 est la fonction a minimiser)
# sans usage a MOA.
# Comme le SPA est totalement entier : arrondi (ceil) sur la valeur de epsilon

function ϵConstraintSPAsingle(
    C::Array{Int,2},
    A::Array{Int,2}, 
    tailleSampling::Int64
    )

    nbctr, nbvar = size(A)
    start = time()

    sumNbFrac = 0.0
    cardSN = 2

    # optimize SPA f1 ---------------------------------------------------------
    mSPA = Model(GLPK.Optimizer)
    @variable(mSPA, 0.0 <= x[1:nbvar] <= 1.0 )   
    @objective(mSPA, Min, sum((C[1,i])*x[i] for i in 1:nbvar))
    @constraint(mSPA, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    optimize!(mSPA)
    xf1RL = value.(x)
    z1f1RL, z2f1RL = evaluerSolution(xf1RL, C)
    verbose ? @printf("  z1_LP: [ %8.2f , %8.2f ] ", z1f1RL, z2f1RL) : nothing

    nbUns,nbFrac = examineVectorVariables( value.(x) )
    sumNbFrac += nbFrac
    verbose ? print("nbOne=$nbUns  nbFrac=$nbFrac  \n") : nothing
    #print("x=", [convert(Int64,value(mSPA[:x][j]; result = i)) for j in 1:nbvar],"\n")

    # optimize SPA f2 ---------------------------------------------------------
    @objective(mSPA, Min, sum((C[2,i])*x[i] for i in 1:nbvar))
    optimize!(mSPA)
    xf2RL = value.(x)
    z1f2RL, z2f2RL = evaluerSolution(xf2RL, C)
    verbose ? @printf("  z2_LP: [ %8.2f , %8.2f ] ", z1f2RL, z2f2RL) : nothing

    nbUns,nbFrac = examineVectorVariables( value.(x) )    
    sumNbFrac += nbFrac
    verbose ? print("nbOne=$nbUns  nbFrac=$nbFrac  \n") : nothing
    #print("x=", [convert(Int64,value(mSPA[:x][j]; result = i)) for j in 1:nbvar],"\n")

    # optimize SPA f1 avec f2(x)≤ϵ  -------------------------------------------
    @objective(mSPA, Min, sum((C[1,i])*x[i] for i in 1:nbvar))
    @constraint(mSPA, epscst, sum((C[2,i])*x[i] for i in 1:nbvar) <= 0.0)

    minf1RL, maxf2RL = z1f1RL, z2f1RL
    maxf1RL, minf2RL = z1f2RL, z2f2RL

    pasSample2 = max(1, floor(Int, (maxf2RL - minf2RL) / (tailleSampling-3))) # pas de l'echantillonage sur z2
    
    SN1 = (Number)[z1f1RL]; SN2 = (Number)[z2f1RL]    
    j1 = 2

    while ceil(Int,maxf2RL) - (j1-1) * pasSample2 > ceil(Int,minf2RL)

            # minimise f1 avec ϵ-contrainte sur f2 -----------------------------
            verbose ? @printf("  z1 %2d : ϵ = %8.2f  ", j1, ceil(Int, ceil(Int,maxf2RL) - (j1-1) * pasSample2)) : nothing # echantillonage sur z2

            # calcul d'une solution epsilon-contrainte
            set_normalized_rhs(epscst, ceil(Int, ceil(Int,maxf2RL) - (j1-1) * pasSample2))
            optimize!(mSPA)
            f1RL = objective_value(mSPA)
            xf1RL = value.(x)
            if termination_status(mSPA) != OPTIMAL
                @assert false "status"
            end

            # recalcule la solution au regard des 2 objectifs
            z1f1RLcourant, z2f1RLcourant = evaluerSolution(xf1RL, C)
            verbose ? @printf("[ %8.2f , %8.2f ]   ", z1f1RLcourant, z2f1RLcourant) : nothing
            push!(SN1, z1f1RLcourant)
            push!(SN2, z2f1RLcourant)

            nbUns,nbFrac = examineVectorVariables( value.(x) )
            cardSN +=1; sumNbFrac += nbFrac
            verbose ? print("nbOne=$nbUns  nbFrac=$nbFrac  \n") : nothing
            #print("x=", [convert(Int64,value(m2SPA[:x][j]; result = i)) for j in 1:nbvar],"\n")

            # maj la valeur limite sur l'objectif 2 pour la solution courante
            j1 = j1+1
    end
    push!(SN1, z1f2RL); push!(SN2, z2f2RL)
    elapsedTime = time()-start
    verbose ? println("  nbVar = $nbvar  moyNbFrac = ",sumNbFrac/cardSN, " ⇒ ", round(100*sumNbFrac/cardSN/nbvar, digits=2),"%") : nothing

    println("  Elapsed time: $(round(elapsedTime,digits=3))s \n ")

    return  SN1, SN2
end


# ==============================================================================
# Calcul des generateurs avec une double ϵ-contrainte sans usage a MOA.
# Comme le SPA est totalement entier : arrondi (ceil) sur la valeur de epsilon

function ϵConstraintSPAdouble(
    C::Array{Int,2},
    A::Array{Int,2}, 
    tailleSampling::Int64
    )

    println("  GLPK  doubleEpsilonConstraint  Con")
    nbctr, nbvar = size(A)
    start = time()

    sumNbFrac = 0.0
    cardSN = 2

    # optimize SPA f1 ---------------------------------------------------------
    mSPAf1 = Model(GLPK.Optimizer)
    @variable(mSPAf1, 0.0 <= x[1:nbvar] <= 1.0 )   
    @objective(mSPAf1, Min, sum((C[1,i])*x[i] for i in 1:nbvar))
    @constraint(mSPAf1, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    optimize!(mSPAf1)
    xf1RL = value.(mSPAf1[:x])
    z1f1RL, z2f1RL = evaluerSolution(xf1RL, C)
    verbose ? @printf("  z1_LP: [ %8.2f , %8.2f ]  ", z1f1RL, z2f1RL) : nothing

    nbUns,nbFrac = examineVectorVariables( value.(mSPAf1[:x]) )
    sumNbFrac += nbFrac
    verbose ? print("nbOne=$nbUns  nbFrac=$nbFrac  \n") : nothing
    #print("x=", [convert(Int64,value(mSPAf1[:x][j]; result = i)) for j in 1:nbvar],"\n")

    # optimize SPA f1 avec f2(x)≤ϵ  -------------------------------------------
    @constraint(mSPAf1, epscst1, sum((C[2,i])*x[i] for i in 1:nbvar) <= 0.0)


    # optimize SPA f2 ---------------------------------------------------------
    mSPAf2 = Model(GLPK.Optimizer)
    @variable(mSPAf2, 0.0 <= x[1:nbvar] <= 1.0 )   
    @objective(mSPAf2, Min, sum((C[2,i])*x[i] for i in 1:nbvar))
    @constraint(mSPAf2, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    optimize!(mSPAf2)
    xf2RL = value.(mSPAf2[:x])
    z1f2RL, z2f2RL = evaluerSolution(xf2RL, C)
    verbose ? @printf("  z2_LP: [ %8.2f , %8.2f ]  ", z1f2RL, z2f2RL) : nothing  

    nbUns,nbFrac = examineVectorVariables( value.(mSPAf2[:x]) )
    sumNbFrac += nbFrac
    verbose ? print("nbOne=$nbUns  nbFrac=$nbFrac  \n") : nothing
    #print("x=", [convert(Int64,value(mSPAf2[:x][j]; result = i)) for j in 1:nbvar],"\n")

    # optimize SPA f2 avec f1(x)≤ϵ  -------------------------------------------
    @constraint(mSPAf2, epscst2, sum((C[1,i])*x[i] for i in 1:nbvar) <= 0.0)    


    minf1RL, maxf2RL = z1f1RL, z2f1RL
    maxf1RL, minf2RL = z1f2RL, z2f2RL

    SN1 = (Number)[z1f1RL]; SN2 = (Number)[z2f1RL]    
    graphic ? figure("Generators GM",figsize=(6.5,5)) : nothing
    graphic ? title("double ϵ-constraint") : nothing
    graphic ? xlabel(L"z^1(x)") : nothing
    graphic ? ylabel(L"z^2(x)") : nothing
    graphic ? scatter(z1f1RL,z2f1RL, c="magenta",marker="+") : nothing

    pasSample2 = max(1, floor(Int, (maxf2RL - minf2RL) / (tailleSampling-3))) # pas de l'echantillonage sur z2
    j1 = 2    
    pasSample1 = max(1, floor(Int, (maxf1RL - minf1RL) / (tailleSampling-3))) # pas de l'echantillonage sur z1
    j2 = 2


    alternance = 1
    maxf2RLlimite = maxf2RL
    minf2RLlimite = minf2RL

    moyNbFrac = 0.0
    while (maxf2RLlimite > minf2RLlimite) &&
          (ceil(Int,maxf2RL) - (j1-1) * pasSample2 > ceil(Int,minf2RL)) &&
          (ceil(Int,maxf1RL) - (j2-1) * pasSample1 > ceil(Int,minf1RL))

          #@show ceil(Int,maxf1RL) - (j2-1) * pasSample1, ceil(Int,minf1RL), (ceil(Int,maxf1RL) - (j2-1) * pasSample1 > ceil(Int,minf1RL))

        alternance += 1
        if alternance % 2 == 0

            # minimise f1 avec ϵ-contrainte sur f2 -----------------------------
            verbose ? @printf("  z1 %2d : ϵ = %8.2f  ", j1, ceil(Int, ceil(Int,maxf2RL) - (j1-1) * pasSample2)) : nothing # echantillonage sur z2

            # calcul d'une solution epsilon-contrainte
            set_normalized_rhs(epscst1, ceil(Int, ceil(Int,maxf2RL) - (j1-1) * pasSample2))
            optimize!(mSPAf1)
            f1RL = objective_value(mSPAf1)
            xf1RL = value.(mSPAf1[:x])
            if termination_status(mSPAf1) != OPTIMAL
                @assert false "error: over f1 in ϵConstraintSPAdouble"
            end

            # recalcule la solution au regard des 2 objectifs
            z1f1RLcourant, z2f1RLcourant = evaluerSolution(xf1RL, C)
            verbose ? @printf("[ %8.2f , %8.2f ]  ", z1f1RLcourant, z2f1RLcourant) : nothing
            push!(SN1, z1f1RLcourant)
            push!(SN2, z2f1RLcourant)
            graphic ? scatter(z1f1RLcourant,z2f1RLcourant, c="magenta") : nothing

            nbUns,nbFrac = examineVectorVariables( value.(mSPAf1[:x]) )
            cardSN +=1; sumNbFrac += nbFrac
            verbose ? print("nbOne=$nbUns  nbFrac=$nbFrac  \n") : nothing
            #print("x=", [convert(Int64,value(m2SPA[:x][j]; result = i)) for j in 1:nbvar],"\n")

            # maj la valeur limite sur l'objectif 2 pour la solution courante
            j1 = j1+1
            maxf2RLlimite = z2f1RLcourant

        else

            # minimise f2 avec ϵ-contrainte sur f1 -----------------------------
            verbose ? @printf("  z2 %2d : ϵ = %8.2f  ", j2, ceil(Int, ceil(Int,maxf1RL) - (j2-1) * pasSample1)) : nothing # echantillonage sur z1

            # calcul d'une solution epsilon-contrainte
            set_normalized_rhs(epscst2, ceil(Int, ceil(Int,maxf1RL) - (j2-1) * pasSample1))
            optimize!(mSPAf2)
            f2RL = objective_value(mSPAf2)
            xf2RL = value.(mSPAf2[:x])
            if termination_status(mSPAf2) != OPTIMAL
                @assert false "error: over f2 in ϵConstraintSPAdouble"
            end

            # recalcule la solution au regard des 2 objectifs
            z1f2RLcourant, z2f2RLcourant = evaluerSolution(xf2RL, C)
            verbose ? @printf("[ %8.2f , %8.2f ]  ", z1f2RLcourant, z2f2RLcourant) : nothing
            push!(SN1, z1f2RLcourant)
            push!(SN2, z2f2RLcourant)
            graphic ? scatter(z1f2RLcourant,z2f2RLcourant, c="cyan") : nothing

            nbUns,nbFrac = examineVectorVariables( value.(mSPAf2[:x]) )
            cardSN +=1; sumNbFrac += nbFrac
            verbose ? print("nbOne=$nbUns  nbFrac=$nbFrac  \n") : nothing
            #print("x=", [convert(Int64,value(m2SPA[:x][j]; result = i)) for j in 1:nbvar],"\n")

            # maj la valeur limite sur l'objectif 2 pour la solution courante
            j2 = j2+1
            minf2RLlimite = z2f2RLcourant
        end
    end
    push!(SN1, z1f2RL); push!(SN2, z2f2RL)
    elapsedTime = time()-start
    verbose ? println("  nbVar = $nbvar  moyNbFrac = ",sumNbFrac/cardSN, " ⇒ ", round(100*sumNbFrac/cardSN/nbvar, digits=2),"%") : nothing

    println("  Elapsed time: $(round(elapsedTime,digits=3))s \n ")

    graphic ? scatter(z1f2RL,z2f2RL, c="cyan", marker="+") : nothing

    return  SN1, SN2
end


# ==============================================================================
# Affichage graphique de tous les points calcules

function displayAll(fname,YN, YSN, LBE, LBD, LBEd1, LBEd2)

    # single ϵ-constraint with MOA vs full dichotomy
    fig1 = figure("Generators",figsize=(6.5,5))
    title("sampling with single ϵ-constraint vs full dichotomy")
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    scatter(LBE[1,:], LBE[2,:], c="red", marker="x", s=80, label=L"$LB$ eps") 
    plot(LBD[1,:], LBD[2,:], c="blue", marker="o", linestyle="dotted", label=L"$LB$ dic", markersize=5) 
    legend() 

    # single ϵ-constraint (16) with/without MOA
    #=
    fig2 = figure("Generators",figsize=(6.5,5))
    title("single ϵ-constraint (16) with/without MOA")
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    # lower bound with epsilon-constraint method using MOA
    scatter(LBE[1,:], LBE[2,:], c="red", marker="x", s=80, label=L"$LB$ eps") 
    # lower bound with a single homemade epsilon-constraint method (f1 opt; f2 eps), 0≤x≤1, ceil 
    scatter(SNs1, SNs2, c="red", marker="+", s=140, label=L"$LB$ eps1") 
    legend() 
    =#

    # lower bound sets
    if exact
        fexact = figure("Objective Space Y",figsize=(6.5,5))
        title("Non-dominated points and lower bound sets")
        z1min = minimum(vcat(YN[1,:],LBD[1,:]))
        z1max = maximum(vcat(YN[1,:],LBD[1,:]))
        z2min = minimum(vcat(YN[2,:],LBD[2,:]))
        z2max = maximum(vcat(YN[2,:],LBD[2,:]))
        # all non-dominated points
        scatter(YN[1,:], YN[2,:], c="lime", s=50, label=L"$Y_N$")
        # all non-dominated supported points
        plot(YSN[1,:], YSN[2,:], c="green", mec="lime", marker="o", linestyle="dotted", label=L"$Y_{SN}$", markersize=7) 
    else
        fapprox = figure("Objective Space Y",figsize=(6.5,5))
        title("lower bound sets")
        z1min = minimum(vcat(LBD[1,:]))
        z1max = maximum(vcat(LBD[1,:]))
        z2min = minimum(vcat(LBD[2,:]))
        z2max = maximum(vcat(LBD[2,:]))
    end

    zmin = min(z1min,z2min) * 0.9
    zmax = max(z1max,z2max) * 1.1       
    xlim(zmin,zmax)
    ylim(zmin,zmax)
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")

    # full lower bound with dichotomy method using MOA 
    plot(LBD[1,:], LBD[2,:], c="blue", marker="o", linestyle="dotted", label=L"$LB$ dic", markersize=5) 

    # lower bound with epsilon-constraint method using MOA
    #scatter(LBE[1,:], LBE[2,:], c="red", marker="x", s=80, label=L"$LB$ eps") 

    # lower bound with a single homemade epsilon-constraint method (f1 opt; f2 eps), 0≤x≤1, ceil 
    #scatter(SNs1, SNs2, c="red", marker="+", s=140, label=L"$LB$ eps1") 

    # lower bound with a double homemade epsilon-constraint method, 0≤x≤1, ceil 
    scatter(LBEd1, LBEd2, c="red", marker="+", s=140, label=L"$LB$ eps2 GM") 

    legend()  
    
end

# ==============================================================================
# Affichage graphique des resultats concis

function displayRes(fname,YN, YLD, YLEd1, YLEd2)

    if exact
        fexact = figure("Objective Space Y",figsize=(6.5,5))
        title("Non-dominated points and lower bound sets")
        z1min = minimum(vcat(YN[1,:],YLD[1,:]))
        z1max = maximum(vcat(YN[1,:],YLD[1,:]))
        z2min = minimum(vcat(YN[2,:],YLD[2,:]))
        z2max = maximum(vcat(YN[2,:],YLD[2,:]))
        # all non-dominated points
        scatter(YN[1,:], YN[2,:], c="lime", s=50, label=L"$Y_N$")
    else
        fapprox = figure("Objective Space Y",figsize=(6.5,5))
        title("lower bound sets")
        z1min = minimum(vcat(YLD[1,:]))
        z1max = maximum(vcat(YLD[1,:]))
        z2min = minimum(vcat(YLD[2,:]))
        z2max = maximum(vcat(YLD[2,:]))
    end

    zmin = min(z1min,z2min) * 0.9
    zmax = max(z1max,z2max) * 1.1       
    xlim(zmin,zmax)
    ylim(zmin,zmax)
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")

    # full lower bound with dichotomy method using MOA 
    plot(YLD[1,:], YLD[2,:], c="blue", marker="o", linestyle="dotted", label=L"$LB$ dic", markersize=5)  
    # lower bound with a double homemade epsilon-constraint method, 0≤x≤1, ceil 
    scatter(YLEd1, YLEd2, c="red", marker="+", s=140, label=L"$LB$ eps2 GM") 
    legend()  
    
end


# ==============================================================================
# ==============================================================================
# point d'entree principal

function main(fname::String)


    YRL = []
    YN  = []
    YSN = []
    LBD = []
    LBE = []
    LBEs1 = []; LBEs2 = []
    LBEd1 = []; LBEd2 = []    


    # --------------------------------------------------------------------------  
    # --------------------------------------------------------------------------  
    # with MOA 


    # --------------------------------------------------------------------------     
    println("\n0) instance and characteristics \n")
    println("  instance = $fname") 
   
    m2SPA = load2SPA(fname)
    println("  nbvar    = ",num_variables(m2SPA))
    println("  nbctr    = ",num_constraints(m2SPA, AffExpr, MathOptInterface.EqualTo{Float64}),"\n\n")  

    # -------------------------------------------------------------------------- 
    println("1) compute the range values of objectives with 0≤x≤1")

    YRL, XRL = solve2SPA(m2SPA, :GLPK, :Lexicographic, :Con)
    minf1RL, maxf2RL = YRL[1,1], YRL[2,1]
    maxf1RL, minf2RL = YRL[1,2], YRL[2,2]

    verbose ? @printf("  f1_min=%8.2f ↔ f1_max=%8.2f (Δ=%.2f) \n",minf1RL, maxf1RL, maxf1RL-minf1RL) : nothing
    verbose ? @printf("  f2_min=%8.2f ↔ f2_max=%8.2f (Δ=%.2f) \n\n\n",minf2RL, maxf2RL, maxf2RL-minf2RL) : nothing

    # -------------------------------------------------------------------------- 
    println("2) compute Y_N with ϵ-constraint method and x∈{0,1}")

    if exact
        YN, XE = solve2SPA(m2SPA, :GLPK, :EpsilonConstraint, :Bin)
        sizeYN = size(YN,2)
        nbProbe = 2^ceil(Int,log2(sizeYN))
    else
        nbProbe = 16
    end

    # -------------------------------------------------------------------------- 
    println("3) compute Y_SN with dichotomic method and x∈{0,1}")

    if exact 
        YSN, XSE = solve2SPA(m2SPA, :GLPK, :Dichotomy, :Bin)
    end

    # -------------------------------------------------------------------------- 
    println("4) compute LB(Y_N) with dichotomic method and 0≤x≤1")

    LBD, XLD = solve2SPA(m2SPA, :GLPK, :Dichotomy, :Con)  

    # -------------------------------------------------------------------------- 
    println("5) probe LB(Y_N) with ϵ-constraint method and 0≤x≤1")

    LBE, XLE = solve2SPA(m2SPA, :GLPK, :EpsilonConstraint, :Con, nbPoints=nbProbe)
    

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------  
    # without MOA

    # --------------------------------------------------------------------------    
    # parsing de l'instance numerique 
    C, A = parse2SPA(fname) 

    # --------------------------------------------------------------------------
    println("6) probe LB(Y_N) with ϵ-constraint method (homemade), ceil and 0≤x≤1")

    LBEs1, LBEs2 = ϵConstraintSPAsingle(C, A, nbProbe)
    
    # --------------------------------------------------------------------------
    println("\n7) probe LB(Y_N) with a double ϵ-constraint method (homemade), ceil and 0≤x≤1") 

    nbProbe = 12 
    LBEd1, LBEd2 = ϵConstraintSPAdouble(C, A, nbProbe)


    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Sortie graphique
    graphic ? displayAll(fname,YN, YSN, LBE, LBD, LBEd1, LBEd2) : nothing

end


# ==============================================================================
# Run all algorithms over all instances

function numericalExperiment(target)

    global graphic = false
    global exact   = false # some instances are very long to be solved with glpk
    
    fnames = getfname(target)
    for instance = 1:length(fnames)
        start = time()
        main( string(target,"/",fnames[instance]) )    
        elapsedTime = time()-start
        println("Elapsed time for $(fnames[instance]) : $elapsedTime (s)")
    end

    return nothing
end


# ==============================================================================

target = "SPA/instances" # path for a standard config on macOS

if experiment
    numericalExperiment(target)
else
    #@time main(target*"/bio"*"sppaa02.txt")
    #@time main(target*"/bio"*"sppnw03.txt")
    #@time main(target*"/bio"*"sppnw04.txt")
    #@time main(target*"/bio"*"sppnw10.txt")
    #@time main(target*"/bio"*"sppnw20.txt")
    #@time main(target*"/bio"*"sppnw25.txt")
    #@time main(target*"/bio"*"didactic3.txt")
    @time main(target*"/bio"*"didactic5.txt")
    #@time main(target*"/bio"*"sppnw29.txt")
    #@time main(target*"/bio"*"sppnw19.txt")
    #@time main(target*"/bio"*"sppnw40.txt")
end

nothing
