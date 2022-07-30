# ============================================================================
# Version academique (minimaliste) de l'algorithme primal simplexe (phase 2)
#
# Implementation revisee de
#   Julia programming for Operations Research (Changhyun Kwon)
#
# X. Gandibleux, Juin 2017 -- validee sur Julia 0.5
# ============================================================================

# Le probleme doit etre donne
# 1) en minimisation
# 2) sous sa forme standard
# 3) toutes les variables non negatives

# La version implementee ne comporte pas
# 1) de regle empechant le cyclage (regle de bland)
# 2) de traitements garantissant une stabilite numerique
# Les evolutions envisageables concernent
# 1) la prise en compte de contraintes d'inegalites
# 2) la phase 1 du simplexe pour trouver une solution initiale admissible
# 3) la prise en compte de bornes sur les variables
# 4) le passage a la version revisee du simplexe

# ----------------------------------------------------------------------------

#module ParametricSimplexMethod

  # --------------------------------------------------------------------------
  # package(s) utilise(s) dans le module

  using Combinatorics

  # --------------------------------------------------------------------------
  # methode(s) visible(s) a l'exterieur du module

#  export parametric_simplex_method

  # --------------------------------------------------------------------------
  # type rassemblant toutes les donnees du tableau simplexe

  type SimplexTableau
    I         ::Array{Int64}   # indices des variables de base
    xB        ::Array{Float64} # activite des variables de base
    a00_f1    ::Float64        # valeur de la fonction 1
    a00_f2    ::Float64        # valeur de la fonction 2
    alpha0_f1 ::Array{Float64} # Couts reduits fonction 1
    alpha0_f2 ::Array{Float64} # Couts reduits fonction 2
    Y         ::Array{Float64} # Y = B-1 * T
  end

  # --------------------------------------------------------------------------
  # verifie si tous les elements du vecteur sont non-negatifs

  function isnonnegative(x::Array{Float64})

    return length( x[ x .< 0] ) == 0

  end

  # --------------------------------------------------------------------------
  # calcule une solution de base admissible initiale

  function initial_BFS(T, d)

    # dimension de la matrice des contraintes
    m, n = size(T)

    # combinaison de m elements parmi n
    comb = collect(combinations(1:n, m))
    for i in length(comb):-1:1
      I = comb[i]
      B = T[:, I]
      xB = inv(B) * d
      if isnonnegative(xB)
        # base admissible trouvee
        return I, xB, B
      end
    end

    error("Infeasible")

  end

  # --------------------------------------------------------------------------
  # affichage ecran d'un tableau

  function print_tableau(t::SimplexTableau)

    m, n = size(t.Y)

    hline0 = repeat("-", 6)
    hline1 = repeat("-", 7*n)
    hline2 = repeat("-", 7)
    hline = join([hline0, "+", hline1, "+", hline2])

    println(hline)

    @printf("%6s|", "")
    for j in 1:length(t.alpha0_f1)
      @printf("%6.2f ", t.alpha0_f1[j])
    end
    @printf("| %6.2f\n", t.a00_f1)

    println(hline)

    for i in 1:m
      @printf("x[%2d] |", t.I[i])
      for j in 1:n
        @printf("%6.2f ", t.Y[i,j])
      end
      @printf("| %6.2f\n", t.xB[i])
    end

    println(hline)

  end

  # --------------------------------------------------------------------------
  # affichage ecran d'un tableau pour un probleme bi-objectif

  function print_tableau2(t::SimplexTableau)

    m, n = size(t.Y)

    hline0 = repeat("-", 6)
    hline1 = repeat("-", 7*n)
    hline2 = repeat("-", 7)
    hline = join([hline0, "+", hline1, "+", hline2])

    println(hline)

    @printf("%6s|", "")
    for j in 1:length(t.alpha0_f1)
      @printf("%6.2f ", t.alpha0_f1[j])
    end
    @printf("| %6.2f\n", t.a00_f1)

    @printf("%6s|", "")
    for j in 1:length(t.alpha0_f2)
      @printf("%6.2f ", t.alpha0_f2[j])
    end
    @printf("| %6.2f\n", t.a00_f2)

    println(hline)

    for i in 1:m
      @printf("x[%2d] |", t.I[i])
      for j in 1:n
        @printf("%6.2f ", t.Y[i,j])
      end
      @printf("| %6.2f\n", t.xB[i])
    end

    println(hline)

  end

  # --------------------------------------------------------------------------
  # pivotage

  function pivoting!(t::SimplexTableau)

    m, n = size(t.Y)

    # indices de la variable entrante et sortante
    entering, exiting = pivot_point(t)
    println("Pivoting: entering = x_$entering, exiting = x_$(t.I[exiting])")

    # Pivoting: exiting-row, entering-column
    # updating exiting-row
    pivot = t.Y[exiting, entering]
    t.Y[exiting, :] /= pivot
    t.xB[exiting] /= pivot

    # updating other rows of Y
    for i in setdiff(1:m, exiting)
      pivot = t.Y[i, entering]
      t.Y[i, :] -= pivot * t.Y[exiting, :]
      t.xB[i] -= pivot * t.xB[exiting]
    end

    # updating the row for the reduced costs
    coef = t.alpha0_f1[entering]
    t.alpha0_f1 -= coef * t.Y[exiting, :]'
    t.a00_f1 -= coef * t.xB[exiting]

    # Updating I
    t.I[ find(t.I.==t.I[exiting]) ] = entering

  end

  # --------------------------------------------------------------------------
  # pivotage pour un probleme bi-objectif

  function pivoting_biobj!(t::SimplexTableau)

    m, n = size(t.Y)

    # indices de la variable entrante et sortante
    entering, exiting = pivot_point_biobj(t)
    println("Pivoting: entering = x_$entering, exiting = x_$(t.I[exiting])")

    # Pivoting: exiting-row, entering-column
    # updating exiting-row
    pivot = t.Y[exiting, entering]
    t.Y[exiting, :] /= pivot
    t.xB[exiting] /= pivot

    # updating other rows of Y
    for i in setdiff(1:m, exiting)
      pivot = t.Y[i, entering]
      t.Y[i, :] -= pivot * t.Y[exiting, :]
      t.xB[i] -= pivot * t.xB[exiting]
    end

    # updating the row for the reduced costs
    coef = t.alpha0_f1[entering]
    t.alpha0_f1 -= coef * t.Y[exiting, :]'
    t.a00_f1 -= coef * t.xB[exiting]

    coef = t.alpha0_f2[entering]
    t.alpha0_f2 -= coef * t.Y[exiting, :]'
    t.a00_f2 -= coef * t.xB[exiting]

    # Updating I
    t.I[ find(t.I.==t.I[exiting]) ] = entering

  end

  # --------------------------------------------------------------------------
  # identification du pivot

  function pivot_point(t::SimplexTableau)

    # Finding the entering variable index
    entering = findfirst(t.alpha0_f1 .> 0)
    if entering == 0
      error("Optimal")
    end

    # min ratio test / finding the exiting variable index
    pos_idx = find( t.Y[:, entering] .> 0 )
    if length(pos_idx) == 0
      error("Unbounded")
    end
    exiting = pos_idx[ indmin( t.xB[pos_idx] ./ t.Y[pos_idx, entering] ) ]

    return entering, exiting
  end

  # --------------------------------------------------------------------------
  # identification du pivot pour un probleme bi-objectif

  function pivot_point_biobj(tableau::SimplexTableau)

    # variable entrante
    Jc =[] # indices des variables candidates
    for i =1:length(tableau.alpha0_f2)
      if tableau.alpha0_f2[i] > 0
        push!(Jc, i)
      end
    end
    if length(Jc) == 0
      error("Optimal")
    end

    println("Jc = ", Jc)

    lambda = zeros(1,length(Jc))
    for j in Jc
      println("j=", j, "lambda[j]=", tableau.alpha0_f2[j]/(tableau.alpha0_f2[j] - tableau.alpha0_f1[j]))
      c2 = tableau.alpha0_f2[j]
      c1 = tableau.alpha0_f1[j]
      lambda[j]= c2/(c2 - c1)
    end
    println(lambda)
    IndlambdaPrim= indmax(lambda)
    entering = Jc[IndlambdaPrim]
    println(entering)

    # variable sortante
    pos_idx = find( tableau.Y[:, entering] .> 0 )
    if length(pos_idx) == 0
      error("Unbounded")
    end
    exiting = pos_idx[ indmin( tableau.xB[pos_idx] ./ tableau.Y[pos_idx, entering] ) ]
    println(exiting, " ", tableau.I[exiting])

    return entering, exiting

  end

  # --------------------------------------------------------------------------
  # construction du premier tableau simplexe

  function initialize(c, T, d)

    c = Array{Float64}(c)
    T = Array{Float64}(T)
    d = Array{Float64}(d)

    m, n = size(T)

    # Finding an initial BFS
    I, xB, B = initial_BFS(T,d)

    Y = inv(B) * T
    c_B = c[I]
    a00_f1 = dot(c_B, xB)
    a00_f2 = 0.0

    # alpha0_f1 is a row vector
    alpha0_f1 = zeros(1,n)
    alpha0_f2 = zeros(1,n)
    J= setdiff(1:n, I)
    alpha0_f1[J] = c_B' * inv(B) * T[:,J] - c[J]'

    return SimplexTableau(I, xB, a00_f1, a00_f2, alpha0_f1, alpha0_f2, Y)

  end

  # --------------------------------------------------------------------------
  # test 1 : solution optimale si couts reduits sont negatifs

  function isOptimal(tableau)

    return findfirst( tableau.alpha0_f1 .> 0 ) == 0

  end

  # --------------------------------------------------------------------------
  # test 1 : solution optimale si couts reduits sont negatifs pour un probleme bi-objectif

  function isOptimal_biobj(tableau)

    return findfirst( tableau.alpha0_f2 .> 0 ) == 0

  end

  # --------------------------------------------------------------------------
  # methode exportee du module

  function parametric_simplex_method(c1, c2, T, d)

    tableau = initialize(c1, T, d)
    print_tableau(tableau)

    while !isOptimal(tableau)
      pivoting!(tableau)
      print_tableau(tableau)
    end

    opt_x = zeros(length(c1))
    opt_x[tableau.I] = tableau.xB


    # Calcule les couts reduits de f2 connaissant l'optimum sur f1
    m, n = size(T)
    I = tableau.I # calcule les info pour la solution de base courante
    B = T[:,I]
    c2_B = c2[I]
    tableau.a00_f2 = dot(c2_B, tableau.xB)

    tableau.alpha0_f2 = zeros(1,n)
    J = setdiff(1:n, I)
    tableau.alpha0_f2[J] = c2_B' * inv(B) * T[:,J] - c2[J]'
    println(" ")
    print_tableau2(tableau)

    while !isOptimal_biobj(tableau)
      pivoting_biobj!(tableau)
      print_tableau2(tableau)
    end

#    entering, exiting = pivot_point_biobj(tableau)

    # ----

    return

  end

#end
