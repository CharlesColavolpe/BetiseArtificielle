#########################################
###                                  ####
###         MODULE D'ALGEBRE         ####
###                                  ####
#########################################

# Recuperation des parametres
include("parametres.jl")

# Affichage propre d'une matrice
# entree : matrice (matrice reel) : matrice a afficher
function Imprimer(matrice::Matrix{Float64})
  for i in 1:dimension
    print("( ")
    for j in 1:dimension
      if (matrice[i,j] >= 0) print(' ') end
      print(matrice[i,j],"  ")
    end
    print(")\n")
  end
  return nothing
end

# Renvoie la norme infinie d'un vecteur
# entree : vecteur   (vecteur reel) : le vecteur a mesurer
# sortie : norme_sup (reel)         : la norme infinie du vecteur
function NormeSup(vecteur::Vector{Float64})::Float64
  norme_sup::Float64 = 0.
  for i in 1:dimension
    norme_sup = max(norme_sup,abs(vecteur[i]))
  end
  return norme_sup
end

# Renvoie le produit de la matrice A et d'un vecteur
# entree : x (vecteur reel) : le facteur
# sortie : y (vecteur reel) : le produit de U et x
function ProdAvec(x::Vector{Float64})::Vector{Float64}
  y = Vector{Float64}(undef, dimension)
  for i in 1:dimension
    for j in 1:dimension
      y[i] = y[i] + A[i,j]*x[j]
    end
  end
  return y
end

# Renvoie le produit de la matrice U et d'un vecteur (avec U la partie triangulaire superieure de A)
# nb : la boucle sur les colonnes est optimisee du fait que la matrice U est triangulaire.
# entree : x (vecteur reel) : le facteur
# sortie : y (vecteur reel) : le produit de U et x
function ProdUvec(x::Vector{Float64})::Vector{Float64}
  y = Vector{Float64}(undef, dimension)
  for i in 1:dimension
    for j in i+1:dimension
      y[i] = y[i] + U[i,j]*x[j]
    end
  end
  return y
end

# Renvoie la solution exacte de l'equation Lx = y (avec L la partie triangulaire inferieure de A)
# entree : y (vecteur reel) : le membre de droite de l'equation
# sortie : x (vecteur reel) : la solution de l'equation
function ResTriInf(y::Vector{Float64})::Vector{Float64}
  x = Vector{Float64}(undef, dimension)
  x[1] = y[1]/L[1,1]
  for i in 2:dimension
    x[i] = y[i]
    for j in 1:i-1
      x[i] = x[i] - L[i,j]*x[j]
    end
    x[i] = x[i]/L[i,i]
  end
  return x
end

# Renvoie la solution du systeme Ax = b par la methode de Gauss-Seidel
# entree : - b    (vecteur reel) : le second membre de l'equation
#          - x0   (vecteur reel) : la premiere valeur de l'algorithme
# sortie : - x    (vecteur reel) : solution approchee
#          - iter (entier)       : nombre d'iterations de l'algorithme
function GaussSeidel(x0::Vector{Float64},b::Vector{Float64})
  rhs::Vector{Float64} = b - ProdUvec(x0)
  x::Vector{Float64}   = ResTriInf(rhs)
  erreur::Float64      = NormeSup(ProdAvec(x)-b)
  iter::UInt16         = 1
  while ( (erreur>tolerance) && (iter<iter_max) )
    rhs    = b - ProdUvec(x)
    x      = ResTriInf(rhs)
    erreur = NormeSup(ProdAvec(x)-b)
    iter   = iter + 1
  end
  return (x,iter)
end