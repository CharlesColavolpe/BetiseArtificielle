#########################################
###                                  ####
###        FICHIER PRINCIPAL         ####
###                                  ####
#########################################

# Recuperation des parametres et des fonctions
include("algebre.jl")

b::Vector{Float64}  = [1.,1.,1.]
x0::Vector{Float64} = [1.0,1.5,1.5]

print("Nous cherchons la solution approchee a ",tolerance," de Ax=b,","\n")
print("avec b = ",b," et A = ","\n")
Imprimer(A)

(x,n) = GaussSeidel(x0,b)
print("Pour x0 = ",x0,"\n")
print("Gauss-Seidel converge en ",n," iteration(s)","\n")