#########################################
###                                  ####
###            PARAMETRES            ####
###                                  ####
#########################################

# Parametres
dimension::UInt16  = 3
tolerance::Float64 = 10^-5
iter_max::UInt16   = 1000

# Definition de la matrice A 
A = Matrix{Float64}(undef, dimension, dimension)
A[1,1] =  2.
A[1,2] = -1.
for i in 2:dimension-1
  A[i,i-1] = -1.
  A[i,i]   =  2.
  A[i,i+1] = -1.
end
A[dimension,dimension-1] = -1.
A[dimension,dimension]   =  2.

# Defintion de la partie superieure U et de la partie inferieure L de A
U = Matrix{Float64}(undef, dimension, dimension)
L = Matrix{Float64}(undef, dimension, dimension)
for i in 1:dimension
  for j in i+1:dimension
    U[i,j] = A[i,j]
  end
  for j in 1:i
    L[i,j] = A[i,j]
  end
end