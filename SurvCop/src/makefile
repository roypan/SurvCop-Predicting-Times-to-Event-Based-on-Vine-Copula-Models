#mylib.so:
#	gcc -fpic -c pbnorm.c R_exchmvn.c frank.c
#	gfortran -fpic -c pnorm.f90 qtnorm.f90 BiCop.f90 bvncop.f90
#	gfortran -shared -o mylib.so *.o

survival_nllk: survival_nllk.f90 BiCop.f90 mixed_copula_llk.f90 pnorm.f90 qtnorm.f90 bvnorm.f
	gfortran -o survival_nllk.exe survival_nllk.f90 BiCop.f90 mixed_copula_llk.f90 pnorm.f90 qtnorm.f90 bvnorm.f

# gcc -c frank.c
# gfortran -c BiCop.f90
# gfortan -o BiCop BiCop.o frank.o
# BiCop < BiCop.in