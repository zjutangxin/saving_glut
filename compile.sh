# Bash script to compile the file
echo "Compile the Saving Glut program."

gfortran ./src/mod_param.f90 ./src/mod_global.f90 ./src/mod_routines.f90 \
     ./src/solve_system.f90 ./src/saving_glut_proto_main.f90 \
     ./obj/minpack.o -llapack -latlas -lblas \
     -o sgdebt_main.out -Ofast -march=native -Wall

echo "Compiled"