gg = g++
LAPACK = -L/Users/noldor/physics/lapack-3.1.1/ -llapack
ARPACK = -L/Users/noldor/physics/ARPACK/ -larpack
BLAS = -L/Users/noldor/physics/BLAS/ -lblas
FORTRAN = -lgfortran

goto_tools.o: goto_tools.cpp goto_tools.h
	$(gg) -c goto_tools.cpp

kd.o: kd.cpp kd.h goto_tools.o
	$(gg) -c kd.cpp goto_tools.o

eigen_wrapper.o: eigen_wrapper.cpp eigen_wrapper.h goto_tools.o
	$(gg) -c eigen_wrapper.cpp goto_tools.o \
	$(LAPACK) $(BLAS) $(FORTRAN)

gaussian_process_noisy.o: gaussian_process_noisy.h gaussian_process_noisy.cpp \
goto_tools.o eigen_wrapper.o kd.o goto_tools.o
	$(gg) -c gaussian_process_noisy.cpp goto_tools.o eigen_wrapper.o kd.o \
	$(LAPACK) $(ARPACK) $(BLAS) $(FORTRAN)
	

gpnoise: gpnoisy_test.cpp gaussian_process_noisy.o
	$(gg) -o gpnoise gpnoisy_test.cpp goto_tools.o \
	kd.o eigen_wrapper.o gaussian_process_noisy.o \
        $(LAPACK) $(BLAS) $(FORTRAN)

gpz_mufold: gpz_code.cpp gaussian_process_noisy.o
	$(gg) -o gpz_mufold gpz_code.cpp goto_tools.o \
	kd.o eigen_wrapper.o gaussian_process_noisy.o \
        $(LAPACK) $(BLAS) $(FORTRAN)

testkd: test_kdtree.cpp kd.o goto_tools.o
	$(gg) -o testkd test_kdtree.cpp goto_tools.o kd.o
