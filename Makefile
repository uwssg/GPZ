gg = g++
LAPACK = -L/Users/noldor/physics/lapack-3.1.1/ -llapack
BLAS = -L/Users/noldor/physics/BLAS/ -lblas
FORTRAN = -lgfortran

goto_tools.o: goto_tools.cpp goto_tools.h
	$(gg) -c goto_tools.cpp

kd.o: kd.cpp kd.h goto_tools.o
	$(gg) -c kd.cpp goto_tools.o

eigen_wrapper.o: eigen_wrapper.cpp eigen_wrapper.h goto_tools.o
	$(gg) -c eigen_wrapper.cpp goto_tools.o \
	$(LAPACK) $(BLAS) $(FORTRAN)

gaussian_process_driver.o: gaussian_process_driver.h \
gaussian_process_driver.cpp \
goto_tools.o eigen_wrapper.o kd.o goto_tools.o
	$(gg) -c gaussian_process_driver.cpp goto_tools.o eigen_wrapper.o kd.o \
	$(LAPACK) $(BLAS) $(FORTRAN)
	

gpz: gpz_code.cpp gaussian_process_driver.o
	$(gg) -o gpz gpz_code.cpp goto_tools.o \
	kd.o eigen_wrapper.o gaussian_process_driver.o \
        $(LAPACK) $(BLAS) $(FORTRAN)
