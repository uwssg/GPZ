GPZ is a code which uses Gaussian Processes to find P(z_p), the probability
distribution that a galaxy of unknown redshift is at redshift z_p.  It requires
training and validation sets of galaxies with photometric data, uncertainties in 
that data, and spectroscopic redshifts.  The training and validation sets 
should be distinct (though the code will combine them when calculating the 
P(z_p)).

The code is written exclusively in C++.

If you have any questions, do not hesitate to email me at
scottvalscott@gmail.com

//////////////////////// Compiling the code /////////////

You will need to download the BLAS and LAPACK linear algebra packages.  Compile
them and generate the library files libblas.a and liblapack.a

Edit the enclosed Makefile to point to those libraries.

Once that is done, typing 'make gpz' should generate the executable gpz

This code can be run by typing ./gpz param_file.sav

param_file.sav is a file that sets various user parameters and tells gpz where
to find the data and where to put the output.  params_example.sav is an example
of such a file with some (hopefully) useful comments.

//////////////////////// Input Data ///////////////////

The training and validation sets should be formatted as for ANNz, to wit: in the
case of N photometric filters, the first N columns should be the fluxes or
magnitudes, the N+1 through 2N columns should be the uncertainties in those
fluxes/magnitudes (in the same order).  The 2N+1 column should be the
spectroscopic redshift.

As it is currently written, the list of unknown galaxies (the "test set") should
be formatted the same way, with a nonsense place holder in the place of the
spectroscopic redshift (the code was written for testing on data sets with known
spectroscopy so that we could assess the code's performance).

////////////////////// Output Data ////////////////

In your parameter file, you will tell GPZ to generate an output file.  That file
will contain the following parameters for each of your test galaxies:

chi1 -- the -2 ln[P(model | data)] for the unimodal model (see equation 22 of
our paper)

mu1 -- the peak redshift of the unimodal model

sig1 -- the variance of the unimodal model

chi2 -- the -2 ln[P(model | data)] for the bi-modal model (see equation 22 of
our paper)

mu2a -- the peak redshift of the low-z mode

sig2a -- the variance of the Gaussian making up the low-z mode

na -- the number of nearest neighbors used to generate the low-z mode

mu2b -- the peak redshift of the high-z mode

sig2b -- the variance of the Gaussian making up the high-z mode

nb -- the number of nearest neighbor galaxies used to generate the high-z mode

ntot -- the total number of nearest neighbor galaxies used (ntot=na+nb)

to turn these parameters into P(z_p) use

P(z_p)=exp(-0.5*chi1)*exp[-0.5*(z_p-mu1)^2/sig1]/sqrt(sig1) +
       exp(-0.5*chi2){ (na/ntot)*exp[-0.5*(z_p-mu2a)^2/sig2a]/sqrt(sig2a) +
                       (nb/ntot)*exp[-0.5*(z_p-mu2b)^2/sig2b]/sqrt(sig2b)}

as presented above, P(z_p) is not normalized

////////////////////// The code ///////////////////

goto_tools.cpp is just a file containing some generically useful subroutines
(mostly a sorting subroutine)

eigen_wrapper.cpp wraps the linear algebra codes from LAPACK into C++ (we
use the invert_lapack() subroutine to do matrix inversion and the
get_determinant() subroutine to get matrix determinants)

gaussian_process_driver.cpp actually does the heavy lifting of the Gaussian
Process.  The subroutine get_pdf() solves for P(z_p)

gpz_code.cpp this is the code that contains main(), i.e. it is the code which
interfaces with the user and calls the Gaussian Process
