      module external_mod
!
!     This is a module file to support the external test program.
!
!
!     H. P. Hratchian
!     Department of Chemistry & Chemical Biology
!     Center for Chemical Computation and Theory
!     University of California, Merced
!     hhratchian@ucmerced.edu
!
!
!
!     USE Connections
!
      use mqc_general
      use mqc_gaussian
      use mqc_algebra2
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=60
!
!
      CONTAINS
!
!
!PROCEDURE newtonRaphsonStep
      function newtonRaphsonStep(gradient,hessian) result(deltaX)
!
!     This function carries out a Newton-Raphson step and returns the NR change
!     in coordinates (deltaX).
!
!
      implicit none
      type(mqc_variable),intent(in)::gradient,hessian
      real(kind=real64),dimension(:),allocatable::deltaX
      integer(kind=int64)::i
      real(kind=real64),dimension(:),allocatable::gPrime
      type(mqc_variable)::hessianEVals,hessianEVecs,deltaXPrime

!
!     Allocate deltaX, form the hessian eigenvectors and eigenvalues, and rotate
!     the gradient to the hessian eigenvector basis (gPrime).
!
      Allocate(deltaX(Size(gradient)))
      call hessian%eigen(hessianEVals,hessianEVecs)
      gPrime = MatMul(Transpose(hessianEVecs),gradient)
!
!     Everything we need is now built and all that's left to do is to compute
!     the NR step. This is first done in the hessian eigenvector basis and then
!     rotate back to the coordinate basis.
!
      do i = 1,Size(gradient)
        if(abs(float(hessianEVals%getVal([i]))).gt.MQC_small) then
          deltaX(i) = -gPrime(i)/float(hessianEVals%getVal([i]))
        else
          deltaX(i) = 0.0
        endIf
      endDo
      deltaXPrime = deltaX
      deltaX = MatMul(hessianEVecs,deltaXPrime)
      return
!
      return
      end function newtonRaphsonStep
!
!
      end module external_mod
