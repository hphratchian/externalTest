INCLUDE "external_mod.f03"
      program external
!
!     This program tests external execution of Gaussian.
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
      use external_mod
!
!     Variable Declarations
!
      implicit none
      integer::i,numCmdLineArgs,systemStatus
      integer::nAtoms,nBasis,nBasisUse
      real(kind=real64)::energySCF
      real(kind=real64),dimension(:),allocatable::carts,deltaX
      character(len=512)::commandLine,gaussianRoute,  &
        gaussianMessageFilename,matrixFilename,timeString
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFileIn,  &
        GMatrixFileOut
!hph      type(mqc_variable)::gradient,hessian,deltaX
      type(mqc_variable)::gradient,hessian
!
!     Format Statements
!
 1000 format(/,1x,'Enter External.',/)
 1010 format(1x,'Matrix File: ',A)
 2000 format(/,1x,'Data from the external Gaussian job..')
 2100 format(3x,a,i)
 2150 format(3x,a,f20.6)
 5000 format(1x,'Time String: ',A)
!
!
      call get_command_argument(4,gaussianMessageFilename)
      open(unit=iOut,file=TRIM(gaussianMessageFilename))
      write(IOut,1000)
      call date_and_time(time=timeString)
      write(iOut,5000) TRIM(timeString)
!
!     Get the name of the input matrix file and open it.
!
      numCmdLineArgs = command_argument_count()
      write(iOut,*)' Hrant - numCmdLineArgs = ',numCmdLineArgs
      do i = 1,numCmdLineArgs
        call get_command_argument(i,matrixFilename)
        write(iOut,*)' command: ',TRIM(matrixFilename)
      endDo

!hph      write(IOut,1010) TRIM(matrixFilename)

!
!     Load the matrix file sent to us from the calling Gaussian job.
!
      call GMatrixFileIn%load(matrixFilename)
      nAtoms = GMatrixFileIn%getVal('nAtoms')
      carts = GMatrixFileIn%getAtomCarts()
      call mqc_print(carts,iOut,header='Cartesian Coordinates In:')
!
!     Set up the "input" file for the external Gaussian job.
!
      GMatrixFileOut = GMatrixFileIn
      GMatrixFileOut%cartesians(3) = -0.9
      GMatrixFileOut%cartesians(6) =  0.9
      call GMatrixFileOut%create('001a.faf')
      carts = GMatrixFileOut%getAtomCarts()
      call mqc_print(carts,iOut,header='Cartesian Coordinates Out:')
      call GMatrixFileOut%closeFile()
      call GMatrixFileIn%closeFile()
!
!     Now, do a loop of Gaussian jobs and optimization steps until a minimum is
!     found.
!
      call GMatrixFileIn%closeFile()
      gaussianRoute = '" -X="#p test hf 3-21g geom=allcheck freq"'
      commandLine = 'gdv -im="001a.faf" -om="'//TRIM(matrixFilename)//  &
        TRIM(gaussianRoute)//' < /dev/null > 001a.log'

      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Hrant - commandLine: ',TRIM(commandLine)
      write(iOut,*)
      write(iOut,*)

      call execute_command_line(TRIM(commandLine),exitStat=systemStatus)
      write(iOut,*)' Hrant - systemStatus = ',systemStatus

      if(systemStatus.ne.0) goto 999

!
!     Load the FAF from the job I just ran, pull data from it, and print the
!     data back to the calling Gaussian job output file.
!
      do i = 1,3
        call GMatrixFileOut%load(TRIM(matrixFilename))
        nBasis    = GMatrixFileOut%getVal('nBasis')
        nBasisUse = GMatrixFileOut%getVal('nBasisUse')
        energySCF = GMatrixFileOut%getValReal('scfEnergy')
        call GMatrixFileOut%getArray('nuclear gradient',mqcVarOut=gradient)
        gradient = RESHAPE(gradient,[3*nAtoms])
        write(iOut,*)
        write(iOut,*)
        write(iOut,*)' Hrant - rank of gradient = ',gradient%getRank()
        write(iOut,*)
        write(iOut,*)
        call GMatrixFileOut%getArray('nuclear force constants',mqcVarOut=hessian)
        write(iOut,2000)
        write(iOut,2100) 'nBasis    =',nBasis
        write(iOut,2100) 'nBasisUse =',nBasisUse
        write(iOut,2150) 'SCF Energy=',energySCF
        call gradient%print(iOut,header='gradient')
        call hessian%print(iOut,header='hessian')
        write(iOut,*)
!
!       Take a NR step and update the cartesian coordinates on the FAF.
!
        deltaX = newtonRaphsonStep(gradient,hessian)
        call mqc_print(deltaX,iOut,header='deltaX')
        GMatrixFileOut%cartesians = GMatrixFileOut%cartesians + deltaX
        call mqc_print(GMatrixFileOut%cartesians,iOut,header='Updated Cartesians')
!
!       Close the output FAF.
!
        call GMatrixFileOut%closeFile()
        call GMatrixFileOut%create('001b.faf')
        call GMatrixFileOut%closeFile()
      endDo
!
!     Jump to the end of the program.
!
      goto 999
!
!
  999 continue
      end program external
