      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'


	DIMENSION TIME(2)

C
      !user coding to set up the FORTRAN environment, open files, close files, 
      !calculate user-defined model-independent history information,
      !write history information to external files,
      !recover history information during restart analyses, etc.
      !do not include calls to utility routine XIT
      call MutexInit( 1 )      ! initialize Mutex #1
      call MutexInit( 2 )
      call MutexInit( 3 )
      call MutexInit( 4 )
      call MutexInit( 5 )
      call MutexInit( 6 )
      call MutexInit( 7 )
      call MutexInit( 8 )
      call MutexInit( 9 )
      call MutexInit( 10 )
      call MutexInit( 11 )
      call MutexInit( 12 )
      call MutexInit( 13 )
      call MutexInit( 14 )
      call MutexInit( 15 )

      RETURN
      END
