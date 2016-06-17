      SUBROUTINE UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     1  PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     2  DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     3  PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT,
     4  JPROPS, NJPROP, PERIOD)
*
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (ZERO= 0.D0, HALF= 0.5D0, ONE=1.D0 )
*
      DIMENSION RHS(MLVARX, *), AMATRX(NDOFEL, NDOFEL),
     1  SVARS(NSVARS), ENERGY(8), PROPS(*), COORDS(MCRD, NNODE),
     2  U(NDOFEL), DU(MLVARX, *), V(NDOFEL), A(NDOFEL), TIME(2),
     3  PARAMS(3), JDLTYP(MDLOAD,*), ADLMAG(MDLOAD, *),
     4  DDLMAG(MDLOAD, *), PREDEF(2, NPREDF, NNODE), LFLAGS(*),
     5  JPROPS(*)
*
*     Redirect

      select case(JTYPE)
      case (201)
        CALL K_LINEAR_SPRING(
     1  RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     2  PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     3  DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     4  PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT,
     5  JPROPS, NJPROP, PERIOD)

      case (212)
        call K_BILINEAR_SPRING(
     1  RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     2  PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     3  DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     4  PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT,
     5  JPROPS, NJPROP, PERIOD)

      case default
        write(*,*) 'Unknown user element id of ', JTYPE
        call xit
      end select

      return
      end

      subroutine UEXTERNALDB(LOP,LRESTART,CTIME,DTIME,KSTEP,KINC)
      implicit none
C  Operation
C  0: Start of the analysis
C  1: start of the current analysis increment.
C  2: end of the analysis increment
C  3: end of the analysis
C  4: beginning of the restart analysis
      INTEGER LOP
C  indicate restart condition
C  0:not restart, 1:writing now, 2:write with overwrite
      INTEGER LRESTART
C  Time increment
      REAL*8 DTIME
C Number of step
      INTEGER KSTEP
C Number of increment
      INTEGER KINC
C  Current Time
      REAL*8 CTIME(2)
      call K_UEXTERNALDB(LOP,LRESTART,CTIME,DTIME,KSTEP,KINC)
      end

