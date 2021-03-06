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
        CALL K_LINEAR_SPRING( RHS,
     1  AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
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

