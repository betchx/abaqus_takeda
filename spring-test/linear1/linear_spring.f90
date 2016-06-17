!
!  2 nodes Linear Spring element with 6 DOFs.
!
!  No mass were considered
!
!  Properties:
!    1: spring constant for dof 1 (translation x)
!    2: spring constant for dof 2 (translation y)
!    3: spring constant for dof 3 (translation z)
!    4: spring constant for dof 4 (rotation x)
!    5: spring constant for dof 5 (rotation y)
!    6: spring constant for dof 6 (rotation z)
!

module LinerSpring
  implicit none
  real, parameter :: one=1.0d0, half=0.5d0, zero=0.0d0
contains


  SUBROUTINE abort(msg)
    implicit none
    character(*) msg
    write(*,*) msg
    call xit
  end SUBROUTINE abort


  subroutine updateStaticLP(AMATRX,RHS,ENERGY,SRESID,PROPS,U,DU)
    implicit none
    real, intent(out) :: AMATRX(12,12)
    real, intent(out) :: RHS(12,*)
    real, intent(out) :: ENERGY(8)
    real, intent(out) :: SRESID(12)
    real, intent(in)  :: PROPS(6)
    real, intent(in)  :: U(12)
    real, intent(in)  :: DU(12,*)

    real, dimension(6) :: forces, dforces, dL, ddL

    call updateKMat(AMATRX, PROPS)

    dL = U(7:12) - U(1:6)
    forces = PROPS*dL
    ddL = DU(7:12,1) - DU(1:6,1)
    dforces = PROPS*ddL
    SRESID(1:6) = -dforces(:)
    SRESID(7:12) = dforces(:)
    RHS(:,1) = RHS(:,1)-SRESID
    ENERGY(2) = half*sum(forces*ddL+dforces*dL+dforces*ddL)
  end subroutine updateStaticLP

  subroutine updateStatic(AMATRX,RHS,ENERGY,SRESID,PROPS,U)
    implicit none
    real, intent(out) :: AMATRX(12,12)
    real, intent(out) :: RHS(12,*)
    real, intent(out) :: ENERGY(8)
    real, intent(out) :: SRESID(12)
    real, intent(in)  :: PROPS(6)
    real, intent(in)  :: U(12)

    real, dimension(6) :: forces, dL

    call updateKMat(AMATRX, PROPS)

    dL = U(7:12) - U(1:6)
    forces = PROPS * dL
    SRESID(1:6) = -forces
    SRESID(7:12) = forces
    RHS(:,1) = RHS(:,1)-SRESID
    !No distributed load can be considered
    ENERGY(2) = half * sum(forces*dL)
  end subroutine updateStatic


  subroutine updateDynamic(AMATRX, RHS, ENERGY, SVARS, SRESID,&
      PARAMS, PROPS, U, V, A, DTIME)
    implicit none
    real, intent(out) :: AMATRX(12,12)
    real, intent(out) :: RHS(12,*)
    real, intent(out) :: ENERGY(8)
    real, intent(inout) :: SVARS(24)
    real, intent(out) :: SRESID(12)
    real, intent(in)  :: PARAMS(3)
    real, intent(in)  :: PROPS(6)
    real, intent(in)  :: U(12)
    real, intent(in)  :: V(12)
    real, intent(in)  :: A(12)
    real, intent(in)  :: DTIME

    real, dimension(6) :: forces, dL
    real :: alpha, beta, gamma, dAdU, dVdU, val
    integer i, k

    alpha = PARAMS(1)
    beta = PARAMS(2)
    gamma = PARAMS(3)

    dAdU = ONE/(BETA*DTIME**2)
    dVdU = gamma / (beta*DTIME)

    ! No operation about mass
    do i=1,6
      k= 1 + 6
      val = (ONE+alpha)*PROPS(i)
      AMATRX(i,i) = AMATRX(i,i) + val
      AMATRX(k,k) = AMATRX(k,k) + val
      AMATRX(i,k) = AMATRX(i,k) - val
      AMATRX(k,i) = AMATRX(k,i) - val
    end do

    dL = U(7:12) - U(1:6)
    forces = PROPS*dL
    SRESID(1:6) = -forces
    SRESID(7:12) = forces
    RHS(:,1) = RHS(:,1)-((ONE+alpha)*SRESID-alpha*SVARS(1:12))
    SVARS(13:24) = SVARS(1:12)
    SVARS(1:12) = SRESID
    ENERGY(1) = zero
    ENERGY(2) = HALF*sum(forces*dL)
  end subroutine updateDynamic


  subroutine updateKmat(AMATRX, PROPS)
    implicit none
    real :: AMATRX(12, 12)   !< 結果のマトリックス
    real :: PROPS(6)        !< バネ剛性（6自由度)

    integer :: i, k

    do i=1, 6
      k= i+6
      AMATRX(i, i) =  PROPS(i)
      AMATRX(k, k) =  PROPS(i)
      AMATRX(i, k) = -PROPS(i)
      AMATRX(k, i) = -PROPS(i)
    end do
  end subroutine

  subroutine updateMMat(AMATRX, PROPS)
    implicit none
    real, intent(out) :: AMATRX(12,12)
    real, intent(in)  :: PROPS(6)
    ! no mass.
    AMATRX(:,:) = zero
  end subroutine updateMMat

  subroutine updateResidual(RHS, SVARS, SRESID, PARAMS,PROPS, U)
    implicit none
    real, intent(out) :: RHS(12,*)
    real, intent(in) :: SVARS(2*12)
    real, intent(out) :: SRESID(12)
    real, intent(in) :: PARAMS(3)
    real, intent(in) :: U(12)
    real, intent(in) :: PROPS(*)

    real :: alpha, forces(6)
    alpha = PARAMS(1)
    forces = PROPS(1:6)*(U(7:12)-U(1:6))
    SRESID(1:6) = -forces
    SRESID(7:12) = forces
    RHS(:,1) = RHS(:,1)+half*alpha*(SVARS(1:12)+SVARS(13:24))
  end subroutine

  subroutine calcInitAcc(AMATRX, RHS, SVARS, ENERGY, SRESID, PROPS,U)
    implicit none
    real, intent(out) :: AMATRX(12,12)
    real, intent(out) :: RHS(12,*)
    real, intent(out) :: SVARS(2*12)
    real, intent(out) :: SRESID(12)
    real, intent(out) :: ENERGY(8)
    real, intent(in) :: PROPS(*)
    real, intent(in) :: U(12)

    real forces(6), du(6)
    call updateMMat(AMATRX, PROPS)
    du = (U(7:12)-U(1:6))
    forces = PROPS(1:6)* du
    SRESID(1:6) = -forces
    SRESID(7:12) = forces
    RHS(:,1) = RHS(:,1)-SRESID(:)
    SVARS(1:12) = SRESID(1:12)
    ENERGY(1) = 0.0  ! because of no mass.
    ENERGY(2) = sum(HALF*forces*du)
  end subroutine calcInitAcc

  subroutine outputForStaticLP(RHS, ENERGY, SVARS, SRESID, PROPS, U, DU)
    implicit none
    real, intent(out) :: RHS(12,*)
    real, intent(out) :: SVARS(12)
    real, intent(out) :: SRESID(12)
    real, intent(out) :: ENERGY(8)
    real, intent(in) :: PROPS(*)
    real, intent(in) :: U(12)
    real, intent(in) :: DU(12)

    real, dimension(6) :: forces, dforces, dL, ddL

    dL = U(7:12)-U(1:6)
    ddL = DU(7:12)-DU(1:6)

    forces = PROPS(1:6)*dL
    dforces = PROPS(1:6)*ddL
    SRESID(1:6) = -dforces
    SRESID(7:12) = dforces
    RHS(:,1) = RHS(:,1)-SRESID(:)
    ENERGY(2)=half*sum(forces*ddL+dforces*dL+dforces*ddL)
    SVARS(:) = zero
    SVARS(1:12) = RHS(:,1)
  end subroutine outputForStaticLP

  subroutine outputForFrequency(RHS, SVARS, SRESID, PROPS, DU, NRHS)
    implicit none
    real, intent(out) :: RHS(12,*)
    real, intent(out) :: SVARS(12*2)
    real, intent(out) :: SRESID(12)
    real, intent(in)  :: PROPS(6)
    real, intent(in)  :: DU(12,*)
    integer, intent(in)  :: NRHS

    integer i
    real, dimension(6) :: dforces

    do i = 1, NRHS
      dforces = PROPS*(DU(7:12,i)-DU(1:6,i))
      SRESID(1:6) = -dforces
      SRESID(7:12) = dforces
      RHS(:,i) = RHS(:,i) - SRESID(:)
    end do
    SVARS(:) = zero
    SVARS(1:12) = RHS(1:12,1)
  end subroutine outputForFrequency

end module


!> Interface
subroutine K_LINEAR_SPRING(&
    RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
    PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME, &
    DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, &
    PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, &
    JPROPS, NJPROP, PERIOD)
  use LinerSpring
  implicit none
  ! need to update
  real, intent(out) :: RHS(MLVARX, NRHS)
  real, intent(out) :: AMATRX(NDOFEL, NDOFEL)
  real, intent(inout) :: SVARS(NSVARS)
  real, intent(out) :: ENERGY(8)
  ! updatable
  real, intent(inout) :: PNEWDT
  ! Refference
  real, intent(in)  :: PROPS(NPROPS)  !< Real properties
  real, intent(in)  :: JPROPS(NJPROP) !< Integer properties
  real, intent(in)  :: COORDS(MCRD, NNODE) !< Initial coordinate of nodes
  real, intent(in)  :: U(NDOFEL)   !< Displacement
  real, intent(in)  :: DU(MLVARX, NDOFEL) !<Incrementation of displacement
  real, intent(in)  :: V(NDOFEL)   !< Velocity by time
  real, intent(in)  :: A(NDOFEL)   !< Acceleration
  integer, intent(in)  :: JDLTYP(MDLOAD,*)!< type of distributed load
  real, intent(in)  :: ADLMAG(MDLOAD, *) !< amplitude of distributed load
  real, intent(in)  :: DDLMAG(MDLOAD, *) !< difference of distributed load
  real, intent(in)  :: PREDEF(2, NPREDF, NNODE) !< Pre Defined Field
  !> Array of Parameters
  !>  In case of implicit dynamic
  !>  Params(1): integration operator alpha
  !>  Params(2): integration operator beta
  !>  Params(3): integration operator gamma
  real, intent(in)  :: PARAMS(3)
  !> Array of flags
  !>  1: procesure type
  !>  2: flag of NLGEOM
  !>  3: porpose of function call
  !>  4: flag of Linear Perturbation
  !>  5: estimation method (0: Newton method, 1:extrapolation)
  Integer, intent(in)  :: LFLAGS(5)
  real, intent(in)  :: TIME(2)  !< (1): step time,  (2):Total time
  real, intent(in)  :: DTIME !< time Incrementation
  real, intent(in)  :: PERIOD !< time width of current step
  integer, intent(in)  :: NDOFEL !< total number of DOF in the element
  integer, intent(in)  :: MLVARX !< number of desp. vectors
  integer, intent(in)  :: NRHS   !< number of lod vrt
  integer, intent(in)  :: NSVARS !< number of svars
  integer, intent(in)  :: NPROPS !< number of PROPS
  integer, intent(in)  :: NJPROP !< number of JPROPS
  integer, intent(in)  :: MCRD   !< Number of coordinate of node
  integer, intent(in)  :: NNODE  !< Number of nodes in the element
  integer, intent(in)  :: JTYPE  !< User Element type ID  (xx of type=Uxx)
  integer, intent(in)  :: KSTEP  !< Current step number
  integer, intent(in)  :: KINC   !< Current increment number
  integer, intent(in)  :: JELEM  !< Element ID
  integer, intent(in)  :: NDLOAD !< Number of distributed load
  integer, intent(in)  :: MDLOAD !< Total number of distributed load
  integer, intent(in)  :: NPREDF !< Number of pre-defined field variables
  !!!
  !Local variables
  real, dimension(12) :: SRESID

  ! argument check
  if(JTYPE.ne.201) call abort('UEL bugs about user element type id')
  if(NNODE.ne.2) call abort('U201 must be 2 nodes')
  if(NDOFEL.ne.12) call abort('nodes of U201 must be three dimension')
  if(NDLOAD.gt.0) call abort('Distributed load on spring is not allowed')

  ! clear
  AMATRX(:,:) = ZERO
  SRESID(:) = ZERO
  RHS(:,:) = ZERO

  !selection by flags
  select case(LFLAGS(3))
  case(1) ! Normal incrementation
    select case(LFLAGS(1))
    case (1:2)
      !Static
      if(LFLAGS(4).eq.1) then
        call updateStaticLP(AMATRX,RHS,ENERGY,SRESID,PROPS,U,DU)
      else
        call updateStatic(AMATRX,RHS,ENERGY,SRESID,PROPS,U)
      endif
    case (11:12)
      !Implicit Dynamic
      call updateDynamic(AMATRX, RHS, ENERGY, SVARS, SRESID, &
        PARAMS, PROPS, U, V, A, DTIME)
    case default
      call abort('This Analysis type is not supported yet.')
    end select
  case(2)
    call updateKMat(AMATRX,PROPS)
  case(3)
    ! no dumping
  case(4)
    ! no mass
    call updateMMat(AMATRX, PROPS)
  case(5)
    call updateResidual(RHS, SVARS, SRESID, PARAMS,PROPS, U)
  case(6)
    call calcInitAcc(AMATRX, RHS, SVARS, ENERGY, SRESID, PROPS,U)
  case(100)
    !output
    select case(LFLAGS(1))
    case(1:2)
      call outputForStaticLP(RHS, ENERGY, SVARS, SRESID, PROPS, U, DU)
    case(41)
      call outputForFrequency(RHS, SVARS, SRESID, PROPS, DU, NRHS)
    end select
  end select

end subroutine K_LINEAR_SPRING

