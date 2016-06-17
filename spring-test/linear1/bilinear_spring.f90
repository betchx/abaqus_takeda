!
!  2 nodes Symmetric Bi-Linear Spring element with 6 DOFs.
!
!  No mass were considered
!
!  Properties: (24)
!    1: spring constant for dof 1 (translation x)
!    2: spring constant for dof 2 (translation y)
!    3: spring constant for dof 3 (translation z)
!    4: spring constant for dof 4 (rotation x)
!    5: spring constant for dof 5 (rotation y)
!    6: spring constant for dof 6 (rotation z)
!    7: Dummy (filler)
!    8: Dummy (filler)
!    9: Yield Force for dof 1
!   10: Yield Force for dof 2
!   11: Yield Force for dof 3
!   12: Yield Moment for dof 4
!   13: Yield Moment for dof 5
!   14: Yield Moment for dof 6
!   15: Dummy (filler)
!   16: Dummy (filler)
!   17: Second slope for dof 1
!   18: Second slope for dof 2
!   19: Second slope for dof 3
!   20: Second slope for dof 4
!   21: Second slope for dof 5
!   22: Second slope for dof 6
!   23: Dummy (filler)
!   24: Dummy (filler)

! SVARS
!   * used for intermediate residual
!   1 to 12 SRESID at end
!  13 to 24 SRESID at start


module BilinerSpring
  implicit none
  real, parameter :: one=1.0d0, half=0.5d0, zero=0.0d0
contains

  SUBROUTINE abort(msg)
    use debug, only: terminate_debug
    implicit none
    character(*) msg
    write(*,*) msg
    call terminate_debug()
    call xit
  end SUBROUTINE abort


  subroutine updateStaticLP(AMATRX,RHS,ENERGY,SRESID,PROPS,U,DU)
    use debug
    implicit none
    real, intent(out) :: AMATRX(12,12)
    real, intent(out) :: RHS(12,*)
    real, intent(out) :: ENERGY(8)
    real, intent(out) :: SRESID(12)
    real, intent(in)  :: PROPS(24)
    real, intent(in)  :: U(12)
    real, intent(in)  :: DU(12,*)

    real, dimension(6) :: forces, dforces, dL, ddL

    call updateKMat(AMATRX, PROPS)

    dL = U(7:12) - U(1:6)
    forces = PROPS(1:6)*dL
    ddL = DU(7:12,1) - DU(1:6,1)
    dforces = PROPS(1:6)*ddL
    SRESID(1:6) = -dforces(:)
    SRESID(7:12) = dforces(:)
    RHS(:,1) = RHS(:,1)-SRESID
    ENERGY(2) = half*sum(forces*ddL+dforces*dL+dforces*ddL)
  end subroutine updateStaticLP

  subroutine updateStatic(AMATRX,RHS,ENERGY,SVARS,SRESID,PROPS,U)
    use debug
    implicit none
    real, intent(out) :: AMATRX(12,12)
    real, intent(out) :: RHS(12,*)
    real, intent(out) :: ENERGY(8)
    real, intent(out) :: SRESID(12)
    real, intent(inout)  :: SVARS(*)
    real, intent(in)  :: PROPS(*)
    real, intent(in)  :: U(12)

    real, dimension(6) :: forces, kk, ee
    character(1000):: tag
    integer :: i

    ! Yield check
    call YieldCheck(forces, kk, ee, PROPS, U)
    call updateKMat(AMATRX, KK)

    SRESID(1:6) = -forces
    SRESID(7:12) = forces

!    tag = '<RHS1'
!    do i = 1, 12
!      tag = trim(tag) // attr('i'//to_s(i),RHS(i,1))
!    enddo
!    call trace(trim(tag) // '/>')

    RHS(:,1) = RHS(:,1)-SRESID

!    tag = '<RHS2'
!    do i = 1, 12
!      tag = trim(tag) // attr('i'//to_s(i),RHS(i,1))
!    enddo
!    call trace(trim(tag) // '/>')

    !No distributed load can be considered
    ENERGY(2) = sum(ee)
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
    real, intent(in)  :: PROPS(24)
    real, intent(in)  :: U(12)
    real, intent(in)  :: V(12)
    real, intent(in)  :: A(12)
    real, intent(in)  :: DTIME

    real, dimension(6) :: forces, kk, ee
    real :: alpha, beta, gamma, dAdU, dVdU, val
    integer i, k

    alpha = PARAMS(1)
    beta = PARAMS(2)
    gamma = PARAMS(3)

    dAdU = ONE/(BETA*DTIME**2)
    dVdU = gamma / (beta*DTIME)

    call YieldCheck(forces, kk, ee, PROPS, U)

    ! No operation about mass
    do i=1,6
      k= 1 + 6
      val = (ONE+alpha)*kk(i)
      AMATRX(i,i) = AMATRX(i,i) + val
      AMATRX(k,k) = AMATRX(k,k) + val
      AMATRX(i,k) = AMATRX(i,k) - val
      AMATRX(k,i) = AMATRX(k,i) - val
    end do

    SRESID(1:6) = -forces
    SRESID(7:12) = forces
    RHS(:,1) = RHS(:,1)-((ONE+alpha)*SRESID-alpha*SVARS(1:12))
    SVARS(13:24) = SVARS(1:12)
    SVARS(1:12) = SRESID
    ENERGY(1) = zero
    ENERGY(2) = sum(ee)
  end subroutine updateDynamic

  subroutine YieldCheck(F, KK, ee, PROPS, U)
    use debug
    implicit none
    real, intent(in)  :: PROPS(24)
    real, intent(in)  :: U(6,2)
    real, intent(out) :: F(6),KK(6),ee(6)
    integer :: i
    real :: dy  ! delta_y
    real :: u2, p2
    real, dimension(6) :: k1, yp, k2, dL

!    call enter('YieldCheck')

    k1 = PROPS(1:6)
    yp = PROPS(9:14)
    k2 = PROPS(17:22)

    dL = U(:,2) - U(:,1)
    KK = K1
    F = K1 * dL
    do i=1, 6
      if (abs(F(i)).gt.yp(i)) then
        dy = sign(yp(i)/K1(i), F(i))
        u2 = dL(i) - dy
        p2 = k2(i)*u2
        F(i) = sign(yp(i), F(i)) + p2
        KK(i) = k2(i)
        ee(i) = half*(yp(i)*abs(dy)+p2*u2)+yp(i)*abs(u2)
        if(i.eq.1) &
          call trace('<DOF'//to_s(i)//attr('F',F(i))//attr('KK',kk(i))//attr('ee',ee(i))//&
          attr('u1',u(i,2))//attr('u0',u(i,1))//attr('dL',dL(i))//&
          attr('dy',dy)//attr('u2',u2)// attr('p2',p2)//'/>')
      else
        ee(i) = half * F(i) * dL(i)
        if(i.eq.1) &
          call trace('<DOF'//to_s(i)//attr('F',F(i))//attr('KK',kk(i))//attr('ee',ee(i))//&
          attr('u1',u(i,2))//attr('u0',u(i,1))//attr('dL',dL(i))//&
          '/>')
      end if
    end do
!    call leave('YieldCheck')
  end subroutine YieldCheck

  subroutine checkAndUpdateKmat(AMATRX, PROPS,U)
    implicit none
    real, intent(out) :: AMATRX(12, 12)
    real, intent(in) :: PROPS(*), U(*)

    real, dimension(6) :: f,k,e

    call YieldCheck(f, k, e, PROPS, u)
    call updateKMat(AMATRX, k)

  end subroutine checkAndUpdateKmat

  subroutine updateKmat(AMATRX, KK)
    implicit none
    real :: AMATRX(12, 12)   !< 結果のマトリックス
    real :: KK(6)        !< バネ剛性（6自由度)

    integer :: i, k

    do i=1, 6
    k= i+6
    AMATRX(i, i) =  KK(i)
    AMATRX(k, k) =  KK(i)
    AMATRX(i, k) = -KK(i)
    AMATRX(k, i) = -KK(i)
    end do
  end subroutine

  subroutine updateMMat(AMATRX, PROPS)
    implicit none
    real, intent(out) :: AMATRX(12,12)
    real, intent(in)  :: PROPS(24)
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
    real, intent(in) :: PROPS(24)

    real :: alpha, forces(6)
    real, dimension(6) :: kk, ee
    alpha = PARAMS(1)
    !forces = k1*(U(7:12)-U(1:6))
    call YieldCheck(forces, kk, ee, PROPS, U)
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
    real, intent(in) :: PROPS(24)
    real, intent(in) :: U(12)

    real, dimension(6) :: forces, kk, ee
    call updateMMat(AMATRX, PROPS)
    call YieldCheck(forces, kk, ee, PROPS, U)
    SRESID(1:6) = -forces
    SRESID(7:12) = forces
    RHS(:,1) = RHS(:,1)-SRESID(:)
    SVARS(1:12) = SRESID(1:12)
    ENERGY(1) = 0.0  ! because of no mass.
    ENERGY(2) = sum(ee)
  end subroutine calcInitAcc

  subroutine outputForStaticLP(RHS, ENERGY, SVARS, SRESID, PROPS, U, DU)
    implicit none
    real, intent(out) :: RHS(12,*)
    real, intent(out) :: SVARS(12)
    real, intent(out) :: SRESID(12)
    real, intent(out) :: ENERGY(8)
    real, intent(in) :: PROPS(24)
    real, intent(in) :: U(12)
    real, intent(in) :: DU(12)

    real, dimension(6) :: forces, dforces, dL, ddL, ff, kk, ee

    call YieldCheck(ff, kk, ee, PROPS, U)

    dL = U(7:12)-U(1:6)
    ddL = DU(7:12)-DU(1:6)

    forces = kk*dL
    dforces = kk*ddL
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
    real, intent(in)  :: PROPS(24)
    real, intent(in)  :: DU(12,*)
    integer, intent(in)  :: NRHS

    integer i
    real, dimension(6) :: dforces, kk, ee


    do i = 1, NRHS
    call YieldCheck(dforces, kk, ee, PROPS, DU(:,i))
    SRESID(1:6) = -dforces
    SRESID(7:12) = dforces
    RHS(:,i) = RHS(:,i) - SRESID(:)
    end do
    SVARS(:) = zero
    SVARS(1:12) = RHS(1:12,1)
  end subroutine outputForFrequency

end module BilinerSpring


!> Interface
subroutine K_BILINEAR_SPRING(&
    RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
    PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME, &
    DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, &
    PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, &
    JPROPS, NJPROP, PERIOD)
  use BilinerSpring
  use debug
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
  character(1000) :: tag
  integer :: i

  ! argument check
  if(JTYPE.ne.212) call abort('UEL bugs about user element type id')
  if(NNODE.ne.2) call abort('U202 must be 2 nodes')
  if(NDOFEL.ne.12) call abort('nodes of U202 must be three dimension')
  if(NDLOAD.gt.0) call abort('Distributed load on spring is not allowed')
  if(NPROPS.ne.24) call abort('U212 requires 24 float parameters')

  !tag = '<RHS0'
  !do i = 1, 12
  !  tag = trim(tag) // attr('i'//to_s(i),RHS(i,1))
  !enddo
  !call trace(trim(tag) // '/>')

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
        call enter('updateStaticLP')
        call updateStaticLP(AMATRX,RHS,ENERGY,SRESID,PROPS,U,DU)
        call leave('updateStaticLP')
      else
        !call enter('updateStatic')
        call updateStatic(AMATRX,RHS,ENERGY,SVARS,SRESID,PROPS,U)
        !call leave('updateStatic')
      endif
    case (11:12)
      !Implicit Dynamic
      call enter('updateDynamic',KINC)
      call updateDynamic(AMATRX, RHS, ENERGY, SVARS, SRESID, &
        PARAMS, PROPS, U, V, A, DTIME)
      call leave('updateDynamic',KINC)
    case default
      call abort('This Analysis type is not supported yet.')
    end select
  case(2)
    call enter('checkAndUpdateKMat')
    call checkAndUpdateKMat(AMATRX,PROPS,U)
    call leave('checkAndUpdateKMat')
  case(3)
    ! no dumping
    call enter('updateDumping')
    call leave('updateDumping')
  case(4)
    ! no mass
    call enter('updateMMat')
    call updateMMat(AMATRX, PROPS)
    call leave('updateMMat')
  case(5)
    call enter('updateResidual')
    call updateResidual(RHS, SVARS, SRESID, PARAMS,PROPS, U)
    call leave('updateResidual')
  case(6)
    call enter('calcInitAcc')
    call calcInitAcc(AMATRX, RHS, SVARS, ENERGY, SRESID, PROPS,U)
    call leave('calcInitAcc')
  case(100)
    !output
    select case(LFLAGS(1))
    case(1:2)
      call enter('outputForStaticLP')
      call outputForStaticLP(RHS, ENERGY, SVARS, SRESID, PROPS, U, DU)
      call leave('outputForStaticLP')
    case(41)
      call enter('outputForFrequency')
      call outputForFrequency(RHS, SVARS, SRESID, PROPS, DU, NRHS)
      call leave('outputForFrequency')
    end select
  end select

end subroutine K_BILINEAR_SPRING

