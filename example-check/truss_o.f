      subroutine UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     1  PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     2  DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     3  PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT,
     4  JPROPS, NJPROP, PERIOD)
*
      INCLUDE 'ABA_PARAM.INC'
      parameter (ZERO= 0.D0, HALF= 0.5D0, ONE=1.D0 )
*
      dimension RHS(MLVARX, *), AMATRX(NDOFEL, NDOFEL),
     1  SVARS(NSVARS), ENERGY(8), PROPS(*), COORDS(MCRD, NNODE),
     2  U(NDOFEL), DU(MLVARX, *), V(NDOFEL), A(NDOFEL), TIME(2),
     3  PARAMS(3), JDLTYP(MDLOAD,*), ADLMAG(MDLOAD, *),
     4  DDLMAG(MDLOAD, *), PREDEF(2, NPREDF, NNODE), LFLAGS(*),
     5  JPROPS(*)
*
      dimension SRESID(6)
*
      real*8 Area, E, aNu, Rho, aLen, aK, aM
      integer k1, k2
*    horizontal truss element sample
*
*
      Area = props(1)  !< 断面積
      E    = props(2)  !< ヤング率
      aNu  = props(3)  !< ポアソン比
      rho  = props(4)  !< 密度
*
      aLen = abs(COORDS(1,2) - COORDS(1, 1))     !< 要素長さ
      aK   = Area * E / aLen                     !< 剛性
      aM   = Half * Area * Rho * ALen            !< 節点質量
*
      ! 初期化
      do K1 = 1,NDOFEL
        SRESID(K1) = ZERO
        do KRHS = 1,NRHS, 1
          RHS(K1,KRHS) = ZERO
        end do
        do K2 = 1,NDOFEL
          AMATRX(K2, K1) = ZERO
        end do
      end do
C
      if (LFLAGS(3).eq.1) then
       ! normal incrementation
       if (LFLAGS(1).eq.1 .or. LFLAGS(1).eq.2) then
         !static
         AMATRX(1,1) = AK
         AMATRX(4,4) = AK
         AMATRX(1,4) = -AK
         AMATRX(4,1) = -AK
         if (LFLAGS(4).ne.0) then
           ! 線形摂動ステップ
           force = ak * (u(4) -u(1))
           dforce = ak * (du(4,1) - du(1,1))
           SRESID(1) = -dforce
           SRESID(4) = dforce
           ! 残差の更新
           rhs(1,1) = rhs(1,1) - SRESID(1)
           rhs(4,1) = rhs(4,1) - SRESID(4)
           ! エネルギーの計算だがよくわからない．．．
           ENERGY(2) = half*force*(du(4,1)-du(1,1))
     1               + half*dforce*(U(4)-u(1))
     2               + half*dforce*(du(4,1)-du(1,1))
         else
           ! 一般解析ステップ
           force = aK * (U(4)-U(1))
           sresid(1) = -force
           rhs(1,1) = rhs(1,1) - sresid(1)
           rhs(4,1) = rhs(4,1) - sresid(4)
           do kDload = 1, ndload
             if(JDLTYP(kDload,1).eq.1001) then
               rhs(4,1) = rhs(4, 1) + adlmag(kdload,1)
               ENERGY(8) = energy(8) + (ADLMAG(kDload,1)
     1                   - half*DDLMAG(kdload,1))*du(4,1)
               if (NRHS.eq.2) then
                 ! riks
                 RHS(4,2) = rhs(4,2)+ddlmag(kDload,1)
               end if
             end if
           enddo
           ENERGY(2) = half * force * (u(4)-u(1))
         end if
       elseif (LFLAGS(1).eq.11 .or. LFLAGS(1).eq.12) then
         ! dynamic
         alpha = PARAMS(1)
         beta  = PARAMS(2)
         gamma = PARAMS(3)
         !
         dAdU = one/(beta*dtime**2)
         dVdU = gamma/(beta*dtime)
         !
         do k1 = 1,NDOFEL
           AMATRX(k1,k1) = aM * dAdU
           RHS(K1,1) = rhs(K1,1) - AM*A(K1)
         end do
         !
         AMATRX(1,1) = AMATRX(1,1)+(one+alpha)*aK
         AMATRX(4,4) = AMATRX(4,4)+(one+alpha)*aK
         AMATRX(1,4) = AMATRX(1,4)-(one+alpha)*aK
         AMATRX(4,1) = AMATRX(4,1)-(one+alpha)*aK
         force = aK * (u(4)-U(1))
         sresid(1) = -force
         sresid(4) = force
         rhs(1,1) = rhs(1,1) - ((one+alpha)*sresid(1)-alpha*svars(1))
         rhs(4,1)=rhs(4,1)-((one+alpha)*sresid(4)-alpha*SVARS(4))
         energy(1) = zero
         ! svarsのバックアップとエネルギーの計算を同時に行っている．
         do k1 = 1,Ndofel
           svars(k1+6) = svars(k1)
           svars(k1) = sresid(k1)
           energy(1) = energy(1)+half*v(k1)*am*v(k1)
         end do
         ENERGY(2) = half*force*(U(4)-u(1))
       end if
      elseif (LFLAGS(3).eq.2) then
        !剛性マトリックスのみを定義
        AMATRX(1,1) =  aK
        AMATRX(4,4) =  aK
        AMATRX(1,4) = -aK
        AMATRX(4,1) = -aK
      elseif (LFLAGS(3).eq.3) then
        !減衰マトリックスのみを定義
        ! Zero
      elseif (LFLAGS(3).eq.4) then
        !質量マトリックスのみを定義
        do k1 = 1,NDOFEL, 1
          AMATRX(k1,k1) = aM
        end do
      elseif (LFLAGS(3).eq.5) then
        ! ハーフインクリメントの残差
        alpha = PARAMS(1)
        force = aK * (u(4)-u(1))
        sresid(1) = -force
        sresid(4) = force
        rhs(1,1) = rhs(1,1) - am*a(1)-(one + alpha) * sresid(1) + half * 
     &    alpha * (svars(1) + svars(1+6))
        rhs(4,1) = rhs(4,1) - am*a(4)-(one + alpha) * sresid(4) + half *
     &    alpha * (svars(4) + svars(4+6))
      elseif (LFLAGS(3).eq.6) then
        ! 初期加速度の計算
        do k1 = 1,ndofel
          AMATRX(k1,k1) = aM
        end do
        force = ak * (U(4)-U(1))
        sresid(1) = -force
        sresid(4) = force
        rhs(1,1) = rhs(1,1) - sresid(1)
        rhs(4,1) = rhs(4,1) - sresid(4)
        ENERGY(1) = zero
        do k1 = 1, NDOFEL
          svars(k1) = sresid(k1)
          ENERGY(1) = ENERGY(1) + half*am*v(k1)**2
        end do
        ENERGY(2) = half * force*(u(4)-u(1))
      elseif (LFLAGS(3).eq.100) then
        ! 出力用摂動量の計算
        if (lflags(1).eq.1 .or. LFLAGS(1).eq.2) then
          ! 静的解析
          force = ak*(U(4)-U(1))
          dforce= ak*(du(4,1)-du(1,1))
          sresid(1) = -dforce
          sresid(4) = dforce
          rhs(1,1) = rhs(1,1)-sresid(1)
          rhs(4,1) = rhs(4,1)-sresid(4)
          ENERGY(2) = half*force*(du(4,1)-du(1,1))
     1                + half * dforce*(u(4)-u(1))
     2                + half * dforce*(du(4,1)-du(1,1))
          do kvar = 1, nsvars
            svars(kvar) = zero
          end do
          svars(1) = rhs(1,1)
          svars(4) = rhs(4,1)
        elseif (LFLAGS(1).eq.41) then
          ! 固有値解析
          do krhs = 1, nrhs
            dforce = aK * (du(4,krhs)-du(1,krhs))
            sresid(1) = -dforce
            sresid(4) = dforce
            rhs(1,krhs) = rhs(1,kRHS)-sresid(1)
            rhs(4,krhs) = rhs(4,kRHS)-sresid(4)
          end do
          do kvar = 1, NSVARS
            svars(kVar) = zero
          enddo
          svars(1) = rhs(1,1)
          svars(4) = rhs(4,1)
        endif
      end if
!
      return
      end
