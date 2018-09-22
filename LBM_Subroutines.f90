
SUBROUTINE Equilibrium_h_g_AMR( C, U, V, ni, nj, Gamma, heq )
	USE D2Q9_Mod, ONLY: Wa, W, M0, ex, ey
	IMPLICIT NONE
	REAL(8), INTENT(IN)  :: C, U, V, ni, nj
	REAL(8), INTENT(OUT) :: Gamma(9), heq(9)

	REAL(8) :: U2, eU(9), eN(9), tmp

	tmp = ( 12*(1-C) )/W * M0

	U2 = U*U + V*V

	eU(:) = ex(:) * U  + ey(:) * V
	eN(:) = ex(:) * ni + ey(:) * nj

	Gamma(:) = Wa(:) * ( 1.d0 + eU(:)*(3.d0 + 4.5d0*eU(:)) - 1.5d0*U2 )

	heq(:) = C * ( Gamma(:) + tmp * Wa(:) * eN(:) )

END

!**********************************************************************

SUBROUTINE normal_AMR( DcDx, DcDy, ni, nj )
	IMPLICIT NONE
	REAL(8), INTENT(IN)  :: DcDx, DcDy
	REAL(8), INTENT(OUT) :: ni, nj

	REAL(8) :: tmp

	tmp = DSQRT( DcDx**2.d0 + DcDy**2.d0 + 1.d-32 )

	ni = DcDx / tmp
	nj = DcDy / tmp

END

!**********************************************************************

SUBROUTINE Perfect_Shift_AMR( nx, ny, f )
	USE D2Q9_Mod, ONLY: ex, ey
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny
	REAL(8), INTENT(INOUT) :: f(9, 0:nx+1,0:ny+1)

	INTEGER :: i, X, Y
	REAL(8) :: fnew(2:9, nx,ny)

	DO Y = 1, ny
	DO X = 1, nx
	DO i = 2, 9

		fnew(i,X,Y) = f(i,X-ex(i),Y-ey(i))

	END DO
	END DO
	END DO

	f(2:9, 1:nx,1:ny) = fnew

END

!**********************************************************************

SUBROUTINE Streaming_LaxWendroff_AMR( nx, ny, cfl, f )
	USE D2Q9_Mod, ONLY: ex, ey
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny
	REAL(8), INTENT(IN)    :: cfl
	REAL(8), INTENT(INOUT) :: f(9, 0:nx+1,0:ny+1)

	INTEGER :: i, X, Y
	REAL(8) :: fnew(2:9, nx,ny)

	DO Y = 1, ny
	DO X = 1, nx
	DO i = 2, 9

		fnew(i,X,Y) = cfl/2*(1+cfl) * f(i,X-ex(i),Y-ey(i)) &
					+ (1 - cfl**2)  * f(i,X,Y) &
					- cfl/2*(1-cfl) * f(i,X+ex(i),Y+ey(i))

	END DO
	END DO
	END DO

	f(2:9, 1:nx,1:ny) = fnew

END

!**********************************************************************

SUBROUTINE Chemical_Potential_AMR( nx, ny, dx, C, mu )
	USE D2Q9_Mod, ONLY: Beta, kappa
	IMPLICIT NONE
	INTEGER, INTENT(IN)  :: nx, ny, dx
	REAL(8), INTENT(IN)  :: C(-1:nx+2,-1:ny+2)
	REAL(8), INTENT(OUT) :: mu(nx,ny)

	INTEGER :: i, j
	REAL(8) :: D2C(nx,ny)

	CALL Laplacian_AMR( nx, ny, C, D2C )

	DO j = 1, ny
	DO i = 1, nx
		mu(i,j) = 4.d0*Beta*C(i,j)*(C(i,j)-1.d0)*(C(i,j)-0.5d0) - kappa*D2C(i,j)/dx**2
	END DO
	END DO

END

!**********************************************************************

SUBROUTINE Collision_h_g_AMR( nx, ny, dx, C, P, Ux, Uy, h, g )
	USE D2Q9_Mod, ONLY: Wa, Gy, taul, tauh, Rhol, Rhoh, wc, ex, ey
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny, dx
	REAL(8), INTENT(IN)    :: C(-1:nx+2,-1:ny+2)
	REAL(8), INTENT(IN)    :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
	REAL(8), INTENT(INOUT) :: h(9, 0:nx+1,0:ny+1), g(9, 0:nx+1,0:ny+1)

	INTEGER :: i, j
	REAL(8) :: dRho3, Rho, tau, s9, ni, nj, Fb
	REAL(8) :: eDc(9),  hlp(9), Gamma(9), heq(9), geq(9)
	REAL(8) :: DcDx, DcDy, mu(nx,ny)

	dRho3 = (Rhoh - Rhol)/3.d0

	CALL Chemical_Potential_AMR( nx, ny, dx, C, mu )

	DO j = 1, ny
	DO i = 1, nx

		CALL ED_C_AMR( dx, C(i-1:i+1,j-1:j+1), eDc(:) )

		CALL Gradient_ED_AMR( eDc(:), DcDx, DcDy )

		CALL normal_AMR( DcDx, DcDy, ni, nj )

		CALL Equilibrium_h_g_AMR( C(i,j), Ux(i,j), Uy(i,j), &
									ni, nj, Gamma(:), heq(:) )

		!*******************	COLLISION (h)	*************
		h(:, i,j) = h(:, i,j) * (1.d0-wc) + heq(:)*wc
		!*******************	COLLISION (h)	*************

		!*******************	COLLISION (g)	*************
		Rho = Rhol + C(i,j) * (Rhoh - Rhol)

		Fb = (Rhoh - Rho) * Gy

		hlp(:) = ( (Gamma(:)-Wa(:))*dRho3 + Gamma(:)*mu(i,j) ) * ( eDc(:) &
				- (Ux(i,j)*DcDx + Uy(i,j)*DcDy) ) + (ey(:)-Uy(i,j)) * Gamma(:) * Fb

		geq(:) = P(i,j)*Wa(:) + (Gamma(:)-Wa(:))*Rho/3.d0 - 0.5d0 * hlp(:)

	!== mixed difference
		CALL ED_M_AMR( dx, C(i-2:i+2,j-2:j+2), eDc(:) )

		CALL Gradient_ED_AMR( eDc(:), DcDx, DcDy )

		hlp(:) = ( (Gamma(:)-Wa(:))*dRho3 + Gamma(:)*mu(i,j) ) * ( eDc(:) &
				- (Ux(i,j)*DcDx + Uy(i,j)*DcDy) ) + (ey(:)-Uy(i,j))*Gamma(:) * Fb

	!============	interpolating the relaxation time	==========
	!!	tau = taul + C(i,j) * (tauh - taul) + 0.5d0
		IF( C(i,j) < 0.d0 ) THEN
			tau = taul + 0.5d0
		ELSEIF( C(i,j) > 1.d0 ) THEN
			tau = tauh + 0.5d0
		ELSE
			tau = C(i,j) * (1.d0/tauh - 1.d0/taul) + 1.d0/taul
			tau = 1.d0/tau + 0.5d0
		END IF
		s9 = 1.d0/tau

		!============	BGK/MRT Collision Models	==========
	!	CALL BGK_Collision_AMR( s9, geq(:), g(:, i,j) )
		CALL MRT_Collision_AMR( s9, geq(:), g(:, i,j) )

		g(:, i,j) = g(:, i,j) + hlp(:)
		!*******************	COLLISION (g)	*************

	END DO
	END DO

END

!**********************************************************************

SUBROUTINE Macroscopic_Values_h_AMR( nx, ny, h, C )
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny
	REAL(8), INTENT(IN)    :: h(9, 0:nx+1, 0:ny+1)
	REAL(8), INTENT(INOUT) :: C(  -1:nx+2,-1:ny+2)

	C(1:nx,1:ny) = SUM(h(:, 1:nx,1:ny), 1)

END

!**********************************************************************

SUBROUTINE Macroscopic_Values_g_AMR( nx, ny, dx, g, C, P, Ux, Uy )
	USE D2Q9_Mod, ONLY: Rhol, Rhoh, Gy, ex, ey
	IMPLICIT NONE
	INTEGER, INTENT(IN)  :: nx, ny, dx
	REAL(8), INTENT(IN)  :: g(9, 0:nx+1,0:ny+1), C(-1:nx+2,-1:ny+2)
	REAL(8), INTENT(OUT) :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)

	INTEGER :: i, j
	REAL(8) :: Rho, dRho6, Fx, Fy, Fb
	REAL(8) :: DcDx(nx,ny), DcDy(nx,ny), mu(nx,ny)

	dRho6 = (Rhoh-Rhol)/6.d0

	CALL Gradient_AMR( nx, ny, dx, C, DcDx, DcDy )

	CALL Chemical_Potential_AMR( nx, ny, dx, C, mu )

	DO j = 1, ny
	DO i = 1, nx

		Rho = Rhol + C(i,j) * (Rhoh - Rhol)

		Fb = (Rhoh - Rho) * Gy

		Fx = mu(i,j) * DcDx(i,j)
		Fy = mu(i,j) * DcDy(i,j) + Fb

		Ux(i,j) = ( 3.d0*SUM(g(2:9, i,j)*ex(2:9)) + 0.5d0*Fx )/Rho
		Uy(i,j) = ( 3.d0*SUM(g(2:9, i,j)*ey(2:9)) + 0.5d0*Fy )/Rho

		P (i,j) = SUM(g(:, i,j)) + dRho6*( Ux(i,j)*DcDx(i,j) + Uy(i,j)*DcDy(i,j) )

	END DO
	END DO

END

!**********************************************************************

SUBROUTINE ED_C_AMR( dx, C, EdotDelC )
	USE D2Q9_Mod, ONLY: ex, ey
	IMPLICIT NONE
	INTEGER, INTENT(IN)  :: dx
	REAL(8), INTENT(IN)  :: C(-1:1,-1:1)
	REAL(8), INTENT(OUT) :: EdotDelC(9)

	INTEGER :: i

	DO i = 1, 9
		EdotDelC(i) = C(ex(i),ey(i)) - C(-ex(i),-ey(i))
	END DO
	EdotDelC(:) = EdotDelC(:) * 0.5d0/DBLE(dx)

END

!**********************************************************************

SUBROUTINE ED_M_AMR( dx, C, EdotDelM )
	USE D2Q9_Mod, ONLY: ex, ey
	IMPLICIT NONE
	INTEGER, INTENT(IN)  :: dx
	REAL(8), INTENT(IN)  :: C(-2:2,-2:2)
	REAL(8), INTENT(OUT) :: EdotDelM(9)

	INTEGER :: i

	DO i = 1, 9
		EdotDelM(i) = ( - C(-ex(i),-ey(i)) - 3.d0*C(0,0) &
						+ 5.d0*C(ex(i),ey(i)) - C(2*ex(i),2*ey(i)) )
	END DO
	EdotDelM(:) = EdotDelM(:) * 0.25d0/DBLE(dx)

END

!**********************************************************************

SUBROUTINE Gradient_ED_AMR( ED, Dx, Dy )
	IMPLICIT NONE
	REAL(8),INTENT(IN) :: ED(9)
	REAL(8),INTENT(OUT):: Dx, Dy

	Dx = ( ED(2)-ED(4) )/3.d0 + ( ED(6)-ED(7)-ED(8)+ED(9) )/12.d0
	Dy = ( ED(3)-ED(5) )/3.d0 + ( ED(6)+ED(7)-ED(8)-ED(9) )/12.d0

END

!**********************************************************************

SUBROUTINE Gradient_AMR( nx, ny, dx, C, DcDx, DcDy )
	IMPLICIT NONE
	INTEGER, INTENT(IN)  :: nx, ny, dx
	REAL(8), INTENT(IN)  :: C(-1:nx+2,-1:ny+2)
	REAL(8), INTENT(OUT) :: DcDx(nx,ny), DcDy(nx,ny)

	INTEGER :: i, j

	DO j = 1, ny
	DO i = 1, nx

		DcDx(i,j) = ( (C(i+1,j) - C(i-1,j))/3.d0  + &
					(C(i+1,j+1) + C(i+1,j-1) - C(i-1,j+1) - C(i-1,j-1))/12.d0 )/DBLE(dx)
		DcDy(i,j) = ( (C(i,j+1) - C(i,j-1))/3.d0 + &
					(C(i+1,j+1) + C(i-1,j+1) - C(i+1,j-1) - C(i-1,j-1))/12.d0 )/DBLE(dx)

	END DO
	END DO

END

!**********************************************************************

SUBROUTINE Laplacian_AMR( nx, ny, C, D2C )
	IMPLICIT NONE
	INTEGER, INTENT(IN)  :: nx, ny
	REAL(8), INTENT(IN)  :: C(-1:nx+2,-1:ny+2)
	REAL(8), INTENT(OUT) :: D2C(nx,ny)

	INTEGER :: i, j

	DO j = 1, ny
	DO i = 1, nx

		D2C(i,j) = ( C(i-1,j-1)+C(i+1,j-1)+C(i-1,j+1)+C(i+1,j+1) &
				 +4*(C(i  ,j-1)+C(i-1,j  )+C(i+1,j  )+C(i  ,j+1)) - 20*C(i,j) )/6

	END DO
	END DO

END

!**********************************************************************	

SUBROUTINE BGK_Collision_AMR( s9, geq, g )
	IMPLICIT NONE
	REAL(8), INTENT(IN)    :: s9, geq(9)
	REAL(8), INTENT(INOUT) :: g(9)

	g(:) = g(:) * (1.d0 - s9) + geq(:) * s9

END

!**********************************************************************	

SUBROUTINE MRT_Collision_AMR( s9, geq, g )
	IMPLICIT NONE
	REAL(8), INTENT(IN)    :: s9, geq(9)
	REAL(8), INTENT(INOUT) :: g(9)

	REAL(8) :: tmp(9), S(9)

	S(1:7) = 1.d0
	S(8:9) = s9
	CALL Convert_D2Q9_standard( g(:) - geq(:), tmp(:) )
	CALL Reconvt_D2Q9_standard( tmp(:)*S(:), tmp(:) )
	g(:) = g(:) - tmp(:)

END

!**********************************************************************	

SUBROUTINE Convert_D2Q9_standard( fin, fout )
	IMPLICIT NONE
	REAL(8), INTENT(IN)  :: fin (9)
	REAL(8), INTENT(OUT) :: fout(9)

	fout( 1) = fin(1)+fin(2)+fin(3)+fin(4)+fin(5)+fin(6)+fin(7)+fin(8)+fin(9)
	fout( 2) =-4.d0*fin(1)-fin(2)-fin(3)-fin(4)-fin(5)+2*(fin(6)+fin(7)+fin(8)+fin(9))
	fout( 3) = 4.d0*fin(1)-2.d0*(fin(2)+fin(3)+fin(4)+fin(5))+fin(6)+fin(7)+fin(8)+fin(9)
	fout( 4) =       fin(2)-fin(4) +fin(6)-fin(7)-fin(8)+fin(9)
	fout( 5) =-2.d0*(fin(2)-fin(4))+fin(6)-fin(7)-fin(8)+fin(9)
	fout( 6) =       fin(3)-fin(5) +fin(6)+fin(7)-fin(8)-fin(9)
	fout( 7) =-2.d0*(fin(3)-fin(5))+fin(6)+fin(7)-fin(8)-fin(9)
	fout( 8) =       fin(2)-fin(3)+fin(4)-fin(5)
	fout( 9) =       fin(6)-fin(7)+fin(8)-fin(9)

END

!**********************************************************************	

SUBROUTINE Reconvt_D2Q9_standard( fin, fout )
	IMPLICIT NONE
	REAL(8), INTENT(IN)  :: fin (9)
	REAL(8), INTENT(OUT) :: fout(9)

	fout( 1) = 4.d0*(fin(1)-fin(2)+fin(3))
	fout( 2) = 4.d0* fin(1)-     fin(2)-2.d0*fin(3) + 6.d0*(fin(4)-fin(5)) + 9.d0*fin(8)
	fout( 3) = 4.d0* fin(1)-     fin(2)-2.d0*fin(3) + 6.d0*(fin(6)-fin(7)) - 9.d0*fin(8)
	fout( 4) = 4.d0* fin(1)-     fin(2)-2.d0*fin(3) - 6.d0*(fin(4)-fin(5)) + 9.d0*fin(8)
	fout( 5) = 4.d0* fin(1)-     fin(2)-2.d0*fin(3) - 6.d0*(fin(6)-fin(7)) - 9.d0*fin(8)
	fout( 6) = 4.d0* fin(1)+2.d0*fin(2)+     fin(3) + 6.d0*fin(4)+3.d0*fin(5)+6.d0*fin(6)+3.d0*fin(7) + 9.d0*fin(9)
	fout( 7) = 4.d0* fin(1)+2.d0*fin(2)+     fin(3) - 6.d0*fin(4)-3.d0*fin(5)+6.d0*fin(6)+3.d0*fin(7) - 9.d0*fin(9)
	fout( 8) = 4.d0* fin(1)+2.d0*fin(2)+     fin(3) - 6.d0*fin(4)-3.d0*fin(5)-6.d0*fin(6)-3.d0*fin(7) + 9.d0*fin(9)
	fout( 9) = 4.d0* fin(1)+2.d0*fin(2)+     fin(3) + 6.d0*fin(4)+3.d0*fin(5)-6.d0*fin(6)-3.d0*fin(7) - 9.d0*fin(9)

	fout(:) = fout(:) / 36.d0

END
