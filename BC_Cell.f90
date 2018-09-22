
SUBROUTINE Bounce_Back_3( nx, ny, f )
	USE D2Q9_Mod, ONLY: ex, ey
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny
	REAL(8), INTENT(INOUT) :: f(9, 0:nx+1,0:ny+1)

	INTEGER :: X, Y, i, J, Xn

	Y = 0
	DO X = 0, nx+1
		f(1, X,Y) = f(1, X,Y+1)
		DO i = 2, 9
			Xn = MIN( MAX(0,X+ex(i)*ey(i)), nx+1 )
			J = (i-2) + 4*MOD(i/2,2)
			f(i, X,Y) = f(J, Xn,Y+1)
		END DO
	END DO

END

!******************************************************************************

SUBROUTINE Bounce_Back_4( nx, ny, f )
	USE D2Q9_Mod, ONLY: ex, ey
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny
	REAL(8), INTENT(INOUT) :: f(9, 0:nx+1,0:ny+1)

	INTEGER :: X, Y, i, J, Xn

	Y = ny+1
	DO X = 0, nx+1
		f(1, X,Y) = f(1, X,Y-1)
		DO i = 2, 9
			Xn = MIN( MAX(0,X-ex(i)*ey(i)), nx+1 )
			J = (i-2) + 4*MOD(i/2,2)
			f(i, X,Y) = f(J, Xn,Y-1)
		END DO
	END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C_3( nx, ny, C )
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny
	REAL(8), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

	INTEGER :: X

	DO X = 0, nx+1
		C(X, 0) = C(X,1)
	END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C_4( nx, ny, C )
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny
	REAL(8), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

	INTEGER :: X

	DO X = 0, nx+1
		C(X,ny+1) = C(X,ny  )
	END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C2_3( nx, ny, C )
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny
	REAL(8), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

	INTEGER :: X

	DO X = -1, nx+2
		C(X,-1) = C(X,2)
		C(X, 0) = C(X,1)
	END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C2_4( nx, ny, C )
	IMPLICIT NONE
	INTEGER, INTENT(IN)    :: nx, ny
	REAL(8), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

	INTEGER :: X

	DO X = -1, nx+2
		C(X,ny+1) = C(X,ny  )
		C(X,ny+2) = C(X,ny-1)
	END DO

END
