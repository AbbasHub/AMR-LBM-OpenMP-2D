
SUBROUTINE Fill_Ghost_Cells_C2
	USE parameters, ONLY: Blk, LevelMin, LevelMax, nLeaf, LeafBlks, nx, ny
	IMPLICIT NONE

	INTEGER :: L, I, m, n, Nghbr(8)

	DO L = LevelMin, LevelMax
!$OMP PARALLEL DO private( n, Nghbr )
		DO m = 1, nLeaf
			n = LeafBlks(m)

			IF( L == Blk(n)%Level ) THEN

				DO I = 1, 8
					Nghbr(I) = Blk(n)%nghbr(I)%n
				END DO

				IF( Nghbr(1) > 0 ) CALL Fill_Left_C2  ( n, L, Nghbr(1) )
				IF( Nghbr(2) > 0 ) CALL Fill_Right_C2 ( n, L, Nghbr(2) )
				IF( Nghbr(3) > 0 ) CALL Fill_Bottom_C2( n, L, Nghbr(3) )
				IF( Nghbr(4) > 0 ) CALL Fill_Top_C2   ( n, L, Nghbr(4) )

				IF( Nghbr(5) > 0 ) CALL Fill_Left_Bottom_C2 ( n, L, Nghbr(5) )
				IF( Nghbr(6) > 0 ) CALL Fill_Left_Top_C2    ( n, L, Nghbr(6) )
				IF( Nghbr(7) > 0 ) CALL Fill_Right_Bottom_C2( n, L, Nghbr(7) )
				IF( Nghbr(8) > 0 ) CALL Fill_Right_Top_C2   ( n, L, Nghbr(8) )

			!	IF( Nghbr(1) < 0 ) apply BC ...
			!	IF( Nghbr(2) < 0 ) apply BC ...
				IF( Nghbr(3) < 0 ) CALL Symmetric_C2_3( nx, ny, Blk(n)%C )
				IF( Nghbr(4) < 0 ) CALL Symmetric_C2_4( nx, ny, Blk(n)%C )

			END IF
		END DO
!$OMP END PARALLEL DO
!$OMP Barrier
	END DO

END

!******************************************************************************

SUBROUTINE Fill_Left_C2( n, L, N1 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N1

	INTEGER	:: L_nbr, child, i, j, Y
	REAL(8) :: C(4,2*ny), C2(0:2,0:ny/2+1)

	L_nbr = Blk(N1)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N1)%IsLeaf ) THEN
			DO j=1,ny
				Blk(n)%C(-1:0,j) = Blk(N1)%C(nx-1:nx,j)
			END DO
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N1)%child(1)
			DO J = 1,    ny
			DO I = nx-3, nx
				C(nx-I+1,J   ) = Blk(child+1)%C(I,J)
				C(nx-I+1,J+ny) = Blk(child+3)%C(I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine( 4,2*ny, 2,ny, C, Blk(n)%C(0:-1:-1,1:ny) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		Y = (ny/2) * MOD(Blk(n)%J,2)

		DO j = 0, ny/2 + 1
			C2(2:0:-1,j) = Blk(N1)%C(nx-1:nx+1,j+Y)
		END DO

		CALL Biquadratic_from_Coarse( 3,ny/2+2, 2,ny, C2, Blk(n)%C(0:-1:-1,1:ny) )

	ELSE
		PRINT '(A,4I5)', " Error in SUB_Fill_Left_C2: ", n, N1, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Right_C2( n, L, N2 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N2

	INTEGER	:: L_nbr, child, i, j, Y
	REAL(8) :: C(4,2*ny), C2(0:2,0:ny/2+1)

	L_nbr = Blk(N2)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N2)%IsLeaf ) THEN
			DO j=1,ny
				Blk(n)%C(nx+1:nx+2,j) = Blk(N2)%C(1:2,j)
			END DO
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N2)%child(1)
			DO J = 1, ny
			DO I = 1, 4
				C(I,J   ) = Blk(child  )%C(I,J)
				C(I,J+ny) = Blk(child+2)%C(I,J)
			END DO
			END DO
		!==============================================

			CALL Biquadratic_from_Fine( 4,2*ny, 2,ny, C, Blk(n)%C(nx+1:nx+2,1:ny) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		Y = (ny/2) * MOD(Blk(n)%J,2)

		DO j = 0, ny/2 + 1
			C2(0:2,j) = Blk(N2)%C(0:2,j+Y)
		END DO

		CALL Biquadratic_from_Coarse( 3,ny/2+2, 2,ny, C2, Blk(n)%C(nx+1:nx+2,1:ny) )

	ELSE
		PRINT '(A,4I5)', " Error in SUB_Fill_Right_C2: ", n, N2, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Bottom_C2( n, L, N3 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N3

	INTEGER :: L_nbr, child, i, j, X
	REAL(8) :: C(2*nx,4), C2(0:nx/2+1,0:2)

	L_nbr = Blk(N3)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N3)%IsLeaf ) THEN
			DO i=1,nx
				Blk(n)%C(i,-1:0) = Blk(N3)%C(i,ny-1:ny)
			END DO
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N3)%child(1)
			DO J = ny-3, ny
			DO I = 1,    nx
				C(I   ,ny-J+1) = Blk(child+2)%C(I,J)
				C(I+nx,ny-J+1) = Blk(child+3)%C(I,J)
			END DO
			END DO
		!============================================

			CALL Biquadratic_from_Fine( 2*nx,4, nx,2, C, Blk(n)%C(1:nx,0:-1:-1) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		X = (nx/2) * MOD(Blk(n)%I,2)

		DO i = 0, nx/2 + 1
			C2(i,2:0:-1) = Blk(N3)%C(i+X,ny-1:ny+1)
		END DO

		CALL Biquadratic_from_Coarse( nx/2+2,3, nx,2, C2, Blk(n)%C(1:nx,0:-1:-1) )

	ELSE
		PRINT '(A,4I5)', " Error in SUB_Fill_Bottom_C2: ", n, N3, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Top_C2( n, L, N4 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N4

	INTEGER	:: L_nbr, child, i, j, X
	REAL(8) :: C(2*nx,4), C2(0:nx/2+1,0:2)

	L_nbr = Blk(N4)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N4)%IsLeaf ) THEN
			DO i=1,nx
				Blk(n)%C(i,ny+1:ny+2) = Blk(N4)%C(i,1:2)
			END DO
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N4)%child(1)
			DO J = 1, 4
			DO I = 1, nx
				C(I   ,J) = Blk(child  )%C(I,J)
				C(I+nx,J) = Blk(child+1)%C(I,J)
			END DO
			END DO
		!==============================================

			CALL Biquadratic_from_Fine( 2*nx,4, nx,2, C, Blk(n)%C(1:nx,ny+1:ny+2) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		X = (nx/2) * MOD(Blk(n)%I,2)

		DO i = 0, nx/2 + 1
			C2(i,0:2) = Blk(N4)%C(i+X,0:2)
		END DO

		CALL Biquadratic_from_Coarse( nx/2+2,3, nx,2, C2, Blk(n)%C(1:nx,ny+1:ny+2) )

	ELSE
		PRINT '(A,4I5)', " Error in SUB_Fill_Top_C2: ", n, N4, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Left_Bottom_C2( n, L, N5 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N5

	INTEGER :: L_nbr, child, i, j
	REAL(8) :: C(4,4), C2(0:2,0:2)

	L_nbr = Blk(N5)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N5)%IsLeaf ) THEN
			Blk(n)%C(-1:0,-1:0) = Blk(N5)%C(nx-1:nx,ny-1:ny)
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N5)%child(1)
			DO J = ny-3, ny
			DO I = nx-3, nx
				C(nx-I+1,ny-J+1) = Blk(child+3)%C(I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine( 4,4, 2,2, C, Blk(n)%C(0:-1:-1,0:-1:-1) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		I = nx - (nx/2) * MOD(Blk(n)%I,2) - 1
		J = ny - (ny/2) * MOD(Blk(n)%J,2) - 1

		C2(2:0:-1,2:0:-1) = Blk(N5)%C(I:I+2,J:J+2)

		CALL Biquadratic_from_Coarse( 3,3, 2,2, C2, Blk(n)%C(0:-1:-1,0:-1:-1) )

	ELSE
		PRINT '(A,4I5)', " Error in Fill_Left_Bottom_C2: ", n, N5, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Left_Top_C2( n, L, N6 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N6

	INTEGER :: L_nbr, child, i, j
	REAL(8) :: C(4,4), C2(0:2,0:2)

	L_nbr = Blk(N6)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N6)%IsLeaf ) THEN
			Blk(n)%C(-1:0,ny+1:ny+2) = Blk(N6)%C(nx-1:nx,1:2)
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N6)%child(1)
			DO J = 1,    4
			DO I = nx-3, nx
				C(nx-I+1,J) = Blk(child+1)%C(I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine( 4,4, 2,2, C, Blk(n)%C(0:-1:-1,ny+1:ny+2) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		I = nx   - (nx/2) * MOD(Blk(n)%I,2) - 1
		J = ny/2 - (ny/2) * MOD(Blk(n)%J,2)

		C2(2:0:-1,0:2) = Blk(N6)%C(I:I+2,J:J+2)

		CALL Biquadratic_from_Coarse( 3,3, 2,2, C2, Blk(n)%C(0:-1:-1,ny+1:ny+2) )

	ELSE
		PRINT '(A,4I5)', " Error in Fill_Left_Top_C2: ", n, N6, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Right_Bottom_C2( n, L, N7 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N7

	INTEGER :: L_nbr, child, i, j
	REAL(8) :: C(4,4), C2(0:2,0:2)

	L_nbr = Blk(N7)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N7)%IsLeaf ) THEN
			Blk(n)%C(nx+1:nx+2,-1:0) = Blk(N7)%C(1:2,ny-1:ny)
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N7)%child(1)
			DO J = ny-3, ny
			DO I = 1,    4
				C(I,ny-J+1) = Blk(child+2)%C(I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine( 4,4, 2,2, C, Blk(n)%C(nx+1:nx+2,0:-1:-1) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		I = nx/2 - (nx/2) * MOD(Blk(n)%I,2)
		J = ny   - (ny/2) * MOD(Blk(n)%J,2) - 1

		C2(0:2,2:0:-1) = Blk(N7)%C(I:I+2,J:J+2)

		CALL Biquadratic_from_Coarse( 3,3, 2,2, C2, Blk(n)%C(nx+1:nx+2,0:-1:-1) )

	ELSE
		PRINT '(A,4I5)', " Error in Fill_Right_Bottom_C2: ", n, N7, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Right_Top_C2( n, L, N8 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N8

	INTEGER :: L_nbr, child, i, j
	REAL(8) :: C(4,4), C2(0:2,0:2)

	L_nbr = Blk(N8)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N8)%IsLeaf ) THEN
			Blk(n)%C(nx+1:nx+2,ny+1:ny+2) = Blk(N8)%C(1:2,1:2)
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N8)%child(1)
			DO J = 1, 4
			DO I = 1, 4
				C(I,J) = Blk(child)%C(I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine( 4,4, 2,2, C, Blk(n)%C(nx+1:nx+2,ny+1:ny+2) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		I = nx/2 - (nx/2) * MOD(Blk(n)%I,2)
		J = ny/2 - (ny/2) * MOD(Blk(n)%J,2)

		C2(0:2,0:2) = Blk(N8)%C(I:I+2,J:J+2)

		CALL Biquadratic_from_Coarse( 3,3, 2,2, C2, Blk(n)%C(nx+1:nx+2,ny+1:ny+2) )

	ELSE
		PRINT '(A,4I5)', " Error in Fill_Right_Top_C2: ", n, N8, L, L_nbr
	END IF

END
