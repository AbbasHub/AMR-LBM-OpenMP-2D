
SUBROUTINE Prolong_to_New_Children( n, c1 )
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, c1

	CALL Prolong_Bilinear( n, c1 )

END

!******************************************************************************

SUBROUTINE Prolong_Bilinear( n, c )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: n, c

	REAL(8) :: hnew_F(9, 2*nx,2*ny), Hnew_C(9, nx,ny)
	REAL(8) :: gnew_F(9, 2*nx,2*ny), Gnew_C(9, nx,ny)
	REAL(8) :: C_F(0:2*nx+1,0:2*ny+1), C_C(0:nx+1,0:ny+1)

	Hnew_C = Blk(n)%h(:, 1:nx,1:ny)
	Gnew_C = Blk(n)%g(:, 1:nx,1:ny)
	C_C    = Blk(n)%C(0:nx+1,0:ny+1)

!== newly added (02/26/2016)
	DEALLOCATE( Blk(n)%h, Blk(n)%g, Blk(n)%C, Blk(n)%P, Blk(n)%Ux, Blk(n)%Uy )
!== newly added (02/26/2016)

	CALL Bilinear_Interpolate_Extrapolate( nx, ny, Hnew_C, hnew_F )
	CALL Bilinear_Interpolate_Extrapolate( nx, ny, Gnew_C, gnew_F )
	CALL Bilinear_Interpolation          ( nx, ny, C_C,    C_F    )

!== now fill children data
	CALL Fill_Children_Data( c  , 0 , 0 , gnew_F, hnew_F, C_F )
	CALL Fill_Children_Data( c+1, nx, 0 , gnew_F, hnew_F, C_F )
	CALL Fill_Children_Data( c+2, 0 , ny, gnew_F, hnew_F, C_F )
	CALL Fill_Children_Data( c+3, nx, ny, gnew_F, hnew_F, C_F )

END

!******************************************************************************

SUBROUTINE Bilinear_Interpolate_Extrapolate( nx, ny, fin_C, fout_F )
	IMPLICIT NONE
	INTEGER, INTENT(IN)  :: nx, ny
	REAL(8), INTENT(IN ) :: fin_C (9,   nx,  ny)
	REAL(8), INTENT(OUT) :: fout_F(9, 2*nx,2*ny)

	INTEGER :: i, X, a, b
	REAL(8) :: tmp(9, 2*nx,ny)

	DO i = 1, 2*nx
		X = MAX( 1, MIN(nx-1,i/2) )
		a = 2*i - 4*X + 1
		b = 4 - a

		tmp(:, i,:) = b * fin_C(:, X,:) + a * fin_C(:, X+1,:)

	END DO
	DO i = 1, 2*ny
		X = MAX( 1, MIN(ny-1,i/2) )
		a = 2*i - 4*X + 1
		b = 4 - a

		fout_F(:, :,i) = b * tmp(:, :,X) + a * tmp(:, :,X+1)

	END DO
	fout_F = fout_F / 16.d0

END

!******************************************************************************

SUBROUTINE Bilinear_Interpolation( nx, ny, Cin_C, Cout_F )
	IMPLICIT NONE
	INTEGER, INTENT(IN)  :: nx, ny
	REAL(8), INTENT(IN ) :: Cin_C (0:  nx+1,0:  ny+1)
	REAL(8), INTENT(OUT) :: Cout_F(0:2*nx+1,0:2*ny+1)

	INTEGER :: i, X, a, b
	REAL(8) :: tmp(0:2*nx+1,0:ny+1)

	DO i = 0, 2*nx+1
		X = i/2
		a = 2*i - 4*X + 1
		b = 4 - a

		tmp(i,:) = b * Cin_C(X,:) + a * Cin_C(X+1,:)

	END DO
	DO i = 0, 2*ny+1
		X = i/2
		a = 2*i - 4*X + 1
		b = 4 - a

		Cout_F(:,i) = b * tmp(:,X) + a * tmp(:,X+1)

	END DO
	Cout_F = Cout_F / 16.d0

END

!******************************************************************************

SUBROUTINE Fill_Children_Data( c, I_idx, J_idx, gnew_F, hnew_F, C_F )
	USE parameters, ONLY: Blk, L0, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: c, I_idx, J_idx
	REAL(8), INTENT(IN) :: gnew_F(9, 2*nx,2*ny)
	REAL(8), INTENT(IN) :: hnew_F(9, 2*nx,2*ny)
	REAL(8), INTENT(IN) :: C_F(0:2*nx+1,0:2*ny+1)

	INTEGER :: I, J

	DO J = 1, ny
	DO I = 1, nx
		Blk(c)%g(:, I,J) = gnew_F(:, I+I_idx,J+J_idx)
		Blk(c)%h(:, I,J) = hnew_F(:, I+I_idx,J+J_idx)
	END DO
	END DO
	DO J = 0, ny+1
	DO I = 0, nx+1
		Blk(c)%C(I,J) = C_F(I+I_idx,J+J_idx)
	END DO
	END DO

	CALL Macroscopic_Values_g_AMR( nx, ny, Blk(c)%dx, Blk(c)%g, &
							Blk(c)%C, Blk(c)%P, Blk(c)%Ux, Blk(c)%Uy )

END
