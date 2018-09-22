
SUBROUTINE Restrict_ghnew( n )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n

	CALL Restrict_Bilinear(n)

END

!******************************************************************************

SUBROUTINE Restrict_Bilinear( n )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n

	INTEGER :: child, i, j, i1, i2, j1, j2

	REAL(8), DIMENSION(9, 2*nx,2*ny) :: h, g
	REAL(8), DIMENSION(   2*nx,2*ny) :: C, P, Ux, Uy

!== newly added (02/26/2016)
	ALLOCATE( Blk(n)%h(9, 0:nx+1, 0:ny+1), Blk(n)%g(9, 0:nx+1,0:ny+1) )
	ALLOCATE( Blk(n)%C(  -1:nx+2,-1:ny+2) )
	ALLOCATE( Blk(n)%P(nx,ny), Blk(n)%Ux(nx,ny), Blk(n)%Uy(nx,ny) )
!-- New Correction! (Added 2016-08-30):
	Blk(n)%h  = 0
	Blk(n)%g  = 0
	Blk(n)%C  = 0
	Blk(n)%P  = 0
	Blk(n)%Ux = 0
	Blk(n)%Uy = 0
!== newly added (02/26/2016) - Correction (2016-08-30)

	child = Blk(n)%child(1)

	DO j = 0, 1
		j1 = 1  + ny*j
		j2 = ny + ny*j
		DO i = 0, 1
			i1 = 1  + nx*i
			i2 = nx + nx*i

			h (:, i1:i2,j1:j2) = Blk(child)%h (:, 1:nx,1:ny)
			g (:, i1:i2,j1:j2) = Blk(child)%g (:, 1:nx,1:ny)
			C (   i1:i2,j1:j2) = Blk(child)%C (   1:nx,1:ny)
			P (   i1:i2,j1:j2) = Blk(child)%P (   1:nx,1:ny)
			Ux(   i1:i2,j1:j2) = Blk(child)%Ux(   1:nx,1:ny)
			Uy(   i1:i2,j1:j2) = Blk(child)%Uy(   1:nx,1:ny)

			child = child + 1

		END DO
	END DO

	CALL Bilinear_from_Fine_f( 2*nx,2*ny, nx,ny, h , Blk(n)%h (:, 1:nx,1:ny) )
	CALL Bilinear_from_Fine_f( 2*nx,2*ny, nx,ny, g , Blk(n)%g (:, 1:nx,1:ny) )
	CALL Bilinear_from_Fine  ( 2*nx,2*ny, nx,ny, C , Blk(n)%C (   1:nx,1:ny) )
	CALL Bilinear_from_Fine  ( 2*nx,2*ny, nx,ny, P , Blk(n)%P (   1:nx,1:ny) )
	CALL Bilinear_from_Fine  ( 2*nx,2*ny, nx,ny, Ux, Blk(n)%Ux(   1:nx,1:ny) )
	CALL Bilinear_from_Fine  ( 2*nx,2*ny, nx,ny, Uy, Blk(n)%Uy(   1:nx,1:ny) )

END
