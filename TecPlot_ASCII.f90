
SUBROUTINE TecPlot_CellCenter
	USE parameters, ONLY: Blk, nLeaf, LeafBlks, nx, ny, L0, T
	USE D2Q9_Mod,	ONLY: Sigma, Gy, D
	IMPLICIT NONE

	INTEGER :: I, J, dx, b, m
	INTEGER :: Nn, Ne, Nbs, N, E
	INTEGER, ALLOCATABLE :: connect(:,:), X(:), Y(:)
	REAL(8), ALLOCATABLE :: u(:), v(:), c(:), p(:)

	CHARACTER(Len= 2):: Id
	CHARACTER(Len=32):: FileName
	INTEGER, SAVE	 :: counter = -1

	counter = counter + 1

	WRITE( Id,'(I2.2)' ) counter

	FileName = 'Rising_Bubble_2D_' // Id // '.dat'

	OPEN (2, file = FileName)
	WRITE(2,*) 'Variables = X, Y, Ux, Uy, C, P'

	!----	calculate total number of leaf blocks
	Nbs = 0
	DO m = 1, nLeaf
		Nbs = Nbs + 1
	END DO
	Nn = Nbs * (nx+1) * (ny+1)
	Ne = Nbs *  nx    *  ny

	WRITE(2,20) 'Zone N =', Nn, ', E =', Ne, ', ZONETYPE=FEquadrilateral, DATAPACKING=Block'
	WRITE(2,* ) 'VARLOCATION=([3,4,6]=CELLCENTERED), T = "', T*DSQRT(Gy/D), '"'
20	FORMAT(2(A,I6),A)

	!----	allocate the connectivity matrix for writing output
	ALLOCATE( connect(4,Ne) )
	ALLOCATE( X(Nn), Y(Nn), u(Ne), v(Ne), c(Nn), p(Ne) )

	!----	construct the connectivity matrix
	N = 1
	E = 1
	DO m = 1, nLeaf
		b = LeafBlks(m)

		dx = Blk(b)%dx

		DO J = 1, ny+1
		DO I = 1, nx+1

			X(N) = Blk(b)%I * dx * nx + dx * (I-1)
			Y(N) = Blk(b)%J * dx * ny + dx * (J-1)

			c(N) = 0.25d0 * ( Blk(b)%C(i-1,j-1) + Blk(b)%C(i,j-1) &
							+ Blk(b)%C(i-1,j  ) + Blk(b)%C(i,j  )) - 0.5d0

			IF( I <= nx .AND. J <= ny ) THEN
				u(E) = Blk(b)%Ux(I,J) !/U0
				v(E) = Blk(b)%Uy(I,J) !/U0

				p(E) = Blk(b)%P(I,J) * 0.5d0*D/Sigma

				connect(1,E) = N
				connect(2,E) = N + 1
				connect(3,E) = N + 1 + nx+1
				connect(4,E) = N     + nx+1

				E = E + 1
			END IF

			N = N + 1

		END DO
		END DO
	END DO

	WRITE(2,'(9F10.3)') ( X(i) / FLOAT(L0), i=1,Nn )
	WRITE(2,'(9F10.3)') ( Y(i) / FLOAT(L0), i=1,Nn )
	WRITE(2,'(8E14.6)') ( u(i), i=1,Ne )
	WRITE(2,'(8E14.6)') ( v(i), i=1,Ne )
	WRITE(2,'(9E14.6)') ( c(i), i=1,Nn )
	WRITE(2,'(8E14.6)') ( p(i), i=1,Ne )
	WRITE(2,'(8I7)'   ) ( connect(1:4,i), i=1,Ne )

END
