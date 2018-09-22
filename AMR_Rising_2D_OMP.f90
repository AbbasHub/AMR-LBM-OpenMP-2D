!================================================
!
!		2D Two-Phase AMR-LBM Solver
!	Conservative Phase-Field + Pressure Evolution
!				Rising Bubble
!
!		based on the following paper:
!	A. Fakhari et al., "A mass-conserving lattice Boltzmann method
!	with dynamic grid refinement for immiscible two-phase flows"
!	Journal of Computational Physics 315 (2016) 434-457
!
!	Written by Abbas Fakhari	06/05/2014
!	updated						08/30/2016
!================================================

PROGRAM Rising_Bubble_AMR_2D_OMP
	USE OMP_LIB
	IMPLICIT NONE

	INTEGER(8) :: t1, t2, Rate

	t2 = OMP_GET_MAX_THREADS( )
	CALL OMP_SET_NUM_THREADS(4)	!specify the number of threads
!$OMP PARALLEL
	t1 = OMP_GET_NUM_THREADS( )
!$OMP END PARALLEL
	PRINT '(/A,I5)', 'OpenMP number of cores used in this run :', t1
	PRINT '( A,I5)', 'OpenMP maximum number of available cores:', t2

	CALL SYSTEM_CLOCK( t1, Rate )

	CALL AMR_LBM_Rising_Bubble

	CALL SYSTEM_CLOCK( t2 )
	PRINT '(/A)',' *****************************************'
	PRINT '(A,F12.1)', ' Time Elapsed:', DBLE(t2-t1)/Rate
	PRINT '( A)',' *****************************************'

END PROGRAM Rising_Bubble_AMR_2D_OMP

!******************************************************************************

SUBROUTINE AMR_LBM_Rising_Bubble
	USE parameters, ONLY: minLevels, maxLevels, LevelMin, LevelMax
	IMPLICIT NONE

	INTEGER(8) :: t1, t2, Rate

	CALL SYSTEM_CLOCK( t1, Rate )

	LevelMin = minLevels
	LevelMax = minLevels

	CALL Initialize_Base_Block( 1 )	! lower  part of the domain
	CALL Initialize_Base_Block( 2 )	! middle part of the domain
	CALL Initialize_Base_Block( 3 )	! upper  part of the domain
	CALL Fix_Initial_Neighbors_3Blocks

	CALL Initial_Mesh( LevelMin )

	CALL Initialize_AMR_Rising_Bubble

	DO WHILE( LevelMax < maxLevels )
		CALL Refine_Derefine_Mesh( 0.002d0, 0.001d0 )	!!refinTol & derefTol
		CALL Initialize_AMR_Rising_Bubble
		LevelMax = LevelMax + 1
	END DO

	CALL SYSTEM_CLOCK( t2 )

	CALL Main_Loop

	PRINT '(A,F12.3)', ' t_initialization:', DBLE(t2-t1)/Rate

END

!******************************************************************************

SUBROUTINE Initialize_Base_Block( N )
	USE parameters, ONLY: Blk, Nblks, L0, nx
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: N

	INTEGER :: I

	Nblks = N

	Blk(N)%Level = 0
	Blk(N)%dx    = L0 / nx
	Blk(N)%I     = 0
	Blk(N)%J     = N-1	!N=1: the lower block

	ALLOCATE( Blk(N)%child(1:4) )
	Blk(N)%child(:)= -1			! Children will be determined soon
	Blk(N)%Parent  = -1			! This is the base block and has no parent

	!============ This is the numbering for the neighboring blocks
	!	6	4	8
	!	1		2
	!	5	3	7
	!============ This is the numbering for the neighboring blocks

	DO I = 1, 8
		ALLOCATE( Blk(N)%nghbr(I)%n )
		Blk(N)%nghbr(I)%n = -I	!External Boundaries
	END DO
	Blk(N)%nghbr(1)%n = N	!Periodic
	Blk(N)%nghbr(2)%n = N	!Periodic

	Blk(N)%IsLeaf  = .True.
	Blk(N)%hasLeaf = .False.
	Blk(N)%Refin   = .False.
	Blk(N)%Deref   = .False.

END SUBROUTINE Initialize_Base_Block

!******************************************************************************

SUBROUTINE Fix_Initial_Neighbors_3Blocks
	USE parameters, ONLY: Blk
	IMPLICIT NONE

	Blk(1)%nghbr(4)%n = 2
	Blk(1)%nghbr(6)%n = 2
	Blk(1)%nghbr(8)%n = 2

	Blk(2)%nghbr(3)%n = 1
	Blk(2)%nghbr(4)%n = 3
	Blk(2)%nghbr(5)%n = 1
	Blk(2)%nghbr(6)%n = 3
	Blk(2)%nghbr(7)%n = 1
	Blk(2)%nghbr(8)%n = 3

	Blk(3)%nghbr(3)%n = 2
	Blk(3)%nghbr(5)%n = 2
	Blk(3)%nghbr(7)%n = 2

END

!******************************************************************************

SUBROUTINE Initialize_AMR_Rising_Bubble
	USE parameters, ONLY: Blk, L0, nLeaf, LeafBlks, nx, ny
	USE D2Q9_Mod,	ONLY: Rhol, Rhoh, Wa, W, D, Sigma, Gy, ex, ey
	IMPLICIT NONE

	LOGICAL, SAVE :: IsFirstTime = .True.

	INTEGER :: dx, i, j, X0, Y0, n, m
	REAL(8) :: R, Ri, Xc, Yc, Fb, Rho, dRho3, ni, nj
	REAL(8) :: eDc(9), hlp(9), Gamma(9)
	REAL(8) :: DcDx, DcDy, mu(nx,ny)

	R  = L0/10.d0
	D  = 2.d0*R

	X0 = L0/2
	Y0 = L0/2

	IF( IsFirstTime ) THEN
		CALL Set_Rising_Bubble_Parameters
		IsFirstTime = .False.
	END IF

	dRho3 = (Rhoh - Rhol)/3.d0

!$OMP PARALLEL DO private( n,dx,Xc,Yc,Ri,mu,eDc,DcDx,DcDy,ni,nj,Rho,Fb,Gamma,hlp )
	DO m = 1, nLeaf
		n = LeafBlks(m)

		dx = Blk(n)%dx

		!=============	first initialize the phase-field
		DO j = -1, ny+2
		DO i = -1, nx+2

			Xc = Blk(n)%I * dx * nx + dx * (i-0.5d0)
			Yc = Blk(n)%J * dx * ny + dx * (j-0.5d0)

			Ri = DSQRT((Xc-X0)**2 + (Yc-Y0)**2)

			Blk(n)%C(i,j) = 0.5d0 + 0.5d0*DTANH(2.d0*(Ri-R)/W)

		END DO
		END DO
		!===========	now initialize the hydrodynamic properties
		DO j = 1, ny
		DO i = 1, nx

			Blk(n)%Ux(i,j) = 0
			Blk(n)%Uy(i,j) = 0

			Blk(n)%P (i,j) = -Sigma/R * Blk(n)%C(i,j)

		END DO
		END DO
		!=====

		CALL Chemical_Potential_AMR( nx, ny, dx, Blk(n)%C, mu )

		!===========
		DO j = 1, ny
		DO i = 1, nx

			CALL ED_C_AMR( dx, Blk(n)%C(i-1:i+1,j-1:j+1), eDc(:) )

			CALL Gradient_ED_AMR( eDc(:), DcDx, DcDy )

			CALL normal_AMR( DcDx, DcDy, ni, nj )

			Rho = Rhol + Blk(n)%C(i,j) * (Rhoh - Rhol)

			Fb = 0	!(Rhoh - Rho) * Gy

			CALL Equilibrium_h_g_AMR( Blk(n)%C(i,j), Blk(n)%Ux(i,j), Blk(n)%Uy(i,j), &
									  ni, nj, Gamma(:), Blk(n)%h(:, i,j) )

			hlp(:) = ( (Gamma(:)-Wa(:))*dRho3 + Gamma(:)*mu(i,j) ) * ( eDc(:) &
					- (Blk(n)%Ux(i,j)*DcDx + Blk(n)%Uy(i,j)*DcDy) ) + (ey(:)-Blk(n)%Uy(i,j))*Gamma(:) * Fb

			Blk(n)%g(:, i,j) = Blk(n)%P(i,j)*Wa(:) &
							 + (Gamma(:)-Wa(:))*Rho/3.d0 - 0.5d0*hlp(:)

		END DO
		END DO
		!=====

	END DO
!$OMP END PARALLEL DO

END

!******************************************************************************

SUBROUTINE Set_Rising_Bubble_Parameters
	USE D2Q9_Mod, ONLY: W, Sigma, Beta, kappa, taul, tauh, Wa, Rhol, Rhoh, Gy, Bo, Ar, D, M0
	IMPLICIT NONE

	REAL(8) :: muh, mul

	Rhol = 0.001d0
	Rhoh = 1

	Bo = 40
	Ar = 40

	Gy = 0.000005d0

	Sigma = Gy * (Rhoh-Rhol) * D**2 / Bo
	muh   = DSQRT(Gy*Rhoh*(Rhoh-Rhol)*D**3) / Ar

	mul   = muh/100

	tauh = 3.d0*(muh/Rhoh)
	taul = 3.d0*(mul/Rhol)

	W     =  5.0d0
	Beta  = 12.0d0 * Sigma/W
	kappa =  1.5d0 * Sigma*W

	PRINT '(/A,4F12.3)', 'muh, mul, Ar, Pe =', muh, mul, Ar, DSQRT(Gy*D)*W/M0

END

!******************************************************************************

SUBROUTINE Main_Loop
	USE parameters, ONLY: Blk, L0, T, nLeaf, LeafBlks, nx, ny
	USE D2Q9_Mod,	ONLY: W, M0, Bo, Ar, taul, tauh, Rhol, Rhoh, wc, Sigma, Gy, D
	IMPLICIT NONE

	INTEGER :: n, m, Tf, step
	INTEGER(8) :: t1, t2, Rate, t_colide, t_hg, t_stream
	REAL(8) :: t_Macroscopic_h, t_Macroscopic_g, t_collision, t_advection
	REAL(8) :: t_Output, t_Ghost_C, t_Ghost_C2, t_Ghost_hg, t_refine
	REAL(8), PARAMETER :: refinTol = 0.002d0
	REAL(8), PARAMETER :: derefTol = 0.001d0

	!==================================
	Tf   = NINT(8 * DSQRT(D/Gy) /20)*20
	step = Tf/20
	!==================================

	t_Output = 0
	t_Macroscopic_h = 0
	t_Macroscopic_g = 0
	t_collision = 0
	t_advection = 0
	t_Ghost_C  = 0
	t_Ghost_C2 = 0
	t_Ghost_hg = 0
	t_refine = 0

	!==========================================================
	PRINT '(/A/)', '    Tf     Sigma    W    Bo     Ar    M0    tauc   taul   tauh     Ref_Tols'
	PRINT '(I8,F8.4,F6.1,I6,5F7.3,2F8.4)', Tf, Sigma, W, Bo, Ar, M0, 1/wc, taul, tauh, refinTol, derefTol
	PRINT '(/A/)','     T   Nblks  nLeaf  Unused  Rho_min    Rho_max   U_max      Mass     Levels'
	!==========================================================

	CALL SYSTEM_CLOCK( t2, Rate )

	DO T = 0, Tf

		CALL Determine_Min_Max_Levels

		IF( MOD(T,step) == 0 ) THEN
			t1 = t2
			CALL TecPlot_CellCenter		!!TecPlot_CellCenter_Binary_2D
			CALL Output
			CALL SYSTEM_CLOCK( t2 )
			t_Output = t_Output + t2-t1
			IF( T==Tf ) EXIT
		END IF

	!== solve the interface-tracking equation and the pressure evolution equation
!$OMP Barrier
		CALL SYSTEM_CLOCK( t1 )
		CALL Fill_Ghost_Cells_C2
		CALL SYSTEM_CLOCK( t2 )
		t_Ghost_C2 = t_Ghost_C2 + t2-t1

		CALL Phase_Field_LBM_h_g_AMR( t_colide, t_hg, t_stream )
		t_collision = t_collision + t_colide
		t_Ghost_hg  = t_Ghost_hg  + t_hg
		t_advection = t_advection + t_stream

		CALL SYSTEM_CLOCK( t1 )
!$OMP PARALLEL DO private( n )
		DO m = 1, nLeaf
			n = LeafBlks(m)

			CALL Macroscopic_Values_h_AMR( nx, ny, Blk(n)%h, Blk(n)%C )

		END DO
!$OMP END PARALLEL DO
!$OMP Barrier
		CALL SYSTEM_CLOCK( t2 )
		t_Macroscopic_h = t_Macroscopic_h + t2-t1

		t1 = t2
		CALL Fill_Ghost_Cells_C
		CALL SYSTEM_CLOCK( t2 )
		t_Ghost_C = t_Ghost_C + t2-t1

		t1 = t2
!$OMP PARALLEL DO private( n )
		DO m = 1, nLeaf
			n = LeafBlks(m)

			CALL Macroscopic_Values_g_AMR( nx, ny, Blk(n)%dx, Blk(n)%g, &
								Blk(n)%C, Blk(n)%P, Blk(n)%Ux, Blk(n)%Uy )

		END DO
!$OMP END PARALLEL DO
!$OMP Barrier
		CALL SYSTEM_CLOCK( t2 )
		t_Macroscopic_g = t_Macroscopic_g + t2-t1

		t1 = t2
		CALL Refine_Derefine_Mesh( refinTol, derefTol )
		CALL SYSTEM_CLOCK( t2 )
		t_refine = t_refine + t2-t1

	END DO

	PRINT '(/A)', '------------------------------'
	PRINT '(A,F12.3)', ' t_Output:        ', t_Output/Rate
	PRINT '(A,F12.3)', ' t_Macroscopic_h: ', t_Macroscopic_h/Rate
	PRINT '(A,F12.3)', ' t_Macroscopic_g: ', t_Macroscopic_g/Rate
	PRINT '(A,F12.3)', ' t_collision:     ', t_collision/Rate
	PRINT '(A,F12.3)', ' t_advection:     ', t_advection/Rate
	PRINT '(A,F12.3)', ' t_Ghost_C :      ', t_Ghost_C /Rate
	PRINT '(A,F12.3)', ' t_Ghost_C2:      ', t_Ghost_C2/Rate
	PRINT '(A,F12.3)', ' t_Ghost_hg:      ', t_Ghost_hg/Rate
	PRINT '(A,F12.3)', ' t_refine:        ', t_refine/Rate

END SUBROUTINE Main_Loop

!******************************************************************************

SUBROUTINE Phase_Field_LBM_h_g_AMR( t_colide, t_hg, t_stream )
	USE parameters, ONLY: Blk, L0, nLeaf, LeafBlks, nx, ny
	IMPLICIT NONE
	INTEGER(8), INTENT(OUT) :: t_colide, t_hg, t_stream

	INTEGER(8) :: t1, t2
	INTEGER    :: m, n
	REAL(8)    :: cfl

	CALL SYSTEM_CLOCK( t1 )
!$OMP PARALLEL DO private( n )
	DO m = 1, nLeaf
		n = LeafBlks(m)

		CALL Collision_h_g_AMR( nx, ny, Blk(n)%dx, Blk(n)%C, Blk(n)%P, &
								Blk(n)%Ux, Blk(n)%Uy, Blk(n)%h, Blk(n)%g )

	END DO
!$OMP END PARALLEL DO
!$OMP Barrier
	CALL SYSTEM_CLOCK( t2 )
	t_colide = t2-t1

	t1 = t2
	CALL Fill_Ghost_Cells_hg
	CALL SYSTEM_CLOCK( t2 )
	t_hg = t2-t1

	t1 = t2
!$OMP PARALLEL DO private( n, cfl )
	DO m = 1, nLeaf
		n = LeafBlks(m)

		cfl = 1.d0/Blk(n)%dx

		IF( cfl > 0.9 ) THEN
			CALL Perfect_Shift_AMR( nx, ny, Blk(n)%h )
			CALL Perfect_Shift_AMR( nx, ny, Blk(n)%g )
		ELSE
			CALL Streaming_LaxWendroff_AMR( nx, ny, cfl, Blk(n)%h )
			CALL Streaming_LaxWendroff_AMR( nx, ny, cfl, Blk(n)%g )
		END IF

	END DO
!$OMP END PARALLEL DO
!$OMP Barrier
	CALL SYSTEM_CLOCK( t2 )
	t_stream = t2-t1

END

!******************************************************************************

SUBROUTINE Determine_Min_Max_Levels
	USE parameters, ONLY: Blk, LevelMin, LevelMax, nLeaf, LeafBlks
	IMPLICIT NONE

	INTEGER	:: L, n, m

	LevelMin =  100
	LevelMax = -100
!$OMP PARALLEL DO private( n, L ) reduction(min:LevelMin) reduction(max:LevelMax)
	DO m = 1, nLeaf
		n = LeafBlks(m)

		L = Blk(n)%Level
		LevelMin = MIN( LevelMin, L )
		LevelMax = MAX( LevelMax, L )
	END DO
!$OMP END PARALLEL DO

END SUBROUTINE Determine_Min_Max_Levels

!******************************************************************************

SUBROUTINE Output
	USE parameters, ONLY: Blk, Nblks, LevelMin, LevelMax, L0, nLeaf, LeafBlks, T, KK, nx, ny
	USE D2Q9_Mod,	ONLY: Rhol, Rhoh, Gy, D
	IMPLICIT NONE

	INTEGER :: m, n, dx, i, j
	REAL(8) :: Mass, Rho, Rho_min, Rho_max, U_max

	Rho_min =  10
	Rho_max = -10
	Mass  = 0
	U_max = 0
!$OMP PARALLEL DO private( n,dx,Rho ) reduction(min:Rho_min) reduction(max:Rho_max,U_max) reduction(+:Mass)
	DO m = 1, nLeaf
		n = LeafBlks(m)

		dx = Blk(n)%dx
		DO j = 1, ny
		DO i = 1, nx
			Rho = Rhol + Blk(n)%C(i,j) * (Rhoh - Rhol)
			Rho_min = MIN( Rho_min, Rho )
			Rho_max = MAX( Rho_max, Rho )
			Mass    = Mass + Rho * dx**2
		END DO
		END DO

		U_max = MAX( U_max, MAXVAL(ABS(Blk(n)%Uy)) )
	END DO
!$OMP END PARALLEL DO
!$OMP Barrier
	WRITE(*,10), T*DSQRT(Gy/D), Nblks, nLeaf, KK, Rho_min, Rho_max, U_max, Mass, LevelMin, LevelMax
10	FORMAT(F7.2,2I7,I6,2F11.5,F9.5,F12.1,2I4)

END SUBROUTINE Output

!******************************************************************************

SUBROUTINE Initial_Mesh( nLevels )
	USE parameters, ONLY: Blk, Nblks, LeafBlks
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: nLevels

	INTEGER :: L, n

	ALLOCATE( LeafBlks(Nblks) )

	DO L = 1, nLevels
		DO n = 1, Nblks		! Nblks will change inside Sub:Add_Child_Blocks
			IF( Blk(n)%Level == L-1 ) THEN
				CALL Add_Child_Blocks( n, Nblks+1 )
			END IF
		END DO
	END DO

	CALL Array_of_Leaf_Blocks

END SUBROUTINE Initial_Mesh

!******************************************************************************

SUBROUTINE Add_Child_Blocks( par, FreeBlk )
	USE parameters, ONLY: Blk, Nblks, MaxBlks, nx, ny
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: par, FreeBlk

	INTEGER	:: N, grandpar, I, J

	N = FreeBlk - 1

	IF( N+4 > MaxBlks ) THEN
		CALL Output
		PRINT '(/A/)', ' number of blocks exceeded MaxBlks! ... Program Terminated.'
		STOP
	END IF

	DO i = 1, 4
		Blk(par)%child(i) = N + i
	END DO
	Blk(par)%IsLeaf  = .False.
	Blk(par)%hasLeaf = .True.

	grandpar = Blk(par)%Parent
	IF( grandpar > 0 ) Blk(grandpar)%hasLeaf = .False.

	DO J = 0, 1
	DO I = 0, 1

		N = N + 1

		Blk(N)%Level = Blk(par)%Level + 1
		Blk(N)%dx    = Blk(par)%dx / 2

		Blk(N)%I = 2 * Blk(par)%I + I
		Blk(N)%J = 2 * Blk(par)%J + J

		Blk(N)%Parent = par

		ALLOCATE( Blk(N)%child(1:4) )
		Blk(N)%child(:) = N

		CALL Find_Neighbors( par, N, I, J )

		Blk(N)%IsLeaf  = .True.
		Blk(N)%hasLeaf = .False.
		Blk(N)%Refin   = .False.
		Blk(N)%Deref   = .False.

		ALLOCATE( Blk(N)%h(9, 0:nx+1, 0:ny+1) )
		ALLOCATE( Blk(N)%g(9, 0:nx+1, 0:ny+1) )
		ALLOCATE( Blk(N)%C(  -1:nx+2,-1:ny+2) )
		ALLOCATE( Blk(N)%P (nx,ny) )
		ALLOCATE( Blk(N)%Ux(nx,ny) )
		ALLOCATE( Blk(N)%Uy(nx,ny) )

		Blk(N)%h  = 0
		Blk(N)%g  = 0
		Blk(N)%C  = 0
		Blk(N)%P  = 0
		Blk(N)%Ux = 0
		Blk(N)%Uy = 0

	END DO
	END DO

	Nblks = MAX(Nblks,N)	! because unallocated blocks might have been used

END SUBROUTINE Add_Child_Blocks

!******************************************************************************

SUBROUTINE Array_of_Leaf_Blocks
	USE parameters, ONLY: Blk, Nblks, LeafBlks, nLeaf
	IMPLICIT NONE

	INTEGER :: n, i
	INTEGER :: tmpArray(Nblks)

	DEALLOCATE( LeafBlks )

	i = 0
	DO n = 1, Nblks
		IF( Blk(n)%IsLeaf ) THEN
			i = i + 1
			tmpArray(i) = n
		END IF
	END DO

	nLeaf = i
	ALLOCATE( LeafBlks(i) )

	LeafBlks(:) = tmpArray(1:i)

END
