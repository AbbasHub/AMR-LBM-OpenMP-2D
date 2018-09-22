
SUBROUTINE Refine_Derefine_Mesh( refinTol, derefTol )
	IMPLICIT NONE
	REAL(8), INTENT(IN) :: refinTol, derefTol

	CALL Flag_for_refinement  ( refinTol )
	CALL Flag_for_derefinement( derefTol )

	CALL Derefine
	CALL Refine

	CALL Array_of_Leaf_Blocks	!Moved here (02/21/2016)

END SUBROUTINE Refine_Derefine_Mesh

!******************************************************************************

SUBROUTINE Flag_for_refinement( refinTol )
	USE parameters, ONLY: Blk, Nblks, minLevels, maxLevels
	IMPLICIT NONE
	REAL(8), INTENT(IN) :: refinTol

	INTEGER :: n
	REAL(8) :: err

!$OMP PARALLEL DO private( err )
	DO n = 1, Nblks
		Blk(n)%Refin = .False.

		IF( Blk(n)%IsLeaf ) THEN
			IF( Blk(n)%Level < minLevels ) THEN
				Blk(n)%Refin = .True.
				CYCLE
			ELSEIF( Blk(n)%Level < maxLevels ) THEN

				CALL Err_Interface( Blk(n)%dx, Blk(n)%C, err )
				IF( err > refinTol ) THEN
					Blk(n)%Refin = .True.
					CYCLE
				END IF

			END IF
		END IF

	END DO
!$OMP END PARALLEL DO

	DO n = 1, Nblks
		IF( Blk(n)%Refin ) THEN
			CALL Enforce_1_Level_Jump( n )
		END IF
	END DO

END SUBROUTINE Flag_for_refinement

!******************************************************************************

SUBROUTINE Enforce_1_Level_Jump( n )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n

	INTEGER :: I

	DO I = 1, 8
		CALL Check_Level_Jump( n, I )
	END DO

END

!******************************************************************************

RECURSIVE SUBROUTINE Check_Level_Jump( n, I )
	USE parameters, ONLY: Blk
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n, I

	INTEGER :: b

	b = Blk(n)%nghbr( I )%n

	IF( b > 0 ) THEN
		IF( .NOT. Blk(b)%Refin .AND. Blk(b)%Level == Blk(n)%Level-1 ) THEN
			IF( Blk(b)%IsLeaf ) THEN
				Blk(b)%Refin = .True.
				CALL Enforce_1_Level_Jump( b )
			END IF
		END IF
	END IF

END SUBROUTINE Check_Level_Jump

!******************************************************************************

SUBROUTINE Flag_for_derefinement( derefTol )
	USE parameters, ONLY: Blk, Nblks, LevelMin, LevelMax, minLevels
	IMPLICIT NONE
	REAL(8), INTENT(IN) :: derefTol

	INTEGER	:: n, L, c1, c4, c
	REAL(8) :: err

!$OMP PARALLEL DO private( c1, c4, err )
	DO n = 1, Nblks
		Blk(n)%Deref = .False.

		IF( Blk(n)%hasLeaf ) THEN

			IF( Blk(n)%Level < minLevels ) CYCLE

			c1 = Blk(n)%child(1)
			c4 = c1 + 3

			IF( ANY(Blk(c1:c4)%Refin) ) CYCLE

		!== derefinement based on the children
			Blk(n)%Deref = .True.

			DO c = c1, c4
				CALL Err_Interface( Blk(c)%dx, Blk(c)%C, err )
				IF( err > derefTol )THEN
					Blk(n)%Deref = .False.
					EXIT
				END IF
			END DO
		!== derefinement based on the children

		END IF

	END DO
!$OMP END PARALLEL DO

	DO L = LevelMax-1, LevelMin-1, -1
		DO n = 1, Nblks
			IF( Blk(n)%Level == L ) THEN
				IF( Blk(n)%Deref ) THEN
					CALL enforce_1_Level_jump_deref_2D( n )
				END IF
			END IF
		END DO
	END DO

END SUBROUTINE Flag_for_derefinement

!******************************************************************************

SUBROUTINE enforce_1_Level_jump_deref_2D( n )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n

	CALL Check_Level_Jump_Deref_2D( n,1, 2, 4 )	! Left
	CALL Check_Level_Jump_Deref_2D( n,2, 1, 3 )	! Right
	CALL Check_Level_Jump_Deref_2D( n,3, 3, 4 )	! Bottom
	CALL Check_Level_Jump_Deref_2D( n,4, 1, 2 )	! Top

	CALL Check_Level_Jump_Deref_2D( n,5, 4, 4 )	! SW (Left -Bottom)
	CALL Check_Level_Jump_Deref_2D( n,6, 2, 2 )	! NW (Left -Top)
	CALL Check_Level_Jump_Deref_2D( n,7, 3, 3 )	! SE (Right-Bottom)
	CALL Check_Level_Jump_Deref_2D( n,8, 1, 1 )	! NE (Right-Top)

END SUBROUTINE enforce_1_Level_jump_deref_2D

!******************************************************************************

SUBROUTINE Check_Level_Jump_Deref_2D( n, Neighbor, c1, c2 )
	USE parameters, ONLY: Blk
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n, Neighbor, c1, c2

	INTEGER :: b, c, I

	IF( .NOT. Blk(n)%Deref ) RETURN

	b = Blk(n)%nghbr( Neighbor )%n

	IF( b > 0 ) THEN
		IF( .NOT. Blk(b)%IsLeaf ) THEN
			c = Blk(b)%child(c1)
			DO I = 1, 2	!I = c1, c2, (c2-c1)/2+1
				IF( Blk(c)%Refin .OR. (.Not.Blk(c)%IsLeaf .AND. .Not.Blk(c)%Deref) ) THEN
					Blk(n)%Deref = .False.
					RETURN
				END IF
				c = Blk(b)%child(c2)
			END DO
		END IF
	END IF

END SUBROUTINE Check_Level_Jump_Deref_2D

!******************************************************************************

SUBROUTINE Derefine
	USE parameters, ONLY: Blk, Nblks
	IMPLICIT NONE

	INTEGER :: n

	DO n = 1, Nblks
		IF( Blk(n)%Deref ) THEN
			CALL Restrict_ghnew( n )
			CALL Remove_Children( n )
		END IF
	END DO

END SUBROUTINE Derefine

!******************************************************************************

SUBROUTINE Remove_Children( par )
	USE parameters, ONLY: Blk, KK, Root, emptyBlk
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: par

	INTEGER :: c1, c4, grandpa, s1, s4, c

	c1 = Blk(par)%child(1)
	c4 = c1 + 3
	Blk(c1:c4)%Level  = -1
	Blk(c1:c4)%IsLeaf = .False.

	Blk(par)%child(:)= par
	Blk(par)%IsLeaf  = .True.
	Blk(par)%hasLeaf = .False.

	grandpa = Blk(par)%Parent
	IF( grandpa > 0 ) THEN
		s1 = Blk(grandpa)%child(1)	! first sibling
		s4 = Blk(grandpa)%child(4)	! last  sibling
		IF( ALL(Blk(s1:s4)%IsLeaf) ) THEN
			Blk(grandpa)%hasLeaf = .True.
		END IF
	END IF

	DO c = c1, c4
		DEALLOCATE( Blk(c)%h, Blk(c)%g, Blk(c)%C, Blk(c)%P, Blk(c)%Ux, Blk(c)%Uy )
		DEALLOCATE( Blk(c)%child )
	END DO

	IF( .NOT. ASSOCIATED(Root) ) THEN
		ALLOCATE( Root )
		emptyBlk => Root
	ELSE
		ALLOCATE( emptyBlk%next )
		emptyBlk => emptyBlk%next
	END IF
	emptyBlk%N    =  c1
	emptyBlk%next => NULL()

	KK = KK + 4

END SUBROUTINE Remove_Children

!******************************************************************************

SUBROUTINE Refine
	USE parameters, ONLY: Blk, Nblks, LevelMin, LevelMax, KK, Root
	IMPLICIT NONE

	INTEGER	:: n, L, c1

	DO L = LevelMin, LevelMax
		DO n = 1, Nblks
			IF( Blk(n)%Level == L )THEN
				IF( Blk(n)%Refin ) THEN

					IF( ASSOCIATED(Root) )THEN
						c1 = Root%N
						CALL Add_Child_Blocks( n, c1 )
						KK = KK - 4
						IF( ASSOCIATED(Root%next) )THEN
							Root => Root%next
						ELSE
							NULLIFY(Root)
						END IF
					ELSE
						c1 = Nblks+1
						CALL Add_Child_Blocks( n, c1 )
					END IF

					CALL Prolong_to_New_Children( n, c1 )

				END IF
			END IF
		END DO
	END DO

END SUBROUTINE Refine

!******************************************************************************

SUBROUTINE Err_vorticity( dx, U, V, err )
	USE parameters, ONLY: nx, ny
	IMPLICIT NONE

	INTEGER, INTENT(IN ) :: dx
	REAL(8), INTENT(IN ) :: U(nx,ny), V(nx,ny)
	REAL(8), INTENT(OUT) :: err

	INTEGER :: I, J
	REAL(8) :: w_z

	err = 0

	DO J = 2, ny-1
	DO I = 2, nx-1

		w_z = ( (V(I+1,J) - V(I-1,J)) - (U(I,J+1) - U(I,J-1)) )**2.d0

		err = MAX( err, w_z )

	END DO
	END DO

	err = SQRT(err)/(2.d0*dx)

END

!******************************************************************************

SUBROUTINE Err_Interface( dx, C, err )
	USE parameters, ONLY: nx, ny
	IMPLICIT NONE

	INTEGER, INTENT(IN ) :: dx
	REAL(8), INTENT(IN ) :: C(-1:nx+2,-1:ny+2)
	REAL(8), INTENT(OUT) :: err

	INTEGER :: I, J
	REAL(8) :: DC

	err = 0

	DO J = 1, ny
	DO I = 1, nx

		DC = (C(I+1,J) - C(I-1,J))**2.d0 + (C(I,J+1) - C(I,J-1))**2.d0

		err = MAX( err, DC )

	END DO
	END DO

	err = SQRT(err)/(2.d0*dx)

END
