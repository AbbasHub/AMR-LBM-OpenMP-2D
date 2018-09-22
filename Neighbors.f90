
SUBROUTINE Find_Neighbors( par, n, I, J )
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: par, n, I, J

	INTEGER :: child

	child = 1 + I + 2*J

	CALL Neighbor( par, n, 1+I, child+1-2*I )	! Left   & Right neighbors
	CALL Neighbor( par, n, 3+J, child+2-4*J )	! Bottom & Top   neighbors

	CALL Diagonal_Neighbor( par, n, child )

END

!******************************************************************************

SUBROUTINE Neighbor( p, n, n2, c )
	USE parameters, ONLY: Blk
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: p, n, n2, c	! n2: neighbor inside the adjacent block

	INTEGER :: n1	! neighbor within the current block

	IF( Blk(p)%nghbr(n2)%n < 0 ) THEN
		Blk(n)%nghbr(n2)%n => Blk(p)%nghbr(n2)%n
	ELSE
		Blk(n)%nghbr(n2)%n => Blk( Blk(p)%nghbr(n2)%n )%child(c)
	END IF

	n1 = n2+1 - 2*MOD(n2-1,2)	! opposite neighbor

	Blk(n)%nghbr(n1)%n => Blk(p)%child(c)

END

!******************************************************************************

SUBROUTINE Diagonal_Neighbor( p, n, c )
	USE parameters, ONLY: Blk
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: p, n, c

	INTEGER, TARGET :: SW = -5
	INTEGER, TARGET :: NW = -6
	INTEGER, TARGET :: SE = -7
	INTEGER, TARGET :: NE = -8

	IF( c == 1 ) THEN
		IF( Blk(p)%nghbr(5)%n < 0 ) THEN
			Blk(n)%nghbr(5)%n => Blk(p)%nghbr(5)%n
		ELSE
			Blk(n)%nghbr(5)%n => Blk(Blk(p)%nghbr(5)%n)%child(4)
		END IF
		IF( Blk(p)%nghbr(1)%n < 0 ) THEN
			Blk(n)%nghbr(6)%n => NW	!Blk(p)%nghbr(1)%n
		ELSE
			Blk(n)%nghbr(6)%n => Blk(Blk(p)%nghbr(1)%n)%child(4)
		END IF
		IF( Blk(p)%nghbr(3)%n < 0 ) THEN
			Blk(n)%nghbr(7)%n => SE	!Blk(p)%nghbr(3)%n
		ELSE
			Blk(n)%nghbr(7)%n => Blk(Blk(p)%nghbr(3)%n)%child(4)
		END IF
		Blk(n)%nghbr(8)%n => Blk(p)%child(4)
	ELSEIF( c == 2 ) THEN
		IF( Blk(p)%nghbr(3)%n < 0 ) THEN
			Blk(n)%nghbr(5)%n => SW	!Blk(p)%nghbr(3)%n
		ELSE
			Blk(n)%nghbr(5)%n => Blk(Blk(p)%nghbr(3)%n)%child(3)
		END IF
		Blk(n)%nghbr(6)%n => Blk(p)%child(3)
		IF( Blk(p)%nghbr(7)%n < 0 ) THEN
			Blk(n)%nghbr(7)%n => Blk(p)%nghbr(7)%n
		ELSE
			Blk(n)%nghbr(7)%n => Blk(Blk(p)%nghbr(7)%n)%child(3)
		END IF
		IF( Blk(p)%nghbr(2)%n < 0 ) THEN
			Blk(n)%nghbr(8)%n => NE	!Blk(p)%nghbr(2)%n
		ELSE
			Blk(n)%nghbr(8)%n => Blk(Blk(p)%nghbr(2)%n)%child(3)
		END IF
	ELSEIF( c == 3 ) THEN
		IF( Blk(p)%nghbr(1)%n < 0 ) THEN
			Blk(n)%nghbr(5)%n => SW	!Blk(p)%nghbr(1)%n
		ELSE
			Blk(n)%nghbr(5)%n => Blk(Blk(p)%nghbr(1)%n)%child(2)
		END IF
		IF( Blk(p)%nghbr(6)%n < 0 ) THEN
			Blk(n)%nghbr(6)%n => Blk(p)%nghbr(6)%n
		ELSE
			Blk(n)%nghbr(6)%n => Blk(Blk(p)%nghbr(6)%n)%child(2)
		END IF
		Blk(n)%nghbr(7)%n => Blk(p)%child(2)
		IF( Blk(p)%nghbr(4)%n < 0 ) THEN
			Blk(n)%nghbr(8)%n => NE	!Blk(p)%nghbr(4)%n
		ELSE
			Blk(n)%nghbr(8)%n => Blk(Blk(p)%nghbr(4)%n)%child(2)
		END IF
	ELSEIF( c == 4 ) THEN
		Blk(n)%nghbr(5)%n => Blk(p)%child(1)
		IF( Blk(p)%nghbr(4)%n < 0 ) THEN
			Blk(n)%nghbr(6)%n => NW	!Blk(p)%nghbr(4)%n
		ELSE
			Blk(n)%nghbr(6)%n => Blk(Blk(p)%nghbr(4)%n)%child(1)
		END IF
		IF( Blk(p)%nghbr(2)%n < 0 ) THEN
			Blk(n)%nghbr(7)%n => SE	!Blk(p)%nghbr(2)%n
		ELSE
			Blk(n)%nghbr(7)%n => Blk(Blk(p)%nghbr(2)%n)%child(1)
		END IF
		IF( Blk(p)%nghbr(8)%n < 0 ) THEN
			Blk(n)%nghbr(8)%n => Blk(p)%nghbr(8)%n
		ELSE
			Blk(n)%nghbr(8)%n => Blk(Blk(p)%nghbr(8)%n)%child(1)
		END IF
	END IF

END
