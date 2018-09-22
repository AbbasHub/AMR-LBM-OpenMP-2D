
MODULE parameters

	IMPLICIT NONE

	INTEGER :: T, Nblks, nLeaf

	INTEGER :: KK = 0

	INTEGER :: LevelMin, LevelMax

	INTEGER,PARAMETER :: nx = 4 *2
	INTEGER,PARAMETER :: ny = 4 *2

	INTEGER,PARAMETER :: L0    = 256 !* 2	!1, 2, 4, 8
	INTEGER,PARAMETER :: maxLevels = NINT( LOG(DBLE(L0/nx)) / LOG(2.0) )
	INTEGER,PARAMETER :: minLevels = 1	!maxLevels

	INTEGER,PARAMETER :: MaxBlks = 6 * (L0/nx)**2

	INTEGER,ALLOCATABLE :: LeafBlks(:)

	TYPE Ptr
		INTEGER, POINTER :: n
	END TYPE Ptr
	TYPE BLOCK
		INTEGER :: Level, dx
		INTEGER :: I, J
		INTEGER :: Parent
		INTEGER,POINTER :: child(:)
		TYPE(Ptr)		:: nghbr(1:8)
		LOGICAL :: IsLeaf
		LOGICAL :: hasLeaf
		LOGICAL :: Refin
		LOGICAL :: Deref
		REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: h, g
		REAL(8), DIMENSION(:,:),   ALLOCATABLE :: C, P, Ux, Uy
	END TYPE BLOCK
	TYPE(BLOCK) :: Blk(MaxBlks)

	TYPE Linked_Free_Blocks
		INTEGER :: N
		TYPE(Linked_Free_Blocks), POINTER :: next
	END TYPE Linked_Free_Blocks
	TYPE(Linked_Free_Blocks), POINTER :: Root, emptyBlk

END MODULE parameters
