
MODULE D2Q9_Mod

	IMPLICIT NONE

	INTEGER,PARAMETER :: ex(9) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
	INTEGER,PARAMETER :: ey(9) = [0, 0, 1, 0,-1, 1, 1,-1,-1]
	REAL(8),PARAMETER :: Wa(9) = [16,4, 4, 4, 4, 1, 1, 1, 1] / 36.d0

	REAL(8),PARAMETER :: M0 = 0.03d0	! 0.05d0
	REAL(8),PARAMETER :: wc = 1.d0/(0.5d0 + M0*3.d0)	!=1/tau_c

	REAL(8) :: W, Sigma, Beta, kappa, taul, tauh, Rhol, Rhoh
	REAL(8) :: Gy, D, Ar

	INTEGER :: Bo

END MODULE D2Q9_Mod
