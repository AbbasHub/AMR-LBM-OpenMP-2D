
SUBROUTINE Biquadratic_from_Fine( i2, j2, X2, Y2, Cin_F, Cout_C )

	IMPLICIT NONE

	INTEGER, INTENT(IN)  :: i2, j2, X2, Y2
	REAL(8), INTENT(IN ) :: Cin_F (i2,j2)
	REAL(8), INTENT(OUT) :: Cout_C(X2,Y2)

	INTEGER :: X, J, i, a, b, c
	REAL(8) :: tmp(i2,Y2)

	b = 6
	DO X = 1, Y2
		J = MOD(X+1,2)
		i = 2*X-1 - J
		a = 3 - 4*J
		c = 2 - a

		tmp(:,X) = a * Cin_F(:,i) + b * Cin_F(:,i+1) + c * Cin_F(:,i+2)

	END DO
	DO X = 1, X2
		J = MOD(X+1,2)
		i = 2*X-1 - J
		a = 3 - 4*J
		c = 2 - a

		Cout_C(X,:) = a * tmp(i,:) + b * tmp(i+1,:) + c * tmp(i+2,:)

	END DO
	Cout_C = Cout_C / 64.d0

END

!******************************************************************************

SUBROUTINE Biquadratic_from_Fine_f( i2, j2, X2, Y2, fin_F, fout_C )

	IMPLICIT NONE

	INTEGER, INTENT(IN)  :: i2, j2, X2, Y2
	REAL(8), INTENT(IN ) :: fin_F (9, i2,j2)
	REAL(8), INTENT(OUT) :: fout_C(9, X2,Y2)

	INTEGER :: X, J, i, a, b, c
	REAL(8) :: tmp(9, i2,Y2)

	b = 6
	DO X = 1, Y2
		J = MOD(X+1,2)
		i = 2*X-1 - J
		a = 3 - 4*J
		c = 2 - a

		tmp(:,:,X) = a * fin_F(:,:,i) + b * fin_F(:,:,i+1) + c * fin_F(:,:,i+2)

	END DO
	DO X = 1, X2
		J = MOD(X+1,2)
		i = 2*X-1 - J
		a = 3 - 4*J
		c = 2 - a

		fout_C(:,X,:) = a * tmp(:,i,:) + b * tmp(:,i+1,:) + c * tmp(:,i+2,:)

	END DO
	fout_C = fout_C / 64.d0

END

!******************************************************************************

SUBROUTINE Biquadratic_from_Coarse( I2, J2, x2, y2, Cin_C, Cout_F )

	IMPLICIT NONE

	INTEGER, INTENT(IN)  :: I2, J2, x2, y2
	REAL(8), INTENT(IN ) :: Cin_C (I2,J2)
	REAL(8), INTENT(OUT) :: Cout_F(x2,y2)

	INTEGER :: x, I, a, b, c
	REAL(8) :: tmp(I2,y2)

	b = 30
	DO x = 1, y2
		I = (x+1)/2
		a = 5 - 8*MOD(x+1,2)
		c = 2-a

		tmp(:,x) = a * Cin_C(:,I) + b * Cin_C(:,I+1) + c * Cin_C(:,I+2)

	END DO
	DO x = 1, x2
		I = (x+1)/2
		a = 5 - 8*MOD(x+1,2)
		c = 2-a

		Cout_F(x,:) = a * tmp(I,:) + b * tmp(I+1,:) + c * tmp(I+2,:)

	END DO
	Cout_F = Cout_F / 1024.d0

END

!******************************************************************************

SUBROUTINE Biquadratic_from_Coarse_f( I2, J2, x2, y2, fin_C, fout_F )

	IMPLICIT NONE

	INTEGER, INTENT(IN)  :: I2, J2, x2, y2
	REAL(8), INTENT(IN ) :: fin_C (9, I2,J2)
	REAL(8), INTENT(OUT) :: fout_F(9, x2,y2)

	INTEGER :: x, I, a, b, c
	REAL(8) :: tmp(9, I2,y2)

	b = 30
	DO x = 1, y2
		I = (x+1)/2
		a = 5 - 8*MOD(x+1,2)
		c = 2-a

		tmp(:,:,x) = a * fin_C(:,:,I) + b * fin_C(:,:,I+1) + c * fin_C(:,:,I+2)

	END DO
	DO x = 1, x2
		I = (x+1)/2
		a = 5 - 8*MOD(x+1,2)
		c = 2-a

		fout_F(:,x,:) = a * tmp(:,I,:) + b * tmp(:,I+1,:) + c * tmp(:,I+2,:)

	END DO
	fout_F = fout_F / 1024.d0

END

!******************************************************************************

SUBROUTINE Bilinear_from_Fine( i2, j2, X2, Y2, Cin_F, Cout_C )

	IMPLICIT NONE

	INTEGER, INTENT(IN)  :: i2, j2, X2, Y2
	REAL(8), INTENT(IN ) :: Cin_F (i2,j2)
	REAL(8), INTENT(OUT) :: Cout_C(X2,Y2)

	INTEGER :: X, i
	REAL(8) :: tmp(i2,Y2)

	DO X = 1, Y2
		i = 2*X-1

		tmp(:,X) = 0.5d0 * ( Cin_F(:,i) + Cin_F(:,i+1) )

	END DO
	DO X = 1, X2
		i = 2*X-1

		Cout_C(X,:) = 0.5d0 * ( tmp(i,:) + tmp(i+1,:) )

	END DO

END

!******************************************************************************

SUBROUTINE Bilinear_from_Fine_f( i2, j2, X2, Y2, hin_F, hout_C )

	IMPLICIT NONE

	INTEGER, INTENT(IN)  :: i2, j2, X2, Y2
	REAL(8), INTENT(IN ) :: hin_F (9, i2,j2)
	REAL(8), INTENT(OUT) :: hout_C(9, X2,Y2)

	INTEGER :: X, i
	REAL(8) :: tmp(9, i2,Y2)

	DO X = 1, Y2
		i = 2*X-1

		tmp(:, :,X) = 0.5d0 * ( hin_F(:, :,i) + hin_F(:, :,i+1) )

	END DO
	DO X = 1, X2
		i = 2*X-1

		hout_C(:, X,:) = 0.5d0 * ( tmp(:, i,:) + tmp(:, i+1,:) )

	END DO

END
