SUBROUTINE matmul_4x4(A, B, C)
     IMPLICIT NONE
     REAL(8), INTENT(IN) :: A(4,4), B(4,4)
     REAL(8), INTENT(OUT) :: C(4,4)
     INTEGER :: i, j, k

     ! Initialize the result matrix to zero
     C = 0.0D0

     ! Perform matrix multiplication
     DO i = 1, 4
        DO j = 1, 4
           DO k = 1, 4
              C(i,j) = C(i,j) + A(i,k) * B(k,j)
           END DO
        END DO
    END DO
END SUBROUTINE matmul_4x4

FUNCTION matvecmul_4x4(A, B) RESULT(C)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(4,4), B(4)
    REAL(8) :: C(4)
    INTEGER :: i, k

    ! Initialize result vector C to zero
    C = 0.0D0

    ! Perform matrix-vector multiplication C = A * B
    DO i = 1, 4
        DO k = 1, 4
            C(i) = C(i) + A(i, k) * B(k)
        END DO
    END DO
END FUNCTION matvecmul_4x4

