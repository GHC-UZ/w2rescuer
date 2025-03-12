FUNCTION matmul_2x2(A, B) RESULT(C)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(2,2), B(2,2)
    REAL(8) :: C(2,2)
    INTEGER :: i, j, k

    ! Initialize result matrix C to zero
    C = 0.0D0

    ! Perform matrix multiplication C = A * B
    DO i = 1, 2
        DO j = 1, 2
            DO k = 1, 2
                C(i, j) = C(i, j) + A(i, k) * B(k, j)
            END DO
        END DO
    END DO
END FUNCTION matmul_2x2

FUNCTION matvecmul_2x2(A, B) RESULT(C)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(2,2), B(2)
    REAL(8) :: C(2)
    INTEGER :: i, k

    ! Initialize result vector C to zero
    C = 0.0D0

    ! Perform matrix-vector multiplication C = A * B
    DO i = 1, 2
        DO k = 1, 2
            C(i) = C(i) + A(i, k) * B(k)
        END DO
    END DO
END FUNCTION matvecmul_2x2

