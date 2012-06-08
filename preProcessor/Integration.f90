MODULE Integration

IMPLICIT NONE

CONTAINS
!-- SUBROUTINE ComputeNDS
    SUBROUTINE ComputeNDS(Mesh,c,iCase,NDS)  
    USE MMesh
    IMPLICIT NONE
    TYPE(TMesh) :: Mesh
    INTEGER :: c,iCase
    REAL,DIMENSION(*) :: NDS
    INTEGER :: i
    SELECT CASE (iCase)
    CASE (1)
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                NDS(i)=-Mesh%N(1,i)*Mesh%A(i)
            ELSE
                NDS(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                NDS(i+Mesh%Npanels)=NDS(i)
            END IF
        END DO
    CASE (2)
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                NDS(i)=-Mesh%N(2,i)*Mesh%A(i)
            ELSE
                NDS(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                NDS(i+Mesh%Npanels)=-NDS(i)
            END IF
        END DO
    CASE (3)
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                NDS(i)=-Mesh%N(3,i)*Mesh%A(i)
            ELSE
                NDS(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                NDS(i+Mesh%Npanels)=NDS(i)
            END IF
        END DO
    CASE (4)
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                NDS(i)=(-(Mesh%XM(3,i)-Mesh%CG(3,c))*Mesh%N(2,i)+(Mesh%XM(2,i)-Mesh%CG(2,c))*Mesh%N(3,i))*-1.*Mesh%A(i)
            ELSE
                NDS(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                NDS(i+Mesh%Npanels)=-NDS(i)
            END IF
        END DO
    CASE (5)
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                NDS(i)=(-(Mesh%XM(1,i)-Mesh%CG(1,c))*Mesh%N(3,i)+(Mesh%XM(3,i)-Mesh%CG(3,c))*Mesh%N(1,i))*-1.*Mesh%A(i)
            ELSE
                NDS(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                NDS(i+Mesh%Npanels)=NDS(i)
            END IF
        END DO
    CASE (6)
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                NDS(i)=(-(Mesh%XM(2,i)-Mesh%CG(2,c))*Mesh%N(1,i)+(Mesh%XM(1,i)-Mesh%CG(1,c))*Mesh%N(2,i))*-1.*Mesh%A(i)
            ELSE
                NDS(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                NDS(i+Mesh%Npanels)=-NDS(i)
            END IF
        END DO
    CASE (7)
        WRITE(*,*) 'Error: force case 7 not implemented yet'
        STOP
    CASE DEFAULT
        WRITE(*,*) 'Error: unknown radiation case'
        STOP
    END SELECT
    END SUBROUTINE
END MODULE Integration