!--------------------------------------------------------------------------------------
! See licence and contributors list in the main directory.
!--------------------------------------------------------------------------------------

MODULE MFace
  ! A single face from a mesh.

  USE MMesh

  IMPLICIT NONE

  ! A single face extracted from the mesh
  TYPE TFace
    REAL, DIMENSION(3, 5) :: X    ! Vertices coordinates, the 5th is the same as the 1st one.
    REAL, DIMENSION(3)    :: XM   ! Center of mass
    REAL, DIMENSION(3)    :: N    ! Normal vector
    REAL                  :: A    ! Area
    REAL                  :: tDis ! Maximal radius of the panel
  END TYPE TFace

CONTAINS

  SUBROUTINE New_Face_Extracted_From_Mesh(Mesh, I, Face)
    TYPE(TMesh), INTENT(IN)  :: Mesh
    INTEGER,     INTENT(IN)  :: i ! The index of the face in the mesh
    TYPE(TFace), INTENT(OUT) :: Face

    ! Local variables
    REAL, DIMENSION(3)       :: M0

    Face%X(1:3, 1:4) = Mesh%X(1:3, Mesh%P(1:4, i))
    Face%X(1:3, 5)   = Mesh%X(1:3, Mesh%P(1, i))

    Face%XM(1:3)     = Mesh%XM(1:3, i)
    Face%N(1:3)      = Mesh%N(1:3, i)

    Face%A           = Mesh%A(i)

    ! Compute max radius
    M0(1:3)   = SUM(Face%X(1:3, 1:4), DIM=2)/4 ! Average vertex
    Face%tDis = MAX(                   &
      NORM2(Face%X(1:3, 1) - M0(1:3)), &
      NORM2(Face%X(1:3, 2) - M0(1:3)), &
      NORM2(Face%X(1:3, 3) - M0(1:3)), &
      NORM2(Face%X(1:3, 4) - M0(1:3))  &
      )
    ! TODO: For a better efficiency, compute it only once and store it in the Mesh type.
  END SUBROUTINE New_Face_Extracted_From_Mesh


  ! SUBROUTINE Reflect_Face_Around_XZ_Plane(Face)
  !   ! Change y coordinate into -y
  !   TYPE(TFace), INTENT(INOUT) :: Face

  !   Face%X(2, :) = -Face%X(2, :)
  !   Face%XM(2)   = -Face%XM(2)
  !   Face%N(2)    = -Face%N(2)
  !   Face%X(:, 1:5) = Face%X(:, 5:1:-1) ! Invert order of vertices to be coherent with normal vector.
  ! END SUBROUTINE Reflect_Face_Around_XZ_Plane


  ! SUBROUTINE Reflect_Face_Around_XY_Plane(Face, z0)
  !   ! Change z coordinate into z0 - z
  !   TYPE(TFace), INTENT(INOUT) :: Face
  !   REAL, INTENT(IN)           :: z0

  !   Face%X(3, :) = z0-Face%X(3, :)
  !   Face%XM(3)   = z0-Face%XM(3)
  !   Face%N(3)    =   -Face%N(3)
  !   Face%X(:, 1:5) = Face%X(:, 5:1:-1) ! Invert order of vertices to be coherent with normal vector.
  ! END SUBROUTINE Reflect_Face_Around_XY_Plane

END MODULE MFace
