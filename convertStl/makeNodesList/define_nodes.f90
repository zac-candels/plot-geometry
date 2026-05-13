!This program reads in a list of vertices and faces defining a solid surface, then defines a cartesian mesh of points for use in an LBM simulation, whereby each point is
!identified as being inside or outside the solid surface.

!------------------------------------------------------------------------------------------------------------------
!MODULES
!------------------------------------------------------------------------------------------------------------------
MODULE GLOBALS
!Global variables
  IMPLICIT NONE
  INTEGER :: GRIDX,GRIDY,GRIDZ,N_VERTS,N_FACES
  DOUBLE PRECISION :: SCALE_FACTOR, Z_FACTOR
  INTEGER, ALLOCATABLE :: NODES(:,:,:), FACES(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: VERTS(:,:)
  
END MODULE GLOBALS
!------------------------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------------------------
!MAIN
!------------------------------------------------------------------------------------------------------------------
PROGRAM DEFINE_NODES
USE GLOBALS
IMPLICIT NONE

!Define the conversion factor between the units in verts.in and lattice units
SCALE_FACTOR=20.0
!Define the domain height relative to the surface texture height
Z_FACTOR=2.0

CALL READ_IN()

CALL PROCESS_SURFACE()

CALL FIND_INTERIOR()

CALL WRITE_OUT()


END PROGRAM DEFINE_NODES
!------------------------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------------------------
!SUBROUTINES
!------------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------------
SUBROUTINE READ_IN()
!Read in the vertex and face files
USE GLOBALS
IMPLICIT NONE
INTEGER :: IO, J1


!Get number of vertices
N_VERTS=0
OPEN (1, FILE = 'verts.in')
DO
  READ(1,*,IOSTAT=IO)
  IF (IO/=0) EXIT
  N_VERTS=N_VERTS+1
END DO
CLOSE (1)

ALLOCATE(VERTS(N_VERTS,3))

!Read vertex list
OPEN(1,FILE='verts.in')
DO J1=1,N_VERTS
	READ(1,*) VERTS(J1,:)
ENDDO
CLOSE(1)

!Get number of faces
N_FACES=0
OPEN (1, FILE = 'faces.in')
DO
  READ(1,*,IOSTAT=IO)
  IF (IO/=0) EXIT
  N_FACES=N_FACES+1
END DO
CLOSE (1)

ALLOCATE(FACES(N_FACES,3))

!Read face list
OPEN(1,FILE='faces.in')
DO J1=1,N_FACES
	READ(1,*) FACES(J1,:)
ENDDO
CLOSE(1)

!If face index counting starts at zero, shift this to 1
IF (MINVAL(FACES)==0) FACES=FACES+1

END SUBROUTINE READ_IN
!------------------------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------------------------
SUBROUTINE PROCESS_SURFACE()
!Scale and adjust the surface to create the NODES array
USE GLOBALS
IMPLICIT NONE
DOUBLE PRECISION :: XMIN, YMIN, ZMIN
INTEGER :: J1

!Scale the domain
VERTS(:,1)=VERTS(:,1)*SCALE_FACTOR
VERTS(:,2)=VERTS(:,2)*SCALE_FACTOR
VERTS(:,3)=VERTS(:,3)*SCALE_FACTOR

!Translate the vertex coordinates so that the minimum lies at 1,1,1
XMIN=MINVAL(VERTS(:,1))
YMIN=MINVAL(VERTS(:,2))
ZMIN=MINVAL(VERTS(:,3))

VERTS(:,1)=VERTS(:,1)-XMIN+1
VERTS(:,2)=VERTS(:,2)-YMIN+1
VERTS(:,3)=VERTS(:,3)-ZMIN+1

!Calculate the domain sizes
GRIDX=MAXVAL(VERTS(:,1))
GRIDY=MAXVAL(VERTS(:,2))
GRIDZ=MAXVAL(VERTS(:,3))*Z_FACTOR

PRINT *, 'Domain size, GRIDX, GRIDY, GRIDZ = ', GRIDX, GRIDY, GRIDZ

!Allocate nodes array
ALLOCATE(NODES(GRIDX,GRIDY,GRIDZ))
NODES(:,:,:)=1

END SUBROUTINE PROCESS_SURFACE
!------------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------------
SUBROUTINE FIND_INTERIOR()
!Find which nodes are inside or outside the solid surface
USE GLOBALS
IMPLICIT NONE
INTEGER :: J1,J2,J3,JF, MINX, MINY, MINZ, MAXX, MAXY, MAXZ
DOUBLE PRECISION :: VE1(3), VE2(3), VE3(3), VE4(3), VN1(3), NORM(3), A0, A12, A13, A23, AREA_CONV


AREA_CONV=1.0D-9

DO JF=1,N_FACES
	!Take one of the faces
	VE1=VERTS(FACES(JF,1),:)
	VE2=VERTS(FACES(JF,2),:)
	VE3=VERTS(FACES(JF,3),:)
					
	!Find the range of nodes that may be in the prism
	MINX=CEILING(MINVAL([VE1(1),VE2(1),VE3(1)]))
	MINY=CEILING(MINVAL([VE1(2),VE2(2),VE3(2)]))
	MINZ=CEILING(MINVAL([VE1(3),VE2(3),VE3(3)]))

	MAXX=FLOOR(MAXVAL([VE1(1),VE2(1),VE3(1)]))
	MAXY=FLOOR(MAXVAL([VE1(2),VE2(2),VE3(2)]))
	MAXZ=FLOOR(MAXVAL([VE1(3),VE2(3),VE3(3)]))
	
	IF (MINX<=0) MINX=1
	IF (MINY<=0) MINY=1
	IF (MINZ<=0) MINZ=1
	IF (MAXX>GRIDX) MAXX=GRIDX
	IF (MAXY>GRIDY) MAXY=GRIDY
	IF (MAXZ>GRIDZ) MAXZ=GRIDZ
	
	!Area of triangular face projected onto the z=0 plane
	CALL AREA_PROJ(A0,VE3-VE1,VE2-VE1)		
	!What this loop does is find whether a point is underneath the lowest-most corner of a triangular face. To do this, all
	!we need to do is find whether the point, when projected onto the z=0 plane, is inside or outside the triangular face that has
	!also been projected onto the z=0 plane. To do this:
	!(1) Take a point (that you want to know if it's inside or outside the shape) and project it onto the z=0 plane (called V4). 
	!(2) You then examine 3 triangles in the z=0 plane, formed by V4 and the each pair of the 3 projected points of the triangular face. 
	!You calculate the areas of these 3 triangles. A neat fact is that if the point V4 lies within the 3 points of the projected face, the
	!area of the sum of the 3 triangles equals the area of the projected face, A0 (within some numerical tolerance AREA_CONV). Otherwise, if
	!V4 is outside the projected face, the area of the sum of the 3 triangles > A0.
	!(3) Now, comes the next neat trick. NODES has a value of 1 if it is outside the solid, or -1 if it is inside the solid. If you detect
	!that a point lies directly below a face, you flip the sign of NODES. Why flip the sign? Well imagine the solid is a sphere in the middle
	!of your simulation domain. All points inside the sphere have exactly 1 face above them. All points outside the sphere have either exactly 0
	!faces above them or 2 faces above them. For a more general shape, all points that have an odd number of faces above them are inside the solid,
	!and all points that have an even number of faces above them are outside the solid.
	DO J1=MINX,MAXX
		DO J2=MINY,MAXY
			VE4=[J1,J2,0]
			CALL AREA_PROJ(A12,VE1-VE4,VE2-VE4)
			CALL AREA_PROJ(A13,VE1-VE4,VE3-VE4)
			CALL AREA_PROJ(A23,VE2-VE4,VE3-VE4)

			IF ((A12+A13+A23<=A0*(1+AREA_CONV))) THEN
				NODES(J1,J2,1:MINZ)=-1*NODES(J1,J2,1:MINZ)
			ENDIF
		ENDDO
	ENDDO
		
		
	!This loop is the more general version of the one above, and can be used for all points (not only those below the lowest vertex of
	!the face. It works identically to the process above with 2 additional steps:
	!(1) Get the normal vector of the face, and make sure it is facing in up (in the +z direction)
	!(2) Get the vector joining any vertex of the face (vertex 1 VE1 is used below) to the point of interest VE4. If the dot product
	!of the vector VE1-VE4 and the upwards-facing face normal is negative, we know the point lies above the face. If the dot product
	!of the vector VE1-VE4 and the upwards-facing face normal is positive, we know the point lies below the face. If the latter
	!condition is satisfied, AND, the projected point lies within the projected face, then we know the point lies beneath
	!the face, and so we flip the value of NODES.
	DO J1=MINX,MAXX
		DO J2=MINY,MAXY
			DO J3=MINZ+1,MAXZ
				VE4=[J1,J2,J3]
				CALL AREA_PROJ(A12,VE1-VE4,VE2-VE4)
				CALL AREA_PROJ(A13,VE1-VE4,VE3-VE4)
				CALL AREA_PROJ(A23,VE2-VE4,VE3-VE4)

				NORM(:)=0
				CALL CROSS_PROD(NORM,VE3-VE1,VE2-VE1)

				IF (DOT_PRODUCT(NORM,[0,0,1])<=0) THEN
					NORM=-NORM
				ENDIF
				
				VN1=VE1-VE4	
			
				IF ( (A12+A13+A23<=A0*(1+AREA_CONV)) .AND. (DOT_PRODUCT(NORM,VN1)>=0) ) THEN
					NODES(J1,J2,J3)=-1*NODES(J1,J2,J3)
				ENDIF
			ENDDO
		ENDDO
	ENDDO
	

ENDDO	


END SUBROUTINE FIND_INTERIOR
!------------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------------
SUBROUTINE WRITE_OUT()
!Write the output files
USE GLOBALS
IMPLICIT NONE
INTEGER :: J1,J2,J3


!Write out the scaled vertices
OPEN(1,FILE='verts_scaled.out')
DO J1=1,N_VERTS
	WRITE(1,*) VERTS(J1,:)
ENDDO
CLOSE(1)

OPEN(1,FILE='faces_scaled.out')
DO J1=1,N_FACES
	WRITE(1,*) FACES(J1,:)
ENDDO
CLOSE(1)

!Write out the nodes array
OPEN(1,FILE='nodes.out')
DO J1=1,GRIDX
	DO J2=1,GRIDY
		DO J3=1,GRIDZ
			WRITE(1,*) NODES(J1,J2,J3)
		ENDDO
	ENDDO
ENDDO


END SUBROUTINE WRITE_OUT
!------------------------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------------------
SUBROUTINE CROSS_PROD(VEC_OUT,R1,R2)
!Cross product
IMPLICIT NONE
DOUBLE PRECISION :: R1(3), R2(3), VEC_OUT(3)
	
VEC_OUT(1)= R1(2)*R2(3)-R1(3)*R2(2)
VEC_OUT(2)= R1(3)*R2(1)-R1(1)*R2(3)
VEC_OUT(3)= R1(1)*R2(2)-R1(2)*R2(1)

END SUBROUTINE CROSS_PROD
!-----------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------
SUBROUTINE AREA_PROJ(AREA_OUT,R1,R2)
!Projected area
IMPLICIT NONE
DOUBLE PRECISION :: R1(3), R2(3), A1(3), A2(3), AREA_OUT, CPROD(3)
A1=R1
A2=R2
A1(3)=0
A2(3)=0

CPROD(:)=0.0
CALL CROSS_PROD(CPROD,A1,A2)
AREA_OUT=ABS(0.5*NORM2(CPROD))
END SUBROUTINE AREA_PROJ
!-----------------------------------------------------------------------------------------------------------------





