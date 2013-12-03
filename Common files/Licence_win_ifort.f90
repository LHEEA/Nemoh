!
!	Vérification de la licence
!
	SUBROUTINE Licence
	
	USE iflport

	IMPLICIT NONE

!	Clef Nom PC
	INTEGER,PARAMETER :: Nchar=36
	CHARACTER(LEN=Nchar) :: NomPC
	CHARACTER(LEN=Nchar),PARAMETER :: Clef='COMPUTERNAME=JEANCHRISTOPHE        '
	CHARACTER(LEN=Nchar),PARAMETER :: Code='The0Key1Is2The3Name4Of5The6Computer'	
	INTEGER,DIMENSION(Nchar) :: Rclef,Rcode
!	Clef Date
	CHARACTER(LEN=8) :: NomDateDebut,Nomcourant
	CHARACTER(LEN=8),PARAMETER :: CodeDate='af0"@lba'
	INTEGER,DIMENSION(8) :: Rdate,Rcourant
	INTEGER*4,DIMENSION(3) :: DateDebut,DateCourante,DateFin
!	
	INTEGER i,j,k
	LOGICAL existence
	INTEGER rsystem

!    WRITE(*,*) ' '
!    WRITE(*,*) '-------------------------------------------------------------------------------'
!    WRITE(*,*) ' Ce logiciel est la propriete de l Ecole Centrale de Nantes et du CNRS.'
!    WRITE(*,*) ' Cette version a fait l objet d un accord de licence pour son exploitation par'
!    WRITE(*,*) ' EDF RD dans le cadre decrit dans l accord du 21 janvier 2009 entre EDF RD et'
!    WRITE(*,*) ' l Ecole Centrale de Nantes. Seule son utilisation dans ce cadre est autorisee.'
!    WRITE(*,*) '-------------------------------------------------------------------------------'
!    WRITE(*,*) ' '
!	Lecture du nom PC
	rsystem=SYSTEM('set computername > key2')
    OPEN(10,FILE='key2')
    READ(10,*) NomPC
    CLOSE(10)
	IF (NomPC.NE.Clef) THEN
!	    WRITE(*,*) 'Licence non valide'
		DO i=1,Nchar
			Rclef(i)=ichar(NomPC(i:i))+ichar(Code(i:i))
			NomPC(i:i)=char(Rclef(i))
		END DO
		OPEN(10,FILE='key2')
		WRITE(10,*) NomPC
		CLOSE(10)
		rsystem=SYSTEM('del key2')
!        STOP
    END IF 
	DO i=1,Nchar
		Rclef(i)=ichar(NomPC(i:i))+ichar(Code(i:i))
		NomPC(i:i)=char(Rclef(i))
	END DO
    OPEN(10,FILE='key2')
    WRITE(10,*) NomPC
    CLOSE(10)
	rsystem=SYSTEM('del key2')
!	Lecture de la date
	CALL DATE_AND_TIME(NomCourant)
	DateCourante(1)=ichar(Nomcourant(8:8))-48+10*(ichar(Nomcourant(7:7))-48)
	DateCourante(2)=ichar(Nomcourant(6:6))-48+10*(ichar(Nomcourant(5:5))-48)
	DateCourante(3)=ichar(Nomcourant(4:4))-48+10*(ichar(Nomcourant(3:3))-48)&
+100*(ichar(Nomcourant(2:2))-48)+1000*(ichar(Nomcourant(1:1))-48)
!	Lecture de la date de derniere utilisation
	OPEN(10,FILE='key')
	READ(10,*) NomDateDebut
	CLOSE(10)
	DO i=1,8
		Rdate(i)=ichar(NomDateDebut(i:i))-ichar(Codedate(i:i))-48
	END DO
	DateDebut(1)=10*Rdate(1)+Rdate(2)
	DateDebut(2)=10*Rdate(3)+Rdate(4)
	DateDebut(3)=1000*Rdate(5)+100*Rdate(6)+10*Rdate(7)+Rdate(8)
	IF ((DateCourante(1)+DateCourante(2)*31+DateCourante(3)*365*31).LT.(DateDebut(1)+DateDebut(2)*31+DateDebut(3)*365*31)) THEN
		WRITE(*,*) DateCourante
		WRITE(*,*) DateDebut
		WRITE(*,*) 'Licence non valide'
		STOP
	ELSE
		Rdate(1)=DateCourante(1)/10
		Rdate(2)=DateCourante(1)-10*Rdate(1)
		Rdate(3)=DateCourante(2)/10
		Rdate(4)=DateCourante(2)-10*Rdate(3)
		Rdate(5)=DateCourante(3)/1000
		Rdate(6)=(DateCourante(3)-1000*Rdate(5))/100
		Rdate(7)=(DateCourante(3)-1000*Rdate(5)-100*Rdate(6))/10
		Rdate(8)=DateCourante(3)-1000*Rdate(5)-100*Rdate(6)-10*Rdate(7)
		DO i=1,8
			Rdate(i)=Rdate(i)+48+ichar(Codedate(i:i))
			NomDatedebut(i:i)=char(Rdate(i))
		END DO
		OPEN(10,FILE='key')
		WRITE(10,*) NomDateDebut
		CLOSE(10)
	END IF
!	Date d expiration de la licence
	Datefin(1)=31
	Datefin(2)=1
	Datefin(3)=2014
!	Test		
	IF ((Datefin(1)+Datefin(2)*31+Datefin(3)*365*31).LT.(DateCourante(1)+DateCourante(2)*31+DateCourante(3)*365*31)) THEN
		WRITE(*,*) 'Licence non valide'
		STOP
	END IF

	RETURN

	END