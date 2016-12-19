c This code reads in the Euler angles for each integration point 
c 3D micrstructure data is from DREAM.3D
c implements anisotropic single crystal elasticity (cubic systems only)
c 
c CODE WRITTEN BY AADITYA LAKSHMANAN and ANAND KANJARLA

      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS,NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT,LAYER,
     4 KSPT, KSTEP, KINC)
	
	
	INCLUDE 'ABA_PARAM.INC'
			
	INTEGER,PARAMETER :: TYP=SELECTED_REAL_KIND(P=15)
	 	
	REAL(KIND=TYP) :: STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 ISO(3,3),AMB,U,ROTMAT(3,3),C11,C12,C44,EUL1,EUL2,EUL3,
     5 EULEAN(10),ADDST=0.000,ADDE=0.000
	CHARACTER :: CMNAME(80),CH,COM
	INTEGER :: I,J,NPT,NOEL,NUMEL,NUMGI,LAB,K1,K2,K3,K4,NO(20),
     1 REM,STAT,NSTATV,KAL
	C11=PROPS(1)
	C12=PROPS(2)
	C44=PROPS(3)
	ISO(1:3,1:3)=0.000
	ISO(1,1)=1.000
	ISO(2,2)=1.000
	ISO(3,3)=1.000
c	Here is dream3d	
c     CONVERTING EULER ANGLES FROM DEGREES TO RADIANS
	
      EUL1=STATEV(1)*4*ATAN(1.d0)/180.
	EUL2=STATEV(2)*4*ATAN(1.d0)/180.
	EUL3=STATEV(3)*4*ATAN(1.d0)/180.
	CALL EULCONV(EUL1,EUL2,EUL3,C11,C12,C44,DDSDDE)
c	
c     CALCULATING THE STRESS 	

      print *, DSTRAN
      print *, '32222222222222222222222222222222222222'
!      print *, 'DDSDDE'
	DO I=1,NTENS
	DO J=1,NTENS
	STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
	END DO
      END DO
      
      print *, STRESS
	RETURN
      
	END SUBROUTINE UMAT
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
      subroutine EULCONV(e1,e2,e3,c11,c12,c44,DDSDDE)
      implicit none
	
	real*8 :: Rmat(3,3),DDSDDE(6,6)
      real*8 :: C_crys(3,3,3,3),C_sam(3,3,3,3)
      real*8 :: delta2(3,3),delta4(3,3,3,3)
	real*8 :: e1,e2,e3,C11,C12,C44
	integer :: i,j,k,l,n,m,p,q,r,s

      
c     euler angles to rotation matrix      
	Rmat(1,1)=cos(e1)*cos(e3)-sin(e1)*sin(e3)*cos(e2)
	Rmat(1,2)=sin(e1)*cos(e3)+cos(e1)*sin(e3)*cos(e2)
	Rmat(1,3)=sin(e3)*sin(e2)
	Rmat(2,1)=-1*cos(e1)*sin(e3)-sin(e1)*cos(e3)*cos(e2)
	Rmat(2,2)=-1*sin(e1)*sin(e3)+cos(e1)*cos(e3)*cos(e2)
	Rmat(2,3)=cos(e3)*sin(e2)
	Rmat(3,1)=sin(e1)*sin(e2)
	Rmat(3,2)=-1*cos(e1)*sin(e2)
	Rmat(3,3)=cos(e2)
c     following tranpose is done because  euler angles maps sample
c     frame to crystal frame. Material properties are usually in 
c     the crsytal frame. 
      Rmat=TRANSPOSE(Rmat)
	delta2(1:3,1:3)=0
	delta4(1:3,1:3,1:3,1:3)=0
	delta2(1,1)=1
	delta2(2,2)=1
	delta2(3,3)=1
	delta4(1,1,1,1)=1
	delta4(2,2,2,2)=1
	delta4(3,3,3,3)=1
	
c constructing the full 4th order stifness tensor in crystal frame
c from the c11,c12 and c44
c this is valis only for cubic cases
	
      do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	C_crys(i,j,k,l)=c12*delta2(i,j)*delta2(k,l)
	C_crys(i,j,k,l)=C_crys(i,j,k,l)+c44*delta2(i,k)*delta2(j,l)
	C_crys(i,j,k,l)=C_crys(i,j,k,l)+c44*delta2(i,l)*delta2(j,k)
	C_crys(i,j,k,l)=C_crys(i,j,k,l)+(c11-c12-2*c44)*delta4(i,j,k,l)
	end do
	end do
	end do
      end do
      
c Obtaining the stifness in sample frame C_sam from C_crys and Rotation matrix 
c simple fourth order tensor transformation  
	
      
	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	 C_sam(i,j,k,l) = 0.
       do p =1,3
       do q=1,3
       do r=1,3
       do s= 1,3
       C_sam(i,j,k,l)=C_sam(i,j,k,l)+Rmat(p,i)*Rmat(q,j)*Rmat(r,k) 
     1 *Rmat(s,l)*C_crys(p,q,r,s)
       end do
       end do
       end do
       end do
      end do
      end do
      end do
      end do
          
c writing the C_sam to the DDSDDE JACOBIAN (6x6) REQUIRED BY ABAQUS	

      DDSDDE(1,1)=C_sam(1,1,1,1)
	DDSDDE(1,2)=C_sam(1,1,2,2)
	DDSDDE(1,3)=C_sam(1,1,3,3)
	DDSDDE(1,4)=C_sam(1,1,1,2)
	DDSDDE(1,5)=C_sam(1,1,1,3)
	DDSDDE(1,6)=C_sam(1,1,2,3)
	DDSDDE(2,2)=C_sam(2,2,2,2)
	DDSDDE(2,3)=C_sam(2,2,3,3)
	DDSDDE(2,4)=C_sam(2,2,1,2)
	DDSDDE(2,5)=C_sam(2,2,1,3)
	DDSDDE(2,6)=C_sam(2,2,2,3)
	DDSDDE(3,3)=C_sam(3,3,3,3)
	DDSDDE(3,4)=C_sam(3,3,1,2)
	DDSDDE(3,5)=C_sam(3,3,1,3)
	DDSDDE(3,6)=C_sam(3,3,2,3)
	DDSDDE(4,4)=C_sam(1,2,1,2)
	DDSDDE(4,5)=C_sam(1,2,1,3)
	DDSDDE(4,6)=C_sam(1,2,2,3)
	DDSDDE(5,5)=C_sam(1,3,1,3)
	DDSDDE(5,6)=C_sam(1,3,2,3)
	DDSDDE(6,6)=C_sam(2,3,2,3)
		
	end subroutine EULCONV
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	
c 	This is for SDVINI
c
c   initialization of State Dependent Variables  
c
	SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
C
       
      INCLUDE 'ABA_PARAM.INC'
      
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
	
       OPEN(12,FILE='/work/btech/mm12b035/Results/Test1/output1.txt')
       READ(12,*)STATEV(1),STATEV(2),STATEV(3)
	RETURN 
	END
