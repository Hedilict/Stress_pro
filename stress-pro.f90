!�˳����Ƕ�ԭ����ƽ��Ӧ������ĸĽ���ϸ��������5x5��
!��������㵯����ѧƽ��Ӧ�����⣨����ȣ�����F90�ϵ���ͨ��
!����������STRESS���ĸ��ӳ���ELSTMX(KK)��MODIFY��DCMPBD��SLVBD
!��������Ҫ�����������롢���������������ļ��㡢��װ����նȾ����Լ����㵥Ԫ��
!��Ӧ����Ӧ��ֵ
!�ӳ���ELSTMX(KK)����նȾ���ϵ�������ɵ�KK��Ԫ�ĸնȾ���K
!�ӳ���MODIFY�����غɽڵ���غ�ֵ��λ�Ʊ߽�ڵ㴦��λ��ֵ��������նȾ���λ������
!�ͽڵ�����������޸�
!�ӳ���DEMPBD�ø�˹��Ԫ�����ԳƵȴ����ܸնȾ�����ϰ벿�ֻ�Ϊ��������
!�ӳ���SLVBD���ӳ���DEMPBD�õ��������Ǿ�����лش����㣬�õ��ڵ�λ������DD
!
	PROGRAM STRESS
	COMMON /ELMATX/ESM(6,6),X(3),Y(3),D(3,3)
	COMMON / GRAD/ B(3,6),AR2
	COMMON /MTL/ EM,PR,TH
	COMMON /AV/A(8500),JGF,JGSM,NP,NBW,JEND
	COMMON /PRO/IRATION,IROWNUM,ICOLNUM,NN,NE
	DIMENSION NS(6),U(6),STRA(3),STRE(3)
	DATA IN/60/,IO/61/
!
!����Ԫ�Լ���ɵ�Ԫ�Ľڵ����鰴�ն�̬������䷽��
	REAL,ALLOCATABLE::XC(:),YC(:)
	INTEGER,ALLOCATABLE::NEL(:,:)
!	INTEGER :: RATION
!	INTEGER :: ROWNUM
!	INTEGER :: COLNUM
    INTEGER::A1
    INTEGER::B1
    INTEGER::C1
    INTEGER::A0
    INTEGER::B0
    INTEGER::C0
    INTEGER::FLAG
!
!ESM(6,6)�D�D��Ԫ�նȾ���X(3),Y(3)�D�D��Ԫ�ڵ����꣬D(3,3)�������ʾ���
!B(3,6)�D�D���ξ���AR2�D�D����������Ķ�����NP�D�D���ɶ�������NBW�D�D�������
!A(8500)�D�D�洢�ڵ�λ���������ڵ�������������նȾ��������A
!JGF��JENDΪ������Ԫ
!JGF=NP-�ڵ�λ������������A�е�λ�ã�JGSM=JGF+NP-�ڵ�������������A�е�λ��
!JEND=JGSM+NP*NBW�D�D�նȾ���������A�е�λ�ã�����A�ܳ���
!NS(6)�D�Dһ����Ԫ�ڵ�����ɶȱ�����飬U(6)�D�Dһ����Ԫ�Ľڵ����ɶ�
!IN/60/,IO/61/һ����������ļ��豸�ţ�STRA(3),STRE(3)�D�D�洢��Ԫ��Ӧ�䡢Ӧ��
!
!---------------------------���������-------------------------------------
! TITLE�D�D�洢�������ݱ���
! NN�D�D�ڵ�������NE�D�D��Ԫ������EM�D�D����ģ����PR�D�D���ɱȣ�TH�D�D��ĺ�ȣ�N�D�D��Ԫ���
! NEL(I)�D�D�ڵ��X�����꣬YC(I)�D�D�ڵ��Y������
! NEL(N,I)�D�D��ɵ�N�������εĵ�I���ڵ�ı��(I=1,2,3)
!
	OPEN(60,FILE='INPUT.DAT',STATUS='UNKNOWN')
	OPEN(61,FILE='OUTPUT.DAT',STATUS='UNKNOWN')
	READ(IN,1) TITLE
1	FORMAT(900A)
	WRITE(IO,1) TITLE
	WRITE(IO,*)
! 
! �������ģ�͵Ľڵ�����NN�͵�Ԫ����NE
!	READ(IN,*)NN
!	READ(IN,*)NE
	IRATION=5
	NN=(1+2*IRATION+1)*(2*IRATION+1)/2
	NE=((2*(2*IRATION)-1)+1)*(2*IRATION)/2 
	WRITE(IO,2)NN,NE
2	FORMAT(/,'��������:',/,'�ڵ���=',I3,5X,'��Ԫ��=',I4,/)
	NP=2*NN
	ALLOCATE(NEL(1:NE,1:3),XC(NN),YC(NN))
! 
! ������ϵ�����ģ��EM�����ɱ�PR��ƽ����TH���ڵ�����XC(I),YC(I)����ɵ�Ԫ��
! �ڵ�NEL(N,I)
! ��ɵ�Ԫ�Ľڵ�ı�Ŷ�����ʱ��˳������
! 
	READ(IN,*)EM
	READ(IN,*)PR
	READ(IN,*)TH
!����Ľڵ�������ó����Զ����룬��Ҫ�������������������ඥ�������
!����һ��SETCOOR����ʵ��
    XLENGTH=2./(2*IRATION)
	YLENGTH=2./(2*IRATION)
	IROWNUM=2*IRATION+1
	ICOLNUM=2*IRATION+1
	DO 456 I=1,IROWNUM
		DO 345 J=1,I
			XC(I*(I-1)/2+J)=(J-1)*XLENGTH  
			YC(I*(I-1)/2+J)=2.0-(I-1)*YLENGTH
345		CONTINUE
456	CONTINUE
! 
! ������ϵ����ʺͼ���ģ���������ݱ��ڼ��ʱ����
! 
	WRITE(IO,3)
3	FORMAT('���ϳ���Ϊ:')
	WRITE(IO,4)EM,PR,TH
4	FORMAT('����ģ��Ϊ:',E12.5,5X,'���ɱ�Ϊ:',E12.5,5X,'���Ϊ:',E12.5)
	WRITE(IO,*)
	DO 6 I=1,NN
	WRITE(IO,5)I,XC(I),YC(I)
5	FORMAT('�ڵ�',I3,'����Ϊ:X=',F8.3,3X,'Y=',F8.3)
6	CONTINUE
	WRITE(IO,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!���NEL�����ֵ��NEL����洢��Ԫ����ڵ�Ķ�Ӧ��ϵ
        A1=1
	    B1=2
	    C1=3
        A0=1
        B0=2
        C0=3
	    NEL(1,1)=A1
	    NEL(1,2)=B1
	    NEL(1,3)=C1
	    M=2
	    FLAG=1
        I=1
	    GOTO 91
90		IF(MOD(FLAG,2).EQ.1) THEN
			B1=C1
			C1=A1+1
			NEL(M,1)=A1
			NEL(M,2)=B1
			NEL(M,3)=C1
			FLAG=FLAG+1
			M=M+1
            GOTO 91
		ENDIF
		IF(MOD(FLAG,2).EQ.0) THEN
			A1=C1
			C1=B1+1
			NEL(M,1)=A1
			NEL(M,2)=B1
			NEL(M,3)=C1
			M=M+1
			FLAG=FLAG+1
            GOTO 91
		ENDIF
91		IF(C1.LT.((I+1)**2+(I+1))/2) THEN
			GOTO 90
		ENDIF
		IF(C1.GE.((I+1)**2+(I+1))/2) THEN 
			IF((I+1).LT.IROWNUM) THEN
                A0=B0
			    B0=B0+I+1
			    C0=B0+1
                A1=A0
			    B1=B0
			    C1=C0
				NEL(M,1)=A1
				NEL(M,2)=B1
				NEL(M,3)=C1
				FLAG=1
                I=I+1
                M=M+1
				GOTO 91
			ENDIF
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	WRITE(IO,*)
! 
! ---------------------------���㿪ʼ-------------------------------
! �����������B=MAX(D+1)*F
! D��һ����Ԫ���ڵ����֮������ֵ��F��һ���ڵ�����ɶ���
! 

	INBW=0
	NBW=0
	DO 20 KK=1,NE
	DO 25 I=1,3
25	NS(I)=NEL(KK,I)
	DO 21 I=1,2
		IJ=I+1
	DO 21 J=IJ,3
		NB=IABS(NS(I)-NS(J))    !�ڵ��֮��ľ���ֵ
		IF(NB.LE.NBW) GOTO 21
		INBW=KK
		NBW=NB
21 CONTINUE
20 CONTINUE
	NBW=(NBW+1)*2   !ƽ������ڵ����ɶ�F=2,NBW��ʱΪ�������
!
!����A�����ݵİ���:A(1,2,3....NP|NP+1.....2NP|2NP+1.....JEND)
! 					�ڵ�λ������  �ڵ�������   ����նȾ���
!��ʼ������A
	JGF=NP     !JGF=2*NN
	JGSM=JGF+NP  !JGSM=4*NN
	JEND=JGSM+NP*NBW     !JEND=4*NN+2*NN*12
! NPΪK��KNN��N��NBWΪ�������
!  	|XXXX|
! 	 |XXXX|
! 	  |XXXX|
! 	   |XXX*|
! 	    |XX**|
! 	     |X***|
!�õȴ����ά���鷽���洢�ĸնȾ����а���*��ʾ��λ�ã�����ʵ�ʳ�����δʹ��
! 
	JL=JEND-JGF
	DO 24 I=1,JEND
24	A(I)=0.0
	GOTO 30
! 
! ���ɲ������ʾ���
! EΪ����ģ��
! VΪ���ɱ�
! 
30	R=EM/(1.-PR**2)
	D(1,1)=R
	D(2,2)=D(1,1)
	D(3,3)=R*(1.-PR)/2.
	D(1,2)=PR*R
	D(2,1)=D(1,2)
	D(1,3)=0.0
	D(3,1)=0.0
	D(2,3)=0.0
	D(3,2)=0.0
! 
! ��Ԫ����ѭ���Ŀ�ʼ
! 
	KK=1
! KKΪ��Ԫ��
! �ڵ����ɶȵ����ɣ��ڵ�����ľֲ��仯
! 
32	DO 31 I=1,3
	J=NEL(KK,I)   !Ԫ��KK�������ڵ���
	NS(2*I-1)=J*2-1
	NS(2*I)=J*2   !Ԫ��KK�������ɶ����ܸնȾ����еķ��̺�
	X(I)=XC(J)
31	Y(I)=YC(J)  !Ԫ��KK�Ľڵ�����ֵ
! 
! �����ӳ���ELSTMX���㵥Ԫ�նȾ���ESM(6,6)�����
! 
	CALL ELSTMX(KK)
! 
! ��Ԫ�նȾ�����װ������նȾ���
	DO 33 I=1,6
		II=NS(I)
		DO 34 J=1,6
			JJ=NS(J)+1-II
			IF(JJ.LE.0)GOTO 34
			J1=JGSM+(JJ-1)*NP+II-(JJ-1)*(JJ-2)/2
			A(J1)=A(J1)+ESM(I,J)
34		CONTINUE
33	CONTINUE
	KK=KK+1
	IF(KK.LE.NE)GOTO 32
! 
! �����ӳ���MODIFY�����غɽڵ㴦���غ�ֵ��λ�Ʊ߽�ڵ㴦��λ��ֵ��������նȾ���
! λ������ͽڵ������������Ӧ���޸�
	CALL MODIFY
		WRITE(IO,35)
35 FORMAT(/,'������Ϊ:',/)
		WRITE(IO,36)NBW
36 FORMAT('�������Ϊ',I3)
		WRITE(IO,37)JEND
37 FORMAT('�������СΪ:',I4)
!
!�����ӳ���DEMPBD���ø�˹��Ԫ�����ԳƵȴ����ܸնȾ���Ϊ�����Ǿ���
	CALL DCMPBD
!
! �����ӳ���SLVBD���ӳ���DEMPBD�õ��������Ǿ�����лش����㣬�õ��ڵ�λ������DD
	CALL SLVBD
! 
! ---------------------------��������ڵ��λ������----------------------------
! 
	WRITE(IO,*)
	DO 45 I=1,NP/2
		WRITE(IO,43)I,A(2*I-1),A(2*I)
43 	FORMAT('�ڵ��',I3,5X,'X�����λ��UX=',E12.5,5X,'Y�����λ��UY=',E12.5)
45 	CONTINUE
! 
! ---------------------------���Ӽ���-----------------------------------------
! ����ڵ㴦��Ӧ��
! 
	DO 96 KK=1,NE
! ���ɽڵ�����ɶ�
! �ڵ�����ֲ��仯
	DO 51 I=1,3
		J=NEL(KK,I)
		NS(2*I-1)=2*J-1
		NS(2*I)=2*J
		X(I)=XC(J)
51		Y(I)=YC(J)
! 
! ��Ԫ�ڵ�λ�ƾֲ��仯
! 
65 DO 73 I=1,6,2
	NS1=NS(I)
	NS2=NS(I+1)
	U(I)=A(NS1)   !�ڵ�X�����λ����
73	U(I+1)=A(NS2)	!�ڵ�Y�����λ����
!
!���㵥ԪӦ��
! 
	CALL ELSTMX(KK)
	DO 52 I=1,3
		STRA(I)=0.0
		DO 52 K=1,6
52			STRA(I)=STRA(I)+B(I,K)*U(K)/AR2
!
! ���㵥ԪӦ��ֵ
! 
	DO 58 I=1,3
		STRE(I)=0
		DO 58 K=1,3
58	STRE(I)=STRE(I)+D(I,K)*(STRA(K))
! 
! ������Ӧ��S1��S2��TM
! 
	AA=(STRE(1)+STRE(2))/2.
	AB=SQRT((ABS(STRE(1)-STRE(2))/2.)**2+STRE(3)**2)
	S1=AA+AB
	S2=AA-AB
	TM=AB
! 
! ������Ӧ��������X��ļн�
! 
	IF(ABS(STRE(1)-STRE(2)).LT.0.0001) GOTO 93
	AC=ATAN2(2*STRE(3),(STRE(1)-STRE(2)))
	THM=(180/3.1415926*AC)/2
	GOTO 94
93	THM=90
! 
! ------------------------�����ԪӦ���Ӧ���ļ�����--------------------------
! 
94	WRITE(IO,57)KK
57	FORMAT(/,'��Ԫ',I4)
	WRITE(IO,95)STRA(1),STRA(2),STRA(3)
95	FORMAT('EPTOX=',E12.5,2X,'EPTOY=',E12.5,2X,'EPTOXY=',E12.5)
	WRITE(IO,97)STRE(1),STRE(2),STRE(3)
97	FORMAT('SX=',E12.5,5X,'SY=',E12.5,5X,'SXY=',E12.5)
	WRITE(IO,98)S1,S2,TM
98	FORMAT('S1=',E12.5,5X,'S2=',E12.5,5X,'TMAX=',E12.5)
	WRITE(IO,99)THM
99	FORMAT(44X,'ANGEL=',F8.2,'��')
96	CONTINUE
	CLOSE(60)
	CLOSE(61)
	STOP
	END

! 
! ����ն�ϵ�����󣬼����KK��Ԫ�ĸնȾ���K
! 
	SUBROUTINE ELSTMX(KK)
		COMMON/MTL/EM,PR,TH
		COMMON/GRAD/B(3,6),AR2
		COMMON/ELMATX/ESM(6,6),X(3),Y(3),D(3,3)
		DIMENSION C(6,3)
		DATA IO/61/ IN/60/
! 
! �������B
!
	DO 20 I=1,3
	DO 20 J=1,6
20	B(I,J)=0.0
	B(1,1)=Y(2)-Y(3)
	B(1,3)=Y(3)-Y(1)
	B(1,5)=Y(1)-Y(2)
	B(2,2)=X(3)-X(2)
	B(2,4)=X(1)-X(3)
	B(2,6)=X(2)-X(1)
	B(3,1)=B(2,2)
	B(3,2)=B(1,1)
	B(3,3)=B(2,4)
	B(3,4)=B(1,3)
	B(3,5)=B(2,6)
	B(3,6)=B(1,5)
	AR2=X(2)*Y(3)+X(3)*Y(1)+X(1)*Y(2)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3)
! 
! AR2=2*DLT
! 
! �������C=(B3*6)*(D3*3)
	DO 22 I=1,6
	DO 22 J=1,3
	C(I,J)=0.0
	DO 22 K=1,3
22	C(I,J)=C(I,J)+B(K,I)*D(K,J)
! 
! �������ESM=DTDB=CB
! 
	DO 27 I=1,6
	DO 27 J=1,6
	SUM=0.0
	DO 28 K=1,3
28	SUM=SUM+C(I,K)*B(K,J)
		ESM(I,J)=SUM*TH/(2.*AR2)
27	CONTINUE
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!�����غɽڵ㴦���غ�ֵ��λ�Ʊ߽�ڵ㴦��λ��ֵ��������նȾ���λ������ͽڵ���
! ��������޸�
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE MODIFY
	 	COMMON/AV/A(8500),JGF,JGSM,NP,NBW,JEND
	 	COMMON /PRO/IRATION,IROWNUM,ICOLNUM,NN,NE
	 	DATA IN/60/,IO/61/
	 	INTEGER :: ROW
		INTEGER :: COL
		INTEGER :: STA
!
!����ڵ�ļ����غɣ��ŵ�A�����нڵ�����������Ӧλ��
!IBΪʩ�����غɽڵ�����ɶ�(X,Y),BVΪ�����ɶ��ϵ��غ�ֵ
202	READ(IN,*)IB	
	IF(IB.LE.0) THEN
    ROW=0
	COL=0
    GOTO 208
    ENDIF
	READ(IN,*)BV
	IF(MOD(IB,2).EQ.1) GOTO 204
	WRITE(IO,203) IB/2,BV
203	FORMAT('�ڵ�',I3,'�غ�Ϊ:PY=',F11.3)
	GOTO 206
204	WRITE(IO,205) IB/2+1,BV
205	FORMAT('�ڵ�',I3,'�غ�Ϊ:PX=',F8.3)
206		A(JGF+IB)=A(JGF+IB)+BV
		GOTO 202

!
!����λ�Ʊ߽�ڵ㴦��λ��ֵ���ŵ�A�����нڵ�λ�����������λ��
!IBΪ�ڵ�λ�Ƶ����ɶȣ�BVΪλ��ֵ
!IBΪ������ʱ����x��ż������y
	
208	IF(ROW.LT.IROWNUM) GOTO 11
!�����һ���ŵ�y����Լ��
	ROW=IROWNUM
	COL=COL+1
	STA=ROW*(ROW-1)/2.0+COL
	IB=2*STA
	GOTO 12
!����һ����x����Լ��Ϊ0
11	COL=1
	ROW=ROW+1
	STA=ROW*(ROW-1)/2.0+COL
	IB=2*STA-1
    COL=0
12	IF(COL.GT.ICOLNUM) RETURN
    BV=0.0 
    IF(MOD(IB,2).EQ.1) GOTO 214
    WRITE(IO,213) (IB+1)/2,BV
213	FORMAT('�ڵ�',I3,'λ��Լ��Ϊ:V=',F8.3)
	GOTO 209
214	WRITE(IO,215)(IB+1)/2,BV
215	FORMAT('�ڵ�',I3,'λ��Լ��Ϊ:U=',F8.3)
209		K=IB-1
		DO 211 J=2,NBW
			M=IB+J-1
			IF(M.GT.NP) GOTO 210
			IJ=JGSM+(J-1)*NP+IB-(J-1)*(J-2)/2
			A(JGF+M)=A(JGF+M)-A(IJ)*BV
			A(IJ)=0.0
210	IF(K.LE.0) GOTO 211
		KJ=JGSM+(J-1)*NP+K-(J-1)*(J-2)/2
		A(JGF+k)=A(JGF+K)-A(KJ)*BV
		A(KJ)=0.0
		K=K-1
211	CONTINUE
		A(JGF+IB)=A(JGSM+IB)*BV
221	CONTINUE
		GOTO 208
		END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     GUASS  ��Ԫ�ӳ���
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE DCMPBD
		COMMON /AV/A(8500),JGF,JGSM,NP,NBW,JEND
!
!�ø�˹��Ԫ�����ԳƵȴ���նȾ���Ϊ�����Ǿ���
!
		NP1=NP-1
		DO 226 I=1,NP1
		MJ=I+NBW-1
		IF(MJ.GT.NP) MJ=NP
		NJ=I+1
		MK=NBW
		IF((NP-I+1).LT.NBW) MK=NP-I+1
		ND=0
		DO 225 J=NJ,MJ
		MK=MK-1
		ND=ND+1
		NL=ND+1
		DO 225 K=1,MK
		NK=ND+k
		JK=JGSM+(K-1)*NP+J-(K-1)*(K-2)/2
		INL=JGSM+(NL-1)*NP+I-(NL-1)*(NL-2)/2
		INK=JGSM+(NK-1)*NP+I-(NK-1)*(NK-2)/2
		II=JGSM+I
225		A(JK)=A(JK)-A(INL)*A(INK)/A(II)
226		CONTINUE
			RETURN
			END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		�ش�����ӳ���
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SLVBD
		COMMON /AV/A(8500),JGF,JGSM,NP,NBW,JEND
		DATA IN/60/,IO/61/
		NP1=NP-1
!
! ���ӳ���DEMPBD�õ��������Ǿ�����лش����㣬�õ��ڵ�λ������DD
!
		DO 250 I=1,NP1
		MJ=I+NBW-1
		IF(MJ.GT.NP)MJ=NP
		NJ=I+1
		L=1
		DO 250 J=NJ,MJ
		L=L+1
		IL=JGSM+(L-1)*NP+I-(L-1)*(L-2)/2
250		A(JGF+J)=A(JGF+J)-A(IL)*A(JGF+I)/A(JGSM+I)
!
!���ڵ����ɶȵĻش����㣬�������ϵ���
!
		A(NP)=A(JGF+NP)/A(JGSM+NP)
		DO 252 K=1,NP1
		I=NP-K
		MJ=NBW
		IF((I+NBW-1).GT.NP) MJ=NP-I+1
		SUM=0.0
		DO 251 J=2,MJ
		N=I+J-1
		IJ=JGSM+(J-1)*NP+I-(J-1)*(J-2)/2
251	SUM=SUM+A(IJ)*A(N)
252	A(I)=(A(JGF+I)-SUM)/A(JGSM+I)
		RETURN
	END

    
