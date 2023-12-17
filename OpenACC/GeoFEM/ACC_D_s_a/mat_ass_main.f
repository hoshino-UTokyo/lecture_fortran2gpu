!C
!C***
!C*** MAT_ASS_MAIN
!C***
!C
      subroutine MAT_ASS_MAIN 
      use hpcmw_all
      implicit REAL*8 (A-H,O-Z)
      integer(kind=kint), dimension(  8) :: nodLOCAL
      real(8),allocatable,dimension(:,:,:,:) :: aa
      real(8),allocatable,dimension(:,:,:,:,:) :: PNX2,PNY2,PNZ2
      real(8),allocatable,dimension(:,:,:,:) :: DETJ2
      integer(kind=kint),allocatable,dimension(:,:,:) :: kk2,IDlu2
      

      call MPI_Barrier (SOLVER_COMM, ierr)
      S1_time= MPI_WTIME()
!C
!C +------------------+
!C | ELEMENT Coloring |
!C +------------------+
!C===
      allocate (ELMCOLORindex(0:NP))
      allocate (ELMCOLORitem (ICELTOT))
      if (allocated (IWKX)) deallocate (IWKX)
      allocate (IWKX(0:NP,3))

      IWKX= 0
      icou= 0
      do icol= 1, NP
        do i= 1, NP
          IWKX(i,1)= 0
        enddo
        do icel= 1, ICELTOT
          if (IWKX(icel,2).eq.0) then
            in1= ICELNOD(icel,1)
            in2= ICELNOD(icel,2)
            in3= ICELNOD(icel,3)
            in4= ICELNOD(icel,4)
            in5= ICELNOD(icel,5)
            in6= ICELNOD(icel,6)
            in7= ICELNOD(icel,7)
            in8= ICELNOD(icel,8)

            ip1= IWKX(in1,1)
            ip2= IWKX(in2,1)
            ip3= IWKX(in3,1)
            ip4= IWKX(in4,1)
            ip5= IWKX(in5,1)
            ip6= IWKX(in6,1)
            ip7= IWKX(in7,1)
            ip8= IWKX(in8,1)

            isum= ip1 + ip2 + ip3 + ip4 + ip5 + ip6 + ip7 + ip8
            if (isum.eq.0) then
              icou= icou + 1
              IWKX(icol,3)= icou
              IWKX(icel,2)= icol
              ELMCOLORitem(icou)= icel

              IWKX(in1,1)= 1
              IWKX(in2,1)= 1
              IWKX(in3,1)= 1
              IWKX(in4,1)= 1
              IWKX(in5,1)= 1
              IWKX(in6,1)= 1
              IWKX(in7,1)= 1
              IWKX(in8,1)= 1
              if (icou.eq.ICELTOT) goto 100            
            endif
          endif
        enddo
      enddo

 100  continue
      ELMCOLORtot= icol
      IWKX(0          ,3)= 0
      IWKX(ELMCOLORtot,3)= ICELTOT

      do icol= 0, ELMCOLORtot
        ELMCOLORindex(icol)= IWKX(icol,3)
      enddo

!      write (*,'(a,2i8)') '### Number of Element Colors', 
!     &                     my_rank, ELMCOLORtot
      deallocate (IWKX)
      X1_time= MPI_WTIME()
!C===

      WEI(1)= +1.0000000000D+00
      WEI(2)= +1.0000000000D+00

      POS(1)= -0.5773502692D+00
      POS(2)= +0.5773502692D+00

      if (allocated(AL)) deallocate (AL)
      if (allocated(AU)) deallocate (AU)
      if (allocated(D )) deallocate (D )
      if (allocated(B )) deallocate (B )
      if (allocated(X )) deallocate (X )
      if (allocated(WW)) deallocate (WW)
      if (allocated(ALUG)) deallocate (ALUG)
      if (allocated(WS)) deallocate (WR)

      allocate (AL(9*NPL), AU(9*NPU), ALUG(9*NP))
      allocate (D(9*NP), B(3*NP), X(3*NP), WW(3*NP,4))
      allocate (WS(3*NP), WR(3*NP))
      allocate (aa(9,8,8,ICELTOT))
      allocate (DETJ2(2,2,2,ICELTOT))
      allocate (PNX2(2,2,2,8,ICELTOT))
      allocate (PNY2(2,2,2,8,ICELTOT))
      allocate (PNZ2(2,2,2,8,ICELTOT))
      allocate (IDlu2(8,8,ICELTOT))
      allocate (kk2(8,8,ICELTOT))
!$acc enter data
!$acc& create(AL,AU,ALUG,D,B,X,WW,WS,WR)
!$acc& create(aa,DETJ2,PNX2,PNY2,PNZ2,IDlu2,kk2)
!$acc& copyin(indexL,itemL,indexU,itemU)
!$acc& copyin(ICELNOD,OLDtoNEW,ELMCOLORindex,ELMCOLORitem)

      if (FTflag.eq.1) then
c$$$          do jv= 1, NHYP
c$$$            jS= IVECT(jv-1) + 1
c$$$            jE= IVECT(jv)
c$$$!$omp parallel do private(ip,j, jS, jE)
c$$$            do ip= 1, PEsmpTOT
c$$$              jS= STACKmc(ip-1,jv) + 1
c$$$              jE= STACKmc(ip  ,jv)
c$$$              do j= jS, jE

!$acc kernels
!$acc loop independent         
         do j=1,N
              X(3*j-2)= 0.d0
              X(3*j-1)= 0.d0
              X(3*j  )= 0.d0
              B(3*j-2)= 0.d0
              B(3*j-1)= 0.d0
              B(3*j  )= 0.d0
              D(9*j-8)= 0.d0
              D(9*j-7)= 0.d0
              D(9*j-6)= 0.d0
              D(9*j-5)= 0.d0
              D(9*j-4)= 0.d0
              D(9*j-3)= 0.d0
              D(9*j-2)= 0.d0
              D(9*j-1)= 0.d0
              D(9*j  )= 0.d0
              ALUG(9*j  )= 0.d0
              ALUG(9*j-8)= 0.d0
              ALUG(9*j-7)= 0.d0
              ALUG(9*j-6)= 0.d0
              ALUG(9*j-5)= 0.d0
              ALUG(9*j-4)= 0.d0
              ALUG(9*j-3)= 0.d0
              ALUG(9*j-2)= 0.d0
              ALUG(9*j-1)= 0.d0
              ALUG(9*j  )= 0.d0
c$$$              WW(3*j-2,1)= 0.d0
c$$$              WW(3*j-1,1)= 0.d0
c$$$              WW(3*j  ,1)= 0.d0
c$$$              WW(3*j-2,2)= 0.d0
c$$$              WW(3*j-1,2)= 0.d0
c$$$              WW(3*j  ,2)= 0.d0
c$$$              WW(3*j-2,3)= 0.d0
c$$$              WW(3*j-1,3)= 0.d0
c$$$              WW(3*j  ,3)= 0.d0
c$$$              WW(3*j-2,4)= 0.d0
c$$$              WW(3*j-1,4)= 0.d0
c$$$              WW(3*j  ,4)= 0.d0
           enddo
!$acc end kernels
c$$$            enddo
c$$$          enddo
       else
!$acc kernels
          B(:)= 0.d0
          X(:)= 0.d0
!$acc end kernels
!$acc kernels
          D(:)= 0.d0
          ALUG(:)= 0.d0
!$acc end kernels
      endif

!$acc kernels
      WW(:,:)= 0.d0
!$acc end kernels

      if (FTflag.eq.1) then
c$$$          do jv= 1, NHYP
c$$$            jS= IVECT(jv-1) + 1
c$$$            jE= IVECT(jv)
c$$$!$omp parallel do private(jS,jE,j,jsL,jeL,jsU,jeU,k)
c$$$            do ip= 1, PEsmpTOT
c$$$              jS= STACKmc(ip-1,jv) + 1
c$$$              jE= STACKmc(ip  ,jv)
c$$$  do j = jS, jE

!$acc kernels
!$acc loop independent         
         do j = 1,N 
              jsL= indexL(j-1)+1
              jeL= indexL(j)
              do k= jsL, jeL
                AL(9*k-8)= 0.d0
                AL(9*k-7)= 0.d0
                AL(9*k-6)= 0.d0
                AL(9*k-5)= 0.d0
                AL(9*k-4)= 0.d0
                AL(9*k-3)= 0.d0
                AL(9*k-2)= 0.d0
                AL(9*k-1)= 0.d0
                AL(9*k  )= 0.d0
              enddo

              jsU= indexU(j-1)+1
              jeU= indexU(j)
              do k= jsU, jeU
                AU(9*k-8)= 0.d0
                AU(9*k-7)= 0.d0
                AU(9*k-6)= 0.d0
                AU(9*k-5)= 0.d0
                AU(9*k-4)= 0.d0
                AU(9*k-3)= 0.d0
                AU(9*k-2)= 0.d0
                AU(9*k-1)= 0.d0
                AU(9*k  )= 0.d0
              enddo
            enddo
!$acc end kernels

c$$$  enddo
c$$$          enddo
       else
!$acc kernels
          AL(:)= 0.d0
          AU(:)= 0.d0
!$acc end kernels
      endif

!C
!C-- INIT.
!C     PNQ   - 1st-order derivative of shape function by QSI
!C     PNE   - 1st-order derivative of shape function by ETA
!C     PNT   - 1st-order derivative of shape function by ZET
!C
      valA=            POISSON /      (1.d0-POISSON)
      valB= (1.d0-2.d0*POISSON)/(2.d0*(1.d0-POISSON))
      valX= ELAST*(1.d0-POISSON)/((1.d0+POISSON)*(1.d0-2.d0*POISSON))

      valA= valA * valX
      valB= valB * valX

      do kp= 1, 2
      do jp= 1, 2
      do ip= 1, 2
        QP1= 1.d0 + POS(ip)
        QM1= 1.d0 - POS(ip)
        EP1= 1.d0 + POS(jp)
        EM1= 1.d0 - POS(jp)
        TP1= 1.d0 + POS(kp)
        TM1= 1.d0 - POS(kp)
        SHAPE(ip,jp,kp,1)= O8th * QM1 * EM1 * TM1
        SHAPE(ip,jp,kp,2)= O8th * QP1 * EM1 * TM1
        SHAPE(ip,jp,kp,3)= O8th * QP1 * EP1 * TM1
        SHAPE(ip,jp,kp,4)= O8th * QM1 * EP1 * TM1
        SHAPE(ip,jp,kp,5)= O8th * QM1 * EM1 * TP1
        SHAPE(ip,jp,kp,6)= O8th * QP1 * EM1 * TP1
        SHAPE(ip,jp,kp,7)= O8th * QP1 * EP1 * TP1
        SHAPE(ip,jp,kp,8)= O8th * QM1 * EP1 * TP1
        PNQ(jp,kp,1)= - O8th * EM1 * TM1
        PNQ(jp,kp,2)= + O8th * EM1 * TM1
        PNQ(jp,kp,3)= + O8th * EP1 * TM1
        PNQ(jp,kp,4)= - O8th * EP1 * TM1
        PNQ(jp,kp,5)= - O8th * EM1 * TP1
        PNQ(jp,kp,6)= + O8th * EM1 * TP1
        PNQ(jp,kp,7)= + O8th * EP1 * TP1
        PNQ(jp,kp,8)= - O8th * EP1 * TP1
        PNE(ip,kp,1)= - O8th * QM1 * TM1
        PNE(ip,kp,2)= - O8th * QP1 * TM1
        PNE(ip,kp,3)= + O8th * QP1 * TM1
        PNE(ip,kp,4)= + O8th * QM1 * TM1
        PNE(ip,kp,5)= - O8th * QM1 * TP1
        PNE(ip,kp,6)= - O8th * QP1 * TP1
        PNE(ip,kp,7)= + O8th * QP1 * TP1
        PNE(ip,kp,8)= + O8th * QM1 * TP1
        PNT(ip,jp,1)= - O8th * QM1 * EM1
        PNT(ip,jp,2)= - O8th * QP1 * EM1
        PNT(ip,jp,3)= - O8th * QP1 * EP1
        PNT(ip,jp,4)= - O8th * QM1 * EP1
        PNT(ip,jp,5)= + O8th * QM1 * EM1
        PNT(ip,jp,6)= + O8th * QP1 * EM1
        PNT(ip,jp,7)= + O8th * QP1 * EP1
        PNT(ip,jp,8)= + O8th * QM1 * EP1
      enddo
      enddo
      enddo

      call MPI_Barrier (SOLVER_COMM, ierr)
      S2_time= MPI_WTIME()

!      do icol= 1, ELMCOLORtot
c$$$!$omp parallel do private (icel0,icel,in1,in2,in3,in4,in5,in6,in7,in8)  &
c$$$!$omp&            private (in10,in20,in30,in40,in50,in60,in70,in80)     &
c$$$!$omp&            private (nodLOCAL,ie,je,ip,jp,kk,iiS,iiE,iDlu,k)      &
c$$$!$omp&            private (PNXi,PNYi,PNZi,PNXj,PNYj,PNZj,a11,a12)       &
c$$$!$omp&            private (a13,a21,a22,a23,a31,a32,a33,ipn,jpn,kpn,coef)&
c$$$!$omp&            private (DETJ,PNX,PNY,PNZ,X1,X2,X3,X4,X5,X6,X7,X8)    &
c$$$!$omp&            private (Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)                     &
c$$$!$omp&            private (Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8)
c$$$      do icel0= ELMCOLORindex(icol-1)+1, ELMCOLORindex(icol)

!$acc parallel 
!$acc loop independent collapse(4)
      do icel0=ELMCOLORindex(0)+1,ELMCOLORindex(ELMCOLORtot) 

      do kp= 1, 2
      do jp= 1, 2
      do ip= 1, 2

        X1= 0.d0
        Y1= 0.d0
        Z1= 0.d0

        X2= DX
        Y2= 0.d0
        Z2= 0.d0

        X3= DX
        Y3= DY
        Z3= 0.d0

        X4= 0.d0
        Y4= DY
        Z4= 0.d0

        X5= 0.d0
        Y5= 0.d0
        Z5= DZ

        X6= DX
        Y6= 0.d0
        Z6= DZ

        X7= DX
        Y7= DY
        Z7= DZ

        X8= 0.d0
        Y8= DY
        Z8= DZ

!C     
!C==   DETERMINANT of the JACOBIAN
        dXdQ =                                                          &
     &           + PNQ(jp,kp,1) * X1 + PNQ(jp,kp,2) * X2                &
     &           + PNQ(jp,kp,3) * X3 + PNQ(jp,kp,4) * X4                &
     &           + PNQ(jp,kp,5) * X5 + PNQ(jp,kp,6) * X6                &
     &           + PNQ(jp,kp,7) * X7 + PNQ(jp,kp,8) * X8                
        dYdQ =                                                          &
     &           + PNQ(jp,kp,1) * Y1 + PNQ(jp,kp,2) * Y2                &
     &           + PNQ(jp,kp,3) * Y3 + PNQ(jp,kp,4) * Y4                &
     &           + PNQ(jp,kp,5) * Y5 + PNQ(jp,kp,6) * Y6                &
     &           + PNQ(jp,kp,7) * Y7 + PNQ(jp,kp,8) * Y8                
        dZdQ =                                                          &
     &           + PNQ(jp,kp,1) * Z1 + PNQ(jp,kp,2) * Z2                &
     &           + PNQ(jp,kp,3) * Z3 + PNQ(jp,kp,4) * Z4                &
     &           + PNQ(jp,kp,5) * Z5 + PNQ(jp,kp,6) * Z6                &
     &           + PNQ(jp,kp,7) * Z7 + PNQ(jp,kp,8) * Z8                
        dXdE =                                                          &
     &           + PNE(ip,kp,1) * X1 + PNE(ip,kp,2) * X2                &
     &           + PNE(ip,kp,3) * X3 + PNE(ip,kp,4) * X4                &
     &           + PNE(ip,kp,5) * X5 + PNE(ip,kp,6) * X6                &
     &           + PNE(ip,kp,7) * X7 + PNE(ip,kp,8) * X8
        dYdE =                                                          &
     &           + PNE(ip,kp,1) * Y1 + PNE(ip,kp,2) * Y2                &
     &           + PNE(ip,kp,3) * Y3 + PNE(ip,kp,4) * Y4                &
     &           + PNE(ip,kp,5) * Y5 + PNE(ip,kp,6) * Y6                &
     &           + PNE(ip,kp,7) * Y7 + PNE(ip,kp,8) * Y8
        dZdE =                                                          &
     &           + PNE(ip,kp,1) * Z1 + PNE(ip,kp,2) * Z2                &
     &           + PNE(ip,kp,3) * Z3 + PNE(ip,kp,4) * Z4                &
     &           + PNE(ip,kp,5) * Z5 + PNE(ip,kp,6) * Z6                &
     &           + PNE(ip,kp,7) * Z7 + PNE(ip,kp,8) * Z8
        dXdT =                                                          &
     &           + PNT(ip,jp,1) * X1 + PNT(ip,jp,2) * X2                &
     &           + PNT(ip,jp,3) * X3 + PNT(ip,jp,4) * X4                &
     &           + PNT(ip,jp,5) * X5 + PNT(ip,jp,6) * X6                &
     &           + PNT(ip,jp,7) * X7 + PNT(ip,jp,8) * X8
        dYdT =                                                          &
     &           + PNT(ip,jp,1) * Y1 + PNT(ip,jp,2) * Y2                &
     &           + PNT(ip,jp,3) * Y3 + PNT(ip,jp,4) * Y4                &
     &           + PNT(ip,jp,5) * Y5 + PNT(ip,jp,6) * Y6                &
     &           + PNT(ip,jp,7) * Y7 + PNT(ip,jp,8) * Y8
        dZdT =                                                          &
     &           + PNT(ip,jp,1) * Z1 + PNT(ip,jp,2) * Z2                &
     &           + PNT(ip,jp,3) * Z3 + PNT(ip,jp,4) * Z4                &
     &           + PNT(ip,jp,5) * Z5 + PNT(ip,jp,6) * Z6                &
     &           + PNT(ip,jp,7) * Z7 + PNT(ip,jp,8) * Z8

        DETJ2(ip,jp,kp,icel0)= dXdQ*(dYdE*dZdT-dZdE*dYdT) +             &
     &                  dYdQ*(dZdE*dXdT-dXdE*dZdT) +                    &
     &                  dZdQ*(dXdE*dYdT-dYdE*dXdT)

!C
!C==   INVERSE JACOBIAN
        coef= 1.d0 / DETJ2(ip,jp,kp,icel0)
        a11= coef * ( dYdE*dZdT - dZdE*dYdT )
        a12= coef * ( dZdQ*dYdT - dYdQ*dZdT )
        a13= coef * ( dYdQ*dZdE - dZdQ*dYdE )

        a21= coef * ( dZdE*dXdT - dXdE*dZdT )
        a22= coef * ( dXdQ*dZdT - dZdQ*dXdT )
        a23= coef * ( dZdQ*dXdE - dXdQ*dZdE )

        a31= coef * ( dXdE*dYdT - dYdE*dXdT )
        a32= coef * ( dYdQ*dXdT - dXdQ*dYdT )
        a33= coef * ( dXdQ*dYdE - dYdQ*dXdE )

        DETJ2(ip,jp,kp,icel0)= dabs(DETJ2(ip,jp,kp,icel0))

!C
!C== set the dNi/dX, dNi/dY & dNi/dZ components
        PNX2(ip,jp,kp,1,icel0)= a11*PNQ(jp,kp,1) + a12*PNE(ip,kp,1) +          &
     &                   a13*PNT(ip,jp,1)
        PNX2(ip,jp,kp,2,icel0)= a11*PNQ(jp,kp,2) + a12*PNE(ip,kp,2) +          &
     &                   a13*PNT(ip,jp,2)
        PNX2(ip,jp,kp,3,icel0)= a11*PNQ(jp,kp,3) + a12*PNE(ip,kp,3) +          &
     &                   a13*PNT(ip,jp,3)
        PNX2(ip,jp,kp,4,icel0)= a11*PNQ(jp,kp,4) + a12*PNE(ip,kp,4) +          &
     &                   a13*PNT(ip,jp,4)
        PNX2(ip,jp,kp,5,icel0)= a11*PNQ(jp,kp,5) + a12*PNE(ip,kp,5) +          &
     &                   a13*PNT(ip,jp,5)
        PNX2(ip,jp,kp,6,icel0)= a11*PNQ(jp,kp,6) + a12*PNE(ip,kp,6) +          &
     &                   a13*PNT(ip,jp,6)
        PNX2(ip,jp,kp,7,icel0)= a11*PNQ(jp,kp,7) + a12*PNE(ip,kp,7) +          &
     &                   a13*PNT(ip,jp,7)
        PNX2(ip,jp,kp,8,icel0)= a11*PNQ(jp,kp,8) + a12*PNE(ip,kp,8) +          &
     &                   a13*PNT(ip,jp,8)

        PNY2(ip,jp,kp,1,icel0)= a21*PNQ(jp,kp,1) + a22*PNE(ip,kp,1) +          &
     &                   a23*PNT(ip,jp,1)
        PNY2(ip,jp,kp,2,icel0)= a21*PNQ(jp,kp,2) + a22*PNE(ip,kp,2) +          &
     &                   a23*PNT(ip,jp,2)
        PNY2(ip,jp,kp,3,icel0)= a21*PNQ(jp,kp,3) + a22*PNE(ip,kp,3) +          &
     &                   a23*PNT(ip,jp,3)
        PNY2(ip,jp,kp,4,icel0)= a21*PNQ(jp,kp,4) + a22*PNE(ip,kp,4) +          &
     &                   a23*PNT(ip,jp,4)
        PNY2(ip,jp,kp,5,icel0)= a21*PNQ(jp,kp,5) + a22*PNE(ip,kp,5) +          &
     &                   a23*PNT(ip,jp,5)
        PNY2(ip,jp,kp,6,icel0)= a21*PNQ(jp,kp,6) + a22*PNE(ip,kp,6) +          &
     &                   a23*PNT(ip,jp,6)
        PNY2(ip,jp,kp,7,icel0)= a21*PNQ(jp,kp,7) + a22*PNE(ip,kp,7) +          &
     &                   a23*PNT(ip,jp,7)
        PNY2(ip,jp,kp,8,icel0)= a21*PNQ(jp,kp,8) + a22*PNE(ip,kp,8) +          &
     &                   a23*PNT(ip,jp,8)

        PNZ2(ip,jp,kp,1,icel0)= a31*PNQ(jp,kp,1) + a32*PNE(ip,kp,1) +          &
     &                   a33*PNT(ip,jp,1)
        PNZ2(ip,jp,kp,2,icel0)= a31*PNQ(jp,kp,2) + a32*PNE(ip,kp,2) +          &
     &                   a33*PNT(ip,jp,2)
        PNZ2(ip,jp,kp,3,icel0)= a31*PNQ(jp,kp,3) + a32*PNE(ip,kp,3) +          &
     &                   a33*PNT(ip,jp,3)
        PNZ2(ip,jp,kp,4,icel0)= a31*PNQ(jp,kp,4) + a32*PNE(ip,kp,4) +          &
     &                   a33*PNT(ip,jp,4)
        PNZ2(ip,jp,kp,5,icel0)= a31*PNQ(jp,kp,5) + a32*PNE(ip,kp,5) +          &
     &                   a33*PNT(ip,jp,5)
        PNZ2(ip,jp,kp,6,icel0)= a31*PNQ(jp,kp,6) + a32*PNE(ip,kp,6) +          &
     &                   a33*PNT(ip,jp,6)
        PNZ2(ip,jp,kp,7,icel0)= a31*PNQ(jp,kp,7) + a32*PNE(ip,kp,7) +          &
     &                   a33*PNT(ip,jp,7)
        PNZ2(ip,jp,kp,8,icel0)= a31*PNQ(jp,kp,8) + a32*PNE(ip,kp,8) +          &
     &                   a33*PNT(ip,jp,8)

      enddo
      enddo
      enddo
      enddo
!$acc end parallel 
        
!$acc parallel vector_length(64)
!$acc loop independent gang
      do icel0=ELMCOLORindex(0)+1,ELMCOLORindex(ELMCOLORtot) 

!C
!C== CONSTRUCT the GLOBAL MATRIX
!$acc loop independent collapse(2) vector
         do ie=1,8
         do je=1,8
            
            PNXi= 0.d0
            PNYi= 0.d0
            PNZi= 0.d0
            PNXj= 0.d0
            PNYj= 0.d0
            PNZj= 0.d0
            
            a11= 0.d0
            a21= 0.d0
            a31= 0.d0
            a12= 0.d0
            a22= 0.d0
            a32= 0.d0
            a13= 0.d0
            a23= 0.d0
            a33= 0.d0
!$acc loop seq
            do kpn= 1, 2
!$acc loop seq
            do jpn= 1, 2
!$acc loop seq
            do ipn= 1, 2
             coef= dabs(DETJ2(ipn,jpn,kpn,icel0))                       &
     &                    *     WEI(ipn)*WEI(jpn)*WEI(kpn)

             PNXi= PNX2(ipn,jpn,kpn,ie,icel0)
             PNYi= PNY2(ipn,jpn,kpn,ie,icel0)
             PNZi= PNZ2(ipn,jpn,kpn,ie,icel0)
             
             PNXj= PNX2(ipn,jpn,kpn,je,icel0)
             PNYj= PNY2(ipn,jpn,kpn,je,icel0)
             PNZj= PNZ2(ipn,jpn,kpn,je,icel0)
             
             a11= a11 + (valX*PNXi*PNXj+valB*(PNYi*PNYj+PNZi*PNZj))*coef
             a22= a22 + (valX*PNYi*PNYj+valB*(PNZi*PNZj+PNXi*PNXj))*coef
             a33= a33 + (valX*PNZi*PNZj+valB*(PNXi*PNXj+PNYi*PNYj))*coef
             
             a12= a12 + (valA*PNXi*PNYj + valB*PNXj*PNYi)*coef
             a13= a13 + (valA*PNXi*PNZj + valB*PNXj*PNZi)*coef
             a21= a21 + (valA*PNYi*PNXj + valB*PNYj*PNXi)*coef
             a23= a23 + (valA*PNYi*PNZj + valB*PNYj*PNZi)*coef
             a31= a31 + (valA*PNZi*PNXj + valB*PNZj*PNXi)*coef
             a32= a32 + (valA*PNZi*PNYj + valB*PNZj*PNYi)*coef
           enddo
           enddo
           enddo
           aa(1,je,ie,icel0) = a11
           aa(2,je,ie,icel0) = a12
           aa(3,je,ie,icel0) = a13
           aa(4,je,ie,icel0) = a21
           aa(5,je,ie,icel0) = a22
           aa(6,je,ie,icel0) = a23
           aa(7,je,ie,icel0) = a31
           aa(8,je,ie,icel0) = a32
           aa(9,je,ie,icel0) = a33
        enddo
        enddo
      enddo
!$acc end parallel 

!$acc parallel vector_length(64)
!$acc loop independent gang
      do icel0=ELMCOLORindex(0)+1,ELMCOLORindex(ELMCOLORtot) 
!C
!C== CONSTRUCT the GLOBAL MATRIX
!$acc loop independent vector collapse(2)
!$acc& private(ip,jp,IDlu,kk,iiS,iiE,k)
         do ie=1,8
         do je=1,8
            ip = OLDtoNEW(ICELNOD(ELMCOLORitem(icel0),ie))
            if (ip.le.N) then
               jp = OLDtoNEW(ICELNOD(ELMCOLORitem(icel0),je))
               
               IDlu= 0
               if (ip.eq.jp) IDlu= 0
               
               kk= 0
               iiS= indexU(ip-1) + 1
               iiE= indexU(ip  )
               do k= iiS, iiE
                  if ( itemU(k).eq.jp ) then
                     kk  = k
                     IDlu= 1
                     exit
                  endif
               enddo
               
               if (kk.eq.0) then
                  iiS= indexL(ip-1) + 1
                  iiE= indexL(ip  )
                  
                  do k= iiS, iiE
                     if ( itemL(k).eq.jp) then
                        kk= k
                        IDlu= -1
                     endif
                  enddo
               endif
            endif
            kk2(je,ie,icel0) = kk
            IDlu2(je,ie,icel0) = IDlu
         enddo
         enddo
      enddo
!$acc end parallel 

!$acc parallel vector_length(192)
!$acc loop independent gang
      do icel0=ELMCOLORindex(0)+1,ELMCOLORindex(ELMCOLORtot) 
!$acc loop independent vector collapse(3) private(IDlu,kk,ip)
         do ie=1,8
         do je=1,8
         do ii=1,9
            IDlu=IDlu2(je,ie,icel0)
!     write (*,'(10i8)') icol,icel0,icel,ip,jp,IDlu,kk
            if (IDlu.eq.1) then
               kk=kk2(je,ie,icel0)-1
               !$acc atomic 
               AU(9*kk+ii)= AU(9*kk+ii) + aa(ii,je,ie,icel0)
               !$acc end atomic 
            else if (IDlu.eq.-1) then
               kk=kk2(je,ie,icel0)-1
               !$acc atomic 
               AL(9*kk+ii)= AL(9*kk+ii) + aa(ii,je,ie,icel0)
               !$acc end atomic 
            else if (IDlu.eq.0) then
               ip = OLDtoNEW(ICELNOD(ELMCOLORitem(icel0),ie)) -1
               !$acc atomic 
               D(9*ip+ii)= D(9*ip+ii) + aa(ii,je,ie,icel0)
               !$acc end atomic 
            endif
         enddo
         enddo
         enddo
      enddo
!$acc end parallel 
!      enddo

      S3_time= MPI_WTIME()
      if (my_rank.eq.0) write (*,'(1pe16.6,a)') X1_time-S1_time,
     &                         '   coloring of elements'
      if (my_rank.eq.0) write (*,'(1pe16.6,a)') S3_time-S2_time,
     &                         '   matrix assembling'
      if (my_rank.eq.0) write (*,'(1pe16.6,a)')
     &     (2528.d0+37376.d0+1728.d0)
     &     * (ELMCOLORindex(ELMCOLORtot)-ELMCOLORindex(0)+1)
     &     / 1000000000.d0 / (S3_time-S2_time),
     &                         '   GFlops'
!$acc exit data delete(aa,DETJ2,PNX2,PNY2,PNZ2,IDlu2,kk2)
      
      return
      end
