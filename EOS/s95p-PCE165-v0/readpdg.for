*     Reading the particle data
*
      I=1
      OPEN(30,FILE=respath(1:IFEOS2)//resfile,STATUS='old')

 159  CONTINUE

C----------Original version-----------------------
C        READ(30,*,END=166) numbe(I),name(I),mass(I),b,degen(I),
C     &               baryon(I),strange(I),d1,d2,d3,charge(I),decay(I),
C     &               npi(I), nK(I), neta(I)
C-------------------------------------------------

CSHEN------Change for pdg.dat---------------------
         READ(30,*,END=166) numbe(I),name(I),mass(I),b,degen(I),
     &                baryon(I),strange(I),d1,d2,d3,charge(I),decay(I)
         write(*,2000) I,numbe(I),name(I),mass(I),degen(I),decay(I),
     &      baryon(I)
CSHEN---------------------------------------------


C        WRITE(*,2000) I,numbe(I),name(I),mass(I),degen(I),decay(I),
C     &     baryon(I), npi(I), nK(I), neta(I)

*        WRITE(*,2005) numbe(I),name(I),mass(I),b,degen(I),
*     &               baryon(I),strange(I),d1,d2,d3,charge(I),decay(I),
*     &               npi(I), nK(I), neta(I)

        DO 165 J=1,decay(I)
          READ(30,*) decaychannel(I,J,1),decaychannel(I,J,2),
     &                branchratio(I,J),(decaychannel(I,J,K),K=4,8)
C          WRITE(*,2006) (decaychannel(I,J,K),K=1,8)
 165    CONTINUE

* Here a little loop as there are no antibaryons in the table
        IF (baryon(I).NE.0) THEN
          I=I+1
          numbe(I)  =-numbe(I-1)
          name(I)   = 'Anti-'//name(I-1) 
          mass(I)   = mass(I-1)
          degen(I)  = degen(I-1)
          baryon(I) =-baryon(I-1)
          strange(I)=-strange(I-1)
          charge(I) =-charge(I-1)
          decay(I)  = decay(I-1)
C          npi(I)    = npi(I-1)
C          nK(I)     = nK(I-1)
C          neta(I)   = neta(I-1)
C          WRITE(*,2000) I,-numbe(I),name(I),mass(I),degen(I),decay(I),
C     &     baryon(I), npi(I), nK(I), neta(I)

          write(*,2000) I, -numbe(I),name(I),mass(I),degen(I),decay(I),
     &      baryon(I)    !CSHEN: for pdg.dat
          do 167 J=1,decay(I)
            decaychannel(I,J,1) = -decaychannel(I-1,J,1)
            decaychannel(I,J,2) = decaychannel(I-1,J,2)
            branchratio(I,J) = branchratio(I-1,J)
            decaychannel(I,J,4) = decaychannel(I-1,J,4)
            decaychannel(I,J,5) = decaychannel(I-1,J,5)
            decaychannel(I,J,6) = decaychannel(I-1,J,6)
            decaychannel(I,J,7) = decaychannel(I-1,J,7)
            decaychannel(I,J,8) = decaychannel(I-1,J,8)
C            write(*,2006) (decaychannel(I,J,K),K=1,8)
167       continue
            
*          WRITE(*,2005) numbe(I),name(I),mass(I),b,degen(I),
*     &         baryon(I),strange(I),d1,d2,d3,charge(I),decay(I),
*     &         npi(I), nK(I), neta(I)

        END IF
*
        I=I+1
        GOTO 159
 166    CONTINUE

        nres=I-1 
*        WRITE(*,*) nres
      CLOSE(30)

 2000 FORMAT(I9,I9,A25,F10.6,I5,I5,I3,3F10.6)
 2005 FORMAT(I9,A25,F10.6,F10.6,I5,I5,I5,I5,I5,I5,I5,I5,F10.6,F10.6,
     &     F10.6)
 2006 FORMAT(I9,I9,F10.6,I5,I5,I5,I5,I5)

CSHEN=================================================================
C======read the chemical potential at freezeout surface for EOS6
         IF (IEOS.eq.6) then
           OPEN(20,FILE='results/decdat_mu.dat',STATUS='OLD')
           do I = 1, datalength
             READ(20,*) (XMUfreeze6(I,J),J=1,ICOLNUMeos6)
           end do
           CLOSE(20)

           do 130 I = 1, datalength  !H.C.  
           do 130 J = 1,maxpar
               CXMUpar(I,J) = 0.0d0    !H.C. initialization for XMUpar
130        continue
C---------------------------------------------------------
C------for stable particles-------------------------------
           do 1310 I = 1, datalength  !* H.C.
           do 131 J = 2,9
              CXMUpar(I,J) = XMUfreeze6(I,J) !* H.C.
131        continue
               CXMUpar(I,18)  = XMUfreeze6(I,10) !* H.C.
               CXMUpar(I,19)  = XMUfreeze6(I,11)
               CXMUpar(I,20)  = XMUfreeze6(I,12)
               CXMUpar(I,21)  = XMUfreeze6(I,13)
               CXMUpar(I,22)  = XMUfreeze6(I,14)
               CXMUpar(I,27)  = XMUfreeze6(I,15)
               CXMUpar(I,28)  = XMUfreeze6(I,16)
               CXMUpar(I,29)  = XMUfreeze6(I,17)
               CXMUpar(I,31)  = XMUfreeze6(I,18)
               CXMUpar(I,32)  = XMUfreeze6(I,19)
               CXMUpar(I,33)  = XMUfreeze6(I,20)
               CXMUpar(I,34)  = XMUfreeze6(I,21)
               CXMUpar(I,35)  = XMUfreeze6(I,22)
               CXMUpar(I,36)  = XMUfreeze6(I,23)
               CXMUpar(I,61)  = XMUfreeze6(I,24)
               CXMUpar(I,62)  = XMUfreeze6(I,25)
               CXMUpar(I,63)  = XMUfreeze6(I,26)
               CXMUpar(I,64)  = XMUfreeze6(I,27)
               CXMUpar(I,171) = XMUfreeze6(I,28)
               CXMUpar(I,172) = XMUfreeze6(I,29)
 1310      continue      
C----end--------------------------------------------------
       write(*,*) 'add chemcial potential for stable particles succeed!'
       write(*,*) 'calculate resonance particles chemical potential'


      do 1320 KK= 1, datalength !* H.C.
         do 132 I = 2, nres
            if (CXMUpar(KK,I) .eq. 0.0) then !*H.C
               do 133 J = 1, decay(I)
                  do 134 K = 1, abs(decaychannel(I,J,2))
                     do 135 L = 1, nres
                        if (decaychannel(I,J,K+3) .eq. numbe(L)) then
                           CXMUpar(KK,I) = CXMUpar(KK,I) 
    &                           + branchratio(I,J)*CXMUpar(KK,L) !* H.C.
                            goto 1335
                         endif
                        if (L .eq. nres) then
                           write(*,*) 'warning: can not find particle ',
    &                           name(I),decaychannel(I,J,K+3)
                        endif
135                  continue
1335                 continue
134               continue
133            continue 
            endif
132      continue
1320  continue 
      
      do 136 I = 1, nres
        write(*,*) name(I), CXMUpar(1,I)
136   continue

      write(*,*)'The chemical pontential of hadron resonance added! :-)'
      
      endif
