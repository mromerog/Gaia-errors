c-----------------------------------------------------------------------
c  This is 'const_pot.h'
c  Constants related to galactic potential

      REAL*8 Mb,Md,Mh,K
      
c general values      
      PARAMETER (R0 = 8.5d0, !solar radius in kpc
     +           Usol = 10.00d0, ! solar velocity(Dehnen i Binney 1998) 
     +           Vsol = 5.25d0,
     +           Wsol = 7.17d0)
     
     
c values of the potential
      PARAMETER (ZMsol = 1.989d33,! Msol in gm
     +           abu = 0.387d0,! en Kpc
     +           ad = 5.3178d0,
     +           ah = 12.0d0,
     +           bd = 0.25d0,
     +           Mb = 1.406d10,!in Msol
     +           Md = 8.561d10,!in Msol
     +           Mh = 1.071d11)!in Msol
     
                
c other variables:
       PARAMETER(K=Gte*ZMsol/pc/1.d13,
     +           abuabu=abu*abu,
     +           bdbd=bd*bd,
     +           rKMb=K*Mb,
     +           q9=-3.d0*rKMb,
     +           rKMd=K*Md,
     +           q8=-3.d0*rKMd,
     +           rKMh=K*Mh)



