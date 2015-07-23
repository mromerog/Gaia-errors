c----------------------------------------------------------------------
c Example of main program that calls the routine Gaia-errors to some 
c input data (provided by the user).
c The input file must contain 9 columns with the following:
c 1-3 cols: x,y,z (in kpc)
c 4-6 cols: vx,vy,vz (in km/s)
c in galactocentric coordinates
c col7: Absolute magnitude in V (mag)
c col8: (V-I) intrinsic colour (mag)
c col9: Absorption in V (mag)
c The code performs the change of coordinates from galactocentric to
c heliocentric and computes the V apparent magnitude and observed (V-I)
c colour, this is the input necessary for the Gaia-errors subroutine.
c The output is transformed back to galactocentric coordinates and it writes
c input and output in the output file (provided by the user).
c
c To compile the code, just type in the command line:
c $ make main_Gaiaerrors
c It uses the gfortran compiler.
c
c It also needs, apart from the usual .h files, the cons_pop.h, that includes
c a specific Teff, logg and [Fe/H] for a population.
c 
c July 2015
c M. Romero-Gomez (ICCUB-IEEC)
c----------------------------------------------------------------------



        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION A(6),AE(6),AO(6),p(2),po(2),pe(2),ap(4),apo(4),ape(4)
        REAL*8 month
        CHARACTER*1 afact
        CHARACTER*64 ARXO,ARXI
        INCLUDE 'const_math.h' 
        INCLUDE 'const_ast.h'
        INCLUDE 'const_pop.h'
        INCLUDE 'const_pot.h'

c Initializations:
        idum=-15350

    
         WRITE(*,*)'GIVE INPUT FILE: '
         READ(*,*)ARXI
         WRITE(*,*)ARXI
         WRITE(*,*)'GIVE OUTPUT FILE: '
         READ(*,*)ARXO
         WRITE(*,*)ARXO
         OPEN(NCO,FILE=ARXO)
         WRITE(NCO,110)'#Av','V','VI','G','GRVS','sigPi/Pi','X','Y','Z',
     &                 'VX','VY','VZ','Xobs','Yobs','Zobs','VXobs',
     &                 'VYobs','VZobs'
 110     FORMAT(18(A8,1X))

         NCI1=27
         OPEN(NCI1,FILE=ARXI)
C READ PARTICLES GIVEN IN INERTIAL COORDINATES 
C POSITIONS IN kpc, VELOCITIES IN km/s IN GALACTOCENTRIC COORDINATES
C ABSOLUTE MAGNITUDE IN V (mag), 
C (V-I)o COLOR OF THE STAR (mag)
C Absorption in V (mag)
 10      READ(NCI1,*,END=11)XIN,YIN,ZIN,VXIN,VYIN,VZIN,aMv,vi_o,Av
C THE CODE SUPOSES SUN ON X-NEGATIVE AXIS.
C IF THE INPUT VARIABLES DO NOT HAVE THE RIGHT DISTRIBUTION, UNCOMMENT
C THE NEXT LINES AND PERFORM THE APPROPRIATE ROTATION:
C PUT THE LOCATION OF THE SUN-GC LINE WITH A CERTAIN ANGLE WITH RESPECT TO
C THE SEMI-MAJOR AXIS OF THE BAR (IF NECESSARY). 
C ANGLE OF THE BAR. THETAS=0 IS ON THE X-POSITIVE AXIS. HERE I ROTATE THE
C BAR SO THAT THE SUN IS ON THE X-NEGATIVE AXIS AND 20 DEG FROM THE SEMI-MAJOR
C AXIS OF THE BAR.
C        thetas=160.d0*pi/180.d0
c         call rot_horari(thetas,xin,yin,x,y)
c         call rot_horari(thetas,vxin,vyin,vx,vy)
         x=XIN
         y=yin
         z=zin
         vx=VXIN 
         vy=VYIN
         vz=vzin

         R=DSQRT(X*X+Y*Y)
c Transform the position of the star from galactocientric to heliocentric 
c coordinates (R0 DEFINED IN const_pot.h)
         dist=dsqrt((x+R0)*(x+R0)+y*y+z*z)   !in kpc
         distpc=dist*1000.d0                 !in pc
         xpi=1.d0/dist                       !parallax in mas
         CALL carte_to_equatorial(X,Y,Z,ALPHA,DELTA)
         a(1)=alpha
         a(2)=delta
         a(3)=xpi
         
C Transforming velocities from galactocentric to galactic helicocentric 
c coordinates and to equatorial heliocentric coordinates
         call Carte_to_UVWH(VX,VY,VZ,UH,VH,WH)
         
         CALL UVWH_to_equatorials(UH,VH,WH,A)
         

c  V apparent magnitude of the source
         V=aMv+5.d0*log10(distpc)-5.d0+Av

c Compute the observed (V-I) of the source, Cardelli et al (1989)
         vi=(1.d0-0.479d0)*Av+vi_o

c Assign atmospheric parameters to the source, defined in const_pop.h
         ap(1)=Teff    
         ap(2)=xlogg  
         ap(3)=FeH    !No radial metallicity gradient is assumed here
         ap(4)=Av     !we hare assuming A0=Av

c Introduce Gaia errors
         jflag=-1                          !weighted errors
c       jflag= 1                          !mean errors
         month=60.d0                          !5 years mission
         CAfactor=1.2d0
         call Gaia_errors(month,CAfactor,jflag,V,VI,a,ao,ae,p,po,pe,
     &                    ap,apo,ape,GRVS,idum)

C transform velocities from equatorial heliocentric to 
c galactocentric 
         CALL equatorials_to_UVWH(ao,UHobs,VHobs,WHobs)
         CALL UVWH_TO_Carte(UHobs,VHobs,WHobs,Vxobs,VYobs,VZobs)


c transform heliocentric to galactocentric positions 
         distobs=1.d0/ao(3)
         call equatorial_to_carte(distobs,ao(1),ao(2),xobs,yobs,zobs)

c WRITE IN OUTPUT FILE
         if (p(1).le.20.d0)then
            WRITE(NCO,100)Av,V,VI,p(1),GRVS,AE(3)/A(3),X,Y,Z,VX,VY,VZ,
     &              Xobs,Yobs,Zobs,VXobs,VYobs,VZobs
  100    FORMAT(18(F18.10,1X))
         ELSE
             BIG=1000000.00000000
             WRITE(NCO,100)Av,V,VI,p(1),GRVS,AE(3)/A(3),X,Y,Z,VX,VY,VZ,
     &              Xobs,Yobs,Zobs,BIG,BIG,BIG
         ENDIF
         GOTO 10
 11      CLOSE(NCI1)
         CLOSE(NCO)

         END




