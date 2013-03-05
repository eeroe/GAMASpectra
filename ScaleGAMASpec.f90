program SimpleMag

!Scales GAMA spectra to match (App.Matched Petrosian) catalogue magnitudes (with a linear scaling in magnitudes). Writes output to ascii files and to FITS (but they don't seem to be working on all platforms - fine for GAIA?). 
IMPLICIT NONE
integer, parameter :: magtype = 2 ! Set to 1 for InputCat-Petro mag scalings, 2 for appmatchedphotom-Petro mag scalings, 3 for appmatchedphotom-kron mag scalings
character(len=200) :: dataset='GAMA_4000randomGroupedGals_GroupInfo'
character(len=200) :: workdir='/disk1/ee/roag/GAMA/NewSpec/Selec2_4000'
character(len=200) :: shortworkdir='Selec2_4000'
character(len=200) :: datdir='/disk1/ee/roag/GAMA/NewSpec/GAMA_GroupedGal_Fits'
character(len=200) :: infilteg = '/disk1/ee/roag/GAMA/SDSSgfilt.dat' !! Location of filter transmissions
character(len=200) :: infilter = '/disk1/ee/roag/GAMA/SDSSrfilt2.dat'
character(len=200) :: infiltei = '/disk1/ee/roag/GAMA/SDSSifilt.dat'
character(len=200) :: outsummaryfile
real, parameter :: l1 = 4686.d0, l2 = 6165.d0, l3=7481.d0 !Effective wavelengths of filters
!real, parameter :: l1 = 4430.d0, l2 = 6005.d0, l3=7230.d0

character(len=200) :: indat,magresults, subdir, inloc
character(len=8) :: key1, key2
character(len=48) :: comm1, comm2
character(len=17) ::val1, val2
real*4 :: gdif = 0.1, rdif=0.1, idif=0.1, magdifg, magdifr, magdifi
real*4 :: inimagg, inimagr, inimagi, inimagdifg, inimagdifr, inimagdifi
integer :: success=0, col, it, ncol, GroupID, GroupID12, GroupID2, Nfof
real*4, dimension(:), allocatable :: sp, fi, wav, magdiff
real*4 :: maggc, magrc, magic, medval
 integer :: ios=0, co=0, colin=1, i=0, k=1, coun=0, segcount=0, counf=0, j=0, n1, n2
real*4 :: cang=2.99792458e18 ! speed of light in Angstroms/s
real :: f1, f2, f3, nf1, nf2, nf3, a, b, c
integer :: CATAID, InterCenCATAID
real*8 :: IterCenRA, IterCenDEC,MassA, MassAfunc, LumB, LumBfunc
!! Median Filter parameters
integer, parameter :: MedFil = 1, MedFilMag=0 ! To calulate magnitudes after median filtering
integer, parameter :: NumBin=301, MedBin=151 !Width of median filter = NumBin*disp
real*4, dimension (NumBin) :: medarray
integer, dimension (NumBin) :: orderarray
real*4, dimension (:,:), allocatable :: interpfiltg, interpfiltr, interpfilti
real*4, dimension (2,5000) :: filting, filtinr, filtini
real*4 :: filtming, filtminr, filtmini, filtmaxg, filtmaxr, filtmaxi
real*4, dimension (:), allocatable ::  filtg, filtr, filti
integer :: coung=0, counr=0, couni=0, speccou=0, totsucc=0, totitt=0

character(len=700) :: outfilespec, outfilemag, infilename, outfitsspec, band, FILENAME, outheaderfile, magty
character :: cCRVAL*80='CRVAL1', cCRPIX*80='CRPIX1', cCD*80='CD1_1', cWMIN*80='WMIN', cWMAX*80='WMAX'
character :: cSPECID*80='SPECID', SPECID*40, comm*48, cN1*80='NAXIS1', cN2*80='NAXIS2',newkeyword*8

integer :: speccount=0,refpoint, nrow, ibitpix=-32, iosg=0, iosr=0, iosi=0, errorcount=0, unsuc=0, IterCenCATAID
real*4 :: magg, magr, magi,gcen, rcen, icen, refval, disp, val,WMIN, WMAX
real :: RA,DEC,PETROMAG_G,PETROMAG_R,PETROMAG_I !! Petrosian Magnitudes from InputCatAv06 in InputCat
real :: MAG_PETRO_U,MAG_PETRO_G,MAG_PETRO_R,MAG_PETRO_I,MAG_PETRO_Z,MAG_PETRO_Y,MAG_PETRO_J,MAG_PETRO_H,MAG_PETRO_K !! Petrosian Magnitudes from ApMatchedCatv03 in ApMatchedPhotom)
real ::  MAG_AUTO_U,MAG_AUTO_G,MAG_AUTO_R,MAG_AUTO_I,MAG_AUTO_Z,MAG_AUTO_Y,MAG_AUTO_J,MAG_AUTO_H,MAG_AUTO_K !! Kron Magnitudes from ApMatchedCatv03 in ApMatchedPhotom)
real ::A_u,A_g,A_r,A_i,A_z,A_y,A_j,A_h,A_k !! Galactic Extinction Corrections (from GalacticExtinctionv02 of InputCat)
real :: G_mag_scale, R_mag_scale, I_mag_scale, SN, IterCenZ, Zfof
real*4, allocatable :: arrayspec(:,:), fitsspec(:,:)

      integer status,unit,readwrite,blocksize, morekeys
      integer group,firstpix,nbuffer,npixels
      real datamin,datamax,nullval,buffer(100)
      logical anynull

INTERFACE 
  FUNCTION freadreal(infilename, keyword)
       character*80, INTENT(IN) :: infilename, keyword
  END FUNCTION freadreal
END INTERFACE

INTERFACE 
  FUNCTION freaddble(infilename, keyword)
       character*80, INTENT(IN) :: infilename, keyword
  END FUNCTION freaddble
END INTERFACE
INTERFACE 

  FUNCTION freadint(infilename, keyword) !! doesn't work yet?
       character*80, INTENT(IN) :: infilename, keyword
  END FUNCTION freadint
END INTERFACE

INTERFACE 
  FUNCTION getaxes(infilename)
       character*80, INTENT(IN) :: infilename
  END FUNCTION getaxes
END INTERFACE


write(6,*);write(6,*);write(6,*)
open(25,file=infilteg)
open(26,file=infilter)
open(27,file=infiltei)

if (magtype.eq.1) then 
magty='IC_Petro'
subdir = 'InputCat_Petro'
inloc = 'InputCatAv06'
write(6,*) 'Scaling flux to match InputCat_Petro mags'
else if (magtype.eq.2) then 
magty=''!'AM_Petro'
subdir = 'AppMatch_Petro'
inloc='ApMatchedCatv03'
write(6,*) 'Scaling flux to match AppMatch_Petro mags'
else if (magtype.eq.3) then 
magty='AM_Kron'
subdir = 'AppMatch_Kron'
inloc='ApMatchedCatv03'
write(6,*) 'Scaling flux to match AppMatch_Kron mags'
end if



 !! Read in input SDSS filter transmission data

	coung=0; counr=0; couni=0; iosg=0; iosr=0; iosi=0
	do WHILE (iosg.ge.0)
		coung=coung+1
		read(25, *,iostat=iosg) filting(1,coung), filting(2,coung)
	end do

	do WHILE (iosr.ge.0)
		counr=counr+1
		read(26, *,iostat=iosr) filtinr(1,counr), filtinr(2,counr)
	end do

	do WHILE (iosi.ge.0)
		couni=couni+1
		read(27, *,iostat=iosi) filtini(1,couni), filtini(2,couni)
	end do

filtming = (filting(1, 1)); filtmaxg=(filting(1, coung-1))
filtminr = (filtinr(1, 1)); filtmaxr=(filtinr(1, counr-1))
filtmini = (filtini(1, 1)); filtmaxi=(filtini(1, couni-1))

outheaderfile = ''// trim(adjustl(workdir)) //'/header.dat'
indat = ''// trim(adjustl(workdir)) //'/'//trim(adjustl(dataset))//'.dat'
magresults = ''// trim(adjustl(workdir)) //'/mag.dat'
open(76, file=magresults)
write(76,*) 'SPECID		','inimagG		','magGcalc	','Catalogue_MAG_G	','inimagR	','magRcalc	','Catalogue_MAG_R	'&
&,'inimagI	','magIcalc	','Catalogue_MAG_I	','success	','A_g	','A_r	','A_i'
write(*,*);write(*,*) 'Number of G band input filter transmission values = ', coung
gcen=(filtmaxg+filtming)/2.0;rcen=(filtmaxr+filtminr)/2.0;icen=(filtmaxi+filtmini)/2.0
write(6,*) 'gcen = ', gcen, ' rcen = ', rcen, ' icen = ', icen
write(6,*) 'Effective wavelengths of filters used (SDSS DR7) : ', l1, l2, l3
close(25);close(26);close(27)

open(56, file=outheaderfile)
write(56,*) 'SPECID		','IN_MAG_G		','IN_MAG_R	','IN_MAG_I	','CatMag_G	','CatMag_R	','CatMag_I	'&
&,'NEWMAG_G	','NEWMAG_R	','NEWMAG_I	', '     S/N'
open(66, file=indat)
read(66, *) 
ios = 0
do WHILE (ios.ge.0)
700	read(66, *, iostat=ios) SPECID,CATAID,RA,DEC,MAG_PETRO_U,MAG_PETRO_G,MAG_PETRO_R,&
	&MAG_PETRO_I,MAG_PETRO_Z,MAG_PETRO_Y,MAG_PETRO_J,MAG_PETRO_H,MAG_PETRO_K,&
	&A_u,A_g,A_r,A_i,A_z,A_y,A_j,A_h,A_k, SN,GroupID, Nfof, IterCenCATAID, IterCenRA, IterCenDEC, IterCenZ, &
	&Zfof, MassA, MassAfunc, LumB, LumBfunc
!write(6,*) magty

!write(6,*)  MAG_PETRO_G,MAG_PETRO_R,MAG_PETRO_I

if (magtype.eq.1) then
G_mag_scale=PETROMAG_G; R_mag_scale=PETROMAG_R; I_mag_scale=PETROMAG_I

else if (magtype.eq.2) then 
G_mag_scale=MAG_PETRO_G; R_mag_scale=MAG_PETRO_R; I_mag_scale=MAG_PETRO_I
else if (magtype.eq.3) then 
G_mag_scale=MAG_AUTO_G; R_mag_scale=MAG_AUTO_R; I_mag_scale=MAG_AUTO_I
end if
write(6,*)  G_mag_scale,R_mag_scale,I_mag_scale

infilename=''// trim(adjustl(datdir)) //'/'// trim(adjustl(SPECID)) // ".fit"

status=0
  !call ftgiou(unit,status)
unit=4
 readwrite=0

	call ftopen(unit,infilename,readwrite,blocksize,status)

   call ftclos(unit, status)
   call ftfiou(unit, status)
!     check for any error, and if so print out error messages
    if (status .gt. 0) call printerror(status)
    if (status .gt. 0) 	errorcount=errorcount+1
    if (status .gt. 0) go to 700

speccou = speccou+1
write(6,*) speccou, '  :  ', infilename

!!! Read in data from Fits Header


n1=freadreal(infilename, cN1)
n2=freadreal(infilename, cN2)

refval=freadreal(infilename, cCRVAL)
refpoint=freadreal(infilename, cCRPIX)
disp=freadreal(infilename, cCD)
WMIN=freadreal(infilename, cWMIN)
WMAX=freadreal(infilename, cWMAX)

nrow=n1; ncol=n2

write(6,*) 'Minimum spectrum wavelength = ' ,(refval+(dble(1.0-refpoint)*disp)), WMIN
write(6,*) 'Maximum spectrum wavelength = ', (refval+(dble(nrow-refpoint)*disp)), WMAX

allocate (sp(nrow));allocate (fi(nrow));allocate (wav(nrow)); allocate(magdiff(nrow))
allocate (filtg(nrow));allocate (filtr(nrow));allocate (filti(nrow))
allocate (interpfiltg(2,nrow));allocate (interpfiltr(2,nrow));allocate (interpfilti(2,nrow))


allocate(arrayspec(n1,n2)); allocate(fitsspec( n1, n2+2))


!! Read in spectra from fits file

status=0
  !call ftgiou(unit,status)
unit=3
 readwrite=0

	call ftopen(unit,infilename,readwrite,blocksize,status)
!     initialize variables
      npixels=nrow*ncol
      group=1
      firstpix=1
      nullval=-999
      datamin=1.0E30
      datamax=-1.0E30

  call ftgpve(unit,group,firstpix,npixels,nullval,arrayspec,anynull,status)

   call ftclos(unit, status)
   call ftfiou(unit, status)
!     check for any error, and if so print out error messages
    if (status .gt. 0)call printerror(status)

!write(6,*) 'First 10 values of flux'
!do i = 1, 10
!write(6,*) arrayspec(i,1)
!end do

do i = 1, nrow
sp(i) = arrayspec(i, 1)
end do

write(6,*) 'Checking for bad values (=-999.000 ??)'
do i  = 1, nrow   !! Remove values which are rediculously low (why are they there?)
		
	if (arrayspec(i, 1).ne.(arrayspec(i,1))) PAUSE 'Error : NaN?'
	if (sp(i).lt.(-100.0)) then
!	write(6,*) 'sp < -100, ', sp(i), i
		if (i.ne.1) then
!			write(6,*) i, (refval+((i-refpoint)*disp)), arrayspec(i,1)
			sp(i)=sp(i-1)
			!	write(6,*) 'sp fixed ?? ', sp(i), i

		else if (i.eq.1) then
!			write(6,*) ' first element error 				!!!!'
			do j = 1, 40
				if ((sp(i+j)).gt.(-1e2)) then
					sp(i)=sp(i+j) 
					write(6,*) sp(i)
					exit
				end if
			end do
		!	write(6,*) i, (refval+((i-refpoint)*disp)), arrayspec(i,1), sp(i)
			if (sp(i).lt.-1e2) then
			unsuc=unsuc+1
			deallocate (sp);deallocate (fi);deallocate (wav); deallocate(magdiff)
			deallocate (filtg);deallocate (filtr);deallocate (filti)
			deallocate (interpfiltg);deallocate (interpfiltr);deallocate (interpfilti)
			deallocate(fitsspec); deallocate(arrayspec)
			go to 700!PAUSE 'ERROR: Flux value not set'
			end if
		end if
	!else
	!sp(i) = arrayspec(i,1)
	end if
end do

!do i = 1, 10
!write(6,*) i, (refval+((i-refpoint)*disp)), sp(i)
!end do

do j = 1, n2
	do i  = 1, nrow
		fitsspec(i,j)=arrayspec(i, j)
	end do
end do


		do i  = 1, nrow
		fitsspec(i,n2+1)=(refval+((i-refpoint)*disp)) !! Wavelengths
		wav(i)=(refval+((i-refpoint)*disp))
		!sp (i) = fitsspec(i,1)
		end do



!!!!!!!!!!!!!! INTERPOLATE FILTERS !!!!!!!!!!!!!!!!!!!!


	filtg(1:nrow)=filt (wav, nrow, filtming, filtmaxg, coung, filting)
do i = 1, nrow
interpfiltg(1, i)=wav(i); interpfiltg(2,i)=filtg(i)
end do
	filtr(1:nrow)=filt (wav, nrow, filtminr, filtmaxr, counr, filtinr)
do i = 1, nrow
interpfiltr(1, i)=wav(i); interpfiltr(2,i)=filtr(i)
end do

	filti(1:nrow)=filt (wav, nrow, filtmini, filtmaxi, couni, filtini)
do i = 1, nrow
interpfilti(1, i)=wav(i); interpfilti(2,i)=filti(i)
end do

!!!!!!!!!!!!!!!!!!! CALCULATE MAGNITUDES !!!!!!!!!!!!!!!!!!!!


maggc= getmag(sp, interpfiltg(2, 1:nrow), wav, disp, nrow)
write(6,*) 'Calculated G magnitude = ', maggc, &
&'Mag check = ', ABmagcheck(sp, interpfiltg(2, 1:nrow), wav, disp, nrow),  ' Cat. Mag - A_g = ', G_mag_scale - A_g
	
magrc= getmag(sp, interpfiltr(2, 1:nrow), wav, disp, nrow)
write(6,*) 'Calculated R magnitude = ', magrc, &
&'Mag check = ', ABmagcheck(sp, interpfiltr(2, 1:nrow), wav, disp, nrow),  ' Cat. Mag - A_r = ', R_mag_scale - A_r

magic= getmag(sp, interpfilti(2, 1:nrow), wav, disp, nrow)
write(6,*) 'Calculated I magnitude = ', magic, &
&'Mag check = ', ABmagcheck(sp, interpfilti(2, 1:nrow), wav, disp, nrow),  ' Cat. Mag - A_i = ', I_mag_scale - A_i


inimagg = maggc; inimagr = magrc; inimagi = magic
inimagdifg = maggc - G_mag_scale+A_g !! Differences between calculated magnitudes and extinction corrected catalogue mags
inimagdifr = magrc - R_mag_scale+A_r
inimagdifi = magic - I_mag_scale+A_i

do i =1, nrow !! Linearly interpolate magnitude differences between effective wavelengths of filters

	if (wav(i).le.l2) then
		magdiff(i)=inimagdifg+(wav(i)-l1)*((inimagdifr-inimagdifg)/(l2-l1))
	else 
		magdiff(i)=inimagdifi+(wav(i)-l3)*((inimagdifi-inimagdifr)/(l3-l2))
	end if


	fitsspec(i,n2+2) = fitsspec(i, 1)*(10.0**(0.4*magdiff(i))) !! Correct flux to match catalogue mags
	success=0
	
	sp(i) = sp(i)*(10.0**(0.4*magdiff(i)))

end do

!! Get new magnitudes from now scaled spectra
maggc= getmag(sp, interpfiltg(2, 1:nrow), wav, disp, nrow)
magrc= getmag(sp, interpfiltr(2, 1:nrow), wav, disp, nrow)
magic= getmag(sp, interpfilti(2, 1:nrow), wav, disp, nrow)


magdifg = ABS(maggc - G_mag_scale+A_g)
magdifr = ABS(magrc - R_mag_scale+A_r)
magdifi = ABS(magic - I_mag_scale+A_i)
if (magdifg.lt.gdif.and.magdifr.lt.rdif.and.magdifi.lt.idif) success = 1

if (success.eq.1) totsucc = totsucc+1

write(6,*) 'Calculated G magnitude after scaling = ', maggc, ' (extinction corrected) CAT mag = ', G_mag_scale-A_g
write(6,*) 'Calculated R magnitude after scaling = ', magrc, ' (extinction corrected) CAT mag = ', R_mag_scale-A_r
write(6,*) 'Calculated I magnitude after scaling = ', magic, ' (extinction corrected) CAT mag = ', I_mag_scale-A_i

if (success.eq.1) then
!outfitsspec=''// trim(adjustl(workdir)) //'/'
outfitsspec=''//trim(adjustl(shortworkdir))//'/Fits/'&
&// trim(adjustl(SPECID)) //'_Scaled.fits'
outfilespec=''// trim(adjustl(shortworkdir)) //'/Data/'&
&// trim(adjustl(SPECID)) // "_Scaled.dat"

else if (success.eq.0) then
unsuc=unsuc+1
outfitsspec=''//trim(adjustl(shortworkdir))//'/Fits/'&
&// trim(adjustl(SPECID)) //'_FAILED.fits'
outfilespec=''// trim(adjustl(shortworkdir)) //'/Data/'&
&// trim(adjustl(SPECID)) //"_FAILED.dat"

end if
write(6,*) 'success? = ', success, ' scaled spectra written to ', outfitsspec



call delete_file(outfitsspec, status) !Delete fits file if it already exists
call writeimage(fitsspec,n1,n2+2,outfitsspec,ibitpix)
call fcopyheader(infilename, outfitsspec)

!Write new results to fits header
key1='ROW6    '; key2='ROW7    '
val1='Wavelength'; val2='Scaled Spectrum'
comm1='Calculated Wavelengths'
comm2='Linearly scaled to match catalogue magnitudes'
call fwritestring(outfitsspec, key1,val1, comm1)
call fwritestring(outfitsspec, key2, val2,comm2)
comm1='G band magnitude calculated from input spectrum'
call fwritereal(outfitsspec, 'IN_MAG_G', inimagg, comm1)
comm1='R band magnitude calculated from input spectrum'
call fwritereal(outfitsspec, 'IN_MAG_R', inimagr, comm1)
comm1='I band magnitude calculated from input spectrum'
call fwritereal(outfitsspec, 'IN_MAG_I', inimagi, comm1)
comm1='Extinction corrected G band magnitude from '//trim(adjustl(inloc))//''
call fwritereal(outfitsspec, 'Cat_Mag_G', G_mag_scale-A_g, comm1)
comm1='Extinction corrected R band magnitude from '//trim(adjustl(inloc))//''
call fwritereal(outfitsspec, 'Cat_Mag_R', R_mag_scale-A_r, comm1)
comm1='Extinction corrected I band magnitude from '//trim(adjustl(inloc))//''
call fwritereal(outfitsspec, 'Cat_Mag_I', I_mag_scale-A_i, comm1)
comm1='New G band magnitude from scaled spectra'
call fwritereal(outfitsspec, 'NEWMAG_G', maggc, comm1)
comm1='New R band magnitude from scaled spectra'
call fwritereal(outfitsspec, 'NEWMAG_R', magrc, comm1)
comm1='New I band magnitude from scaled spectra'
call fwritereal(outfitsspec, 'NEWMAG_I', magic, comm1)

write(56, *) trim(adjustl(SPECID)), inimagg, inimagr, inimagi, G_mag_scale,R_mag_scale,I_mag_scale,maggc,magrc,magic, SN

write(76,*) trim(adjustl(SPECID)), inimagg, maggc, G_mag_scale, inimagr, magrc, R_mag_scale,&
&inimagi, magic, I_mag_scale, success, A_g,A_r,A_i
open(34, file=outfilespec)
do i = 1, nrow
write(34,*) fitsspec(i, 1), fitsspec(i, 2),fitsspec(i, 3), fitsspec(i, 4), fitsspec(i, 5),&
& fitsspec(i, 6),  fitsspec(i, 7), (10.d0**(0.4d0*magdiff(i))) , magdiff(i)
sp(i)=0.d0; fi(i)=0.d0; wav(i)=0.d0; filtg(i)=0.d0; filtr(i)=0.d0; filti(i)=0.d0
end do
close(34)
deallocate (sp);deallocate (fi);deallocate (wav); deallocate(magdiff)
deallocate (filtg);deallocate (filtr);deallocate (filti)
deallocate (interpfiltg);deallocate (interpfiltr);deallocate (interpfilti)
deallocate(fitsspec); deallocate(arrayspec)

end do
close(56)
close(76)
write(6,*) 'magnitudes written to file : ', magresults
write(6,*) 'spectra written to file : ', outfilespec
write(6,*) totsucc, ' of ', speccou,' spectra successfully scaled (',dble(totsucc)*100.0/dble(speccou),')%'
write(6,*) unsuc, ' of ', speccou,' spectra unsuccessfully scaled (',dble(unsuc)*100.0/dble(speccou),')%'
write(6,*) 'Error opening ', errorcount, ' fits files'


outsummaryfile=''// trim(adjustl(workdir)) //'/ResultsSummary.dat'
open(17, file=outsummaryfile)

write(17,*) 'header data written to file : ', outheaderfile
write(17,*) totsucc, ' of ', speccou,' spectra successfully scaled (',dble(totsucc)*100.0/dble(speccou),')%'
write(17,*) unsuc, ' of ', speccou,' spectra unsuccessfully scaled (',dble(unsuc)*100.0/dble(speccou),')%'
write(17,*) 'Error opening ', errorcount, ' fits files'

close(17)

CONTAINS

 FUNCTION filt (wav, N, filtmin, filtmax, coun, filtin)
      IMPLICIT NONE
      INTEGER j, N, m, coun
	real*4 :: wspec, filtmin, filtmax, lowoff, highoff
	real*4, dimension (N) :: wav, filt, intfilt
	real*4, dimension(2,coun) :: filtin


do j = 1, N 
wspec = wav(j)
		if (wspec.ge.filtmin.and.wspec.lt.filtmax) then
		!if (wspec.ge.100.and.wspec.lt.10000) then

		do m=1,coun-1
			if (wspec.ge.filtin(1,m).and.wspec.lt.filtin(1,m+1)) then

				lowoff=(wspec-filtin(1,m))
				highoff=(filtin(1,m+1)-wspec)
				intfilt(j)=((lowoff*filtin(2,m+1))/(lowoff+highoff))+&
				&((highoff*filtin(2,m))/(lowoff+highoff))
		
			end if
		end do
		else 
				intfilt(j)=0.d0

		end if

	end do
				do j = 1, N
				filt(j)=intfilt(j)
				end do

      RETURN 
END FUNCTION filt



real*4 FUNCTION getmag (spe, fi, wav, disp, N)
      IMPLICIT NONE

      INTEGER N
      real*4 :: spe(1:N), fi(1:N), wav (1:N)
      real*4 ::getmag, freq, dfreq ,disp, Ospec=0.d0, stspec=1.0e-6, Ostandard=0.d0
	real*4 :: cang=2.99792458e18 
	integer :: i
      
      Ospec = 0.d0		!! Makes a difference if you set them to 0 here instead of
      Ostandard = 0.d0  !! just in the definitions above
      
	DO i=1,N
	freq=cang/wav(i)
	dfreq=freq*freq*disp/cang
	Ospec=Ospec+(spe(i)*disp*fi(i)/freq)
	Ostandard=Ostandard+(stspec*dfreq*fi(i)/freq)

	end do
!write(6,*) (Ospec/Ostandard), log10(Ospec/Ostandard)
	getmag=8.9d0-2.5d0*log10(Ospec/Ostandard)

      RETURN 
END FUNCTION getmag

real*4 FUNCTION ABmagcheck (spe, fi, wav, disp, N)
      IMPLICIT NONE

      INTEGER N
      real*4 :: spe(1:N), fi(1:N), wav (1:N)
      real*4 ::ABmagcheck, freq, dfreq ,disp, Ospec=0.d0, stspec=3631.0e-6, Ostandard=0.d0
	real*4 :: cang=2.99792458e18 , hp=6.626068e-27
	integer :: i
      
      Ospec = 0.d0		!! Makes a difference if you set them to 0 here instead of
      Ostandard = 0.d0  	!! just in the definitions above
      
	DO i=1,N
	freq=cang/wav(i)
	dfreq=freq*freq*disp/cang
	!Ospec=Ospec+(spe(i)*disp*fi(i)*(1e-17)/freq)
	Ostandard=Ostandard+(stspec*dfreq*fi(i)/freq)
	Ospec=Ospec+(spe(i)*disp*fi(i)/freq)
	end do

	ABmagcheck=-2.5d0*log10(Ospec/Ostandard)

      RETURN 
END FUNCTION ABmagcheck

end program SimpleMag


include '/disk1/ee/roag/Fortran/FitsSubroutines.f90'

include '/disk1/ee/roag/Fortran/cookbook.f90'
!include '/disk1/ee/roag/Fortran/GetFitsSpecInfo.f90'
!include '/disk1/ee/roag/GAMA/NewSpec/TestSpecCode/fitsio.f90'
