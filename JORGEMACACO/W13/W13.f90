program colSphN

	implicit none
	
	double precision, dimension (:,:), allocatable :: Sa, S, Sp, V
	double precision, dimension (:,:), allocatable :: F
	double precision, dimension (:,:), allocatable :: Magnet
	double precision, dimension (:), allocatable :: D
	double precision, dimension (:), allocatable :: Center, R, MoI, M
	double precision, dimension (:), allocatable :: ThetaP, Theta, ThetaA, W, Tor, Itheta
	
	double precision :: h, K, t, CsiW, Room, Circumradius, P
	integer :: l, i, j, N , dimen, TTW, alpha, beta, cont
	
		N = 2			! Especificar número de esferas
		dimen = 2		! Dimensões calculadas/analisadas
		Room = 30.d0
		P = dacos(-1.d0)
		
	allocate(Sa(dimen,N), S(dimen,N), Sp(dimen,N), V(dimen,N), D(dimen))
	allocate(ThetaP(N),Theta(N),ThetaA(N), W(N), Itheta(2))
	
	allocate(F(dimen,N), Tor(N))
	
	allocate(Center(dimen))
	
	allocate(R(N), MoI(N), M(N), Magnet(dimen,N))
	
		W = 0.d0
		V = 0.d0
		
		h = 10.d0**(-4.d0)		! Passo de tempo	 
		TTW = 0 				! ttwrite que ainda não entendi direito
		
		M = 0.01d0				! Massa das esferas
		R = 1.d0				! Raio das esferas
		t = 0.d0				! Tempo
		Itheta(2) = 0.d0
		
	do alpha=1,N
		MoI(alpha) = 0.4d0*M(alpha)*R(alpha)*R(alpha)
	end do

	open(unit=129, file='div13.dat', status='unknown')
	open(unit=120, file='conv13.dat', status='unknown')

	do i = 1, 360
	
		S(1,2) = 10.d0
		S(2,2) = 15.d0
		V = 0.d0
		W = 0.d0
		
		Itheta(1) = 0.d0
		Itheta(2) = Itheta(2) + (P/180.d0)
		
		Magnet(1,2) = 10.d0*cos(Itheta(2))
		Magnet(2,2) = 10.d0*sin(Itheta(2))
		
		do j = 1, 360
		
			t = 0.d0
			l = 0
			Sa(1,1) = 20.d0
			Sa(2,1) = 15.d0
			V = 0.d0
			W = 0.d0
			
			W(1) = 13.d0*(P/10.d0)
	
			S(:,1) = Sa(:,1)
			ThetaA(1) = Itheta(1)
			Theta(1) = ThetaA(1) + h*W(1)	
			Magnet(1,1) = cos(Theta(1))
			Magnet(2,1) = sin(Theta(1))
			
			do while (l.eq.0) 
				
				! Wall contact check: -----------------------------------------------------------------
					if (S(1,1).lt.R(1)) then
						l = -1
					else if ((S(1,1)+R(1)).gt.Room) then
						l = -1
					end if
					if (S(2,1).lt.R(1)) then
						l = -1	
					else if ((S(2,1)+R(1)).gt.Room) then
						l = -1	
					end if
				! -------------------------------------------------------------------------------------
				
				! Magnetic interaction: ---------------------------------------------------------------
					D(:) = S(:,2) - S(:,1)
					call Mag (D,Magnet(:,1),Magnet(:,2),F(:,1),Tor(1))
					if (NORM2(D).lt.(R(1)+R(2))) then
						l = 1
					end if
				! -------------------------------------------------------------------------------------
			
				! Verlet integration: -----------------------------------------------------------------
					Sp(:,1) = 2.d0*S(:,1) - Sa(:,1) + (((h*h)/M(1))*F(:,1))
					ThetaP(1) = 2.d0*Theta(1) - ThetaA(1) + (((h*h)/MoI(1))*Tor(1))
					Magnet(1,1) = cos(Theta(1))
					Magnet(2,1) = sin(Theta(1))
					
					V(:,1) = (Sp(:,1) - Sa(:,1))/(2.d0*h)
					W(1) = (ThetaP(1) - ThetaA(1))/(2.d0*h)
					
					Sa(:,1) = S(:,1)
					S(:,1) = Sp(:,1)
					ThetaA(1) = Theta(1)
					Theta(1) = ThetaP(1)
					
					Tor = 0.d0
					F = 0.d0
					t = t+h
				! -------------------------------------------------------------------------------------
			end do
			
			if (l.eq.(-1)) then
				write(129,*) 13.d0*(P/10.d0), t, i, j
			else 
				write(120,*) 13.d0*(P/10.d0), t, i, j
			end if
			
			Itheta(1) = Itheta(1) + (P/180.d0)
			
		end do
	end do
	
	close (unit=129)
	close (unit=120)
	
	deallocate(Sa, S, Sp, V)
	deallocate(ThetaA, Theta, ThetaP, W)
	deallocate(F, Tor)
	deallocate(R, M, MoI, Magnet)
	deallocate(Center)

contains

	subroutine Mag ( D, MagnetA, MagnetB, FA, TorA)
		
		double precision, dimension(2), intent (in) :: D
		double precision, dimension(2), intent (in) :: MagnetA, MagnetB
		
		double precision, intent (inout) :: TorA
		double precision, dimension(2), intent (inout) :: FA
		
		double precision, dimension(2) :: NORMAL, Fmag, B
		double precision :: DotA, DotB, DotC, Dif
		
		Dif = Norm2(D)
		NORMAL = -D/Dif
		
		DotA = Dot_product(MagnetA(:), Normal(:))
		DotB = Dot_product(MagnetB(:), Normal(:))
		DotC = Dot_Product(MagnetA(:), MagnetB(:))
		
		Fmag(:) = 3.d0*(DotA*MagnetB(:) + DotB*MagnetA(:) + (DotC - 5.d0*DotA*DotB)*Normal(:))/(dif*dif*dif*dif)
		B(:) = (3.d0*DotB*Normal(:) - MagnetB(:))/(Dif*Dif*Dif)
		
		FA(:) = FA(:) + Fmag(:)
		TorA = TorA + (MagnetA(1)*B(2) - MagnetA(2)*B(1))
		
	end subroutine

	subroutine CIPoligon (Center, S, Vertices, Circumradius)
		
		integer, intent (in) :: Vertices
		double precision, intent (in) :: Circumradius
		double precision, dimension(2), intent (in) :: Center
		double precision, dimension(2,Vertices), intent (inout) :: S
		
		integer :: k
		double precision :: pi, Theta
		
		pi = 3.141596d0
		Theta = 2.d0*pi/Vertices
		
		do k = 0, (Vertices - 1)
			S(1,k+1) = Center(1) + Circumradius*cos(k*Theta)
			S(2,k+1) = Center(2) + Circumradius*sin(k*Theta)
		end do
		
	end subroutine

	subroutine wall (alpha, cord, ori, CsiW, V, F, R, W, Tor)
	
		implicit none
		integer, intent(in) :: alpha, cord
		double precision :: Fn, Fs, K, mi, GamaN, GamaS, Vs
		double precision, intent(in) ::  CsiW, ori
		double precision, dimension(N), intent(inout) :: Tor, W, R
		double precision, dimension(2,N), intent(in) :: V
		double precision, dimension(2,N), intent(inout) :: F

			K = 10.d0**(5.d0)		
			GamaN = 50.d0			
			GamaS = 20.d0
			Mi = 0.5d0				

			Fn = - K*CsiW - GamaN*V(cord,alpha)
			Vs = V(3-cord,alpha) + ori*R(alpha)*W(alpha)
			if((abs(GamaS*Vs)).lt.(abs(Mi*Fn))) then
				Fs = - GamaS*Vs
			else if (abs(Vs).gt.(10**(-6))) then
				Fs = - Mi*Vs*abs(Fn/Vs)
			else 
				Fs = 0.d0
			end if
			
			F(cord,alpha) = F(cord,alpha) + Fn
			F(3-cord,alpha) = F(3-cord,alpha) + Fs
			Tor(alpha) = Tor(alpha) + ori*R(alpha)*Fs
			
	end subroutine
	
	subroutine ContForThirdL (D, alpha, beta, V, W, R, F, Tor)
	
		implicit none
		integer, intent (in) :: alpha, beta
		double precision :: K, GamaN, GamaS, Mi
		double precision, dimension(2) :: NORMAL, TANGENCIAL, Vr, Vn, Vs, Fn, Fs, CSI
		double precision, dimension(2), intent (in) :: D
		double precision, dimension(2,N), intent (in) :: V
		double precision, dimension(2,N), intent (inout) :: F
		double precision, dimension(N), intent (in) :: R, W
		double precision, dimension(N), intent (inout) :: Tor
		
		K = 10.d0**(6.d0)		! Constante elástica
		GamaN = 50.d0			! Coeficiente de dissipação (Normal)
		GamaS = 45.d0			! Coeficiente de dissipação (Shear)
		Mi = 0.5d0				! Coeficiente de atrito dinâmico
		
		NORMAL = D/(NORM2(D))	
		TANGENCIAL(1) = -NORMAL(2)	
		TANGENCIAL(2) = NORMAL(1)
		
		CSI = ( R(alpha) + R(beta) - NORM2(D) )*NORMAL
		
		Vr(:) = V(:,alpha) - V(:,beta)
		Vn = (Dot_Product(Vr,NORMAL))*NORMAL
		Vs = ((Dot_Product(Vr,TANGENCIAL)) + W(alpha)*R(alpha) + W(beta)*R(beta))*TANGENCIAL
		
		Fn = - K*CSI - GamaN*Vn
		if	((GamaS*NORM2(Vs)).lt.(Mi*NORM2(Fn))) then
			Fs = - GamaS*Vs
		else
			Fs = - (Mi*NORM2(Fn)/NORM2(Vs))*Vs
		end if
		
		F(:,alpha) = F(:,alpha) + Fn + Fs
		F(:,beta)  = F(:,beta) - Fn - Fs
		Tor(alpha) = Tor(alpha) + R(alpha)*Dot_Product(Fs,TANGENCIAL)
		Tor(beta) = Tor(beta) + R(beta)*Dot_Product(Fs,TANGENCIAL)
	
	end subroutine

	subroutine salva_eps(ttwrite,parts,part_raio,part_pos_x,part_pos_y,part_ang)
	
	  implicit none

	  integer, intent (in) :: ttwrite, parts
	  
	  double precision, dimension (1:N), intent (in) :: part_raio
	  double precision, dimension (1:N), intent (in) :: part_pos_x
	  double precision, dimension (1:N), intent (in) :: part_pos_y
	  double precision, dimension (1:N), intent (in) :: part_ang

	  ! integer, intent (in) :: arq !Se diferente de 0, arquivo de CI

	  integer :: i
	  
	  character (len = 15) :: filename1
	 
	  !Criação do arquivo de saída - número entre 000 e 999!!!
	 
		filename1="outputs.000.eps"
		
		filename1(11:11)=CHAR(48+mod(ttwrite,100)-10*(mod(ttwrite,100)/10))
		filename1(10:10)=CHAR(48+(mod(ttwrite,100)/10))
		filename1(9:9)=CHAR(48+ttwrite/100)
	  

	  open(unit=210,file=filename1,status='unknown')

	  !Escreve o cabeçalho inicial do eps
	  write(210,90)
	  write(210,91)
	  write(210,92) filename1
	  write(210,93)
	  write(210,94)
	  write(210,95)
	  write(210,96) 0, 0, 200, 200
	  write(210,97)
	  write(210,98)
	  write(210,99)
	  write(210,100)
	  write(210,101)
	  write(210,102)
	  write(210,103)
	  !write particles without orientation
	  do i=1,parts
		 write(210,104) int(10*part_pos_x(i)), int(10*part_pos_y(i)), int(10*part_raio(i))
		 write(210,*) int(10*part_pos_x(i)), int(10*part_pos_y(i)), ' moveto' 
		 write(210,*) int(10*part_raio(i)*cos(part_ang(i))), int(10*part_raio(i)*sin(part_ang(i))),' rlineto' 
		 write(210,*) 'stroke'
	  end do
	  write(210,105)
	  write(210,106)
	  write(210,107)
	  close(unit=210)
	  
	90 format('%%!PS-Adobe-3.0 EPSF-3.0')
	91 format('%%Document-Fonts: Times-Roman')
	92 format('%%Title: ', A15)
	93 format('%%Creator: YDS')
	94 format('%%CreationDate: unknown')
	95 format('%%Pages: 1')
	96 format('%%BoundingBox: ',4I5)
	97 format('%%LanguageLevel: 1')
	98 format('%%EndComments')
	99 format('%%BeginProlog')
	100 format('%%EndProlog')
	101 format('0.0000 0.0000 0.0000 setrgbcolor')
	102 format('%% Page:     1    1')
	103 format('save')
	104 format(3I8,'  0   360  arc')
	105 format('restore showpage')
	106 format('%%Trailer')
	107 format('%%EOF')

	end subroutine salva_eps

end program colSphN