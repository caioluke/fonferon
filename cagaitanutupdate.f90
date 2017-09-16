program colisao3d
implicit none
double precision, dimension (:,:), allocatable :: Sold,Snow,S,v     
double precision, dimension (2) :: D
double precision, dimension (:,:,:), allocatable :: Fnormal
double precision, dimension (:,:), allocatable :: FR,Fparede,Fatrito
double precision, dimension (:), allocatable :: R,Inercia,tetaold,tetanow,teta,omega,Torque,ang,m
integer, dimension (:,:,:), allocatable :: Cell
integer, dimension (:,:), allocatable :: marcabola
double precision :: h,K,gamaN,gamaS,Lx,Ly,g,minormal,miparede,Lcell,rmedio
double precision :: Ldown, amplitude
integer :: i,j,n,Nballs,cont,a,b,c,veri,verfim,hori,horfim,penbola,ultbola
integer :: xis,ypsilon,nxis,nyip,hor,ver
real :: start,finish
call cpu_time(start)

 Nballs = 251 !Número de bolinhas

allocate(Sold(Nballs,2),Snow(Nballs,2),S(Nballs,2),v(Nballs,2),m(Nballs),R(Nballs))
allocate(tetaold(Nballs),tetanow(Nballs),teta(Nballs),omega(Nballs),Inercia(Nballs))
allocate(Fnormal(Nballs,Nballs,2),Fparede(Nballs,2),FR(Nballs,2))
allocate(Torque(Nballs),Fatrito(Nballs,2),ang(Nballs))

! As matrizez posição,velocidade e força resultante são estruturadas de forma que a linha i seja a bolinha i 
!e as colunas 1,2,3 as direções x,y,z respectivamente
! A matriz das forças de contato é estruturada sendo a força da bolinha i com a bolinha j localizada na posição Fij
!só que uma para cada dimensão sendo 1,2,3 x,y,z respectivamente
!**UPDATE** - Retirada da coordenada z para poder trabalhar com a rotação
!**************************************************
!Botar as paredes da caixa, só colocar um moveto e rlineto igual da linha do norte
!**************************************************
!Mudar o sinal da velocidade relativa nas paredes de cima e da esquerda, colocar -w*R
!!!!!!Tá com problema nas paredes de cima e da esquerda, ta girando ao contrário


 h = 10.0**(-5.d0)      !Passo de tempo
 rmedio = 0.3d0          !Raio medio das bolinhas
 K = 10.d0**(5.d0)      !Constante de elasticidade
 m = 10.d0**(-3.d0)     !Massa da bolinha
 gamaN = 9.d0           !Coeficiente da força dissipativa (Normal)
 gamaS = 9.d0           !Coeficiente da força dissipativa (Shear)
 Lx = 10.8d0!25.d0  para 1001            !Tamanho da parede do confinamento em x
 Ly = 20.d0!32.d0  para 1001          !Tamanho da parede do confinamento em y
 Ldown = 0.45d0         !Altura da parede de baixo
 g = -10.d0             !Gravidade
 minormal = 0.5d0       !Coeficiente de atrito entre as bolinhas
 miparede = 0.5d0       !Coeficiente de atrito da parede

!Amplitude da vibração da caixa
 amplitude = 0.4d0

!Subrotina da posiçao inicial
call pos_init(S,R,rmedio,Lx,Ly,Ldown,Lcell,nxis,nyip,Nballs,amplitude)
 v = 0.d0
 omega = 0.d0
 teta = 0.d0
 
 m(Nballs) = m(1)*(R(Nballs)/rmedio)**3
 
!!!!!!!!!!!!Para testar a rotação na parede
!Parede de baixo, para direita e para a esquerda check
!Parede da direita, para cima e para baixo check
!Parede da de cima, para direita e para a esquerda
 
allocate(Cell(0:nxis,0:nyip,Nballs))
allocate(marcabola(0:nxis,0:nyip))

!Como o método é de passos múltiplos, calcula-se uma iteração fazendo "s = s0 + v*t" para cada bolinha, analogamente para o ângulo
do i=1,Nballs
	Sold(i,:)= S(i,:)
	Snow(i,:)= S(i,:) + v(i,:)*h
	tetaold(i)= teta(i)
	tetanow(i)= teta(i) + omega(i)*h
end do

!Momento de inércia da ESFERA: (2/5)*mr^2
 Inercia = 2.d0*m*R*R/5.d0

 ang = mod(teta,6.28318530718)
!Para gerar o arquivo da condiçao inicial
 cont = 0
call salva_eps(cont , Lx, Ly, Nballs, R, S(:,1), S(:,2), ang)

!Loop para correr o tempo
do n=1,20*100000
	Cell = -1
	marcabola = -1
	
	do a=1,Nballs
		xis = int(S(a,1)/Lcell)
		ypsilon = int(S(a,2)/Lcell)
		
		select case(marcabola(xis,ypsilon))
		case(-1)
			marcabola(xis,ypsilon) = a
		case default
			Cell(xis,ypsilon,a) = marcabola(xis,ypsilon)
			marcabola(xis,ypsilon) = a
		end select
	end do

	!Forças e Torques recalculados a cada passo, logo, devem ser zerados
	Fnormal = 0.d0 
	Fparede = 0.d0
	Fatrito = 0.d0
	Torque = 0.d0
	FR(:,1) = 0.d0
	FR(:,2) = m(:)*g

	!Para chacoalhar a caixa
	Ldown = 0.45d0 + funcparede(n*h,amplitude)
	
	do b=0,nxis
		do c=0,nyip
			veri = 0
			verfim = 1
			
			hori= 0
			horfim= 1
			
			!If para perguntar se está em contato com alguma parede
				if(b.eq.0) then !Parede 1 (esquerda)
					j = marcabola(b,c)
				
					do while (j.ne.-1) !While para variar as bolas da célula
					
						if((S(j,1)-R(j)).lt.0.d0) then !Força da parede 1 (esquerda)
							call wall_force(Lx,Ly,Ldown,Nballs,j,K,gamaS,gamaN,minormal,S,v,omega,R,Torque,Fatrito,Fparede,FR,1)
						end if
						
						j = Cell(b,c,j)
					end do
					
				else if(b.eq.nxis-1) then !Parede 3 (direita)
					j = marcabola(b,c)
					
					do while (j.ne.-1) !While para variar as bolas da célula 
					
						if((S(j,1)+R(j)).gt.Lx) then !Força da parede 3 (direita)
							call wall_force(Lx,Ly,Ldown,Nballs,j,K,gamaS,gamaN,minormal,S,v,omega,R,Torque,Fatrito,Fparede,FR,3)
						end if
						
						j = Cell(b,c,j)
					end do
					
				else if(b.eq.nxis) then !Parede 3 (direita)
					horfim = 0
					j = marcabola(b,c)
					
					do while (j.ne.-1) !While para variar as bolas da célula	
					
						if((S(j,1)+R(j)).gt.Lx) then !Força da parede 3 (direita)
							call wall_force(Lx,Ly,Ldown,Nballs,j,K,gamaS,gamaN,minormal,S,v,omega,R,Torque,Fatrito,Fparede,FR,3)
						end if
						
						j = Cell(b,c,j)
					end do
				
				end if !Aqui acaba o if das paredes < e >
				
				if(c.eq.0) then !Parede 2 (de baixo)
					j = marcabola(b,c)
					
					do while (j.ne.-1) !While para variar as bolas da célula 
					
						if((S(j,2)-R(j)).lt.Ldown) then !Força da parede 2 (de baixo)
							call wall_force(Lx,Ly,Ldown,Nballs,j,K,gamaS,gamaN,minormal,S,v,omega,R,Torque,Fatrito,Fparede,FR,2)
						end if
						
						j = Cell(b,c,j)
					end do
					
				else if(c.eq.nyip-1) then !Parede 4 (de cima)
					j = marcabola(b,c)
					
					do while (j.ne.-1) !While para variar as bolas da célula 
					
						if((S(j,2)+R(j)).gt.Ly) then !Força da parede 4 (de baixo)
							call wall_force(Lx,Ly,Ldown,Nballs,j,K,gamaS,gamaN,minormal,S,v,omega,R,Torque,Fatrito,Fparede,FR,4)
						end if
						
						j = Cell(b,c,j)
					end do
					
				else if(c.eq.nyip) then !Parede 4 (de cima)
					verfim = 0
					j = marcabola(b,c)
					
					do while (j.ne.-1) !While para variar as bolas da célula 
					
						if((S(j,2)+R(j)).gt.Ly) then !Força da parede 4 (de baixo)
							call wall_force(Lx,Ly,Ldown,Nballs,j,K,gamaS,gamaN,minormal,S,v,omega,R,Torque,Fatrito,Fparede,FR,4)
						end if 
						
						j = Cell(b,c,j)
					end do
					
				end if !Aqui acaba o if das paredes \/ e /\ 
			
			!Última bola só varia dentro da própria célula
			!Penúltima bola muda de célula, fazendo penbola = marcabola
			
			do ver= veri,verfim
				do hor= hori,horfim
				
				ultbola = marcabola(b,c)
				do while (ultbola.ne.-1) !While para mudar a última bola
				
					if((hor.eq.0).and.(ver.eq.0)) then
						penbola = Cell(b,c,ultbola)
					else
						penbola = marcabola(b+hor,c+ver)
					end if
					
					do while (penbola.ne.-1) !While para perguntar entre as bolas de índices menores
						j = ultbola
						i = penbola
						D(:) = S(j,:)- S(i,:)
			
						if(norm2(D).lt.(R(j)+R(i))) then !Força de contato = ...
							call contact_force(Nballs,j,i,K,gamaS,gamaN,minormal,D,v,omega,R,Torque,Fatrito,Fnormal,FR)
						end if !Aqui acaba o if das forças de contato
			
			
						penbola = Cell(b+hor,c+ver,penbola)
					end do !Aqui acaba o while de dentro
		
					ultbola = Cell(b,c,ultbola)
				end do !Aqui acaba o while de fora
				
				end do !Aqui acaba o loop do hor
				
				if(b.ne.0) then
					hori = -1
				end if
				
			end do !Aqui acaba o loop do ver
		end do !Aqui acaba o loop do c
	end do !Aqui acaba o loop do b
		
	do j=1,Nballs
		do i=1,Nballs !Loop para calcular a força resultante de cada bolinha, somando suas coordenadas
			FR(j,:) = FR(j,:) + Fnormal(i,j,:)
		end do
		FR(j,:) = FR(j,:) + Fparede(j,:)
		
		!Metodo de Verlet
		!Sold: S(t-h) // tetaold: teta(t-h)
		!Snow: S(t)   // tetanow: teta(t)
		!S   : S(t+h) // teta   : teta(t+h)
		S(j,:) = 2.d0*Snow(j,:) - Sold(j,:) + (FR(j,:)/m(j))*(h**2.d0)
		v(j,:) = (S(j,:) - Sold(j,:))/(2.d0*h)
		Sold(j,:) = Snow(j,:)
		Snow(j,:) = S(j,:)
		
		teta(j) = 2.d0*tetanow(j) - tetaold(j) + (Torque(j)/Inercia(j))*(h**2.d0)
		omega(j) = (teta(j) - tetaold(j))/(2.d0*h)
		tetaold(j) = tetanow(j)
		tetanow(j) = teta(j)
	end do
	
	ang = mod(teta,6.28318530718)
	if(mod(n,7500).eq.0) then
		cont = cont + 1
		call salva_eps(cont , Lx, Ly, Nballs, R, S(:,1), S(:,2), ang)
	end if
end do !Aqui termina o loop do tempo

call cpu_time(finish)

deallocate(Sold,Snow,S,v,m,R)
deallocate(tetaold,tetanow,teta,omega,Inercia)
deallocate(Fnormal,Fparede,FR)
deallocate(Torque,Fatrito,ang)

contains

function sinal(x)
	implicit none
	double precision :: sinal,x

	if(x.eq.0.d0) then
		sinal = 0.d0
	else
		sinal = x/abs(x)
	end if
end function

function funcparede(x,amp)
 implicit none
 double precision:: funcparede,x,amp,w
 
 !Para calcular a força depois de ter aumentado até onde quiser
 ! Faça (tempofim - tempoinicio)**EXPOENTE = Aumento, ai resolve aplicando o log
 ! EX: tempofim do aumento = 4, tempo para começar a chacoalhar = 2, Aumento = 524 vezes
 !           (4-2)**EXP = 524 ->  2**EXP = 524  -> EXP = 9
 !Ou faz com exponencial rs
 
 !w = 2*pi*f
 w = 18.84955d0

 if((x.gt.1.5d0).and.(x.le.3.5d0)) then !Aumenta amplitude
	funcparede = (amp/1000.d0)*(exp(3.454d0*(x-1.5d0)))*sin(w*(x-1.5d0))
 else if((x.gt.3.5d0).and.(x.le.18.d0)) then !Mantem amplitude
	funcparede = amp*sin(w*(x-1.5d0))
 else if((x.gt.18.d0).and.(x.le.19.5d0)) then !Diminui amplitude ate parar
	funcparede = amp*(exp(-10.82*(x-18.d0)))*sin(w*(x-1.5d0))
 else !Antes de estabilizar
	funcparede = 0.d0
 end if

end function

subroutine contact_force(Nballs,j,i,K,gamaS,gamaN,minormal,D,v,omega,R,Torque,Fatrito,Fnormal,FR)
implicit none
integer, intent (in) :: Nballs,j,i
double precision, intent (in) :: K,gamaS,gamaN,minormal
double precision, dimension (2), intent (in) :: D
double precision, dimension (Nballs,2), intent (in) :: v
double precision, dimension (Nballs), intent (in) :: R,omega
double precision, dimension (Nballs), intent (inout) :: Torque
double precision, dimension (Nballs,2), intent (inout) :: Fatrito,FR
double precision, dimension (Nballs,Nballs,2), intent (inout) :: Fnormal

double precision, dimension (2) :: csi,n_dir,s_dir,vrel
double precision :: vrelnormal,vreltangente,rel,nor,sinalvtan

	n_dir= D/norm2(D)
	s_dir(1) = n_dir(2)
	s_dir(2) = -n_dir(1)
	
	vrel(:) = v(j,:) - v(i,:)
	
	vrelnormal = vrel(1)*n_dir(1) + vrel(2)*n_dir(2)
	vreltangente = vrel(1)*s_dir(1) + vrel(2)*s_dir(2) + omega(i)*R(i) + omega(j)*R(j)
	sinalvtan = sinal(vreltangente)
	
	csi = (R(j)+R(i) - norm2(D))*n_dir
	Fnormal(i,j,:) = K*csi(:) - gamaN*vrelnormal*n_dir(:)
	Fnormal(j,i,:) = -Fnormal(i,j,:)
	
	rel = gamaS*abs(vreltangente)
	nor = minormal*norm2(Fnormal(j,i,:))
	
	if(nor.gt.rel) then
		Fatrito(j,:) = - rel*sinalvtan*s_dir(:)
		Fatrito(i,:) =   rel*sinalvtan*s_dir(:)
		FR(j,:) = FR(j,:) + Fatrito(j,:)
		FR(i,:) = FR(i,:) + Fatrito(i,:)
		Torque(j) = Torque(j) - R(j)*norm2(Fatrito(j,:))*sinalvtan
		Torque(i) = Torque(i) - R(i)*norm2(Fatrito(i,:))*sinalvtan
	else
		Fatrito(j,:) = - nor*sinalvtan*s_dir(:)
		Fatrito(i,:) =   nor*sinalvtan*s_dir(:)
		FR(j,:) = FR(j,:) + Fatrito(j,:)
		FR(i,:) = FR(i,:) + Fatrito(i,:)
		Torque(j) = Torque(j) - R(j)*norm2(Fatrito(j,:))*sinalvtan
		Torque(i) = Torque(i) - R(i)*norm2(Fatrito(i,:))*sinalvtan
	end if

end subroutine

subroutine wall_force(Lx,Ly,Ldown,Nballs,j,K,gamaS,gamaN,minormal,S,v,omega,R,Torque,Fatrito,Fparede,FR,parede)
implicit none
integer, intent (in) :: Nballs,j,parede
double precision, intent (in) :: K,gamaS,gamaN,minormal,Lx,Ly,Ldown
double precision, dimension (Nballs,2), intent (in) :: S,v
double precision, dimension (Nballs), intent (in) :: R,omega
double precision, dimension (Nballs), intent (inout) :: Torque
double precision, dimension (Nballs,2), intent (inout) :: Fatrito,Fparede,FR

double precision, dimension (2) :: csi
double precision :: vrelnormal,vreltangente,rel,nor,sinalvtan

!Nas paredes de cima e da esquerda tem que trocar o sinal do omega*R para a velocidade relativa e do Torque também
!Parede 1 é a da esquerda, 2 de baixo, 3 da direita e 4 de cima
 select case(parede)
 case(1) !Parede da esquerda 
 	csi(1) = 0.d0 -(S(j,1)-R(j))
	Fparede(j,1) = K*csi(1) - gamaN*v(j,1)
	
	vrelnormal = v(j,1)
	vreltangente = v(j,2) - omega(j)*R(j)	
	sinalvtan = sinal(vreltangente)
	
	rel = gamaS*abs(vreltangente)
	nor = miparede*abs(Fparede(j,1))
	if(nor.gt.rel) then
		Fatrito(j,2) = - rel*sinalvtan
		FR(j,2) = FR(j,2) + Fatrito(j,2)
		Torque(j) = Torque(j) + R(j)*abs(Fatrito(j,2))*sinalvtan
	else
		Fatrito(j,2) = - nor*sinalvtan
		FR(j,2) = FR(j,2) + Fatrito(j,2)
		Torque(j) = Torque(j) + R(j)*abs(Fatrito(j,2))*sinalvtan
	end if	
	
 case(2) !Parede de baixo
 	csi(2) = Ldown -(S(j,2)-R(j))
	Fparede(j,2) = K*csi(2) - gamaN*v(j,2)
	
	vrelnormal = v(j,2)
	vreltangente = v(j,1) + omega(j)*R(j)	
	sinalvtan = sinal(vreltangente)
	
	rel = gamaS*abs(vreltangente)
	nor = miparede*abs(Fparede(j,2))
	if(nor.gt.rel) then
		Fatrito(j,1) = - rel*sinalvtan
		FR(j,1) = FR(j,1) + Fatrito(j,1)
		Torque(j) = Torque(j) - R(j)*abs(Fatrito(j,1))*sinalvtan
	else
		Fatrito(j,1) = - nor*sinalvtan
		FR(j,1) = FR(j,1) + Fatrito(j,1)
		Torque(j) = Torque(j) - R(j)*abs(Fatrito(j,1))*sinalvtan
	end if
	
 case(3) !Parede da direita
 	csi(1) = Lx-(S(j,1)+R(j))
	Fparede(j,1) = K*csi(1) - gamaN*v(j,1)
	
	vrelnormal = v(j,1)
	vreltangente = v(j,2) + omega(j)*R(j)	
	sinalvtan = sinal(vreltangente)
	
	rel = gamaS*abs(vreltangente)
	nor = miparede*abs(Fparede(j,1))
	if(nor.gt.rel) then
		Fatrito(j,2) = - rel*sinalvtan
		FR(j,2) = FR(j,2) + Fatrito(j,2)
		Torque(j) = Torque(j) - R(j)*abs(Fatrito(j,2))*sinalvtan
	else
		Fatrito(j,2) = - nor*sinalvtan
		FR(j,2) = FR(j,2) + Fatrito(j,2)
		Torque(j) = Torque(j) - R(j)*abs(Fatrito(j,2))*sinalvtan
	end if	
	
 case(4) !Parede de cima
	csi(2) = Ly-(S(j,2)+R(j))
	Fparede(j,2) = K*csi(2) - gamaN*v(j,2)
	
	vrelnormal = v(j,2)
	vreltangente = v(j,1) - omega(j)*R(j)
	sinalvtan = sinal(vreltangente)

	rel = gamaS*abs(vreltangente)
	nor = miparede*abs(Fparede(j,2))
	if(nor.gt.rel) then
		Fatrito(j,1) = - rel*sinalvtan
		FR(j,1) = FR(j,1) + Fatrito(j,1)
		Torque(j) = Torque(j) + R(j)*abs(Fatrito(j,1))*sinalvtan
	else
		Fatrito(j,1) = - nor*sinalvtan
		FR(j,1) = FR(j,1) + Fatrito(j,1)
		Torque(j) = Torque(j) + R(j)*abs(Fatrito(j,1))*sinalvtan
	end if	
	
 end select
	

end subroutine

subroutine pos_init(S,R,rmedio,Lx,Ly,Ldown,Lcell,nxis,nyip,Nballs,amp)
 implicit none

 double precision, dimension (Nballs,2), intent (inout) :: S
 double precision, dimension (Nballs), intent (inout) :: R
 double precision, intent (in) :: rmedio,Lx,Ly,Ldown,amp
 double precision, intent (inout) :: Lcell
 integer, intent (inout) :: nxis,nyip
 integer, intent (in) :: Nballs
 
 integer :: i,j,initrandom,premier,t,bolanut
 integer, dimension (3) :: timeArray
 double precision :: rnut,prox,rmax,Lepsilon,alt,proxb
 
 call itime(timeArray)
 initrandom = rand (timeArray(1)+timeArray(2)+timeArray(3))

 bolanut = Nballs
 
 rmax = rmedio*(1.1d0)
 
 do t= 1,Nballs-1
 	R(t) = rmedio*(1.d0 + rand(0)*(14.d0/100.d0))
 end do
 
 R(bolanut) = 4*rmedio
 
 Lcell = 2.d0*R(bolanut)
 Lepsilon = 2.d0*rmax
 
 nxis = int(Lx/Lcell)
 nyip = int((Ly+amp)/Lcell)
 
 !!!!!!!!!!!!!!Para a bola grande ficar embaixo
 S(bolanut,1) = Lx/2.d0
 S(bolanut,2) = Ldown + R(bolanut) + Lepsilon*0.1d0

 alt = Ldown + 0.55*Lepsilon 

 j= 1
 !!!!!!!While para colocar do lado da bola grande
 do while ((alt.lt.0.9d0*(S(bolanut,2) + R(bolanut))).and.(j.le.Nballs-1))
	S(j,1) = rand(0)*0.3*rmedio + 0.65*Lepsilon
	S(j,2) = alt
	
	j = j+1
	
	proxb = S(j-1,1) + R(j-1) + (rand(0)*0.3*rmedio + 0.15d0*Lepsilon) + 2.d0*R(j) 
	
	!While para colocar a esquerda da bola grande
 	do while ((proxb.lt.(S(bolanut,1)-R(bolanut)-0.05d0*Lepsilon)).and.(j.le.Nballs-1))
		S(j,1) = proxb - R(j)
		S(j,2) = S(j-1,2)
		
		j = j+1
		proxb = S(j-1,1) + R(j-1) + (rand(0)*rmedio + 0.1d0*Lepsilon) + 2.d0*R(j)
	end do
	
	proxb = S(bolanut,1) + R(bolanut) + (rand(0)*0.3*rmedio + 0.15d0*Lepsilon) + 2.d0*R(j)
	
	!While para colocar a direita da bola grande
	do while ((proxb.lt.(0.98d0*Lx)).and.(j.le.Nballs-1))
		S(j,1) = proxb - R(j)
		S(j,2) = S(j-1,2)
		
		j = j+1
		proxb = S(j-1,1) + R(j-1) + (rand(0)*rmedio + 0.1d0*Lepsilon) + 2.d0*R(j)
	end do
	
	alt = alt + 1.1*Lepsilon
 end do


 i= 0
 !j= 1
 !!!!!!!!!While para colocar as bolas em cima da bola grande
 do while (j.le.Nballs-1)
 	S(j,1) = rand(0)*rmedio + 0.8*Lepsilon
 	S(j,2) = Ldown + Lepsilon*(1.3d0*i + 0.55d0) + S(bolanut,2) + R(bolanut)
	
 	j= j+1
	
 	prox = S(j-1,1) + R(j-1) + (rand(0)*0.5*rmedio + 0.15d0*Lepsilon) + 2.d0*R(j) 
	
 	do while ((prox.lt.(0.98d0*Lx)).and.(j.le.Nballs-1))	
		S(j,1) = prox - R(j)
		S(j,2) = S(j-1,2)
		
		j = j+1
		if(j.ne.Nballs+1) then
			prox = S(j-1,1) + R(j-1) + (rand(0)*0.5*rmedio + 0.15d0*Lepsilon) + 2.d0*R(j)
		end if
	end do

	i= i+1
 end do

 
end subroutine

subroutine salva_eps(ttwrite,Lx,Ly,parts,part_raio,part_pos_x,part_pos_y,part_ang)
  implicit none

  integer, intent (in) :: ttwrite, parts
  
  double precision, dimension (1:parts), intent (in) :: part_raio
  double precision, dimension (1:parts), intent (in) :: part_pos_x
  double precision, dimension (1:parts), intent (in) :: part_pos_y
  double precision, dimension (1:parts), intent (in) :: part_ang
  double precision, intent (in) :: Lx,Ly

  integer :: i
  
  double precision :: scale = 10.d0
  integer :: xinfesq = -4
  integer :: yinfesq = -2
  integer :: xsupdir = 255
  integer :: ysupdir = 309
  
  character (len = 15) :: filename1
 
  !Criação do arquivo de saída - número entre 000 e 999!!!
  if(ttwrite.ne.0) then
     filename1="outputi.000.eps"
     filename1(11:11)=CHAR(48+mod(ttwrite,100)-10*(mod(ttwrite,100)/10))
     filename1(10:10)=CHAR(48+(mod(ttwrite,100)/10))
     filename1(9:9)=CHAR(48+ttwrite/100)
  else
	 filename1="outputi.000.eps"
  end if

  open(unit=210,file=filename1,status='unknown')

  !Escreve o cabeçalho inicial do eps
  write(210,90)
  write(210,91)
  write(210,92) filename1
  write(210,93)
  write(210,94)
  write(210,95)
  write(210,96) xinfesq, yinfesq, xsupdir,ysupdir
  write(210,97)
  write(210,98)
  write(210,99)
  write(210,100)
  write(210,101)
  write(210,102)
  write(210,103)
  !write particles without orientation
  do i=1,parts
     write(210,104) int(scale*part_pos_x(i)), int(scale*part_pos_y(i)), int(scale*part_raio(i))
     write(210,*) int(scale*part_pos_x(i)), int(scale*part_pos_y(i)), ' moveto' 
     write(210,*) int(scale*part_raio(i)*cos(part_ang(i))), int(scale*part_raio(i)*sin(part_ang(i))),' rlineto' 
     write(210,*) 'stroke'
  end do
  !!!!! Colocar as paredes
  write(210,*) 0,int(scale*Ldown),  ' moveto' 
  write(210,*) int(scale*Lx),0,     ' rlineto'
  write(210,*) 0,int(scale*Ly),     ' rlineto'
  write(210,*) -int(scale*Lx),0,    ' rlineto'
  write(210,*) 0,-int(scale*Ly),    ' rlineto'
  write(210,*) 'stroke'
  !!!!!!!!!!!!!!!!!!
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


end program colisao3d