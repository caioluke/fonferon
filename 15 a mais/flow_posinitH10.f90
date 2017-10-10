program flow_posinit
implicit none
double precision, dimension (:,:), allocatable :: Sold,Snow,S,v
double precision, dimension (2) :: D
double precision, dimension (:,:,:), allocatable :: Fnormal
double precision, dimension (:,:), allocatable :: FR,Fparede,Fatrito,Froughness
double precision, dimension (:), allocatable :: R,Inercia,tetaold,tetanow,teta,omega,Torque,ang,m,Srough
integer, dimension (:,:,:), allocatable :: Cell,Cellrough
integer, dimension (:,:), allocatable :: marcabola,marcarough
double precision :: h,K,gamaN1,gamaN2,gamaS,Lx,Ly,g,minormal,miparede,Lcell,rmax,densidade
double precision :: rrough,deltarough,saveme,Hmax,flow_angle1
double precision :: scale,xinfesq,yinfesq,xsupdir,ysupdir
integer :: i,j,n,Nballs,cont,a,b,c,veri,verfim,hori,horfim,penbola,ultbola
integer :: xis,ypsilon,nxis,nyip,hor,ver,verclone,verclonei,verclonefim
integer :: Nroughs,rugoshor,rugosver,rugosi,rugosfim,nyiprough
real :: start,finish,hours,mins,secs
call cpu_time(start)

! As matrizez posição, velocidade e força resultante são estruturadas de forma que a linha i seja a bolinha i 
!e as colunas 1,2 as direções x,y respectivamente
! A matriz das forças de contato é estruturada sendo a força da bolinha i com a bolinha j localizada na posição Fij
!só que uma para cada dimensão sendo 1,2 x,y respectivamente
!**UPDATE** Subrotina das forças de contato e com a parede
!**UPDATE** Rotação certa nas paredes
!**UPDATE** Escala de unidades certa

!******************************IMPORTANTE*****************************************
!- O sistema de unidades é o SI
!- O número de bolas é o máximo para preencher o volume da caixa, exceto as 2 primeiras fileiras
!- As divisões de Lx/2*rmax e Ly/2*rmax, ou seja, Lx/Lcell e Ly/Lcell têm que ser inteiras
!- A divisão de Lx/(2.d0*rrough+deltarough) tem que ser inteira
!- Lcell = 2*rmax, ou seja, nenhum raio pode ser maior que rmax, portanto,
!para flutuar os outros raios subtraia desse valor

!ANTES DE SIMULAR:
!- Ajeitar o scale e a bounding box
!- Verificar se as condições estão sendo satisfeitas
!*********************************************************************************
!Dados que precisam satisfazer condições especiais

 rmax = 1.d0 / 100.d0                         !Maior raio das bolinhas (m)
 Lcell = 2.d0*rmax                            !Tamanho da célula (m)
 rrough = 2.d0*rmax                           !Raio da rugosidade (m)
 deltarough = 0.d0 / 100.d0                   !Espaço entre as rugosidades (m)
 Lx = 100.d0*Lcell                            !Tamanho da parede em x (m)
 Ly = (10.d0 + 2.d0 + 15.d0)*Lcell            !Tamanho da parede em y (m)
!*********************************************************************************
!Scale e bounding box

 if(Lx.lt.Ly) then
	scale = 500.d0/Lx
 else
	scale = 500.d0/Ly
 end if
 xinfesq = 0
 yinfesq = -1
 xsupdir = scale*Lx
 ysupdir = scale*Ly + 1
!*********************************************************************************
!Número de células em x, x tem uma a mais para clonar
 saveme = (Lx/Lcell)                          
 if((saveme - int(saveme)).le.0.1d0) then
	nxis = int(saveme)
 else if((saveme - int(saveme)).gt.0.9d0) then
	nxis = int(saveme)+1
 end if 
 
!Número de células em y
 saveme = (Ly/Lcell)                        
 if((saveme - int(saveme)).le.0.1d0) then
	nyip = int(saveme)-1
 else if((saveme - int(saveme)).gt.0.9d0) then
	nyip = int(saveme)
 end if
 
!Número de células em y para as rugosidades
 saveme = rrough/Lcell                      
 if((saveme - int(saveme)).le.0.1d0) then
	nyiprough = int(saveme)
 else if((saveme - int(saveme)).gt.0.9d0) then
	nyiprough = int(saveme)+1
 end if
 
!Número de rugosidades
 saveme = (Lx/(2.d0*rrough+deltarough))                   
 if((saveme - int(saveme)).le.0.1d0) then
	Nroughs = int(saveme)
 else if((saveme - int(saveme)).gt.0.9d0) then
	Nroughs = int(saveme)+1
 end if
 
 !Número de bolas
 !Tirei um a mais pra poder dar DEZ AAAAA células antes de bater na parede
 Nballs = (nxis)*(nyip-nyiprough-15)
 
!*********************************************************************************
 h = 10.0**(-4.d0)                            !Passo de tempo (s)
 K = 2.d0*10.d0**(5.d0)                       !Constante de elasticidade (N/m)
 densidade = 7860.d0                          !Densidade da bolinha (kg/m³)
 g = -9.8665d0                                !Gravidade (m/s²)
 
 !Coeficiente de restituição = 0.92
 gamaN1 = 3.0506d0                            !Coeficiente de dissipação entre as bolinhas
 gamaN2 = 4.3372d0                            !Coeficiente de dissipação entre a bola e a parede/rugosidade
 gamaS = 5.d0                                 !Coeficiente da força tangente, protege para o fat não dar pau
 minormal = 0.5d0                             !Coeficiente de atrito entre as bolinhas
 miparede = 0.5d0                             !Coeficiente de atrito da parede
!*********************************************************************************

 allocate(Sold(Nballs,2),Snow(Nballs,2),S(Nballs,2),v(Nballs,2),m(Nballs),R(Nballs))
 allocate(tetaold(Nballs),tetanow(Nballs),teta(Nballs),omega(Nballs),Inercia(Nballs))
 allocate(Fnormal(Nballs,Nballs,2),Fparede(Nballs,2),FR(Nballs,2),Froughness(Nballs,2))
 allocate(Torque(Nballs),Fatrito(Nballs,2),ang(Nballs))
 allocate(Cell(-1:nxis,0:nyip,Nballs),marcabola(-1:nxis,0:nyip))
 allocate(Cellrough(0:nxis-1,0:nyiprough,Nroughs),marcarough(0:nxis-1,0:nyiprough),Srough(Nroughs))
 
!Colocar as rugosidades em suas células
 Cellrough = -1
 marcarough = -1
 
 do a=1,Nroughs
	Srough(a) = rrough + (a-1)*1.d0*(2.d0*rrough + deltarough) 
	
	xis = int(Srough(a)/Lcell)
	!Case para saber se já existe alguma bola naquela célula de rugosidade
	!Caso marcarough = -1, esta bola é a primeira / Caso marcarough != -1, já existem bolas na célula	
	select case(marcarough(xis,0))
	case(-1)
		do b=0,nyiprough
			marcarough(xis,b) = a
		end do
	case default
		do b=0,nyiprough
			Cellrough(xis,0,a) = marcarough(xis,0)
			marcarough(xis,0) = a
		end do
	end select
	
 end do

!Subrotina da posiçao inicial
 call pos_init(S,v,teta,omega,ang,R,m,Inercia,densidade,rmax,Lcell,nxis,nyip,nyiprough,Nballs)

!Como o método é de passos múltiplos, calcula-se uma iteração fazendo "s = s0 + v*t" para cada bolinha, analogamente para o ângulo
 do i=1,Nballs
	Sold(i,:)= S(i,:)
	Snow(i,:)= S(i,:) + v(i,:)*h
	tetaold(i)= teta(i)
	tetanow(i)= teta(i) + omega(i)*h
 end do

!Para gerar o arquivo da condiçao inicial
 cont = 0
 call salva_eps(cont,Lx,Ly,Nballs,R,S(:,1),S(:,2),ang,Nroughs,rrough,Srough,scale,xinfesq,yinfesq,xsupdir,ysupdir)
 
!Loop para correr o tempo
do n=1,4*10000
	Cell = -1
	marcabola = -1
	
	do a=1,Nballs
		!A última célula (nxis) e a primeira (-1) servem apenas para calcular as forças NÃO HÁ BOLINHAS LÁ
		!OU SEJA, não tem nenhuma bola na célula -1 e nem na nxis, SÓ CLONES
		!A célula -1 representa a nxis-1, para as forças
		!A célula nxis representa a 0, para as forças
		!Para colocar na célula de verdade só deslocar Lx da posição dela
		!Desloca-se as posições anteriores para poder calcular o processo iterativo corretamente
		
		if(S(a,1).lt.0.d0) then
			Sold(a,1) = Sold(a,1) + Lx
			Snow(a,1) = Snow(a,1) + Lx
			S(a,1) = S(a,1) + Lx
		else if(S(a,1).gt.Lx) then	
			Sold(a,1) = Sold(a,1) - Lx
			Snow(a,1) = Snow(a,1) - Lx
			S(a,1) = S(a,1)	- Lx
		end if
		
		xis = int(S(a,1)/Lcell)	
		ypsilon = int(S(a,2)/Lcell)
		
		!Case para saber se já existe alguma bola naquela célula
		!Caso marcabola = -1, esta bola é a primeira / Caso marcabola != -1, já existem bolas na célula	
		select case(marcabola(xis,ypsilon))
		case(-1)
			marcabola(xis,ypsilon) = a
		case default
			Cell(xis,ypsilon,a) = marcabola(xis,ypsilon)
			marcabola(xis,ypsilon) = a
		end select
		
	end do
	
	!Clonar as bolas de nxis-1 para -1
	marcabola(-1,:) = marcabola(nxis-1,:)
	Cell(-1,:,:) = Cell(nxis-1,:,:)
	
	!Clonar as bolas de 0 para nxis
	marcabola(nxis,:) = marcabola(0,:)
	Cell(nxis,:,:) = Cell(0,:,:)
	
	!Forças e Torques recalculados a cada passo, logo, devem ser zerados
	Fnormal = 0.d0 
	Fparede = 0.d0
	Fatrito = 0.d0
	Froughness = 0.d0
	Torque = 0.d0
	FR(:,1) = 0.d0
	FR(:,2) = m(:)*g

	do c=0,nyip
		do b=0,nxis-1
			veri = 0
			verfim = 1
			
			hori= 0
			horfim= 1
			
			!O 0 tem que perguntar pro -1, a esquerda dele e na diagonal pra esquerda
			if(b.eq.0) then !Portal da esquerda
				verclonei = 0
				if(c.eq.nyip) then
					verclonefim = 0
				else
					verclonefim = 1
				end if
				
				!Essa parte pergunta para a esquerda e diagonal esquerda
				!Pergunta de 0 para -1 (nxis-1), quem está lá CLONADO de nxis-1 (último x do domínio)
				do verclone = verclonei,verclonefim !Do para variar verticalmente	
				
					ultbola = marcabola(0,c)
					do while (ultbola.ne.-1) !While para mudar a última bola					
						penbola = marcabola(-1,c+verclone)
					
						do while (penbola.ne.-1) !While para perguntar entre as bolas de índices menores
							j = ultbola
							i = penbola
							!Ajusta-se a posição da bola que é clonada, deslocando Lx para a esquerda
							D(1) = S(j,1) - (S(i,1) - Lx)
							D(2) = S(j,2) - S(i,2)
					
							if(norm2(D).lt.(R(j)+R(i))) then !Força de contato = ...
								call contact_force(Nballs,j,i,K,gamaS,gamaN1,minormal,D,v,omega,R,Torque,Fatrito,Fnormal,FR)
							end if
		
							penbola = Cell(-1,c+verclone,penbola)
						end do !Aqui acaba o while de dentro
						
						ultbola = Cell(0,c,ultbola)
					end do !Aqui acaba o while de fora
			
				end do !Aqui acaba o loop do verclone
					
			else if(b.eq.nxis-1) then !Portal da direita
				horfim = 0
				
				!Essa parte pergunta a força pra diagonal direita
				!Pergunta de nxis-1 para nxis(0), quem está lá CLONADO de 0
				if(c.ne.nyip) then 
					ultbola = marcabola(nxis-1,c)
					
					do while (ultbola.ne.-1) !While para mudar a última bola
						penbola = marcabola(nxis,c+1)
					
						do while (penbola.ne.-1) !While para perguntar entre as bolas de índices menores
							j = ultbola
							i = penbola
							!Ajusta-se a posição da bola que é clonada, deslocando Lx para a direita
							D(1) = S(j,1) - (S(i,1) + Lx)
							D(2) = S(j,2) - S(i,2)
						
							if(norm2(D).lt.(R(j)+R(i))) then !Força de contato = ...
								call contact_force(Nballs,j,i,K,gamaS,gamaN1,minormal,D,v,omega,R,Torque,Fatrito,Fnormal,FR)
							end if
		
							penbola = Cell(nxis,c+1,penbola)
						end do !Aqui acaba o while de dentro
									ultbola = Cell(nxis-1,c,ultbola)
					end do !Aqui acaba o while de fora
				end if !Aqui acaba o if de perguntar se o y é diferente do último
			
			end if !Aqui acaba o if dos portais < >
			
			!Interação com a parede
			if(c.eq.0) then !Parede 2 (de baixo)		
				j = marcabola(b,c)
				
				do while (j.ne.-1) !While para variar as bolas da célula 				
					if((S(j,2)-R(j)).lt.0.d0) then !Força da parede 2 (de baixo)
						call wall_force(Lx,Ly,Nballs,j,K,gamaS,gamaN2,miparede,S,v,omega,R,Torque,Fatrito,Fparede,FR,2)
					end if					
					j = Cell(b,c,j)
				end do
				
			else if(c.eq.nyip) then !Parede 4 (de cima)
				verfim = 0
				j = marcabola(b,c)
				
				do while (j.ne.-1) !While para variar as bolas da célula 				
					if((S(j,2)+R(j)).gt.Ly) then !Força da parede 4 (de baixo)
						call wall_force(Lx,Ly,Nballs,j,K,gamaS,gamaN2,minormal,S,v,omega,R,Torque,Fatrito,Fparede,FR,4)
					end if 						
					j = Cell(b,c,j)
				end do
				
			end if !Aqui acaba o if das paredes \/ e /\ 
			
			!Para calcular a força com as rugosidades
			if(c.le.(nyiprough)) then
				if(b.eq.0) then
					rugosi = 0
					rugosfim = 1
					
					!Tem que calcular a força com a rugosidade do outro lado, e corrigir a posição
					do rugosver = 0,nyiprough
				
						ultbola = marcabola(b,c)
					
						do while (ultbola.ne.-1) !While para mudar a bola dentro da célula escolhida
					
							penbola = marcarough(nxis-1,rugosver)
				
							do while (penbola.ne.-1) !While para perguntar entre as rugosidades de índices menores
								j = ultbola
								i = penbola
								D(1) = S(j,1) + Lx - Srough(i)
								D(2) = S(j,2)
			
								if(norm2(D).lt.(R(j)+rrough)) then !Força de contato = ...
									call roughness_force(Nballs,j,rrough,K,gamaS,gamaN2,minormal,D,v,omega,R,Torque,Fatrito,Froughness,FR)
								end if 
			
								penbola = Cellrough(nxis-1,rugosver,penbola)
							end do !Aqui acaba de variar as rugosidades
			
							ultbola = Cell(b,c,ultbola)
						end do !Aqui acaba de variar as bolas da célula escolhida				
					end do
					
					
				else if (b.eq.nxis-1) then
					rugosi = -1
					rugosfim = 0
					
					!Tem que calcular a força com a rugosidade do outro lado, e corrigir a posição
					do rugosver = 0,nyiprough
				
						ultbola = marcabola(b,c)
					
						do while (ultbola.ne.-1) !While para mudar a bola dentro da célula escolhida
					
							penbola = marcarough(0,rugosver)
				
							do while (penbola.ne.-1) !While para perguntar entre as rugosidades de índices menores
								j = ultbola
								i = penbola
								D(1) = S(j,1) - Lx - Srough(i)
								D(2) = S(j,2)
			
								if(norm2(D).lt.(R(j)+rrough)) then !Força de contato = ...
									call roughness_force(Nballs,j,rrough,K,gamaS,gamaN2,minormal,D,v,omega,R,Torque,Fatrito,Froughness,FR)
								end if 
			
								penbola = Cellrough(0,rugosver,penbola)
							end do !Aqui acaba de variar as rugosidades
			
							ultbola = Cell(b,c,ultbola)
						end do !Aqui acaba de variar as bolas da célula escolhida				
					end do
					
				else
					rugosi = -1
					rugosfim = 1
				end if
				
				do rugosver = 0,nyiprough
					do rugoshor = rugosi,rugosfim	
				
						ultbola = marcabola(b,c)
					
						do while (ultbola.ne.-1) !While para mudar a bola dentro da célula escolhida
					
							penbola = marcarough(b+rugoshor,rugosver)
				
							do while (penbola.ne.-1) !While para perguntar entre as rugosidades de índices menores
								j = ultbola
								i = penbola
								D(1) = S(j,1) - Srough(i)
								D(2) = S(j,2)
			
								if(norm2(D).lt.(R(j)+rrough)) then !Força de contato = ...
									call roughness_force(Nballs,j,rrough,K,gamaS,gamaN2,minormal,D,v,omega,R,Torque,Fatrito,Froughness,FR)
								end if 
			
								penbola = Cellrough(b+rugoshor,rugosver,penbola)
							end do !Aqui acaba de variar as rugosidades
			
							ultbola = Cell(b,c,ultbola)
						end do !Aqui acaba de variar as bolas da célula escolhida
				
					end do
				end do
			end if
			
			
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
								call contact_force(Nballs,j,i,K,gamaS,gamaN1,minormal,D,v,omega,R,Torque,Fatrito,Fnormal,FR)
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
		end do !Aqui acaba o loop do b
	end do !Aqui acaba o loop do c
		
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
	! if(mod(n,7500).eq.0) then
		! cont = cont + 1
		! call salva_eps(cont,Lx,Ly,Nballs,R,S(:,1),S(:,2),ang,Nroughs,rrough,Srough,scale,xinfesq,yinfesq,xsupdir,ysupdir)
		! write(*,*) cont,n*h
	! end if
	
 end do !Aqui termina o loop do tempo
 
 cont = cont + 1
 call salva_eps(cont,Lx,Ly,Nballs,R,S(:,1),S(:,2),ang,Nroughs,rrough,Srough,scale,xinfesq,yinfesq,xsupdir,ysupdir)
 
 Hmax = 0.d0
 do j=1,Nballs
	if(S(j,2).gt.Hmax) then
		Hmax = S(j,2)
		i = j
	end if
 end do
 Hmax = Hmax + R(i)

open(unit=51,file='initflowH10.dat',status='unknown')
 write(51,*) Nballs,Nroughs,rmax,Lx,Ly,Lcell,nxis,nyip,rrough,deltarough,nyiprough
 write(51,*) h,K,densidade,g,gamaN1,gamaN2,gamaS,minormal,miparede,Hmax
 write(51,*) scale,xinfesq,yinfesq,xsupdir,ysupdir

 do j=1,Nballs
	write(51,*) S(j,1),S(j,2),v(j,1),v(j,2),m(j)
	write(51,*) R(j),teta(j),omega(j),Inercia(j)
 end do
 
close(unit=51)
 
 
 
deallocate(Sold,Snow,S,v,m,R)
deallocate(tetaold,tetanow,teta,omega,Inercia)
deallocate(Fnormal,Fparede,FR)
deallocate(Torque,Fatrito,ang)
deallocate(Cell,marcabola)
deallocate(Cellrough,marcarough,Srough,Froughness)

 call cpu_time(finish)
 hours = int((finish-start)/3600)
 mins = int((finish-start)/60)
 secs = int((finish - start) - hours - mins)

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

subroutine contact_force(Nballs,j,i,K,gamaS,gamaN1,minormal,D,v,omega,R,Torque,Fatrito,Fnormal,FR)
 implicit none
 integer, intent (in) :: Nballs,j,i
 double precision, intent (in) :: K,gamaS,gamaN1,minormal
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
 Fnormal(i,j,:) = K*csi(:) - gamaN1*vrelnormal*n_dir(:)
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

subroutine wall_force(Lx,Ly,Nballs,j,K,gamaS,gamaN2,miparede,S,v,omega,R,Torque,Fatrito,Fparede,FR,parede)
 implicit none
 integer, intent (in) :: Nballs,j,parede
 double precision, intent (in) :: K,gamaS,gamaN2,miparede,Lx,Ly
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
	Fparede(j,1) = K*csi(1) - gamaN2*v(j,1)
	
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
 	csi(2) = 0.d0 -(S(j,2)-R(j))
	Fparede(j,2) = K*csi(2) - gamaN2*v(j,2)
	
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
	Fparede(j,1) = K*csi(1) - gamaN2*v(j,1)
	
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
	Fparede(j,2) = K*csi(2) - gamaN2*v(j,2)
	
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

subroutine roughness_force(Nballs,j,rrough,K,gamaS,gamaN2,miparede,D,v,omega,R,Torque,Fatrito,Froughness,FR)
 implicit none
 integer, intent (in) :: Nballs,j
 double precision, intent (in) :: K,gamaS,gamaN2,miparede,rrough
 double precision, dimension (2), intent (in) :: D
 double precision, dimension (Nballs,2), intent (in) :: v
 double precision, dimension (Nballs), intent (in) :: R,omega
 double precision, dimension (Nballs), intent (inout) :: Torque
 double precision, dimension (Nballs,2), intent (inout) :: Fatrito,FR
 double precision, dimension (Nballs,2), intent (inout) :: Froughness

 double precision, dimension (2) :: csi,n_dir,s_dir,vrel
 double precision :: vrelnormal,vreltangente,rel,nor,sinalvtan

 n_dir= D/norm2(D)
 s_dir(1) = n_dir(2)
 s_dir(2) = -n_dir(1)
 
 vrel(:) = v(j,:)
	
 vrelnormal = vrel(1)*n_dir(1) + vrel(2)*n_dir(2)
 vreltangente = vrel(1)*s_dir(1) + vrel(2)*s_dir(2) + omega(j)*R(j)
 sinalvtan = sinal(vreltangente)
 
 csi = (R(j)+ rrough - norm2(D))*n_dir
 Froughness(j,:) = K*csi(:) - gamaN2*vrelnormal*n_dir(:)
 
 rel = gamaS*abs(vreltangente)
 nor = miparede*norm2(Froughness(j,:))
 
 if(nor.gt.rel) then
	Fatrito(j,:) = - rel*sinalvtan*s_dir(:)
	FR(j,:) = FR(j,:) + Fatrito(j,:) + Froughness(j,:)
	Torque(j) = Torque(j) - R(j)*norm2(Fatrito(j,:))*sinalvtan
 else
	Fatrito(j,:) = - nor*sinalvtan*s_dir(:)
	FR(j,:) = FR(j,:) + Fatrito(j,:) + Froughness(j,:)
	Torque(j) = Torque(j) - R(j)*norm2(Fatrito(j,:))*sinalvtan
 end if

end subroutine

subroutine pos_init(S,v,teta,omega,ang,R,m,Inercia,densidade,rmax,Lcell,nxis,nyip,nyiprough,Nballs)
 implicit none

 double precision, dimension (Nballs,2), intent (inout) :: S,v
 double precision, dimension (Nballs), intent (inout) :: teta,omega,ang,R,m,Inercia
 double precision, intent (in) :: rmax,Lcell,densidade
 integer, intent (in) :: Nballs,nxis,nyip,nyiprough
 
 integer :: i,j,initrandom,bola,xis,yip
 integer, dimension (3) :: timeArray
 double precision :: prox
 
 call itime(timeArray)
 initrandom = rand (timeArray(1)+timeArray(2)+timeArray(3))

 bola = 0
 do yip= nyiprough+1,nyip
	do xis= 0,nxis-1
		bola = bola + 1
		
		if (bola.le.Nballs) then		
			S(bola,1) = Lcell*( ((0.5d0 + 1.0d0*xis)) + 0.03d0*(rand(0)-0.5d0) )
			S(bola,2) = Lcell*( ((0.5d0 + 1.0d0*yip)) + 0.03d0*(rand(0)-0.5d0) )
			R(bola)   = rmax*(0.9d0 + (rand(0)-0.5d0)*(5.d0/100.d0))			
		end if
		
	end do
 end do
 
 m = densidade*4.18879020479*R*R*R ! m = p*(4/3)*pi*R³
 Inercia = 0.4d0*m*R*R !Momento de inércia da ESFERA: (2/5)*mr^2
 v = 0.d0
 omega = 0.d0
 teta = 0.d0
 ang = mod(teta,6.28318530718)
 
end subroutine

subroutine salva_eps(ttwrite,Lx,Ly,parts,part_raio,part_pos_x,part_pos_y,part_ang,nroughs,rough_raio,rough_pos_x &
& ,scale,xinfesq,yinfesq,xsupdir,ysupdir)
 implicit none

 integer, intent (in) :: ttwrite, parts
  
 double precision, dimension (1:parts), intent (in) :: part_raio
 double precision, dimension (1:parts), intent (in) :: part_pos_x
 double precision, dimension (1:parts), intent (in) :: part_pos_y
 double precision, dimension (1:parts), intent (in) :: part_ang
 double precision, intent (in) :: Lx,Ly,scale,xinfesq,yinfesq,xsupdir,ysupdir

 double precision, dimension (1:nroughs), intent (in) :: rough_pos_x
 double precision, intent (in) :: rough_raio
 integer, intent (in) :: nroughs
 
 integer :: i
 double precision :: sclone
  
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
 
 !Escreve o cabeçalho inicial do EPS
 write(210,90)
 write(210,91)
 write(210,92) filename1
 write(210,93)
 write(210,94)
 write(210,95)
 write(210,96) int(xinfesq),int(yinfesq),int(xsupdir),int(ysupdir)
 write(210,97)
 write(210,98)
 write(210,99)
 write(210,100)
 write(210,101) 0.d0,0.d0,0.d0
 write(210,102)
 write(210,103)
  
 !Escreve as bolinhas
 do i=1,parts
 
	!Caso a bolinha esteja atravessando o portal, escrever o clone dela
	if((part_pos_x(i)-part_raio(i)).lt.0.d0) then	
		sclone = part_pos_x(i) + Lx
		write(210,104) (scale*sclone), (scale*part_pos_y(i)), (scale*part_raio(i))
		write(210,105) (scale*sclone), (scale*part_pos_y(i))
		write(210,106) (scale*part_raio(i)*cos(part_ang(i))), (scale*part_raio(i)*sin(part_ang(i)))
		write(210,*) 'stroke'  
	else if(part_pos_x(i)+part_raio(i).gt.Lx) then
		sclone = part_pos_x(i) - Lx
		write(210,104) (scale*sclone), (scale*part_pos_y(i)), (scale*part_raio(i))
		write(210,105) (scale*sclone), (scale*part_pos_y(i))
		write(210,106) (scale*part_raio(i)*cos(part_ang(i))), (scale*part_raio(i)*sin(part_ang(i)))
		write(210,*) 'stroke'	  
	end if
	
	write(210,104) (scale*part_pos_x(i)), (scale*part_pos_y(i)), (scale*part_raio(i))
	write(210,105) (scale*part_pos_x(i)), (scale*part_pos_y(i))
	write(210,106) (scale*part_raio(i)*cos(part_ang(i))), (scale*part_raio(i)*sin(part_ang(i)))
	write(210,*) 'stroke'
 end do

 !Coloca as paredes
 write(210,105) 0.d0,0.d0
 write(210,106) (scale*Lx),0.d0
 write(210,106) 0.d0,(scale*Ly)
 write(210,106) -(scale*Lx),0.d0
 write(210,106) 0.d0,-(scale*Ly)
 write(210,*) 'stroke'
 
 !Colocar as rugosidades
 write(210,101) 0.d0,0.d0,1.d0
 do i=1,nroughs
	write(210,104) (scale*rough_pos_x(i)), 0.d0, (scale*rough_raio)
	write(210,*) 'stroke'
 end do
 
 write(210,107)
 write(210,108)
 write(210,109)
 
90 format('%%!PS-Adobe-3.0 EPSF-3.0')
91 format('%%Document-Fonts: Times-Roman')
92 format('%%Title: ', A15)
93 format('%%Creator: CL')
94 format('%%CreationDate: unknown')
95 format('%%Pages: 1')
96 format('%%BoundingBox: ',4I5)
97 format('%%LanguageLevel: 1')
98 format('%%EndComments')
99 format('%%BeginProlog')
100 format('%%EndProlog')
101 format(3F10.3,' setrgbcolor')
102 format('%% Page:     1    1')
103 format('save')
104 format(3F10.3,'  0   360  arc')
105 format(2F10.3,' moveto')
106 format(2F10.3,' rlineto')
107 format('restore showpage')
108 format('%%Trailer')
109 format('%%EOF')
 
 close(unit=210)
  
end subroutine salva_eps

! subroutine outputfiles(flow_angle1)
! implicit none

! double precision, intent(in) :: flow_angle1

! integer :: j,flowint
! double precision:: saveme
! character (len = 19) :: filename1
!character (len = 17) :: filename2
! character (len = 18) :: filename2

 ! filename1="initflow0teta00.dat"
 !filename2="callprogteta00.sh"
 ! filename2="callprogteta00.bat"

 ! saveme = flow_angle1/0.01745329252d0                
 ! if((saveme - int(saveme)).le.0.1d0) then
	! flowint = int(saveme)
 ! else if((saveme - int(saveme)).gt.0.9d0) then
	! flowint = int(saveme) + 1
 ! end if

 ! filename1(9:9)=CHAR(48+)
 ! filename1(14:14)=CHAR(48+(mod(flowint,100)/10))
 ! filename1(15:15)=CHAR(48+mod(flowint,100)-10*(mod(flowint,100)/10))

 ! filename2(13:13)=CHAR(48+(mod(flowint,100)/10))
 ! filename2(14:14)=CHAR(48+mod(flowint,100)-10*(mod(flowint,100)/10))
 
 ! open(unit=51,file=filename1,status='unknown')
 ! open(unit=52,file=filename2,status='unknown')

 ! write(52,70)
 ! write(52,71)

 ! write(51,*) Nballs,Nroughs,rmax,Lx,Ly,Lcell,nxis,nyip,rrough,deltarough,nyiprough
 ! write(51,*) h,K,densidade,g,gamaN1,gamaN2,gamaS,minormal,miparede,Hmax
 ! write(51,*) scale,xinfesq,yinfesq,xsupdir,ysupdir
 ! write(51,*) flow_angle1
 
 ! do j=1,Nballs
	! write(51,*) S(j,1),S(j,2),v(j,1),v(j,2),m(j)
	! write(51,*) R(j),teta(j),omega(j),Inercia(j)
 ! end do

 ! call execute_command_line (",filename2,")
 
! 70 format('gfortran flow_flow')
!71 format ('./a.out &')
! 71 format ('a.exe')

! close(unit=51)
! close(unit=52)

! end subroutine

end program flow_posinit