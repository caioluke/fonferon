program flow
implicit none
double precision, dimension (:,:), allocatable :: Sold,Snow,S,v
double precision, dimension (2) :: D
double precision, dimension (:,:,:), allocatable :: Fnormal
double precision, dimension (:,:), allocatable :: FR,Fparede,Fatrito,Froughness,phi
double precision, dimension (:), allocatable :: R,Inercia,tetaold,tetanow,teta,omega,Torque,ang,m,Srough,packfrac
integer, dimension (:,:,:), allocatable :: Cell,Cellrough
integer, dimension (:,:), allocatable :: marcabola,marcarough
double precision :: h,K,gamaN1,gamaN2,gamaS,Lx,Ly,g,minormal,miparede,Lcell,rmax,densidade
double precision :: flow_angle,flow_angle1,rrough,deltarough,saveme,Hmax
double precision :: scale,tempogif
integer :: i,j,n,Nballs,cont,a,b,c,veri,verfim,hori,horfim,penbola,ultbola
integer :: xis,ypsilon,nxis,nyip,hor,ver,verclone,verclonei,verclonefim
integer :: Nroughs,rugoshor,rugosver,rugosi,rugosfim,nyiprough,ninit,filecont,ty
real :: start,finish,days,hours,mins
call cpu_time(start)

 cont = 0  !Cont começa em 0 mas se tiver que ler do backup só comentar aqui
 ninit = 1 !O tempo começa em 1, mas caso tenha que ler do backup ele atualiza o valor ali embaixo

 filecont = 0 !Cont para salvar os dados
 
open(unit=31,file='packfrac.dat',status='unknown') 
open(unit=30,file='initflow.dat',status='old')
!open(unit=69,file='backup.dat',status='old')

 read(30,*) Nballs,Nroughs,rmax,Lx,Ly,Lcell,nxis,nyip,rrough,deltarough,nyiprough
 read(30,*) h,K,densidade,g,gamaN1,gamaN2,gamaS,minormal,miparede,Hmax
 read(30,*) scale
 !read(69,*) ninit,cont
 
!close(69)

!*********************************************************************************

 allocate(Sold(Nballs,2),Snow(Nballs,2),S(Nballs,2),v(Nballs,2),m(Nballs),R(Nballs))
 allocate(tetaold(Nballs),tetanow(Nballs),teta(Nballs),omega(Nballs),Inercia(Nballs))
 allocate(Fnormal(Nballs,Nballs,2),Fparede(Nballs,2),FR(Nballs,2),Froughness(Nballs,2))
 allocate(Torque(Nballs),Fatrito(Nballs,2),ang(Nballs))
 allocate(Cell(-1:nxis,0:nyip,Nballs),marcabola(-1:nxis,0:nyip))
 allocate(Cellrough(0:nxis-1,0:nyiprough,Nroughs),marcarough(0:nxis-1,0:nyiprough),Srough(Nroughs))
 allocate(phi(0:nxis-1,0:nyip),packfrac(0:nyip))
 
!******************************IMPORTANTE*****************************************
!- A condição inicial e os parâmetros foram importados do programa 'posinitTESTE.f90'
!- Arquivo 30, ler a condição inicial normal
!- Arquivo 69, ler do backup
!- Ajustar o ângulo antes de simular

 flow_angle1 = (14.d0)*0.01745329252d0

 do j=1,Nballs
	read(30,*) S(j,1),S(j,2),v(j,1),v(j,2),m(j)
	read(30,*) R(j),teta(j),omega(j),Inercia(j)
 end do
!*********************************************************************************

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

!Como o método é de passos múltiplos, calcula-se uma iteração fazendo "s = s0 + v*t" para cada bolinha, analogamente para o ângulo
 do i=1,Nballs
	Sold(i,:)= S(i,:)
	Snow(i,:)= S(i,:) + v(i,:)*h
	tetaold(i)= teta(i)
	tetanow(i)= teta(i) + omega(i)*h
 end do

!Para gerar o arquivo da condiçao inicial
 call salva_eps(cont,Lx,Ly,Nballs,R,S(:,1),S(:,2),ang,Nroughs,rrough,Srough,v,scale,0.d0)

!Loop para correr o tempo
 do n=ninit,80*100000
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
	
	!Para o ângulo de escoamento, inicialmente dá-se um ângulo maior para retirarmos o efeito da "história"
	if(n.le.1000) then
		flow_angle = (28.d0)*0.01745329252d0
	else 
		flow_angle = flow_angle1
	end if
	FR(:,1) = m(:)*g*sin(flow_angle)
	FR(:,2) = m(:)*g*cos(flow_angle)
	
	!Forças de contato
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
	
	!Para escrever o EPS
	if(n.ge.7000000) then
		if(mod(n,10000).eq.0) then
			cont = cont + 1
			tempogif = n*h
			call salva_eps(cont,Lx,Ly,Nballs,R,S(:,1),S(:,2),ang,Nroughs,rrough,Srough,v,scale,tempogif)
		end if
	end if
	
	!Para salvar os dados
	if(n.ge.4000000) then
		if(mod(n,50000).eq.0) then		
			filecont = filecont +1
		
			call salva_phi(filecont,phi,packfrac,S,R,Nballs,nxis,nyip,Lcell,Cell,marcabola)
			do ty=1,nyip
				write(31,*) ty,packfrac(ty),n*h
			end do
		
		end if
	end if
	
	!Para escrever o arquivo de backup
	if(mod(n,10*100000).eq.0) then
		open(unit=69,file='backup.dat',status='unknown')
		
		write(69,*) Nballs,Nroughs,rmax,Lx,Ly,Lcell,nxis,nyip,rrough,deltarough,nyiprough
		write(69,*) h,K,densidade,g,gamaN1,gamaN2,gamaS,minormal,miparede,Hmax
		write(69,*) scale
		write(69,*) n+1,cont

		do j=1,Nballs
			write(69,*) S(j,1),S(j,2),v(j,1),v(j,2),m(j)
			write(69,*) R(j),teta(j),omega(j),Inercia(j)
		end do	
		
		close(69)
	end if
	
 end do !Aqui termina o loop do tempo
 
close(unit=30)
close(unit=31)
 
deallocate(Sold,Snow,S,v,m,R)
deallocate(tetaold,tetanow,teta,omega,Inercia)
deallocate(Fnormal,Fparede,FR)
deallocate(Torque,Fatrito,ang)
deallocate(Cell,marcabola)
deallocate(Cellrough,marcarough,Srough,Froughness)
deallocate(phi,packfrac)

 call cpu_time(finish)
 days = int((finish-start)/(24*60*60))
 hours = mod(int((finish-start)/(60*60)),24)
 mins = mod(int((finish-start)/60),60)
 
open(unit=25,file='dadosflow.txt',status='unknown')
 write(25,*) "Número de bolas :", Nballs
 write(25,*) "Altura no inicio do escoamento (Hmax) :", Hmax
 write(25,*) "Altura que o cara usa no artigo (N*d/L) :", (Nballs*2.d0*rmax)/(Lx)
 write(25,*) "Ângulo de escoamento (graus) :", flow_angle*57.29577951290d0
 write(25,*) "Tempo de simulação (segundos) :", n*h
 write(25,*) "O tempo para rodar foi de ",days," dias ",hours," horas ",mins," minutos"
close(25)
 
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

subroutine salva_eps(ttwrite,Lx,Ly,parts,part_raio,part_pos_x,part_pos_y,part_ang,nroughs,rough_raio,rough_pos_x &
& ,v,scale,tempogif)
 implicit none

 integer, intent (in) :: ttwrite, parts
  
 double precision, dimension (1:parts), intent (in) :: part_raio
 double precision, dimension (1:parts), intent (in) :: part_pos_x
 double precision, dimension (1:parts), intent (in) :: part_pos_y
 double precision, dimension (1:parts), intent (in) :: part_ang
 double precision, intent (in) :: Lx,Ly,scale,tempogif
 
 double precision, dimension (1:parts,2), intent (in) :: v
 
 double precision, dimension (1:nroughs), intent (in) :: rough_pos_x
 double precision, intent (in) :: rough_raio
 integer, intent (in) :: nroughs
 
 double precision, dimension (1:parts,3) :: RGB
 double precision :: sclone,vmax
 integer :: i
  
 character (len = 16) :: filename1
 
 !Criação do arquivo de saída - número entre 000 e 999!!!

 filename1="outputi.0000.eps"
 filename1(12:12)=CHAR(48+mod(ttwrite,10))
 filename1(11:11)=CHAR(48+(mod(ttwrite,100)/10))
 filename1(10:10)=CHAR(48+(mod(ttwrite,1000)/100))
 filename1(9:9)  =CHAR(48+ttwrite/1000)

 open(unit=210,file=filename1,status='unknown')
 
 !Escreve o cabeçalho inicial do EPS
 write(210,90)
 write(210,91)
 write(210,92) filename1
 write(210,93)
 write(210,94)
 write(210,95)
 write(210,96) -int(scale*Ly*0.5d0),0,int(scale*Lx),int(scale*(0.5d0*Lx+0.866d0*Ly))
 write(210,97)
 write(210,98)
 write(210,99)
 write(210,100)
 write(210,102)
 write(210,103)
 write(210,*) '/Times-Roman 24 selectfont'
 write(210,105) -scale*Ly*0.29d0,scale*0.9d0*Ly
 write(210,111) tempogif
 write(210,*) 'stroke'
 write(210,*) '/Times-Roman 9 selectfont'
 write(210,*) '20 rotate'
 
  vmax = 0.d0
 do i=1,parts
	if(vmax.lt.norm2(v(i,:))) then
		vmax = norm2(v(i,:))
	end if
 end do
  
 do i=1,parts
	RGB(i,1) = norm2(v(i,:))/vmax
	RGB(i,2) = 0.d0
	RGB(i,3) = 1.d0 - (norm2(v(i,:))/vmax)
 end do
 
 !Escreve as bolinhas
 do i=1,parts
 
 !**************************************************************
 !                       IMPORTANTE
 !**************************************************************
 !        Se quiser escrever a bola com o número, tem que tirar
 !a linha tanto da bola normal, quanto das bolas que atravessam
 !o portal
 
	!!!Caso a bolinha esteja atravessando o portal, escrever o clone dela
	write(210,101) RGB(i,1),RGB(i,2),RGB(i,3) !Colocar aqui a parada da velocidade
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
	
	!Sem número
	write(210,104) (scale*part_pos_x(i)), (scale*part_pos_y(i)), (scale*part_raio(i))
	write(210,105) (scale*part_pos_x(i)), (scale*part_pos_y(i))
	write(210,106) (scale*part_raio(i)*cos(part_ang(i))), (scale*part_raio(i)*sin(part_ang(i)))
	write(210,*) 'stroke'
	
	!Com número
	! write(210,101) 0.d0,0.d0,1.d0
	! write(210,105) (scale*(part_pos_x(i))),&!-1.d0*part_raio(i)*cos(part_ang(i))) ), &
	! & (scale*(part_pos_y(i)))!-1.d0*part_raio(i)*sin(part_ang(i))) )	
	! write(210,112) i
	! write(210,*) 'stroke'
	
	! write(210,101) 0.d0,0.d0,0.d0
	! write(210,104) (scale*part_pos_x(i)), (scale*part_pos_y(i)), (scale*part_raio(i))
	! write(210,105) (scale*part_pos_x(i)), (scale*part_pos_y(i))
	! write(210,106) (scale*part_raio(i)*cos(part_ang(i))), (scale*part_raio(i)*sin(part_ang(i)))
	! write(210,*) 'stroke'
	
 end do
 
 !Colocar as rugosidades, tem que vir antes de limpar em volta da caixa
 write(210,101) 0.467d0,0.533d0,0.6d0
 do i=1,nroughs
	write(210,104) (scale*rough_pos_x(i)), 0.d0, (scale*rough_raio)
	write(210,*) 'closepath gsave fill grestore'
	write(210,*) 'stroke'
 end do

 
 !Pra poder rotacionar, limpar em volta da caixa, pra tirar as transições das bolinhas
 write(210,101) 1.d0,1.d0,1.d0
 write(210,105) -(scale*Lx), -(scale*Ly)
 write(210,106) 0.d0, 2.d0*(scale*Ly)
 write(210,106) (scale*Lx), 0.d0
 write(210,106) 0.d0, -(scale*Ly)
 write(210,106) (scale*Lx), 0.d0
 write(210,106) 0.d0, (scale*Ly)
 write(210,106) (scale*Lx), 0.d0
 write(210,106) 0.d0, -2.d0*(scale*Ly)
 write(210,*) 'closepath gsave fill grestore'
 write(210,*) 'stroke'
 
 !Coloca as paredes
 write(210,101) 0.467d0,0.533d0,0.6d0
 write(210,105) 0.d0,0.d0
 write(210,106) (scale*Lx),0.d0
 write(210,106) 0.d0,(scale*Ly)
 write(210,106) -(scale*Lx),0.d0
 write(210,106) 0.d0,-(scale*Ly)
 write(210,*) 'stroke'
 
 !Pra testar o cutoff
 ! write(210,101) 0.d0,0.d0,0.d0
 ! do i=1,nxis
	! write(210,105) i*scale*Lcell,0.d0
	! write(210,106) 0.d0,(scale*Ly)
	! write(210,*) 'stroke'
 
	! write(210,105) 0.d0,i*scale*Lcell
	! write(210,106) (scale*Lx),0.d0
	! write(210,*) 'stroke'
 ! end do
 
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
110 format(5F10.3,'  arc')
111 format('(',F7.2,') show' )
112 format('(',I5,') show' )
 
 close(unit=210)
  
end subroutine salva_eps

subroutine salva_phi(ttwrite,phi,packfrac,S,R,Nballs,nxis,nyip,Lcell,Cell,marcabola)
 implicit none
 double precision, dimension (0:nxis-1,0:nyip), intent (inout) :: phi
 double precision, dimension (0:nyip), intent (inout) :: packfrac
 double precision, dimension (Nballs,2), intent (in) :: S
 double precision, dimension (Nballs), intent (in) :: R
 integer, dimension (-1:nxis,0:nyip,Nballs), intent (in) :: Cell
 integer, dimension (-1:nxis,0:nyip), intent (in) :: marcabola
 double precision, intent (in) :: Lcell
 integer, intent (in) :: Nballs,nxis,nyip,ttwrite

 integer :: i,j,highestball,ultbola,penbola,b,c,ver,hor,tx,ty
 double precision :: Hmax,saveme,area
 
 character (len = 11) :: filename1

	filename1="phi.000.dat"
    filename1(7:7)=CHAR(48+mod(ttwrite,100)-10*(mod(ttwrite,100)/10))
    filename1(6:6)=CHAR(48+(mod(ttwrite,100)/10))
    filename1(5:5)=CHAR(48+ttwrite/100)
 
 !**************************************************************
 !                       IMPORTANTE
 !**************************************************************
 ! PHI foi feito pra rrough = 2*rmax ***
 ! Tirei 3 células da esquerda <, 3 da direita > para não ter problemas de contar duas vezes a mesma bola (portal)
 ! Tirei 3 (por isso foi feito p este rrough) células de baixo, para não ter problema com a rugosidade
 !e tirei 2 de cima da maior altura, para não ter problema com pegar células vazias
 
 Hmax = 0.d0
 do j=1,Nballs
	if(S(j,2).gt.Hmax) then
		Hmax = S(j,2)
		i = j
	end if
 end do

 saveme = (S(i,2)/Lcell)

 if((saveme - int(saveme)).le.0.1d0) then
	highestball = int(saveme)-1
 else
	highestball = int(saveme)
 end if
 
 phi = 0.d0
 packfrac = 0.d0
 
 
 do c=3,highestball-2!2,highestball-1  
	do b=3,nxis-4!2,nxis-3         
		
		area = 0.d0
		
		do ver= -2,2	!-1,1
			do hor= -2,2	!-1,1
				
				ultbola = marcabola(b+hor,c+ver)
				
				do while (ultbola.ne.-1) !While para mudar a última bola					
					area = area + 3.14159265359d0*R(ultbola)*R(ultbola)
					ultbola = Cell(b,c,ultbola)						
				end do		
	
			end do !Aqui acaba o loop do hor				
		end do !Aqui acaba o loop do ver
		
		phi(b,c) = area/(25.d0*Lcell*Lcell)!(9.d0*Lcell*Lcell) !(25.d0*Lcell*Lcell)
 
	end do !Aqui acaba o loop do b
 end do !Aqui acaba o loop do c
 
 do j=0,nyip
	do i=0,nxis-1!3,nxis-4!2,nxis-3
		packfrac(j) = packfrac(j) + phi(i,j)
	end do
	
	packfrac(j) = packfrac(j)/(1.d0*(nxis-6))
 end do
 
open(unit=210,file=filename1,status='unknown') 
 do ty=1,nyip
	do tx=1,nxis-1
		write(210,*) tx,ty,phi(tx,ty)
	end do
 end do
 
close(unit=210)
 
end subroutine

end program flow