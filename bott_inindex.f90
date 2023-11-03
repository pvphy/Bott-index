subroutine bottindex(i1,i2,trace)
    use global 
    use array
    implicit none
    integer::i,j,l,k,n,i1,nfill,i2
    doubleprecision::pi,l1,k1,x2,y2,z2,phi1,theta1,det,det1,trace,mu1,mu2
    complex*16::proj,proj_op(2*nos,2*nos),check(2*nos,2*nos),e_theta(2*nos,2*nos),e_phi(2*nos,2*nos),u1(2*nos,2*nos),w1(2*nos,2*nos)
    complex*16::u2(2*nos,2*nos),w2(2*nos,2*nos),prod(2*nos,2*nos),u3(2*nos,2*nos),w3(2*nos,2*nos)
    complex*16::u4(2,2),FindDet,proj_op1(2*nos,2*nos),det2,prod1(2*nos,2*nos),prod2(2*nos,2*nos)

    pi=4.0*atan(1.0)
    proj_op=complex(0.0d0,0.0d0)
    mu1=0.410d0
    mu2=2.40d0
    !nfill=2*nos/6
    nfill=0
    do i=1,2*nos
       if (evl_s(i).le.mu) nfill= nfill+1
    enddo
    !print*,nfill

    do i=1,2*nos    
        do j=1,2*nos 
            proj=cmplx(0.0d0,0.0d0)          
            do l=(i1-1)*nfill+1,(i2-1)*nfill
              
                proj=proj+h(i,l)*conjg(h(j,l))
                
            enddo
            proj_op(i,j)=proj
        enddo
    enddo
     

    do i=1,nos           
        if(tag(i).eq.1)then
            phi(i)=2*pi*x(i)/d
            theta(i)=2*pi*y(i)/d
        endif

        if(tag(i).eq.2)then 
            phi(i)=2*pi*x(i)/d
            theta(i)=2*pi*y(i)/d
        endif      

        if(tag(i).eq.3)then
            phi(i)=2*pi*x(i)/d
            theta(i)=2*pi*y(i)/d 
        endif     
    enddo


    do i=1,nos
        x2=(2+cos(theta(i)))*cos(phi(i))
        y2=(2+cos(theta(i)))*sin(phi(i))
        z2=sin(theta(i))
        write(57,*) phi(i),theta(i),x2,y2,z2
    enddo


    do i=1,nos
        e_theta(2*i-1,2*i-1)=exp(cmplx(0.0d0,theta(i)))
        e_theta(2*i,2*i)=exp(cmplx(0.0d0,theta(i)))
        e_phi(2*i-1,2*i-1)=exp(cmplx(0.0d0,phi(i)))
        e_phi(2*i,2*i)=exp(cmplx(0.0d0,phi(i)))
    enddo
 
    do i=1,2*nos
        write(30,*) ((e_theta(i,j)),j=1,2*nos)
        write(31,*) ((e_phi(i,j)),j=1,2*nos)
    enddo

    u3=matmul(proj_op,e_theta)
    u1=matmul(u3,proj_op)

    w3=matmul(proj_op,e_phi)
    w1=matmul(w3,proj_op)
 
    do i=1,2*nos
        write(34,*) ((u1(i,j)),j=1,2*nos)
        write(35,*) ((w1(i,j)),j=1,2*nos)
    enddo

    do i=1,2*nos
        do j=1,2*nos
            u2(i,j)=conjg(u1(j,i)) 
            w2(i,j)=conjg(w1(j,i)) 
        enddo
    enddo
 

    ! u2=transpose(conjg(u1))
    ! w2=transpose(conjg(w1))

    prod1=matmul(u1,w1)
    prod2=matmul(prod1,u2)
    prod=matmul(prod2,w2)

    
    ! call s_v_d(u1,2*nos,2*nos)

    ! do i=1,2*nos
    !     do j=1,2*nos
    !         write(3000,*) u_1(i,j),vt(i,j)  
    !     enddo
    ! enddo    

    call diagonalization1(prod)



    trace=0.0d0

    do i=1,2*nos
        if(abs(real(alpha(i))).gt.1e-6)then
            trace=trace+aimag(log(alpha(i)))/(2*pi)
        endif
       !print*,alpha(i),log(alpha(i))
    enddo

    
    ! u4(1,1)=cmplx(2,3)
    ! u4(1,2)=cmplx(0,5)
    ! u4(2,1)=cmplx(0,1)
    ! u4(2,2)=3

    ! u4=transpose(conjg(u4))
    !print*,u4(1,2),u4(2,1)
    ! det1=finddet(prod,2*nos)
    ! print*,det1

    
endsubroutine bottindex
