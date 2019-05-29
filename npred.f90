program main

  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer i,lfn,mlL,rep,species,order; character*80 fileName; character*3 tpc
  real*8  freqL1(4),tranL1(4,4),tranL2(16,4),tranL3(64,4),tranL4(256,4)
  real*8  freqN4(64,4),tranN4((147-4)*256,4),freqN1(147,4),tranN1(584,4)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(1,file='freqL.txt')  !##### freqL1
  read(1,*) freqL1; close(1)
  open(1,file='tranL.txt')  !##### tranL1
  do i=1,4; read(1,*) tranL1(i,:); end do; close(1)
  open(1,file='tranL2.txt') !##### tranL2
  do i=1,16; read(1,*) tranL2(i,:); end do; close(1)
  open(1,file='tranL3.txt') !##### tranL3
  do i=1,64; read(1,*) tranL3(i,:); end do; close(1)
  open(1,file='tranL4.txt') !##### tranL4
  do i=1,256; read(1,*) tranL4(i,:); end do; close(1)

  open(1,file='147freqN.txt') !freqN1
  do i=1,147; read(1,*) freqN1(i,:); end do
  close(1)
  open(2,file='147tranN.txt') !tranN1
  do i=1,584; read(2,*) tranN1(i,:); end do
  close(2)

  open(1,file='146-149freqN4.txt') !##### freqN4
  do i=1,64; read(1,*) freqN4(i,:); end do; close(1)
  open(1,file='146-149tranN4.txt') !##### tranN4
  do i=1,(147-4)*256; read(1,*) tranN4(i,:); end do; close(1)

  write(*,'(a)',advance='no') 'Input the fileName of sequence                  :  '; read*,fileName
  write(*,'(a)',advance='no') 'Input the maximum length of linker              :  '; read*,mlL
  write(*,'(a)',advance='no') 'Input the order of selected model (1 or 4)      :  '; read*,order
  if(order/=1.and.order/=4) then; print*,'1 or 4 should be inputed!'; stop; end if
  write(*,'(a)',advance='no') 'Update the linker length distribution? yes or no:  '; read*,tpc
  if(tpc=='yes') then; rep=2; else if(tpc=='no ') then; rep=1; else; print*,'yes or no should be inputed!'; stop; end if
  print*,' '
  print*,'Select the species from the following list:' 
  print*,'1=Human;         2=Mouse;           3=Rat;'
  print*,'4=Zebrafish;     5=D. melanogaster; 6=C. elegans;'
  print*,'7=S. cerevisiae; 8=C. albicans;     9=S. pombe;'
  print*,'10=A. thaliana;  11=Maize.'
  write(*,'(a)',advance='no') 'Input the number of selected species            :  '; read*,species
  print*,'Calculating......'

  lfn=len_trim(fileName)
  if(order==1) then
    call vtbfb(lfn,trim(fileName),freqL1,tranL1,freqN1,tranN1,mlL,rep,species)
  else if(order==4) then
    call vtbfbNL4(lfn,trim(fileName),freqL1,tranL1,tranL2,tranL3,tranL4,freqN4,tranN4,mlL,rep,species)
  end if
  
end
!#################################################################################################
subroutine vtbfb(lfn,filename,freqL,tranL,freqN,tranNN,mlL,rep,species)

  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer   t,d,i,j,k,l,m,n,z,mlL,lfn,rep,nvt,species
  integer   zt,np,zp(100),q,st(100),ed(100),ext,rts(100),rte(100)
  integer   ztt,nsec,x,base; integer,allocatable:: nrcv(:),nbrk(:)
  real*8    tempN,tempL,temp,tp,sc(4)
  real*8    freqL(4),tranL(4,4),freqN(147,4),tranN(146,4,4),tranNN(584,4),Pd(mlL),tmp(mlL)

  character*80 seqName,out1,out2
  character(len=lfn) filename
  character(len=10000000),allocatable:: a(:)

  integer*1,allocatable:: wt(:),w(:),c(:)
  integer,allocatable:: change(:),Nst(:)
  real*8,allocatable:: r(:),hatFN(:),hatFL(:),ra(:),hatAN(:),hatAL(:),rb(:),hatBN(:),hatBL(:)
  real*8,allocatable:: ppEndN(:),ppN(:)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1,146; do j=1,4; tranN(i,j,:)=tranNN(4*(i-1)+j,:); end do; end do

  if(species==1) then !Human
    sc(1)=0.9690; sc(2)=1.050
  else if(species==2) then !Mouse
    sc(1)=0.9540; sc(2)=1.074
  else if(species==3) then !Rat
    sc(1)=0.9280; sc(2)=1.116
  else if(species==4) then !Zebrafish
    sc(1)=1.028; sc(2)=0.9543
  else if(species==5) then !D.melanogaster
    sc(1)=0.9313; sc(2)=1.111
  else if(species==6) then !C.elegans
    sc(1)=1.045; sc(2)=0.9269
  else if(species==7) then !S.cerevisiae
    sc(1)=1.0; sc(2)=1.0
  else if(species==8) then !C.albicans
    sc(1)=1.077; sc(2)=0.8757
  else if(species==9) then !S.pombe
    sc(1)=1.033; sc(2)=0.9462
  else if(species==10) then !A.thaliana
    sc(1)=1.038; sc(2)=0.9384
  else if(species==11) then !Maize
    sc(1)=0.8596; sc(2)=1.226
  else
    print*,'unavailable for this species number, quit!'; stop
  end if
  sc(3)=sc(2); sc(4)=sc(1)
  !------------------ scaling -----------------------
  do k=1,4
    freqL(k)=freqL(k)*sc(k); tranL(:,k)=tranL(:,k)*sc(k)
    freqN(:,k)=freqN(:,k)*sc(k); tranN(:,:,k)=tranN(:,:,k)*sc(k)
  end do
  freqL=freqL/sum(freqL); do j=1,4; tranL(j,:)=tranL(j,:)/sum(tranL(j,:)); end do
  do i=1,147
    freqN(i,:)=freqN(i,:)/sum(freqN(i,:))
    if(i<147) then
      do j=1,4; tranN(i,j,:)=tranN(i,j,:)/sum(tranN(i,j,:)); end do
    end if
  end do
  !-------------------------------------------------
  allocate(a(1))
  open(1,file=filename)
  read(1,'(1x,a80)') seqName
  read(1,'(a)') a(1); m=len_trim(a(1)); k=1; j=m
  read(1,'(a)',iostat=l) a(1)(1:m)
  do while(l==0)
    k=k+1; j=len_trim(a(1)(1:m))
    read(1,'(a)',iostat=l) a(1)(1:m)
  end do
  ztt=(k-1)*m+j
  close(1)
  if(ztt<148) then
    print*,'The sequence length is shorter than 148bp, quit!'; stop
  else if(148<=ztt.and.ztt<200000.and.rep==2) then
    print*,'Updating linker length distribution is not allowed since the sequence is not long enough.'
    print*,'The linker length distribution is still assumed uniform.'
    rep=1
  end if; allocate(wt(0:ztt+1))

  wt(0)=0; nsec=0
  open(1,file=filename)
  read(1,'(1x,a80)') seqName
  do i=1,k
    read(1,'(a)') a(1)(1:m); n=m; if(i==k) n=j
    do l=1,n; t=(i-1)*m+l
      if(a(1)(l:l)=='A'.or.a(1)(l:l)=='a') then
        wt(t)=1
      else if(a(1)(l:l)=='C'.or.a(1)(l:l)=='c') then
        wt(t)=2
      else if(a(1)(l:l)=='G'.or.a(1)(l:l)=='g') then
        wt(t)=3
      else if(a(1)(l:l)=='T'.or.a(1)(l:l)=='t') then
        wt(t)=4
      else
        wt(t)=0
      end if
      if(wt(t-1)==0.and.wt(t)>0) nsec=nsec+1
    end do
  end do
  close(1); deallocate(a); allocate(nrcv(nsec),nbrk(nsec))

  nsec=0; wt(ztt+1)=0
  do i=1,ztt+1
    if(wt(i-1)==0.and.wt(i)>0) then
      nsec=nsec+1; nrcv(nsec)=i
    else if(wt(i-1)>0.and.wt(i)==0) then
      nbrk(nsec)=i
      if(nbrk(nsec)-nrcv(nsec)<148) nsec=nsec-1
    end if
  end do

  out1(1:lfn)=filename; out1(lfn+1:lfn+16)='_Prediction1.txt'
  open(1,file=out1(1:lfn+16),status='replace'); temp=-0.05; t=0
  write(1,*) 'Position  P-start,occup  N/L'
  do i=1,nrcv(1)-1; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1)') i,temp,temp,t; end do
  close(1)
  !out2(1:lfn)=filename; out2(lfn+1:lfn+17)='_viterbiPosi1.txt'
  !open(2,file=out2(1:lfn+17),status='replace'); close(2)

  Pd=1.0/real(mlL); ext=7000
  !!!!!!!!!!!!
  do l=1,rep !
  !!!!!!!!!!!!
  tmp=0.0

  !============================
  do x=1,nsec !start of section
  !============================
  zt=nbrk(x)-nrcv(x); base=nrcv(x)-1

  np=1; do while(zt/np>5000000); np=np+1; end do
  do k=1,np-1; zp(k)=(zt/np)*k; end do; zp(np)=zt

  if(np==1) then
    st(1)=1; ed(1)=zt; rts(1)=0; rte(1)=0
  else
    do q=1,np
      if(q==1) then
        st(q)=1; ed(q)=zp(q)+ext; rts(q)=0; rte(q)=ext
      else if(1<q.and.q<np) then
        st(q)=(zp(q-1)+1)-ext; ed(q)=zp(q)+ext; rts(q)=ext; rte(q)=ext
      else if(q==np) then
        st(q)=(zp(q-1)+1)-ext; ed(q)=zp(q); rts(q)=ext; rte(q)=0
      end if
    end do
  end if

  !----------------------------
  do q=1,np !start of partition
  !----------------------------

  z=ed(q)-st(q)+1
  allocate(w(z),c(z),change(z))
  allocate(r(z),hatFN(0:z),hatFL(z))
  allocate(ra(z),hatAN(0:z),hatAL(z),rb(0:(z-1)),hatBN(0:(z-1)),hatBL(0:z))

  w(1:z)=wt(base+st(q):base+ed(q)); do i=1,z; c(i)=5-w(z-i+1); end do
  !!!!!!!!!!!!!!!!!!!!!!!  viterbi algorithm  !!!!!!!!!!!!!!!!!!!!!!!!

  change=0; hatFN(0)=1.0; hatAN(0)=1.0; hatBL(z)=1.0
  !!!!!!!!!!!!!!!! 1 -- 147 !!!!!!!!!!!!!!

  do t=1,min(mlL,147)
    hatFN(t)=0.0; hatFL(t)=1.0
    r(t)=Pd(t)*freqL(w(1))*freqL(c(z-t+1))
    do i=1,t-1
      r(t)=r(t)*tranL(w(i),w(i+1))*tranL(c(z-t+i),c(z-t+i+1))/r(i)
    end do
  end do

  do t=mlL+1,147; hatFN(t)=0.0; hatFL(t)=0.0; r(t)=1.0/16.0; end do

  if(l==rep) then
  !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  ra(1:147)=r(1:147); hatAN(1:147)=hatFN(1:147); hatAL(1:147)=hatFL(1:147)
  !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB 
  end if
  !!!!!!!!!!!!!!!! 148 -- z !!!!!!!!!!!!!!!
      
  do t=148,z

    tempN=0.0
    if(hatFL(t-147)/=0.0.and.t<z) then
      tempN=hatFL(t-147)*freqN(1,w(t-146))*freqN(1,c(z-t+1))
      do i=1,146
        tempN=tempN*tranN(i,w(t-147+i),w(t-146+i))*tranN(i,c(z-t+i),c(z-t+i+1))/r(t-147+i)
      end do
    end if

    tempL=0.0
    do d=1,min(mlL,t)
      if(d==1) then
        tp=freqL(c(z-t+1))
      else
        tp=tranL(w(t-d+1),w(t-d+2))*tranL(c(z-t+d-1),c(z-t+d))/r(t-d+1)*tp
      end if

      if(hatFN(t-d)/=0.0) then
        temp=hatFN(t-d)*Pd(d)*freqL(w(t-d+1))*tp
        if(temp>tempL) then
          tempL=temp; change(t)=t-d
        end if
      end if 
    end do

    r(t)=tempN+tempL
    if(r(t)==0.0) then
      r(t)=1.0/16.0; hatFN(t)=0.0; hatFL(t)=0.0
    else
      hatFN(t)=tempN/r(t); hatFL(t)=1.0-hatFN(t)
    end if

    if(l==rep) then
    !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    tempN=0.0
    if(hatAL(t-147)/=0.0.and.t<z) then
      tempN=hatAL(t-147)*freqN(1,w(t-146))*freqN(1,c(z-t+1))
      do i=1,146
        tempN=tempN*tranN(i,w(t-147+i),w(t-146+i))*tranN(i,c(z-t+i),c(z-t+i+1))/ra(t-147+i)
      end do
    end if

    tempL=0.0
    do d=1,min(mlL,t)
      if(d==1) then
        tp=freqL(c(z-t+1))
      else
        tp=tranL(w(t-d+1),w(t-d+2))*tranL(c(z-t+d-1),c(z-t+d))/ra(t-d+1)*tp
      end if

      if(hatAN(t-d)/=0.0) then
        temp=hatAN(t-d)*Pd(d)*freqL(w(t-d+1))*tp
        tempL=tempL+temp
      end if 
    end do

    ra(t)=tempN+tempL
    if(ra(t)==0.0) then
      ra(t)=1.0/16.0; hatAN(t)=0.0; hatAL(t)=0.0
    else
      hatAN(t)=tempN/ra(t); hatAL(t)=1.0-hatAN(t)
    end if
    !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
    end if

  end do
  !!!!!!!!!!!!!!! backward - backward - backward !!!!!!!!!!!!

  if(l==rep) then
  !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  do t=z-1,max(z-mlL,z-147),-1
    hatBN(t)=1.0; hatBL(t)=0.0
    rb(t)=Pd(z-t)*freqL(w(t+1))*freqL(c(1))
    do i=1,z-t-1
      rb(t)=rb(t)*tranL(w(t+i),w(t+i+1))*tranL(c(i),c(i+1))/rb(t+i)
    end do
  end do

  do t=z-mlL-1,z-147,-1; hatBN(t)=0.0; hatBL(t)=0.0; rb(t)=1.0/16.0; end do

  do t=z-148,0,-1
    tempL=0.0
    if(hatBN(t+147)>0.0.and.t>0) then
      tempL=freqN(1,w(t+1))*freqN(1,c(z-t-146))*hatBN(t+147)
      do i=1,146
        tempL=tempL*tranN(i,w(t+i),w(t+i+1))*tranN(i,c(z-t-147+i),c(z-t-146+i))/rb(t+i)
      end do
    end if

    tempN=0.0
    do d=1,min(mlL,z-t)
      if(d==1) then
        tp=freqL(w(t+1))
      else
        tp=tp*tranL(w(t+d-1),w(t+d))*tranL(c(z-t-d+1),c(z-t-d+2))/rb(t+d-1)
      end if

      if(hatBL(t+d)>0.0) then
        temp=tp*freqL(c(z-t-d+1))*Pd(d)*hatBL(t+d)
        tempN=tempN+temp
      end if 
    end do

    rb(t)=tempN+tempL
    if(rb(t)==0.0) then
      rb(t)=1.0/16.0; hatBN(t)=0.0; hatBL(t)=0.0
    else
      hatBN(t)=tempN/rb(t); hatBL(t)=1.0-hatBN(t)
    end if
  end do
  !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  end if

  !!!!!!!!!!!!!!!!!!!!!! update Pd !!!!!!!!!!!!!!!!!!!!!!!!!!

  if(l<rep) then
    n=z; m=change(z)
    do while(m>0)
      tmp(n-m)=tmp(n-m)+1.0
      n=m-147; m=change(n)
    end do
  end if

  deallocate(r,hatFN,hatFL,hatAL,hatBL); allocate(ppEndN(z+146),ppN(z),Nst(0:z/147))
  if(l==rep) then

    ppEndN=0.0; ppN=0.0; temp=1.0
    do t=z-1,1,-1; temp=temp*rb(t)/ra(t+1); ppEndN(t)=hatAN(t)*hatBN(t)*temp; end do
    do t=1,z; do i=t,t+146; ppN(t)=ppN(t)+ppEndN(i); end do; end do

    n=z; m=change(z); k=0; w=0
    do while(m>0)
      k=k+1; Nst(k)=m-146; w(m-146:m)=1
      n=m-147; m=change(n)
    end do; nvt=k; Nst(0)=z+1

    open(1,file=out1(1:lfn+16),position='append')
    do t=1+rts(q),z-rte(q); write(1,'(i9,2x,f5.3,3x,f5.3,3x,i1)') base+st(q)-1+t,ppEndN(t+146),ppN(t),w(t); end do
    close(1)

    !open(2,file=out2(1:lfn+17),position='append')
    !do i=nvt,1,-1; if(1+rts(q)<=Nst(i).and.Nst(i)<=z-rte(q)) then 
    !  write(2,'(a,3x,i9,8x,f8.5,8x,i4)') trim(seqName),base+st(q)-1+Nst(i)+73,ppEndN(Nst(i)+146),Nst(i-1)-Nst(i)-147
    !end if; end do
    !close(2)

  end if
  deallocate(w,c,change,ra,hatAN,rb,hatBN,ppEndN,ppN,Nst)

  !-----------------------
  end do !end of partition
  !-----------------------

  if(l==rep) then; open(1,file=out1(1:lfn+16),position='append')
    temp=-0.05; t=0
    if(x<nsec) then
      do i=nbrk(x),nrcv(x+1)-1; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1)') i,temp,temp,t; end do
    else if(x==nsec) then
      if(nbrk(nsec)<ztt+1) then
        do i=nbrk(nsec),ztt; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1)') i,temp,temp,t; end do
      end if
    end if
  close(1); end if

  !=====================
  end do !end of section
  !=====================

  if(l<rep) then
    Pd=0.1
    do k=1,mlL; do i=1,mlL
      if(0.5*((k-i)/1.5)**2<20) Pd(k)=Pd(k)+tmp(i)*exp(-0.5*((k-i)/1.5)**2)
    end do; end do
    temp=sum(Pd); Pd=Pd/temp  

    !out2(1:lfn)=filename; out2(lfn+1:lfn+9)='_Lld1.txt'
    !open(2,file=out2(1:lfn+9))
    !do i=1,mlL; write(2,'(f7.5)') Pd(i); end do
    !close(2); !stop
  end if

  !!!!!!!!
  end do !
  !!!!!!!!
  
  deallocate(wt,nrcv,nbrk)

end subroutine vtbfb
!############################ end of subroutine vtbfb #################################################

subroutine vtbfbNL4(lfn,filename,freqL1,tranL1,TtranL2,TtranL3,TtranL4,TfreqN4,TtranN4,mlL,rep,species)

  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer   t,d,i,j,k,l,m,n,z,mlL,rep,lfn,nvt,species
  integer   zt,np,zp(100),q,st(100),ed(100),ext,rts(100),rte(100)
  integer   ztt,nsec,x,base; integer,allocatable:: nrcv(:),nbrk(:)
  real*8    tempN,tempL,temp,tp,Pd(mlL),tmp(mlL),sc(4)
  real*8    TtranL2(16,4),TtranL3(64,4),TtranL4(256,4),TfreqN4(64,4),TtranN4((147-4)*256,4)
  real*8    freqL1(4),tranL1(4,4),tranL2(4,4,4),tranL3(4,4,4,4),tranL4(4,4,4,4,4)
  real*8    freqN4(4,4,4,4),tranN4(5:147,4,4,4,4,4),freqL4(4,4,4,4)

  character*80 seqName,out1,out2
  character(len=lfn) filename
  character(len=10000000),allocatable:: a(:)

  integer*1,allocatable:: wt(:),w(:),c(:)
  integer,allocatable:: change(:),Nst(:); real*8,allocatable:: ppEndN(:),ppN(:)
  real*8,allocatable:: hatFN(:),hatFL(:),r(:),hatAN(:),hatAL(:),ra(:),hatBN(:),hatBL(:),rb(:)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1,4; do j=1,4; tranL2(i,j,:)=TtranL2((i-1)*4+j,:)
  do k=1,4; tranL3(i,j,k,:)=TtranL3((i-1)*16+(j-1)*4+k,:)
  do l=1,4; tranL4(i,j,k,l,:)=TtranL4((i-1)*64+(j-1)*16+(k-1)*4+l,:)
  end do; end do; end do; end do

  do i=1,4; do j=1,4; do k=1,4
    freqN4(i,j,k,:)=TfreqN4((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do n=5,147; do i=1,4; do j=1,4; do k=1,4; do l=1,4
    tranN4(n,i,j,k,l,:)=TtranN4((n-5)*256+(i-1)*64+(j-1)*16+(k-1)*4+l,:)
  end do; end do; end do; end do; end do

  if(species==1) then !Human
    sc(1)=0.9690; sc(2)=1.050
  else if(species==2) then !Mouse
    sc(1)=0.9540; sc(2)=1.074
  else if(species==3) then !Rat
    sc(1)=0.9280; sc(2)=1.116
  else if(species==4) then !Zebrafish
    sc(1)=1.028; sc(2)=0.9543
  else if(species==5) then !D.melanogaster
    sc(1)=0.9313; sc(2)=1.111
  else if(species==6) then !C.elegans
    sc(1)=1.045; sc(2)=0.9269
  else if(species==7) then !S.cerevisiae
    sc(1)=1.0; sc(2)=1.0
  else if(species==8) then !C.albicans
    sc(1)=1.077; sc(2)=0.8757
  else if(species==9) then !S.pombe
    sc(1)=1.033; sc(2)=0.9462
  else if(species==10) then !A.thaliana
    sc(1)=1.038; sc(2)=0.9384
  else if(species==11) then !Maize
    sc(1)=0.8596; sc(2)=1.226
  else
    print*,'unavailable for this species number, quit!'; stop
  end if
  sc(3)=sc(2); sc(4)=sc(1)
  !--------------- scaling ---------------
  do k=1,4
    freqL1(k)=freqL1(k)*sc(k); tranL1(:,k)=tranL1(:,k)*sc(k); tranL2(:,:,k)=tranL2(:,:,k)*sc(k)
    tranL3(:,:,:,k)=tranL3(:,:,:,k)*sc(k); tranL4(:,:,:,:,k)=tranL4(:,:,:,:,k)*sc(k)
    tranN4(:,:,:,:,:,k)=tranN4(:,:,:,:,:,k)*sc(k)
  end do
  do i=1,4; do j=1,4; do k=1,4; do l=1,4
    freqN4(i,j,k,l)=freqN4(i,j,k,l)*sc(i)*sc(j)*sc(k)*sc(l)
  end do; end do; end do; end do
  freqL1=freqL1/sum(freqL1); freqN4=freqN4/sum(freqN4)
  do i=1,4; tranL1(i,:)=tranL1(i,:)/sum(tranL1(i,:))
  do j=1,4; tranL2(i,j,:)=tranL2(i,j,:)/sum(tranL2(i,j,:))
  do k=1,4; tranL3(i,j,k,:)=tranL3(i,j,k,:)/sum(tranL3(i,j,k,:))
  do l=1,4; tranL4(i,j,k,l,:)=tranL4(i,j,k,l,:)/sum(tranL4(i,j,k,l,:))
  do n=5,147; tranN4(n,i,j,k,l,:)=tranN4(n,i,j,k,l,:)/sum(tranN4(n,i,j,k,l,:))
  end do; end do; end do; end do; end do
  !----------------------------------------

  do i=1,4; do j=1,4; do k=1,4; do l=1,4
    freqL4(i,j,k,l)=freqL1(i)*tranL1(i,j)*tranL2(i,j,k)*tranL3(i,j,k,l)
  end do; end do; end do; end do

  allocate(a(1))
  open(1,file=filename)
  read(1,'(1x,a80)') seqName
  read(1,'(a)') a(1); m=len_trim(a(1)); k=1; j=m
  read(1,'(a)',iostat=l) a(1)(1:m)
  do while(l==0)
    k=k+1; j=len_trim(a(1)(1:m))
    read(1,'(a)',iostat=l) a(1)(1:m)
  end do
  ztt=(k-1)*m+j
  close(1)
  if(ztt<148) then
    print*,'The sequence length is shorter than 148bp, quit!'; stop
  else if(148<=ztt.and.ztt<200000.and.rep==2) then
    print*,'Updating linker length distribution is not allowed since the sequence is not long enough.'
    print*,'The linker length distribution is still assumed uniform.'
    rep=1
  end if; allocate(wt(0:ztt+1))

  wt(0)=0; nsec=0
  open(1,file=filename)
  read(1,'(1x,a80)') seqName
  do i=1,k
    read(1,'(a)') a(1)(1:m); n=m; if(i==k) n=j
    do l=1,n; t=(i-1)*m+l
      if(a(1)(l:l)=='A'.or.a(1)(l:l)=='a') then
        wt(t)=1
      else if(a(1)(l:l)=='C'.or.a(1)(l:l)=='c') then
        wt(t)=2
      else if(a(1)(l:l)=='G'.or.a(1)(l:l)=='g') then
        wt(t)=3
      else if(a(1)(l:l)=='T'.or.a(1)(l:l)=='t') then
        wt(t)=4
      else
        wt(t)=0
      end if
      if(wt(t-1)==0.and.wt(t)>0) nsec=nsec+1
    end do
  end do
  close(1); deallocate(a); allocate(nrcv(nsec),nbrk(nsec))

  nsec=0; wt(ztt+1)=0
  do i=1,ztt+1
    if(wt(i-1)==0.and.wt(i)>0) then
      nsec=nsec+1; nrcv(nsec)=i
    else if(wt(i-1)>0.and.wt(i)==0) then
      nbrk(nsec)=i
      if(nbrk(nsec)-nrcv(nsec)<148) nsec=nsec-1
    end if
  end do

  out1(1:lfn)=filename; out1(lfn+1:lfn+16)='_Prediction4.txt'
  open(1,file=out1(1:lfn+16),status='replace'); temp=-0.05; t=0
  write(1,*) 'Position  P-start,occup  N/L'
  do i=1,nrcv(1)-1; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1)') i,temp,temp,t; end do
  close(1)
  !out2(1:lfn)=filename; out2(lfn+1:lfn+17)='_viterbiPosi4.txt'
  !open(2,file=out2(1:lfn+17),status='replace'); close(2)

  Pd=1.0/real(mlL); ext=7000
  !!!!!!!!!!!!
  do l=1,rep !
  !!!!!!!!!!!!
  tmp=0.0

  !============================
  do x=1,nsec !start of section
  !============================
  zt=nbrk(x)-nrcv(x); base=nrcv(x)-1

  np=1; do while(zt/np>5000000); np=np+1; end do
  do k=1,np-1; zp(k)=(zt/np)*k; end do; zp(np)=zt

  if(np==1) then
    st(1)=1; ed(1)=zt; rts(1)=0; rte(1)=0
  else
    do q=1,np
      if(q==1) then
        st(q)=1; ed(q)=zp(q)+ext; rts(q)=0; rte(q)=ext
      else if(1<q.and.q<np) then
        st(q)=(zp(q-1)+1)-ext; ed(q)=zp(q)+ext; rts(q)=ext; rte(q)=ext
      else if(q==np) then
        st(q)=(zp(q-1)+1)-ext; ed(q)=zp(q); rts(q)=ext; rte(q)=0
      end if
    end do
  end if

  !----------------------------
  do q=1,np !start of partition
  !----------------------------

  z=ed(q)-st(q)+1
  allocate(w(z),c(z),change(z))
  allocate(r(z),hatFN(0:z),hatFL(z))
  allocate(ra(z),hatAN(0:z),hatAL(z),rb(0:(z-1)),hatBN(0:(z-1)),hatBL(0:z))

  w(1:z)=wt(base+st(q):base+ed(q)); do i=1,z; c(i)=5-w(z-i+1); end do
  !!!!!!!!!!!!!!!!!!!!!!!  viterbi algorithm  !!!!!!!!!!!!!!!!!!!!!!!!

  change=0; hatFN(0)=1.0; hatAN(0)=1.0; hatBL(z)=1.0
  !!!!!!!!!!!!!! 1 -- 147 !!!!!!!!!!!!!!

  hatFN(1)=0.0; hatFL(1)=1.0
  r(1)=Pd(1)*freqL1(w(1))*freqL1(c(z))

  hatFN(2)=0.0; hatFL(2)=1.0
  r(2)=Pd(2)*freqL1(w(1))*freqL1(c(z-1))*tranL1(w(1),w(2))*tranL1(c(z-1),c(z))/r(1)

  hatFN(3)=0.0; hatFL(3)=1.0
  r(3)=Pd(3)*freqL1(w(1))*freqL1(c(z-2))*tranL1(w(1),w(2))*tranL1(c(z-2),c(z-1))/r(1)&
  &*tranL2(w(1),w(2),w(3))*tranL2(c(z-2),c(z-1),c(z))/r(2)

  hatFN(4)=0.0; hatFL(4)=1.0
  r(4)=Pd(4)*freqL4(w(1),w(2),w(3),w(4))*freqL4(c(z-3),c(z-2),c(z-1),c(z))/r(1)/r(2)/r(3)

  do t=5,min(mlL,147)
    hatFN(t)=0.0; hatFL(t)=1.0

    r(t)=Pd(t)*freqL4(w(1),w(2),w(3),w(4))*freqL4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))/r(1)/r(2)/r(3)
    do i=5,t
      r(t)=r(t)*tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))*tranL4(c(z-t+i-4),c(z-t+i-3),c(z-t+i-2),c(z-t+i-1),c(z-t+i))/r(i-1)
    end do
  end do

  do t=mlL+1,147; hatFN(t)=0.0; hatFL(t)=0.0; r(t)=1.0/16.0; end do

  if(l==rep) then; hatAN(1:147)=hatFN(1:147); hatAL(1:147)=hatFL(1:147); ra(1:147)=r(1:147); end if
  !!!!!!!!!!!!!!!! 148 -- z !!!!!!!!!!!!!
      
  do t=148,z

    tempN=0.0
    if(hatFL(t-147)/=0.0.and.t<z) then
      tempN=hatFL(t-147)*freqN4(w(t-146),w(t-145),w(t-144),w(t-143))*freqN4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))&
      &/r(t-146)/r(t-145)/r(t-144)
      do i=5,147
        tempN=tempN*tranN4(i,w(t-151+i),w(t-150+i),w(t-149+i),w(t-148+i),w(t-147+i))/r(t-148+i)&
        &*tranN4(i,c(z-t+i-4),c(z-t+i-3),c(z-t+i-2),c(z-t+i-1),c(z-t+i))
      end do
    end if

    tempL=0.0
    do d=1,min(mlL,t)
      temp=0.0

      if(d<=3.and.hatFN(t-d)/=0.0) then
        temp=hatFN(t-d)*Pd(d)*freqL1(w(t-d+1))*freqL1(c(z-t+1))
        if(d>=2) temp=temp*tranL1(w(t-d+1),w(t-d+2))*tranL1(c(z-t+1),c(z-t+2))/r(t-d+1)
        if(d>=3) temp=temp*tranL2(w(t-d+1),w(t-d+2),w(t-d+3))*tranL2(c(z-t+1),c(z-t+2),c(z-t+3))/r(t-d+2)	  
      end if

      if(d>=4) then
        if(d==4) then
          tp=freqL4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))
        else if(d>4) then
          tp=tranL4(w(t-d+1),w(t-d+2),w(t-d+3),w(t-d+4),w(t-d+5))*tranL4(c(z-t+d-4),c(z-t+d-3),c(z-t+d-2),c(z-t+d-1),c(z-t+d))&
          &*tp/r(t-d+4)
        end if
        if(hatFN(t-d)/=0.0) temp=hatFN(t-d)*Pd(d)*freqL4(w(t-d+1),w(t-d+2),w(t-d+3),w(t-d+4))*tp/r(t-d+1)/r(t-d+2)/r(t-d+3)
      end if

      if(temp>tempL) then
        tempL=temp; change(t)=t-d
      end if	
    end do

    r(t)=tempN+tempL
    if(r(t)==0.0) then
      r(t)=1.0/16.0; hatFN(t)=0.0; hatFL(t)=0.0
    else
      hatFN(t)=tempN/r(t); hatFL(t)=1.0-hatFN(t)
    end if

    if(l==rep) then
    !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    tempN=0.0
    if(hatAL(t-147)/=0.0.and.t<z) then
      tempN=hatAL(t-147)*freqN4(w(t-146),w(t-145),w(t-144),w(t-143))*freqN4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))&
      &/ra(t-146)/ra(t-145)/ra(t-144)
      do i=5,147
        tempN=tempN*tranN4(i,w(t-151+i),w(t-150+i),w(t-149+i),w(t-148+i),w(t-147+i))/ra(t-148+i)&
        &*tranN4(i,c(z-t+i-4),c(z-t+i-3),c(z-t+i-2),c(z-t+i-1),c(z-t+i))
      end do
    end if

    tempL=0.0
    do d=1,min(mlL,t)
      temp=0.0

      if(d<=3.and.hatAN(t-d)/=0.0) then
        temp=hatAN(t-d)*Pd(d)*freqL1(w(t-d+1))*freqL1(c(z-t+1))
        if(d>=2) temp=temp*tranL1(w(t-d+1),w(t-d+2))*tranL1(c(z-t+1),c(z-t+2))/ra(t-d+1)
        if(d>=3) temp=temp*tranL2(w(t-d+1),w(t-d+2),w(t-d+3))*tranL2(c(z-t+1),c(z-t+2),c(z-t+3))/ra(t-d+2)	  
      end if

      if(d>=4) then
        if(d==4) then
          tp=freqL4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))
        else if(d>4) then
          tp=tranL4(w(t-d+1),w(t-d+2),w(t-d+3),w(t-d+4),w(t-d+5))*tranL4(c(z-t+d-4),c(z-t+d-3),c(z-t+d-2),c(z-t+d-1),c(z-t+d))*tp&
          &/ra(t-d+4)
        end if
        if(hatAN(t-d)/=0.0) temp=hatAN(t-d)*Pd(d)*freqL4(w(t-d+1),w(t-d+2),w(t-d+3),w(t-d+4))*tp/ra(t-d+1)/ra(t-d+2)/ra(t-d+3)
      end if

      tempL=tempL+temp	
    end do

    ra(t)=tempN+tempL
    if(ra(t)==0.0) then
      ra(t)=1.0/16.0; hatAN(t)=0.0; hatAL(t)=0.0
    else
      hatAN(t)=tempN/ra(t); hatAL(t)=1.0-hatAN(t)
    end if
    !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
    end if

  end do
  !!!!!!!!!!!!!!!!!!!!!!! backward !!!!!!!!!!!!!!!!!!!!!!!!!

  if(l==rep) then
  !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  do t=z-1,max(z-mlL,z-147),-1
    hatBN(t)=1.0; hatBL(t)=0.0
    rb(t)=Pd(z-t)*freqL1(w(t+1))*freqL1(c(1))
    if(z-t>=2) rb(t)=rb(t)*tranL1(w(t+1),w(t+2))*tranL1(c(1),c(2))/rb(t+1)
    if(z-t>=3) rb(t)=rb(t)*tranL2(w(t+1),w(t+2),w(t+3))*tranL2(c(1),c(2),c(3))/rb(t+2)	  
    if(z-t>=4) rb(t)=rb(t)*tranL3(w(t+1),w(t+2),w(t+3),w(t+4))*tranL3(c(1),c(2),c(3),c(4))/rb(t+3)
    do i=5,z-t
      rb(t)=rb(t)*tranL4(w(t+i-4),w(t+i-3),w(t+i-2),w(t+i-1),w(t+i))&
      &*tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))/rb(t+i-1)
    end do
  end do

  do t=z-mlL-1,z-147,-1; hatBN(t)=0.0; hatBL(t)=0.0; rb(t)=1.0/16.0; end do

  do t=z-148,0,-1

    tempL=0.0
    if(hatBN(t+147)>0.0.and.t>0) then
      tempL=hatBN(t+147)*freqN4(w(t+1),w(t+2),w(t+3),w(t+4))*freqN4(c(z-t-146),c(z-t-145),c(z-t-144),c(z-t-143))&
      &/rb(t+1)/rb(t+2)/rb(t+3)
      do i=5,147
        tempL=tempL*tranN4(i,w(t+i-4),w(t+i-3),w(t+i-2),w(t+i-1),w(t+i))/rb(t+i-1)&
        &*tranN4(i,c(z-t-151+i),c(z-t-150+i),c(z-t-149+i),c(z-t-148+i),c(z-t-147+i))
      end do
    end if

    tempN=0.0
    do d=1,min(mlL,z-t)
      temp=0.0

      if(d<=3.and.hatBL(t+d)>0.0) then
        temp=Pd(d)*freqL1(w(t+1))*freqL1(c(z-t-d+1))*hatBL(t+d)
        if(d>=2) temp=temp*tranL1(w(t+1),w(t+2))*tranL1(c(z-t-d+1),c(z-t-d+2))/rb(t+1)
        if(d>=3) temp=temp*tranL2(w(t+1),w(t+2),w(t+3))*tranL2(c(z-t-d+1),c(z-t-d+2),c(z-t-d+3))/rb(t+2)	  
      end if

      if(d>=4) then
        if(d==4) then
          tp=freqL4(w(t+1),w(t+2),w(t+3),w(t+4))/rb(t+1)/rb(t+2)/rb(t+3)
        else if(d>4) then
          tp=tranL4(w(t+d-4),w(t+d-3),w(t+d-2),w(t+d-1),w(t+d))*tranL4(c(z-t-d+1),c(z-t-d+2),c(z-t-d+3),c(z-t-d+4),c(z-t-d+5))&
          &*tp/rb(t+d-1)
        end if
        if(hatBL(t+d)>0.0) temp=tp*freqL4(c(z-t-d+1),c(z-t-d+2),c(z-t-d+3),c(z-t-d+4))*Pd(d)*hatBL(t+d)
      end if

      tempN=tempN+temp
    end do

    rb(t)=tempN+tempL
    if(rb(t)==0.0) then
      rb(t)=1.0/16.0; hatBN(t)=0.0; hatBL(t)=0.0
    else
      hatBN(t)=tempN/rb(t); hatBL(t)=1.0-hatBN(t)
    end if

  end do
  !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  end if
  !^^^^^^^^^^^^^^^^^^^ end of backward ^^^^^^^^^^^^^^^^^^^

  if(l<rep) then
    n=z; m=change(z)
    do while(m>0)
      tmp(n-m)=tmp(n-m)+1.0
      n=m-147; m=change(n)
    end do
  end if

  deallocate(r,hatFN,hatFL,hatAL,hatBL); allocate(ppEndN(z+146),ppN(z),Nst(0:z/147))
  if(l==rep) then

    ppEndN=0.0; ppN=0.0; temp=1.0
    do t=z-1,1,-1; temp=temp*rb(t)/ra(t+1); ppEndN(t)=hatAN(t)*hatBN(t)*temp; end do
    do t=1,z; do i=t,t+146; ppN(t)=ppN(t)+ppEndN(i); end do; end do

    n=z; m=change(z); k=0; w=0
    do while(m>0)
      k=k+1; Nst(k)=m-146; w(m-146:m)=1
      n=m-147; m=change(n)
    end do; nvt=k; Nst(0)=z+1

    open(1,file=out1(1:lfn+16),position='append')
    do t=1+rts(q),z-rte(q); write(1,'(i9,2x,f5.3,3x,f5.3,3x,i1)') base+st(q)-1+t,ppEndN(t+146),ppN(t),w(t); end do
    close(1)

    !open(2,file=out2(1:lfn+17),position='append')
    !do i=nvt,1,-1; if(1+rts(q)<=Nst(i).and.Nst(i)<=z-rte(q)) then 
    !  write(2,'(a,3x,i9,8x,f8.5,8x,i4)') trim(seqName),base+st(q)-1+Nst(i)+73,ppEndN(Nst(i)+146),Nst(i-1)-Nst(i)-147
    !end if; end do
    !close(2)

  end if
  deallocate(w,c,change,ra,hatAN,rb,hatBN,ppEndN,ppN,Nst)

  !-----------------------
  end do !end of partition
  !-----------------------

  if(l==rep) then; open(1,file=out1(1:lfn+16),position='append')
    temp=-0.05; t=0
    if(x<nsec) then
      do i=nbrk(x),nrcv(x+1)-1; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1)') i,temp,temp,t; end do
    else if(x==nsec) then
      if(nbrk(nsec)<ztt+1) then
        do i=nbrk(nsec),ztt; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1)') i,temp,temp,t; end do
      end if
    end if
  close(1); end if

  !=====================
  end do !end of section
  !=====================

  if(l<rep) then
    Pd=0.1
    do k=1,mlL; do i=1,mlL
      if(0.5*((k-i)/1.5)**2<20) Pd(k)=Pd(k)+tmp(i)*exp(-0.5*((k-i)/1.5)**2)
    end do; end do
    temp=sum(Pd); Pd=Pd/temp  

    !out2(1:lfn)=filename; out2(lfn+1:lfn+9)='_Lld4.txt'
    !open(2,file=out2(1:lfn+9))
    !do i=1,mlL; write(2,'(f7.5)') Pd(i); end do
    !close(2); !stop
  end if

  !!!!!!!!
  end do !
  !!!!!!!!

  deallocate(wt,nrcv,nbrk)

end subroutine vtbfbNL4
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^ end of program ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
