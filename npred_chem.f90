program main

  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer i,lfn,mlL,rep,species,order, ind; character*80 fileName; character*3 tpc
  real*8  cfreqL1_h(4),ctranL1_h(4,4),ctranL2_h(16,4),ctranL3_h(64,4),ctranL4_h(256,4),Pdd(500,11),Pd(500)
  real*8  cfreqN4_h(64,4),ctranN4_h((147-4)*256,4),cfreqN1_h(147,4),ctranN1_h(584,4)
  real*8  cfreqL1_m(4),ctranL1_m(4,4),ctranL2_m(16,4),ctranL3_m(64,4),ctranL4_m(256,4)
  real*8  cfreqN4_m(64,4),ctranN4_m((147-4)*256,4),cfreqN1_m(147,4),ctranN1_m(584,4)
  real*8  cfreqL1_p(4),ctranL1_p(4,4),ctranL2_p(16,4),ctranL3_p(64,4),ctranL4_p(256,4)
  real*8  cfreqN4_p(64,4),ctranN4_p((147-4)*256,4),cfreqN1_p(147,4),ctranN1_p(584,4)
  real*8  cfreqL1_y(4),ctranL1_y(4,4),ctranL2_y(16,4),ctranL3_y(64,4),ctranL4_y(256,4)
  real*8  cfreqN4_y(64,4),ctranN4_y((147-4)*256,4),cfreqN1_y(147,4),ctranN1_y(584,4)
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ind=0
  
  !#human
  open(1,file='./h/cfreqL_h.txt')  !##### freqL1
  read(1,*) cfreqL1_h; close(1)
  open(1,file='./h/ctranL_h.txt')  !##### tranL1
  do i=1,4; read(1,*) ctranL1_h(i,:); end do; close(1)
  open(1,file='./h/ctranL2_h.txt') !##### tranL2
  do i=1,16; read(1,*) ctranL2_h(i,:); end do; close(1)
  open(1,file='./h/ctranL3_h.txt') !##### tranL3
  do i=1,64; read(1,*) ctranL3_h(i,:); end do; close(1)
  open(1,file='./h/ctranL4_h.txt') !##### tranL4
  do i=1,256; read(1,*) ctranL4_h(i,:); end do; close(1)

  open(1,file='./h/cfreqN_h.txt') !freqN1
  do i=1,147; read(1,*) cfreqN1_h(i,:); end do
  close(1)
  open(2,file='./h/ctranN_h.txt') !tranN1
  do i=1,584; read(2,*) ctranN1_h(i,:); end do
  close(2)

  open(1,file='./h/cfreqN4_h.txt') !##### freqN4
  do i=1,64; read(1,*) cfreqN4_h(i,:); end do; close(1)
  open(1,file='./h/ctranN4_h.txt') !##### tranN4
  do i=1,(147-4)*256; read(1,*) ctranN4_h(i,:); end do; close(1)
  
  !#mouse
  
  open(1,file='./m/cfreqL_m.txt')  !##### freqL1
  read(1,*) cfreqL1_m; close(1)
  open(1,file='./m/ctranL_m.txt')  !##### tranL1
  do i=1,4; read(1,*) ctranL1_m(i,:); end do; close(1)
  open(1,file='./m/ctranL2_m.txt') !##### tranL2
  do i=1,16; read(1,*) ctranL2_m(i,:); end do; close(1)
  open(1,file='./m/ctranL3_m.txt') !##### tranL3
  do i=1,64; read(1,*) ctranL3_m(i,:); end do; close(1)
  open(1,file='./m/ctranL4_m.txt') !##### tranL4
  do i=1,256; read(1,*) ctranL4_m(i,:); end do; close(1)

  open(1,file='./m/cfreqN_m.txt') !freqN1
  do i=1,147; read(1,*) cfreqN1_m(i,:); end do
  close(1)
  open(2,file='./m/ctranN_m.txt') !tranN1
  do i=1,584; read(2,*) ctranN1_m(i,:); end do
  close(2)

  open(1,file='./m/cfreqN4_m.txt') !##### freqN4
  do i=1,64; read(1,*) cfreqN4_m(i,:); end do; close(1)
  open(1,file='./m/ctranN4_m.txt') !##### tranN4
  do i=1,(147-4)*256; read(1,*) ctranN4_m(i,:); end do; close(1)
  
  !#pombe
  
  open(1,file='./p/cfreqL_p.txt')  !##### freqL1
  read(1,*) cfreqL1_p; close(1)
  open(1,file='./p/ctranL_p.txt')  !##### tranL1
  do i=1,4; read(1,*) ctranL1_p(i,:); end do; close(1)
  open(1,file='./p/ctranL2_p.txt') !##### tranL2
  do i=1,16; read(1,*) ctranL2_p(i,:); end do; close(1)
  open(1,file='./p/ctranL3_p.txt') !##### tranL3
  do i=1,64; read(1,*) ctranL3_p(i,:); end do; close(1)
  open(1,file='./p/ctranL4_p.txt') !##### tranL4
  do i=1,256; read(1,*) ctranL4_p(i,:); end do; close(1)

  open(1,file='./p/cfreqN_p.txt') !freqN1
  do i=1,147; read(1,*) cfreqN1_p(i,:); end do
  close(1)
  open(2,file='./p/ctranN_p.txt') !tranN1
  do i=1,584; read(2,*) ctranN1_p(i,:); end do
  close(2)

  open(1,file='./p/cfreqN4_p.txt') !##### freqN4
  do i=1,64; read(1,*) cfreqN4_p(i,:); end do; close(1)
  open(1,file='./p/ctranN4_p.txt') !##### tranN4
  do i=1,(147-4)*256; read(1,*) ctranN4_p(i,:); end do; close(1)
  

  !#yeast
  
  open(1,file='./y/cfreqL_y.txt')  !##### freqL1
  read(1,*) cfreqL1_y; close(1)
  open(1,file='./y/ctranL_y.txt')  !##### tranL1
  do i=1,4; read(1,*) ctranL1_y(i,:); end do; close(1)
  open(1,file='./y/ctranL2_y.txt') !##### tranL2
  do i=1,16; read(1,*) ctranL2_y(i,:); end do; close(1)
  open(1,file='./y/ctranL3_y.txt') !##### tranL3
  do i=1,64; read(1,*) ctranL3_y(i,:); end do; close(1)
  open(1,file='./y/ctranL4_y.txt') !##### tranL4
  do i=1,256; read(1,*) ctranL4_y(i,:); end do; close(1)

  open(1,file='./y/cfreqN_y.txt') !freqN1
  do i=1,147; read(1,*) cfreqN1_y(i,:); end do
  close(1)
  open(2,file='./y/ctranN_y.txt') !tranN1
  do i=1,584; read(2,*) ctranN1_y(i,:); end do
  close(2)

  open(1,file='./y/cfreqN4_y.txt') !##### freqN4
  do i=1,64; read(1,*) cfreqN4_y(i,:); end do; close(1)
  open(1,file='./y/ctranN4_y.txt') !##### tranN4
  do i=1,(147-4)*256; read(1,*) ctranN4_y(i,:); end do; close(1)
  
  
  !##input linker length distribution

  open(1,file='Pdd.txt') !##### freqN4
  do i=1,500; read(1,*) Pdd(i,:); end do; close(1)
  

  write(*,'(a)',advance='no') 'Input the fileName of sequence                  :  '; read*,fileName
  !maximum linker length prespecified =500
  !!write(*,'(a)',advance='no') 'Input the maximum length of linker              :  '; read*,mlL
  write(*,'(a)',advance='no') 'Input the order of selected model (1 or 4)      :  '; read*,order
  if(order/=1.and.order/=4) then; print*,'1 or 4 should be inputed!'; stop; end if
  
  !###### no linker update, input linker length distribution, rep=1
  !write(*,'(a)',advance='no') 'Update the linker length distribution? yes or no:  '; read*,tpc
  !if(tpc=='yes') then; rep=2; else if(tpc=='no ') then; rep=1; else; print*,'yes or no should be inputed!'; stop; end if
  
  rep=1; mlL=500
  
  print*,' '
  print*,'Select the species from the following list:' 
  print*,'1=Human;         2=Mouse;           3=Rat;'
  print*,'4=Zebrafish;     5=D. melanogaster; 6=C. elegans;'
  print*,'7=S. cerevisiae; 8=C. albicans;     9=S. pombe;'
  print*,'10=A. thaliana;  11=Maize.'
  write(*,'(a)',advance='no') 'Input the number of selected species            :  '; read*,species
  print*,'Calculating......'

  !## input linker length distribtuion
  lfn=len_trim(fileName)
  Pd(1:mlL)=Pdd(1:mlL,species)
  
  if(order==1) then
    if(species==1) then
        call cvtbfb(lfn,trim(fileName),cfreqL1_h,ctranL1_h,cfreqN1_h,ctranN1_h,mlL,rep,species,Pdd,ind)
    else if (species==2) then
        call cvtbfb(lfn,trim(fileName),cfreqL1_m,ctranL1_m,cfreqN1_m,ctranN1_m,mlL,rep,species,Pdd,ind)
    else if (species==7) then
        call cvtbfb(lfn,trim(fileName),cfreqL1_y,ctranL1_y,cfreqN1_y,ctranN1_y,mlL,rep,species,Pdd,ind)
    else if (species==9) then
        call cvtbfb(lfn,trim(fileName),cfreqL1_p,ctranL1_p,cfreqN1_h,ctranN1_p,mlL,rep,species,Pdd,ind)
    else 
        call cvtbfb(lfn,trim(fileName),cfreqL1_y,ctranL1_y,cfreqN1_y,ctranN1_y,mlL,rep,species,Pdd,ind)
    end if
  else if(order==4) then
    if(species==1) then
    call cvtbfbNL4(lfn,trim(fileName),cfreqL1_h,ctranL1_h,ctranL2_h,ctranL3_h,ctranL4_h,cfreqN4_h,ctranN4_h,mlL,rep,species,Pdd,ind)
    else if (species==2) then
    call cvtbfbNL4(lfn,trim(fileName),cfreqL1_m,ctranL1_m,ctranL2_m,ctranL3_m,ctranL4_m,cfreqN4_m,ctranN4_m,mlL,rep,species,Pdd,ind)
    else if (species==7) then
    call cvtbfbNL4(lfn,trim(fileName),cfreqL1_y,ctranL1_y,ctranL2_y,ctranL3_y,ctranL4_y,cfreqN4_y,ctranN4_y,mlL,rep,species,Pdd,ind)
    else if (species==9) then
    call cvtbfbNL4(lfn,trim(fileName),cfreqL1_p,ctranL1_p,ctranL2_p,ctranL3_p,ctranL4_p,cfreqN4_p,ctranN4_p,mlL,rep,species,Pdd,ind)
    else 
    call cvtbfbNL4(lfn,trim(fileName),cfreqL1_y,ctranL1_y,ctranL2_y,ctranL3_y,ctranL4_y,cfreqN4_y,ctranN4_y,mlL,rep,species,Pdd,ind)
    end if
  end if

end


!!######### chemical map-trained prediction

subroutine cvtbfb(lfn,filename,freqL,tranL,freqN,tranNN,mlL,rep,species,Pdd,ind)

  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer   t,d,i,j,k,l,m,n,z,mlL,lfn,rep,nvt,species,ind,tep,lfn1
  integer   zt,np,zp(100),q,st(100),ed(100),ext,rts(100),rte(100),file_name_num(lfn)
  integer   ztt,nsec,x,base; integer,allocatable:: nrcv(:),nbrk(:)
  real*8    tempN,tempL,temp,tp,sc(4),scc(2,0:11)
  real*8    freqL(4),tranL(4,4),freqN(147,4),tranN(146,4,4),tranNN(584,4),Pd(mlL),tmp(mlL),Pdd(mlL,11)

  character*80 seqName,out1; character*2 tpc
  character(len=lfn) filename
  character(len=10000000),allocatable:: a(:)

  integer*1,allocatable:: wt(:),w(:),c(:)
  integer,allocatable:: change(:),Nst(:)
  real*8,allocatable:: r(:),hatFN(:),hatFL(:),ra(:),hatAN(:),hatAL(:),rb(:),hatBN(:),hatBL(:)
  real*8,allocatable:: ppEndN(:),ppN(:),asc(:),aas(:)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do i=1,146; do j=1,4; tranN(i,j,:)=tranNN(4*(i-1)+j,:); end do; end do

  allocate(a(1))
  l=1; open(1,file=filename,status='old',err=6); l=0
6 if(l==1) then
    ind=1; return
  end if

  read(1,'(1x,a80)',iostat=l) seqName
  if(l/=0) then
    ind=2; return
  end if
  read(1,'(a)',iostat=l) a(1)
  if(l/=0) then
    ind=2; return
  end if
  m=len_trim(a(1)); k=1; j=m
  read(1,'(a)',iostat=l) a(1)(1:m)
  do while(l==0)
    k=k+1; j=len_trim(a(1)(1:m))
    read(1,'(a)',iostat=l) a(1)(1:m)
  end do
  ztt=(k-1)*m+j
  close(1)
  if(ztt<148) then
    ind=2; return
  else if(148<=ztt.and.ztt<200000.and.rep>=2) then
    print*,'Updating linker length distribution is not allowed since the sequence is not long enough.'
    print*,'The linker length distribution is still assumed uniform.'
    rep=1
  end if; allocate(wt(0:ztt+1))

  wt(0)=0; nsec=0; scc(1:2,0)=0.0
  open(1,file=filename)
  read(1,'(1x,a80)') seqName
  do i=1,k
    read(1,'(a)') a(1)(1:m); n=m; if(i==k) n=j
    do l=1,n; t=(i-1)*m+l
      if(a(1)(l:l)=='A'.or.a(1)(l:l)=='a') then
        wt(t)=1; scc(1,0)=scc(1,0)+1.0
      else if(a(1)(l:l)=='C'.or.a(1)(l:l)=='c') then
        wt(t)=2; scc(2,0)=scc(2,0)+1.0
      else if(a(1)(l:l)=='G'.or.a(1)(l:l)=='g') then
        wt(t)=3; scc(2,0)=scc(2,0)+1.0
      else if(a(1)(l:l)=='T'.or.a(1)(l:l)=='t') then
        wt(t)=4; scc(1,0)=scc(1,0)+1.0
      else if(a(1)(l:l)=='N'.or.a(1)(l:l)=='n') then
        wt(t)=0
      else
        ind=2; return
      end if
      if(wt(t-1)==0.and.wt(t)>0) nsec=nsec+1
    end do
  end do
  close(1); deallocate(a); allocate(nrcv(nsec),nbrk(nsec))
  temp=sum(scc(1:2,0)); scc(1,0)=scc(1,0)/temp/0.617; scc(2,0)=scc(2,0)/temp/0.383

  nsec=0; wt(ztt+1)=0
  do i=1,ztt+1
    if(wt(i-1)==0.and.wt(i)>0) then
      nsec=nsec+1; nrcv(nsec)=i
    else if(wt(i-1)>0.and.wt(i)==0) then
      nbrk(nsec)=i
      if(nbrk(nsec)-nrcv(nsec)<148) nsec=nsec-1
    end if
  end do

  do i=lfn,1,-1; if(filename(i:i)=='/'.or.filename(i:i)=='\') exit; end do
  lfn1=lfn-i
  out1(1:lfn1)=filename(i+1:lfn); out1(lfn1+1:lfn1+16)='_P1_chem.txt'
  open(1,file=out1(1:lfn1+16),status='replace'); temp=-0.05; t=0; tpc='NA'
  write(1,*) 'Position  P-start Occup  N/L  Affinity'
  do i=1,nrcv(1)-1; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1,4x,a)') i,temp,temp,t,tpc; end do
  close(1)
  !out2(1:lfn)=filename; out2(lfn+1:lfn+17)='_viterbiPosi1.txt'
  !open(2,file=out2(1:lfn+17),status='replace'); close(2)

  !species==1 !Human
    scc(1,1)=1.0; scc(2,1)=1.0
  !species==2 !Mouse
    scc(1,2)=1.0; scc(2,2)=1.0
  !species==3 !Rat
    scc(1,3)=0.9280; scc(2,3)=1.116
  !species==4 !Zebrafish
    scc(1,4)=1.028; scc(2,4)=0.9543
  !species==5 !D.melanogaster
    scc(1,5)=0.9313; scc(2,5)=1.111
  !species==6 !C.elegans
    scc(1,6)=1.045; scc(2,6)=0.9269
  !species==7 !S.cerevisiae
    scc(1,7)=1.0; scc(2,7)=1.0
  !species==8 !C.albicans
    scc(1,8)=1.077; scc(2,8)=0.8757
  !species==9 !S.pombe
    scc(1,9)=1.0; scc(2,9)=1.0
  !species==10 !A.thaliana
    scc(1,10)=1.038; scc(2,10)=0.9384
  !species==11 !Maize
    scc(1,11)=0.8596; scc(2,11)=1.226
  if(species<0.or.species>11) then
    ind=3; return
  end if
  sc(1)=scc(1,species); sc(2)=scc(2,species)
  sc(3)=sc(2); sc(4)=sc(1)
  if(species==0) then
    tep=1; temp=abs(scc(1,1)-scc(1,0))
    do i=2,11; tp=abs(scc(1,i)-scc(1,0))
      if(tp<temp) then; tep=i; temp=tp; end if
    end do
  else
    tep=species
  end if; Pd(1:mlL)=Pdd(1:mlL,tep)
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

  ext=7000
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

  w(1:z)=wt(base+st(q):base+ed(q)); do i=1,z; c(i)=5_1-w(z-i+1); end do
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
      hatFN(t)=tempN/r(t); hatFL(t)=tempL/r(t)
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
      hatAN(t)=tempN/ra(t); hatAL(t)=tempL/ra(t)
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
      hatBN(t)=tempN/rb(t); hatBL(t)=tempL/rb(t)
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

  deallocate(r,hatFN,hatFL,hatAL,hatBL); allocate(ppEndN(z+146),ppN(z),asc(74:z-73),aas(74:z-73),Nst(0:z/147))
  if(l==rep) then

    ppEndN=0.0; ppN=0.0; temp=1.0
    do t=z-1,1,-1; temp=temp*rb(t)/ra(t+1); ppEndN(t)=hatAN(t)*hatBN(t)*temp; end do
    do t=1,z; do i=t,t+146; ppN(t)=ppN(t)+ppEndN(i); end do; end do

    do t=74,z-73
      asc(t)=freqN(1,w(t-73))/freqL(w(t-73))*freqN(1,c(z-t-72))/freqL(c(z-t-72))
      do i=1,146
        asc(t)=asc(t)*tranN(i,w(t-74+i),w(t-73+i))/tranL(w(t-74+i),w(t-73+i))&
        &*tranN(i,c(z-t-73+i),c(z-t-72+i))/tranL(c(z-t-73+i),c(z-t-72+i))
      end do
      asc(t)=log(asc(t))
    end do; i=27
    if(z-73-i>=74+i+1) then
      do t=74,74+i; aas(t)=sum(asc(74:t+i))/(t+i-73.0); end do
      do t=74+i+1,z-73-i; aas(t)=aas(t-1)+(asc(t+i)-asc(t-1-i))/(2*i+1.0); end do
      do t=z-73-i+1,z-73; aas(t)=sum(asc(t-i:z-73))/(z-72.0-t+i); end do
    else
      aas=asc
    end if
    temp=sum(aas(:))/real(z-146); tp=sqrt(sum((aas(:)-temp)**2)/(z-147.0))
    aas=(aas-temp)/tp

    n=z; m=change(z); k=0; w=0
    do while(m>0)
      k=k+1; Nst(k)=m-146; w(m-146:m)=1
      n=m-147; m=change(n)
    end do; nvt=k; Nst(0)=z+1

    open(1,file=out1(1:lfn1+16),position='append')
    do t=1+rts(q),z-rte(q)
      if(t<74.or.z-73<t) then
        write(1,'(i9,2x,f5.3,3x,f5.3,3x,i1,4x,a)') base+st(q)-1+t,ppEndN(t+146),ppN(t),w(t),tpc
      else
        write(1,'(i9,2x,f5.3,3x,f5.3,3x,i1,1x,f8.3)') base+st(q)-1+t,ppEndN(t+146),ppN(t),w(t),aas(t)
      end if
    end do
    close(1)

    !open(2,file=out2(1:lfn+17),position='append')
    !do i=nvt,1,-1; if(1+rts(q)<=Nst(i).and.Nst(i)<=z-rte(q)) then 
    !  write(2,'(a,3x,i9,8x,f8.5,8x,i4)') trim(seqName),base+st(q)-1+Nst(i)+73,ppEndN(Nst(i)+146),Nst(i-1)-Nst(i)-147
    !end if; end do
    !close(2)

  end if
  deallocate(w,c,change,ra,hatAN,rb,hatBN,ppEndN,ppN,asc,aas,Nst)

  !-----------------------
  end do !end of partition
  !-----------------------

  if(l==rep) then; open(1,file=out1(1:lfn1+16),position='append')
    temp=-0.05; t=0
    if(x<nsec) then
      do i=nbrk(x),nrcv(x+1)-1; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1,4x,a)') i,temp,temp,t,tpc; end do
    else if(x==nsec) then
      if(nbrk(nsec)<ztt+1) then
        do i=nbrk(nsec),ztt; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1,4x,a)') i,temp,temp,t,tpc; end do
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

end subroutine cvtbfb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!###4th order chemical map-trained models

subroutine cvtbfbNL4(lfn,filename,freqL1,tranL1,TtranL2,TtranL3,TtranL4,TfreqN4,TtranN4,mlL,rep,species,Pdd,ind)

  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer   t,d,i,j,k,l,m,n,z,mlL,rep,lfn,nvt,species,ind,tep,lfn1
  integer   zt,np,zp(100),q,st(100),ed(100),ext,rts(100),rte(100),file_name_num(lfn)
  integer   ztt,nsec,x,base; integer,allocatable:: nrcv(:),nbrk(:)
  real*8    tempN,tempL,temp,tp,Pd(mlL),tmp(mlL),sc(4),Pdd(mlL,11),scc(2,0:11)
  real*8    TtranL2(16,4),TtranL3(64,4),TtranL4(256,4),TfreqN4(64,4),TtranN4((147-4)*256,4)
  real*8    freqL1(4),tranL1(4,4),tranL2(4,4,4),tranL3(4,4,4,4),tranL4(4,4,4,4,4)
  real*8    freqN4(4,4,4,4),tranN4(5:147,4,4,4,4,4),freqL4(4,4,4,4)

  character*80 seqName,out1; character*2 tpc
  character(len=lfn) filename
  character(len=10000000),allocatable:: a(:)

  integer*1,allocatable:: wt(:),w(:),c(:)
  integer,allocatable:: change(:),Nst(:); real*8,allocatable:: ppEndN(:),ppN(:),asc(:),aas(:)
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

  allocate(a(1))
  l=1; open(1,file=filename,status='old',err=6); l=0
6 if(l==1) then
    ind=1; return
  end if

  read(1,'(1x,a80)',iostat=l) seqName
  if(l/=0) then
    ind=2; return
  end if
  read(1,'(a)',iostat=l) a(1)
  if(l/=0) then
    ind=2; return
  end if
  m=len_trim(a(1)); k=1; j=m
  read(1,'(a)',iostat=l) a(1)(1:m)
  do while(l==0)
    k=k+1; j=len_trim(a(1)(1:m))
    read(1,'(a)',iostat=l) a(1)(1:m)
  end do
  ztt=(k-1)*m+j
  close(1)
  if(ztt<148) then
    ind=2; return
  else if(148<=ztt.and.ztt<200000.and.rep>=2) then
    print*,'Updating linker length distribution is not allowed since the sequence is not long enough.'
    print*,'The linker length distribution is still assumed uniform.'
    rep=1
  end if; allocate(wt(0:ztt+1))

  wt(0)=0; nsec=0; scc(1:2,0)=0.0
  open(1,file=filename)
  read(1,'(1x,a80)') seqName
  do i=1,k
    read(1,'(a)') a(1)(1:m); n=m; if(i==k) n=j
    do l=1,n; t=(i-1)*m+l
      if(a(1)(l:l)=='A'.or.a(1)(l:l)=='a') then
        wt(t)=1; scc(1,0)=scc(1,0)+1.0
      else if(a(1)(l:l)=='C'.or.a(1)(l:l)=='c') then
        wt(t)=2; scc(2,0)=scc(2,0)+1.0
      else if(a(1)(l:l)=='G'.or.a(1)(l:l)=='g') then
        wt(t)=3; scc(2,0)=scc(2,0)+1.0
      else if(a(1)(l:l)=='T'.or.a(1)(l:l)=='t') then
        wt(t)=4; scc(1,0)=scc(1,0)+1.0
      else if(a(1)(l:l)=='N'.or.a(1)(l:l)=='n') then
        wt(t)=0
      else
        ind=2; return
      end if
      if(wt(t-1)==0.and.wt(t)>0) nsec=nsec+1
    end do
  end do
  close(1); deallocate(a); allocate(nrcv(nsec),nbrk(nsec))
  temp=sum(scc(1:2,0)); scc(1,0)=scc(1,0)/temp/0.617; scc(2,0)=scc(2,0)/temp/0.383

  nsec=0; wt(ztt+1)=0
  do i=1,ztt+1
    if(wt(i-1)==0.and.wt(i)>0) then
      nsec=nsec+1; nrcv(nsec)=i
    else if(wt(i-1)>0.and.wt(i)==0) then
      nbrk(nsec)=i
      if(nbrk(nsec)-nrcv(nsec)<148) nsec=nsec-1
    end if
  end do

  do i=lfn,1,-1; if(filename(i:i)=='/'.or.filename(i:i)=='\') exit; end do
  lfn1=lfn-i
  out1(1:lfn1)=filename(i+1:lfn); out1(lfn1+1:lfn1+16)='_P4_chem.txt'
  open(1,file=out1(1:lfn1+16),status='replace'); temp=-0.05; t=0; tpc='NA'
  write(1,*) 'Position  P-start Occup  N/L  Affinity'
  do i=1,nrcv(1)-1; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1,4x,a)') i,temp,temp,t,tpc; end do
  close(1)
  !out2(1:lfn)=filename; out2(lfn+1:lfn+17)='_viterbiPosi4.txt'
  !open(2,file=out2(1:lfn+17),status='replace'); close(2)

  !species==1 !Human
    scc(1,1)=1.0; scc(2,1)=1.0
  !species==2 !Mouse
    scc(1,2)=1.0; scc(2,2)=1.0
  !species==3 !Rat
    scc(1,3)=0.9280; scc(2,3)=1.116
  !species==4 !Zebrafish
    scc(1,4)=1.028; scc(2,4)=0.9543
  !species==5 !D.melanogaster
    scc(1,5)=0.9313; scc(2,5)=1.111
  !species==6 !C.elegans
    scc(1,6)=1.045; scc(2,6)=0.9269
  !species==7 !S.cerevisiae
    scc(1,7)=1.0; scc(2,7)=1.0
  !species==8 !C.albicans
    scc(1,8)=1.077; scc(2,8)=0.8757
  !species==9 !S.pombe
    scc(1,9)=1.0; scc(2,9)=1.0
  !species==10 !A.thaliana
    scc(1,10)=1.038; scc(2,10)=0.9384
  !species==11 !Maize
    scc(1,11)=0.8596; scc(2,11)=1.226
  if(species<0.or.species>11) then
    ind=3; return
  end if
  sc(1)=scc(1,species); sc(2)=scc(2,species)
  sc(3)=sc(2); sc(4)=sc(1)
  if(species==0) then
    tep=1; temp=abs(scc(1,1)-scc(1,0))
    do i=2,11; tp=abs(scc(1,i)-scc(1,0))
      if(tp<temp) then; tep=i; temp=tp; end if
    end do
  else
    tep=species
  end if; Pd(1:mlL)=Pdd(1:mlL,tep)
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

  ext=7000
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

  w(1:z)=wt(base+st(q):base+ed(q)); do i=1,z; c(i)=5_1-w(z-i+1); end do
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
      hatFN(t)=tempN/r(t); hatFL(t)=tempL/r(t)
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
      hatAN(t)=tempN/ra(t); hatAL(t)=tempL/ra(t)
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
      hatBN(t)=tempN/rb(t); hatBL(t)=tempL/rb(t)
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

  deallocate(r,hatFN,hatFL,hatAL,hatBL); allocate(ppEndN(z+146),ppN(z),asc(74:z-73),aas(74:z-73),Nst(0:z/147))
  if(l==rep) then

    ppEndN=0.0; ppN=0.0; temp=1.0
    do t=z-1,1,-1; temp=temp*rb(t)/ra(t+1); ppEndN(t)=hatAN(t)*hatBN(t)*temp; end do
    do t=1,z; do i=t,t+146; ppN(t)=ppN(t)+ppEndN(i); end do; end do

    do t=74,z-73
      asc(t)=freqN4(w(t-73),w(t-72),w(t-71),w(t-70))/freqL4(w(t-73),w(t-72),w(t-71),w(t-70))&
      &*freqN4(c(z-t-72),c(z-t-71),c(z-t-70),c(z-t-69))/freqL4(c(z-t-72),c(z-t-71),c(z-t-70),c(z-t-69))
      do i=5,147
        asc(t)=asc(t)*tranN4(i,w(t-78+i),w(t-77+i),w(t-76+i),w(t-75+i),w(t-74+i))&
                    &/tranL4(w(t-78+i),w(t-77+i),w(t-76+i),w(t-75+i),w(t-74+i))&
                    &*tranN4(i,c(z-t-77+i),c(z-t-76+i),c(z-t-75+i),c(z-t-74+i),c(z-t-73+i))&
                    &/tranL4(c(z-t-77+i),c(z-t-76+i),c(z-t-75+i),c(z-t-74+i),c(z-t-73+i))
      end do
      asc(t)=log(asc(t))
    end do; i=27
    if(z-73-i>=74+i+1) then
      do t=74,74+i; aas(t)=sum(asc(74:t+i))/(t+i-73.0); end do
      do t=74+i+1,z-73-i; aas(t)=aas(t-1)+(asc(t+i)-asc(t-1-i))/(2*i+1.0); end do
      do t=z-73-i+1,z-73; aas(t)=sum(asc(t-i:z-73))/(z-72.0-t+i); end do
    else
      aas=asc
    end if
    temp=sum(aas(:))/real(z-146); tp=sqrt(sum((aas(:)-temp)**2)/(z-147.0))
    aas=(aas-temp)/tp

    n=z; m=change(z); k=0; w=0
    do while(m>0)
      k=k+1; Nst(k)=m-146; w(m-146:m)=1
      n=m-147; m=change(n)
    end do; nvt=k; Nst(0)=z+1

    open(1,file=out1(1:lfn1+16),position='append')
    do t=1+rts(q),z-rte(q)
      if(t<74.or.z-73<t) then
        write(1,'(i9,2x,f5.3,3x,f5.3,3x,i1,4x,a)') base+st(q)-1+t,ppEndN(t+146),ppN(t),w(t),tpc
      else
        write(1,'(i9,2x,f5.3,3x,f5.3,3x,i1,1x,f8.3)') base+st(q)-1+t,ppEndN(t+146),ppN(t),w(t),aas(t)
      end if
    end do
    close(1)

    !open(2,file=out2(1:lfn+17),position='append')
    !do i=nvt,1,-1; if(1+rts(q)<=Nst(i).and.Nst(i)<=z-rte(q)) then 
    !  write(2,'(a,3x,i9,8x,f8.5,8x,i4)') trim(seqName),base+st(q)-1+Nst(i)+73,ppEndN(Nst(i)+146),Nst(i-1)-Nst(i)-147
    !end if; end do
    !close(2)

  end if
  deallocate(w,c,change,ra,hatAN,rb,hatBN,ppEndN,ppN,asc,aas,Nst)

  !-----------------------
  end do !end of partition
  !-----------------------

  if(l==rep) then; open(1,file=out1(1:lfn1+16),position='append')
    temp=-0.05; t=0
    if(x<nsec) then
      do i=nbrk(x),nrcv(x+1)-1; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1,4x,a)') i,temp,temp,t,tpc; end do
    else if(x==nsec) then
      if(nbrk(nsec)<ztt+1) then
        do i=nbrk(nsec),ztt; write(1,'(i9,2x,f5.2,3x,f5.2,3x,i1,4x,a)') i,temp,temp,t,tpc; end do
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

end subroutine cvtbfbNL4

