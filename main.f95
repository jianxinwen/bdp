

program main
    use const
    use formula
    use lsm
    implicit none
    !=====================
    !�����������
    !=====================
    real,allocatable::container_x(:),container_y(:)  !!�ع�����(xi,yi)
    real::lsm_k,lsm_b,vielle_b,vielle_n
    real::min_p,max_p,step_p   !!ָ��ѹǿ���Է�Χ�Ͳ���
    integer::idx=1      !!idx ����ֵ
    real::i,j

    real::DELTA=0.01    !!��⾫��
    real::T
    logical::found
    real::RHO_M         !!�ƽ����ܶ�(kg/m^3)




    RHO_M=CONST_RHO_ox*MUTABLE_ALPHA+CONST_RHO_f*(1-MUTABLE_ALPHA)
    found=.false.

    !=====================
    !���άҲ�﹫ʽ r=b P^n
    !=====================
    write(*,"('��ʼ���άҲ��ȼ�ٹ�ʽ')")
    write(*,"('���ڼ��㡤����')")
    min_p=50.0;max_p=150.0;step_p=10.0
    allocate(container_x(int((max_p-min_p)/step_p)+1))
    allocate(container_y(int((max_p-min_p)/step_p)+1))
    do j=min_p,max_p,step_p
        MUTABLE_P=j
        container_x(idx)=log(j)  !!ln(r)
        do i=800.0,1500.0,0.001
            if(abs(T_s(i)-i)<=DELTA)then
                T=T_s(i)
                container_y(idx)=log(m(t)*1E3/RHO_M)  !!ln(p)
                exit
            end if
        end do
!        write(*,"('x',I3,'=',F6.2,'   ','y',I3,'=',F10.5,'')") idx,container_x(idx),idx,container_y(idx)
        idx=idx+1
    end do
    lsm_k=get_k(container_x,container_y,size(container_x))
    lsm_b=get_b(container_x,container_y,size(container_x))
    vielle_b=exp(lsm_b)
    vielle_n=lsm_k
    write(*,"('άҲ��ȼ�ٹ�ʽ��r=',F8.5,'*P^',F7.5)")vielle_b,vielle_n
    write(*,*)
    write(*,"('==================================')")
    write(*,*)
!    write(*,"('��ʼ����ض������µ����ݣ�����������','P=',F8.3,'atm,',' ��=',F6.3)") MUTABLE_P,MUTABLE_ALPHA

    !=====================
    !�ض������µ�����
    !=====================
    write(*,"('��ʼ����ض������µ�����')")
    write(*,"('������ѹǿ(atm):')")
    read *,MUTABLE_P
    write(*,"('��������������������(��):')")
    read *,MUTABLE_ALPHA
    write(*,"('����������','P=',F8.3,'atm,',' ��=',F6.3)") MUTABLE_P,MUTABLE_ALPHA

    write(*,"('���ڼ��㡤����')")
    write(*,*)
    do i=800.0,1500.0,0.0001
        if(abs(T_s(i)-i)<=DELTA)then
            T=T_s(i)
            found=.true.
            write(*,"('�ҵ����Ͼ���(��=',F5.3,')���¶�')")DELTA
            write(*,"('ȼ���¶ȣ�',F10.5,' K')") T
            write(*,"('�ƽ���ȼ�٣�',F10.5,' mm/s')") m(t)*1E3/RHO_M
            write(*,"('�������')")
            exit
        end if
    end do
    write(*,*)
    if(found .neqv. .true.) write(*,"('δ�ҵ����Ͼ��ȵ��¶ȣ���������Ȼ��ߵ�������')")










end
