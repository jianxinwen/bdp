

program main
    use const
    use formula
    use lsm
    implicit none
    !=====================
    !申明所需变量
    !=====================
    real,allocatable::container_x(:),container_y(:)  !!回归数据(xi,yi)
    real::lsm_k,lsm_b,vielle_b,vielle_n
    real::min_p,max_p,step_p   !!指定压强测试范围和步长
    integer::idx=1      !!idx 坐标值
    real::i,j

    real::DELTA=0.01    !!求解精度
    real::T
    logical::found
    real::RHO_M         !!推进剂密度(kg/m^3)




    RHO_M=CONST_RHO_ox*MUTABLE_ALPHA+CONST_RHO_f*(1-MUTABLE_ALPHA)
    found=.false.

    !=====================
    !求解维也里公式 r=b P^n
    !=====================
    write(*,"('开始求解维也里燃速公式')")
    write(*,"('正在计算···')")
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
    write(*,"('维也里燃速公式：r=',F8.5,'*P^',F7.5)")vielle_b,vielle_n
    write(*,*)
    write(*,"('==================================')")
    write(*,*)
!    write(*,"('开始求解特定条件下的数据，计算条件：','P=',F8.3,'atm,',' α=',F6.3)") MUTABLE_P,MUTABLE_ALPHA

    !=====================
    !特定条件下的数据
    !=====================
    write(*,"('开始求解特定条件下的数据')")
    write(*,"('请输入压强(atm):')")
    read *,MUTABLE_P
    write(*,"('请输入氧化剂质量分数(α):')")
    read *,MUTABLE_ALPHA
    write(*,"('计算条件：','P=',F8.3,'atm,',' α=',F6.3)") MUTABLE_P,MUTABLE_ALPHA

    write(*,"('正在计算···')")
    write(*,*)
    do i=800.0,1500.0,0.0001
        if(abs(T_s(i)-i)<=DELTA)then
            T=T_s(i)
            found=.true.
            write(*,"('找到符合精度(δ=',F5.3,')的温度')")DELTA
            write(*,"('燃面温度：',F10.5,' K')") T
            write(*,"('推进剂燃速：',F10.5,' mm/s')") m(t)*1E3/RHO_M
            write(*,"('计算完成')")
            exit
        end if
    end do
    write(*,*)
    if(found .neqv. .true.) write(*,"('未找到符合精度的温度，请调整精度或者迭代步长')")










end
