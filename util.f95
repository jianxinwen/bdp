!!����ģ��
!!��λͳһΪK,kg,m,J,s
module const
    implicit none
    !=====================
    !ʵ��ָ���Ĳ���
    !=====================
    real::MUTABLE_P=70.0         !ѹǿ(atm)
    real::MUTABLE_ALPHA=0.85     !�ƽ���������������������
    !=====================
    !��������
    !=====================
    real,parameter::CONST_R_0=8.315              !����״̬����(J/(mol k))
    real,parameter::CONST_CAL_TO_R0_RATIO=4.18   !1cal=4.18J
    !=====================
    !�ƽ�������
    !=====================
    real,parameter::CONST_T_f=2545.0            !�ƽ������Ȼ����¶�(K)
    real,parameter::CONST_T_0=300.0             !(K)
    real,parameter::CONST_AVERAGE_M=2.62E-2     !ƽ��Ħ������(kg/mol)
    real,parameter::CONST_D_0=2.0E-5            !AP���ӳ�ʼֱ��(m)
    real,parameter::CONST_RHO_f=1.21E3          !ȼ���ܶ�(kg/m^3)
    real,parameter::CONST_RHO_ox=1.95E3         !�������ܶ�(kg/m^3)
    !=====================
    !�������
    !=====================
    real,parameter::CONST_DELTA_PF=1.5                      !PF�淴Ӧ����
    real,parameter::CONST_k_PF=3.0E4                        !PF�淴Ӧ���ʳ���(kg/(m^3 s atm^��PF))
    real,parameter::CONST_D_a=1.6E-5                        !�ο��¶�Ta�Ͳο�ѹǿPa�µ���ɢϵ�� (m^2/s(���³�ѹ))
    real,parameter::CONST_LAMDA=0.3*CONST_CAL_TO_R0_RATIO   !�ȴ���ϵ��(J/(m s K))
    real,parameter::CONST_C_s=3.0E2*CONST_CAL_TO_R0_RATIO   !�ƽ�������͹����ƽ������ (J/(kg K))
    real,parameter::CONST_A_fh=0.3                          !ƽ������߶�����
    real,parameter::CONST_E_PF=1.5E4*CONST_CAL_TO_R0_RATIO  !(J/mol)
    !=====================
    !AP����
    !=====================
    real,parameter::CONST_T_AP=1400.0                       !AP��Ԫ�ƽ����ľ���ȼ���¶�(K)
    real,parameter::CONST_DELTA_AP=1.8                      !AP�淴Ӧ����
    real,parameter::CONST_E_ox=2.2E4*CONST_CAL_TO_R0_RATIO  !(J/mol)
    real,parameter::CONST_A_ox=3.0E6                        !(kg/(m^2 s))
    real,parameter::CONST_Q_L=-1.2E5*CONST_CAL_TO_R0_RATIO  !�����������෴Ӧ��(J/kg)
    real,parameter::CONST_E_AP=3.0E4*CONST_CAL_TO_R0_RATIO  !(J/mol)
    real,parameter::CONST_k_AP=1.12E3                       !AP�淴Ӧ���ʳ��� (kg/(m^3 s atm^��AP))
    !=====================
    !ճ�ϼ�����
    !=====================
    real,parameter::CONST_E_f=1.5E4*CONST_CAL_TO_R0_RATIO   !(J/mol)
    real,parameter::CONST_A_f=2.5E4                         !(kg/(m^2 s))
    real,parameter::CONST_Q_f=5.0E4*CONST_CAL_TO_R0_RATIO   !ճ�ϼ������෴Ӧ��(J/kg)
end module


!���㹫ʽģ��
module formula
    use const
    implicit none

    contains
    !!�����������ֽ��ٶ�
    real function m_ox(T_s_)
        real::T_s_
        m_ox=CONST_A_ox*exp(-CONST_E_ox/(CONST_R_0*T_s_))
        return
    end function

    !!ճ�ϼ������ֽ��ٶ�
    real function m_f(T_s_)
        real::T_s_
        m_f=CONST_A_ox*exp(-CONST_E_f/(CONST_R_0*T_s_))
        return
    end function

    !!��� r_ox,r_f
    real function r_ox(T_s_)
        real::T_s_
        r_ox=m_ox(T_s_)/CONST_RHO_ox
    end function
    real function r_f(T_s_)
        real::T_s_
        r_f=m_f(T_s_)/CONST_RHO_f
    end function

    !!��� h/D0
    real function h_div_D0(T_s_)
        real::T_s_
        real::m=0.8,n=0.72,C_ign=1.90E12      !nҲ����ȡֵ0.75,C_ign��λ��atm/m^5
        real::t_ign
        t_ign=C_ign*(CONST_D_0**(m+1.0))/(MUTABLE_P**n)
        h_div_D0=(1.0/2.0)*(1.0+1.0/sqrt(3.0))*(1.0-r_ox(T_s_)/r_f(T_s_))+r_ox(T_s_)*t_ign/CONST_D_0
    end function

    !!�ƽ����������������������
    real function zeta()
        zeta=1.0/(1.0+((1.0-MUTABLE_ALPHA)*CONST_RHO_ox)/(MUTABLE_ALPHA*CONST_RHO_f))
    end function

    !!��� Sox/S0
    real function Sox_div_S0(T_s_)
        real::T_s_
        real::hd2,zt
        real::numerator,denominator
        hd2=h_div_D0(T_s_)**2
        zt=zeta()
        numerator=zt*(6.0*hd2+1.0)
        denominator=1.0+6.0*zt*hd2
        Sox_div_S0=numerator/denominator
    end function

    !!��� m ,�ƽ�������ȼ��
    real function m(T_s_)
        real::T_s_
        m=m_ox(T_s_)*Sox_div_S0(T_s_)/MUTABLE_ALPHA
    end function


    !!���AP�淴Ӧ���� x_AP
    real function x_AP(T_s_)
        real::T_s_
         x_AP=m_ox(T_s_)/(CONST_k_AP*(MUTABLE_P**CONST_DELTA_AP))
    end function

    !!�����淴Ӧ���� x_Pch
    real function x_Pch(T_s_)
        real::T_s_
         x_Pch=m(T_s_)/(CONST_k_PF*(MUTABLE_P**CONST_DELTA_PF))
    end function

    !!��� x_Fd
    real function x_Fd(T_s_)
        real::T_s_
        real::R,T
        R=CONST_R_0/CONST_AVERAGE_M
        T=CONST_T_f
        x_Fd=(R*m_ox(T_s_)*(CONST_D_0**2))/(4*CONST_D_a*(T**0.75))
    end function

    !!��� x_Pd
    real function x_Pd(T_s_)
        real::T_s_
        real::R,T
        R=CONST_R_0/CONST_AVERAGE_M
        T=CONST_T_AP
        x_Pd=(R*m_ox(T_s_)*(CONST_D_0**2))/(4*CONST_D_a*(T**0.75))
    end function

    !!��� average_x_Pd
    real function average_x_Pd(T_s_)
        real::T_s_
        average_x_Pd=CONST_A_fh*x_Pd(T_s_)

    end function
    !!��� average_x_Fd
    real function average_x_Fd(T_s_)
        real::T_s_
        average_x_Fd=CONST_A_fh*x_Fd(T_s_)
    end function




    !!��� XI_AP
    real function XI_AP(T_s_)
        real::T_s_
        XI_AP=m_ox(T_s_)*CONST_C_s*x_AP(T_s_)/CONST_LAMDA
    end function
    !!��� XI_FF
    real function XI_FF(T_s_)
        real::T_s_,x_FF
        x_FF=x_AP(T_s_)+average_x_Fd(T_s_)
        XI_FF=m_ox(T_s_)*CONST_C_s*x_FF/CONST_LAMDA
    end function
    !!��� PHI_PF
    real function PHI_PF(T_s_)
        real::T_s_,x_PF
        x_PF=x_Pd(T_s_)+x_Pch(T_s_)
        PHI_PF=m(T_s_)*CONST_C_s*x_PF/CONST_LAMDA
    end function

    !!��� Q_AP
    real function Q_AP()
        Q_AP=CONST_C_s*(CONST_T_AP-CONST_T_0)+CONST_Q_L
    end function
    !!��� Q_FF
    real function Q_FF()
        Q_FF=CONST_C_s*(CONST_T_f-CONST_T_0)-MUTABLE_ALPHA*CONST_C_s*(CONST_T_AP-CONST_T_0)+(1.0-MUTABLE_ALPHA)*CONST_Q_f
    end function
    !!��� Q_PF
    real function Q_PF()
        Q_PF=CONST_C_s*(CONST_T_f-CONST_T_0)+MUTABLE_ALPHA*CONST_Q_L+(1.0-MUTABLE_ALPHA)*CONST_Q_f
    end function
    !!��� BETA_F
    real function BETA_F(T_s_)
        real::T_s_
        BETA_F=(x_AP(T_s_)-x_Pch(T_s_))/x_Pd(T_s_)
    end function

    !!��� T_s
    real function T_s(T_s_)
        real::T_s_
        T_s=CONST_T_0-MUTABLE_ALPHA*CONST_Q_L/CONST_C_s-&
        (1.0-MUTABLE_ALPHA)*CONST_Q_f/CONST_C_s+BETA_F(T_s_)*&
        (Q_PF()/CONST_C_s)*exp(-PHI_PF(T_s_))+&
        (1.0-BETA_F(T_s_))*&
        (MUTABLE_ALPHA/CONST_C_s)*(Q_AP()*exp(-XI_AP(T_s_))+Q_FF()*exp(-XI_FF(T_s_)))
    end function




end module

!!��С���˷�
module lsm
    implicit none


    contains
    !!���ϵ��k
    real function get_k(x,y,n)
        integer::n
        real::x(n),y(n),sum_x,sum_y,average_x,average_y,num=0,den=0
        integer::i,j
        sum_x=sum(x);sum_y=sum(y)
        average_x=sum_x/n;average_y=sum_y/n
        do i=1,n
            num=num+(x(i)-average_x)*(y(i)-average_y)
            den=den+(x(i)-average_x)**2
        end do
        get_k=num/den
    end function
    !!���ϵ��b
    real function get_b(x,y,n)
        integer::n
        real::x(n),y(n),sum_x,sum_y,average_x,average_y
        integer::i,j
        sum_x=sum(x);sum_y=sum(y)
        average_x=sum_x/n;average_y=sum_y/n
        get_b=average_y-get_k(x,y,n)*average_x
    end function

end module



