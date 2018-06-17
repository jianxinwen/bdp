!!参数模块
!!单位统一为K,kg,m,J,s
module const
    implicit none
    !=====================
    !实验指定的参数
    !=====================
    real::MUTABLE_P=70.0         !压强(atm)
    real::MUTABLE_ALPHA=0.85     !推进剂中氧化剂的质量分数
    !=====================
    !其他参数
    !=====================
    real,parameter::CONST_R_0=8.315              !气体状态常数(J/(mol k))
    real,parameter::CONST_CAL_TO_R0_RATIO=4.18   !1cal=4.18J
    !=====================
    !推进剂参数
    !=====================
    real,parameter::CONST_T_f=2545.0            !推进剂绝热火焰温度(K)
    real,parameter::CONST_T_0=300.0             !(K)
    real,parameter::CONST_AVERAGE_M=2.62E-2     !平均摩尔质量(kg/mol)
    real,parameter::CONST_D_0=2.0E-5            !AP粒子初始直径(m)
    real,parameter::CONST_RHO_f=1.21E3          !燃料密度(kg/m^3)
    real,parameter::CONST_RHO_ox=1.95E3         !氧化剂密度(kg/m^3)
    !=====================
    !火焰参数
    !=====================
    real,parameter::CONST_DELTA_PF=1.5                      !PF焰反应级数
    real,parameter::CONST_k_PF=3.0E4                        !PF焰反应速率常数(kg/(m^3 s atm^δPF))
    real,parameter::CONST_D_a=1.6E-5                        !参考温度Ta和参考压强Pa下的扩散系数 (m^2/s(常温常压))
    real,parameter::CONST_LAMDA=0.3*CONST_CAL_TO_R0_RATIO   !热传导系数(J/(m s K))
    real,parameter::CONST_C_s=3.0E2*CONST_CAL_TO_R0_RATIO   !推进剂气相和固相的平均比热 (J/(kg K))
    real,parameter::CONST_A_fh=0.3                          !平均火焰高度因子
    real,parameter::CONST_E_PF=1.5E4*CONST_CAL_TO_R0_RATIO  !(J/mol)
    !=====================
    !AP参数
    !=====================
    real,parameter::CONST_T_AP=1400.0                       !AP单元推进剂的绝热燃烧温度(K)
    real,parameter::CONST_DELTA_AP=1.8                      !AP焰反应级数
    real,parameter::CONST_E_ox=2.2E4*CONST_CAL_TO_R0_RATIO  !(J/mol)
    real,parameter::CONST_A_ox=3.0E6                        !(kg/(m^2 s))
    real,parameter::CONST_Q_L=-1.2E5*CONST_CAL_TO_R0_RATIO  !氧化剂凝聚相反应热(J/kg)
    real,parameter::CONST_E_AP=3.0E4*CONST_CAL_TO_R0_RATIO  !(J/mol)
    real,parameter::CONST_k_AP=1.12E3                       !AP焰反应速率常数 (kg/(m^3 s atm^δAP))
    !=====================
    !粘合剂参数
    !=====================
    real,parameter::CONST_E_f=1.5E4*CONST_CAL_TO_R0_RATIO   !(J/mol)
    real,parameter::CONST_A_f=2.5E4                         !(kg/(m^2 s))
    real,parameter::CONST_Q_f=5.0E4*CONST_CAL_TO_R0_RATIO   !粘合剂凝聚相反应热(J/kg)
end module


!计算公式模块
module formula
    use const
    implicit none

    contains
    !!氧化剂质量分解速度
    real function m_ox(T_s_)
        real::T_s_
        m_ox=CONST_A_ox*exp(-CONST_E_ox/(CONST_R_0*T_s_))
        return
    end function

    !!粘合剂质量分解速度
    real function m_f(T_s_)
        real::T_s_
        m_f=CONST_A_ox*exp(-CONST_E_f/(CONST_R_0*T_s_))
        return
    end function

    !!求解 r_ox,r_f
    real function r_ox(T_s_)
        real::T_s_
        r_ox=m_ox(T_s_)/CONST_RHO_ox
    end function
    real function r_f(T_s_)
        real::T_s_
        r_f=m_f(T_s_)/CONST_RHO_f
    end function

    !!求解 h/D0
    real function h_div_D0(T_s_)
        real::T_s_
        real::m=0.8,n=0.72,C_ign=1.90E12      !n也可以取值0.75,C_ign单位：atm/m^5
        real::t_ign
        t_ign=C_ign*(CONST_D_0**(m+1.0))/(MUTABLE_P**n)
        h_div_D0=(1.0/2.0)*(1.0+1.0/sqrt(3.0))*(1.0-r_ox(T_s_)/r_f(T_s_))+r_ox(T_s_)*t_ign/CONST_D_0
    end function

    !!推进剂中氧化剂的体积分数ζ
    real function zeta()
        zeta=1.0/(1.0+((1.0-MUTABLE_ALPHA)*CONST_RHO_ox)/(MUTABLE_ALPHA*CONST_RHO_f))
    end function

    !!求解 Sox/S0
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

    !!求解 m ,推进剂质量燃速
    real function m(T_s_)
        real::T_s_
        m=m_ox(T_s_)*Sox_div_S0(T_s_)/MUTABLE_ALPHA
    end function


    !!求解AP焰反应距离 x_AP
    real function x_AP(T_s_)
        real::T_s_
         x_AP=m_ox(T_s_)/(CONST_k_AP*(MUTABLE_P**CONST_DELTA_AP))
    end function

    !!求解初焰反应距离 x_Pch
    real function x_Pch(T_s_)
        real::T_s_
         x_Pch=m(T_s_)/(CONST_k_PF*(MUTABLE_P**CONST_DELTA_PF))
    end function

    !!求解 x_Fd
    real function x_Fd(T_s_)
        real::T_s_
        real::R,T
        R=CONST_R_0/CONST_AVERAGE_M
        T=CONST_T_f
        x_Fd=(R*m_ox(T_s_)*(CONST_D_0**2))/(4*CONST_D_a*(T**0.75))
    end function

    !!求解 x_Pd
    real function x_Pd(T_s_)
        real::T_s_
        real::R,T
        R=CONST_R_0/CONST_AVERAGE_M
        T=CONST_T_AP
        x_Pd=(R*m_ox(T_s_)*(CONST_D_0**2))/(4*CONST_D_a*(T**0.75))
    end function

    !!求解 average_x_Pd
    real function average_x_Pd(T_s_)
        real::T_s_
        average_x_Pd=CONST_A_fh*x_Pd(T_s_)

    end function
    !!求解 average_x_Fd
    real function average_x_Fd(T_s_)
        real::T_s_
        average_x_Fd=CONST_A_fh*x_Fd(T_s_)
    end function




    !!求解 XI_AP
    real function XI_AP(T_s_)
        real::T_s_
        XI_AP=m_ox(T_s_)*CONST_C_s*x_AP(T_s_)/CONST_LAMDA
    end function
    !!求解 XI_FF
    real function XI_FF(T_s_)
        real::T_s_,x_FF
        x_FF=x_AP(T_s_)+average_x_Fd(T_s_)
        XI_FF=m_ox(T_s_)*CONST_C_s*x_FF/CONST_LAMDA
    end function
    !!求解 PHI_PF
    real function PHI_PF(T_s_)
        real::T_s_,x_PF
        x_PF=x_Pd(T_s_)+x_Pch(T_s_)
        PHI_PF=m(T_s_)*CONST_C_s*x_PF/CONST_LAMDA
    end function

    !!求解 Q_AP
    real function Q_AP()
        Q_AP=CONST_C_s*(CONST_T_AP-CONST_T_0)+CONST_Q_L
    end function
    !!求解 Q_FF
    real function Q_FF()
        Q_FF=CONST_C_s*(CONST_T_f-CONST_T_0)-MUTABLE_ALPHA*CONST_C_s*(CONST_T_AP-CONST_T_0)+(1.0-MUTABLE_ALPHA)*CONST_Q_f
    end function
    !!求解 Q_PF
    real function Q_PF()
        Q_PF=CONST_C_s*(CONST_T_f-CONST_T_0)+MUTABLE_ALPHA*CONST_Q_L+(1.0-MUTABLE_ALPHA)*CONST_Q_f
    end function
    !!求解 BETA_F
    real function BETA_F(T_s_)
        real::T_s_
        BETA_F=(x_AP(T_s_)-x_Pch(T_s_))/x_Pd(T_s_)
    end function

    !!求解 T_s
    real function T_s(T_s_)
        real::T_s_
        T_s=CONST_T_0-MUTABLE_ALPHA*CONST_Q_L/CONST_C_s-&
        (1.0-MUTABLE_ALPHA)*CONST_Q_f/CONST_C_s+BETA_F(T_s_)*&
        (Q_PF()/CONST_C_s)*exp(-PHI_PF(T_s_))+&
        (1.0-BETA_F(T_s_))*&
        (MUTABLE_ALPHA/CONST_C_s)*(Q_AP()*exp(-XI_AP(T_s_))+Q_FF()*exp(-XI_FF(T_s_)))
    end function




end module

!!最小二乘法
module lsm
    implicit none


    contains
    !!求解系数k
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
    !!求解系数b
    real function get_b(x,y,n)
        integer::n
        real::x(n),y(n),sum_x,sum_y,average_x,average_y
        integer::i,j
        sum_x=sum(x);sum_y=sum(y)
        average_x=sum_x/n;average_y=sum_y/n
        get_b=average_y-get_k(x,y,n)*average_x
    end function

end module



