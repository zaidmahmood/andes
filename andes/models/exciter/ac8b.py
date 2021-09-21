from andes.core.param import NumParam
from andes.core.var import Algeb

from andes.core.service import ConstService, VarService, FlagValue

from andes.core.block import LagAntiWindup, LeadLag, Washout, Lag, HVGate
from andes.core.block import LessThan, LVGate, AntiWindup, IntegratorAntiWindup
from andes.core.block import Integrator, PIDAWHardLimit, Piecewise, PIDController

from andes.models.exciter.excbase import ExcBase, ExcBaseData
from andes.models.exciter.saturation import ExcQuadSat


class AC8BData(ExcBaseData):
    def __init__(self):
        ExcBaseData.__init__(self)
        self.TR = NumParam(info='Sensing time constant',
                           tex_name='T_R',
                           default=0.01,
                           unit='p.u.',
                           )

        self.kP = NumParam(info='PID proportional coeff.',
                           tex_name='k_P',
                           default=10,
                           vrange=(10, 500),
                           )
        self.kI = NumParam(info='PID integrative coeff.',
                           tex_name='k_I',
                           default=10,
                           vrange=(10, 500),
                           )
        self.kD = NumParam(info='PID direvative coeff.',
                           tex_name='k_D',
                           default=10,
                           vrange=(10, 500),
                           )
        self.Td = NumParam(info='PID direvative time constant.',
                           tex_name='T_d',
                           default=0.2,
                           vrange=(0, 0.5),
                           )

        self.VPMAX = NumParam(info='PID maximum limit',
                              tex_name='V_{PMAX}',
                              default=999,
                              unit='p.u.')
        self.VPMIN = NumParam(info='PID minimum limit',
                              tex_name='V_{PMIN}',
                              default=-999,
                              unit='p.u.')

        self.VRMAX = NumParam(info='Maximum excitation limit',
                              tex_name='V_{RMAX}',
                              default=7.3,
                              unit='p.u.',
                              vrange=(1, 10))
        self.VRMIN = NumParam(info='Minimum excitation limit',
                              tex_name='V_{RMIN}',
                              default=1,
                              unit='p.u.',
                              vrange=(-1, 1.5))

        # TODO: check default value for VFEMAX
        self.VFEMAX = NumParam(info='Maximum VFE',
                               tex_name=r'V_{FEMAX}',
                               default=999,
                               unit='p.u.')

        # TODO: check default value for VEMIN
        self.VEMIN = NumParam(info='Minimum excitation output',
                              tex_name=r'V_{EMIN}',
                              default=-999,
                              unit='p.u.')

        self.TA = NumParam(info='Lag time constant in anti-windup lag',
                           tex_name='T_A',
                           default=0.04,
                           unit='p.u.',
                           )
        self.KA = NumParam(info='Gain in anti-windup lag TF',
                           tex_name='K_A',
                           default=40,
                           unit='p.u.',
                           )
        self.TE = NumParam(info='Exciter integrator time constant',
                           tex_name='T_E',
                           default=0.8,
                           unit='p.u.',
                           )

        self.E1 = NumParam(info='First saturation point',
                           tex_name='E_1',
                           default=0.,
                           unit='p.u.',
                           )
        self.SE1 = NumParam(info='Value at first saturation point',
                            tex_name='S_{E1}',
                            default=0.,
                            unit='p.u.',
                            )
        self.E2 = NumParam(info='Second saturation point',
                           tex_name='E_2',
                           default=1.,
                           unit='p.u.',
                           )
        self.SE2 = NumParam(info='Value at second saturation point',
                            tex_name='S_{E2}',
                            default=1.,
                            unit='p.u.',
                            )

        self.KE = NumParam(info='Gain added to saturation',
                           tex_name='K_E',
                           default=1,
                           unit='p.u.',
                           )
        self.KD = NumParam(default=0,
                           info='Ifd feedback gain',
                           tex_name='K_C',
                           vrange=(0, 1),
                           )

        self.KC = NumParam(default=0.1,
                           info='Rectifier loading factor proportional to commutating reactance',
                           tex_name='K_C',
                           vrange=(0, 1),
                           )


class AC8BModel(ExcBase):
    def __init__(self, system, config):
        ExcBase.__init__(self, system, config)

        # Assume FEX is in (0, 0.433) at initial
        self.VE0 = ConstService(info='Initial VE',
                                 tex_name=r'V_{E0}',
                                 v_str='vf0 + 0.577 * KC * XadIfd')

        self.VR0 = ConstService(info='Initial VR',
                                tex_name=r'V_{R0}',
                                v_str='VFE0 + VE0')

        self.vref0 = ConstService(info='Initial reference voltage input',
                                  tex_name=r'V_{ref0}',
                                  v_str='v')

        # control block begin
        self.LG = Lag(self.v, T=self.TR, K=1,
                      info='Voltage transducer',
                      )

        self.UEL = Algeb(info='Interface var for under exc. limiter',
                         tex_name='U_{EL}',
                         v_str='0',
                         e_str='0 - UEL'
                         )
        self.OEL = Algeb(info='Interface var for over exc. limiter',
                         tex_name='O_{EL}',
                         v_str='0',
                         e_str='0 - OEL'
                         )
        self.Vs = Algeb(info='Voltage compensation from PSS',
                        tex_name='V_{s}',
                        v_str='0',
                        e_str='0 - Vs'
                        )

        self.vref = Algeb(info='Reference voltage input',
                          tex_name='V_{ref}',
                          unit='p.u.',
                          v_str='vref0',
                          e_str='vref0 - vref'
                          )

        self.vi = Algeb(info='Total input voltages',
                        tex_name='V_i',
                        unit='p.u.',
                        e_str='ue * (-LG_y + vref + UEL + OEL + Vs - vi)',
                        v_str='-v + vref',
                        diag_eps=True,
                        )

        self.PIDAW = PIDAWHardLimit(u=self.vi, kp=self.kP, ki=self.kI,
                                    kd=self.kD, Td=self.Td, x0='VR0 / KA',
                                    aw_lower=self.VPMIN, aw_upper=self.VPMAX,
                                    lower=self.VPMIN, upper=self.VPMAX,
                                    tex_name='PID', info='PID', name='PIDAW',
                                    )

        self.LA = LagAntiWindup(u=self.PIDAW_y,
                                T=self.TA,
                                K=self.KA,
                                upper=self.VRMAX,
                                lower=self.VRMIN,
                                info=r'V_{R}, Anti-windup lag',
                                )

        self.VEMAX = VarService(info='Maximum excitation output',
                                tex_name=r'V_{EMAX}',
                                v_str='safe_div(VFEMAX - KD * XadIfd, KE + Se)')

        # --- debug
        self.se0=ConstService(v_str='Indicator(INT_y > SAT_A) * SAT_B * (INT_y - SAT_A) ** 2')
        # --- debug end

        # LA_y is VR
        # TODO: check max and min
        self.INT = IntegratorAntiWindup(u='ue * (LA_y - VFE)',
                                        T=self.TE,
                                        K=1,
                                        y0=self.VE0,
                                        lower=self.VEMIN,
                                        upper=self.VEMAX,
                                        info=r'V_{E}, Integrator Anti-windup',
                                        )

        self.SAT = ExcQuadSat(self.E1, self.SE1, self.E2, self.SE2,
                              info='Field voltage saturation',
                              )

        self.SL = LessThan(u=self.INT_y, bound=self.SAT_A, equal=False, enable=True, cache=False)

        # SL_z0 indicates saturation
        self.Se = Algeb(tex_name=r"V_{out}*S_e(|V_{out}|)", info='saturation output',
                        v_str='Indicator(INT_y > SAT_A) * SAT_B * (INT_y - SAT_A) ** 2',
                        e_str='ue * (SL_z0 * (INT_y - SAT_A) ** 2 * SAT_B - Se)',
                        diag_eps=True,
                        )


        self.VFE0 = ConstService(info='Initial VFE', tex_name=r'V_{FE0}',
                                 v_str='INT_y * KE + Se + XadIfd * KD',
                                 )

        # INT_y is VE
        self.VFE = Algeb(info='Combined saturation feedback',
                         tex_name='V_{FE}',
                         unit='p.u.',
                         v_str='INT_y * KE + Se + XadIfd * KD',
                         e_str='ue * (INT_y * KE + Se + XadIfd * KD - VFE)',
                         diag_eps=True,
                         )

        self.IN = Algeb(tex_name='I_N',
                        info='Input to FEX',
                        v_str='safe_div(KC * XadIfd, INT_y)',
                        e_str='ue * (KC * XadIfd - INT_y * IN)',
                        diag_eps=True,
                        )

        # TODO: Check funs
        # Copy from ESST3A
        self.FEX = Piecewise(u=self.IN,
                             points=(0, 0.433, 0.75, 1),
                             funs=('1', '1 - 0.577*IN', 'sqrt(0.75 - IN ** 2)', '1.732*(1 - IN)', 0),
                             info='Piecewise function FEX',
                             )

        self.vout.e_str = 'FEX_y * INT_y - vout'


class AC8B(AC8BData, AC8BModel):
    """
    Exciter AC8B model.

    Reference:

    [1] PowerWorld, Exciter AC8B, [Online],

    [2] NEPLAN, Exciters Models, [Online],

    Available:

    https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20AC8B.htm

    https://www.neplan.ch/wp-content/uploads/2015/08/Nep_EXCITERS1.pdf
    """
    def __init__(self, system, config):
        AC8BData.__init__(self)
        AC8BModel.__init__(self, system, config)
