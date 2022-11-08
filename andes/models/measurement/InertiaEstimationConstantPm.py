"""
Inertia estimation model based on swing equation with constant Pm
 
"""
from andes.core import ConstService, NumParam, ModelData, Model, IdxParam, ExtState, State, ExtAlgeb, ExtParam, Algeb
from andes.core.block import  Lag, Piecewise, Washout, Gain, Integrator, LagFreeze
from andes.core.discrete import Limiter, Derivative


class InertiaEstimationConstantPm(ModelData, Model):
    """
    Estimates inertia of a device. Outputs estimation in pu value.
    """ 

    def __init__(self, system, config):
        ModelData.__init__(self)
        Model.__init__(self, system, config)
        self.flags.update({'tds': True})
        #parameters
        self.syn = IdxParam(model='SynGen',
                            info='Synchronous generator idx',
                            mandatory=True,
                            unique=True,
                            )
        self.gov = IdxParam(model='TurbineGov',
                            info='tgov idx',
                            mandatory=True,
                            unique=True,
                            )
        self.epsilon = NumParam(default=0.000001,
                           info="tolerance",
                           unit="p.u.",
                           tex_name=r'\epsilon',
                           )
        self.negepsilon = ConstService(v_str = '-1 * epsilon')
        self.pe_tol = NumParam(default=0.000001,
                           info="tolerance of pe_dot",
                           unit="p.u."
                           )
        self.neg_pe_tol = ConstService(v_str = '-1 * pe_tol')

        self.Tm = NumParam(default=0.01,
                           info="Time Constant",
                           unit="sec",
                           tex_name='T_m',
                           )
        self.Tpe = NumParam(default=0.01,
                           info="Time Constant",
                           unit="sec",
                           tex_name='T_{pe}',
                           )
        self.Tsignal = NumParam(default=0.01,
                           info="Time Constant",
                           unit="sec",
                           tex_name='T_m',
                           )
        self.iTm = ConstService(v_str = "1/Tm",
                           tex_name='1/Tm',
                           )
        self.Two = NumParam(default=0.0001,
                           info="washout time const",
                           unit="sec",
                           tex_name='T_wo',
                           )
        self.Tf = NumParam(default=0.0001,
                           info="filter time const",
                           unit="sec",
                           tex_name='T_f',
                           )
        self.Kp = NumParam(default=50,
                           info="proportional constant",
                           unit="p.u.",
                           tex_name='K_p',
                           )
        self.Ki = NumParam(default=1,
                           info="integral constant",
                           unit="p.u.",
                           tex_name='K_i'
                           )
        self.damping = ExtParam(src='D',
                                model='SynGen',
                                indexer=self.syn,
                                tex_name = 'damping',
                                export = True
                                )
        self.Mg = ExtParam(src='M',
                            model='SynGen',
                            indexer=self.syn,
                            tex_name = 'generator inertia',
                            export = True
                            )
        self.ug = ExtParam(src='u',
                           model='SynGen',
                           indexer=self.syn,
                           tex_name = 'generator connectivity',
                           export = True
                           )
        #variables
        self.omega = ExtState(src='omega',
                                model='SynGen',
                                indexer=self.syn,
                                tex_name = r'\dot \omega',
                                export = True
                                )
        self.tm = ExtAlgeb(src='tm',
                          model='SynGen',
                          indexer=self.syn,
                          tex_name = 'mechanical torque of generator',
                          export = True
                          )
        self.te = ExtAlgeb(src='te',
                          model='SynGen',
                          indexer=self.syn,
                          tex_name = 'electrical torque of generator',
                          export = True
                          )          
        self.omega_dot = Algeb(tex_name = r'\dot{\omega}', info = r'\dot{\omega}',
                              v_str = '0', 
                              e_str = 'ug * (-1 * damping * (omega - 1) - te + tm) / Mg - omega_dot'
                              )
        
        self.Pe = ExtAlgeb(src='Pe',
                           model='SynGen',
                           indexer=self.syn,
                           tex_name = 'Pe',
                           export = True
                           )
###Treatment of Pm to reset Pm for multidisturbance cases
        self.Peint = Algeb(v_str = 'Pe',
                            e_str = 'Pe - Peint')
        self.Pm0 = ConstService(v_str='Pe', info='initial Pe',
                               tex_name='P_{m0}',
                               )
        self.Pedot = Washout(u = self.Peint, T = 1, K = 1 )
        self.Pedotlimit = Piecewise(u = self.Pedot_y, points= ['neg_pe_tol','pe_tol'], funs= [1, 0, 1]
                               )

        self.Pm = LagFreeze(u = self.Pe, freeze= self.Pedotlimit_y, T = self.Tpe, K = 1)
########################################################### 
# 

        self.windowswitch = NumParam(default = 1,
                                    info = 'Switch to control the window of ON period')
        self.pdiff = Algeb(v_str = '(Pe - Pm_y)',
                           e_str = '(Pe - Pm_y) - pdiff',
                           tex_name= 'P_{diff}'
                           )
        self.PmTest = Algeb(v_str = 'Pm_y',
                            e_str = 'Pm_y - PmTest')
        #washout to find omegadoubledot
        self.omegawashout = NumParam(default=0.1,
                           info="Time Constant for omega washout",
                           unit="sec",
                           tex_name='washout constant',
                           )
        self.omegadoubledot = Washout(u = self.omega_dot, K = 1,
                            T = self.omegawashout)   
        #self.signal = Lag(u = self.omegadoubledot_y,
        #                  K = 1, T = self.Tsignal)
        ## Liu Blocks ############################################################################################################
        self.k_omega = Gain(u = "omega_dot - omegadot_star_y", 
                            K = self.Kp
                            )        
        self.omegadot_star = Integrator(u = self.k_omega_y, T = 1, K = self.Ki, 
                                        y0 = '0', check_init = False
                                        )
        self.omegadoubledotliu = Lag(u = "k_omega_y - omegadoubledotliu_y",
                                  K = 1, T = self.Tf
                                  )        
        self.signal = Lag(u = self.omegadoubledotliu_y,
                          K = 0.0001, T = self.Tsignal)
        #########################################################################################################################
        #main blocks
        self.peak = Piecewise(u = self.signal_y, points= ['negepsilon', 'epsilon'], funs= [1, 0, 1], 
                               name = 'peak')
        
        self.sign = Piecewise(u = self.omega_dot, points= ['negepsilon', 'epsilon'], funs= [1, 0, -1], 
                               name = 'sign')

        self.M_star = State(v_str = 'windowswitch * sign_y * ( M_star * omega_dot + (Pe - Pm_y))',
                            e_str = 'windowswitch * sign_y * ( M_star * omega_dot + (Pe - Pm_y))',
                            t_const= self.Tm,
                            info = "Estimated Inertia",
                            tex_name= 'M^{*}'

                            )