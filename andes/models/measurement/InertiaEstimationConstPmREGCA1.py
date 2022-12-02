"""
Inertia estimation model based on swing equation with constant Pm
 
"""
from andes.core.service import ExtService, PostInitService
from andes.core import ConstService, NumParam, ModelData, Model, IdxParam, ExtState, State, ExtAlgeb, ExtParam, Algeb
from andes.core.block import  DeadBand1, Piecewise, Gain, Integrator, Lag, Washout


class InertiaEstimationConstPmREGCA1(ModelData, Model):
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
        self.res = IdxParam(model='REGCA1',
                            info='RES idx',
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
        self.Tm = NumParam(default=0.01,
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
        self.Mg = ExtService(src='M',
                            model='SynGen',
                            indexer=self.syn,
                            tex_name = 'generator inertia'
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
                           model='REGCA1',
                           indexer=self.res,
                           tex_name = 'Pe',
                           export = True
                           )
        self.PeAlgeb = Algeb(v_str = 'Pe',
                           e_str = 'Pe-PeAlgeb'
                           )

        #self.Pm = PostInitService(info='Initial Pe',
        #                     tex_name='P_{e0}', v_str='PeAlgeb' )
        self.Pm = ConstService(v_str='Pe', info='initial Pe',
                               tex_name='P_{m}',
                               )
        
        self.pdiff = Algeb(v_str = '0',
                           e_str = '(Pe - Pm) - pdiff',
                           tex_name = 'P_{diff}'
                           )
        #PmTest = Pm
        self.PmTest = Algeb(v_str = 'Pe',
                            e_str = 'Pm - PmTest')
        self.windowswitch = NumParam(default = 1,
                                    info = 'Switch to control the window of ON period')
        ## Liu Blocks ############################################################################################################
        
        self.Tsignal = NumParam(default=0.01,
                           info="Time Constant",
                           unit="sec",
                           tex_name='T_m',
                           )
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
                          K = 0.001, T = self.Tsignal)
        #########################################################################################################################
        #main blocks
        #self.condition = Algeb(v_str = 'M_star * omega_dot + pdiff',
        #                    e_str = 'M_star * omega_dot + pdiff - condition')

        #self.piececondition = Piecewise(u = self.condition , points= [-0.0000001, +0.00000010], funs= [1, 0, -1],    
        #                       name = 'piececondition')
        self.peak = Piecewise(u = self.signal_y, points= ['negepsilon', 'epsilon'], funs= [1, 0, 1], 
                               name = 'peak')
        
        self.sign = Piecewise(u = self.omega_dot, points= ['negepsilon', 'epsilon'], funs= [1, 0, -1], 
                               name = 'sign') 
        #self.pdiffswitch = Piecewise(u = self.pdiff, points= [-.0001, 0.0001], funs= [1, 0, 1], 
        #                       ) 
 
        self.M_star = State(v_str = '0',
                            e_str = 'windowswitch  * sign_y * (M_star*omega_dot+(Pe-Pm))',
                            t_const= self.Tm,
                            info = "Estimated Inertia",
                            tex_name= 'M^{*}'
                            )
        self.T_mlag = NumParam(default=0.1,
                           info="Time Constant for M Lag Filter",
                           unit="sec",
                           tex_name='Lag time constant',
                           )
        self.D_mlag = NumParam(default=1,
                           info="D Constant for M Lag Filter",
                           tex_name='Lag D constant',
                           )
        self.K_mlag = NumParam(default=1,
                           info="Gain Constant for M Lag Filter",
                           tex_name='Lag Gain constant',
                           )
        self.M_lag = Lag(u = self.M_star, K = self.K_mlag, D = self.D_mlag,
                            T = self.T_mlag)  