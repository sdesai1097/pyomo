
'''
Problem from Example 3.5.1: Water Treatment Network
(Generalized Disjunctive Programming: A Framework for Formulation and
 Alternative Algorithms for MINLP Optimization - Grossmann and Ruiz, 2012)

This example corresponds to a synthesis problem of a distributed
wastewater multicomponent network. Given a set of process liquid
streams with known composition, a set of technologies for the removal
of pollutants, and a set of mixers and splitters, the objective is to
find the interconnections of the technologies and their flowrates to 
meet the specified discharge composition of pollutant at minimum
total cost. Discrete choices involve deciding what equipment to use for 
each treatment unit.

Link to Article: 
https://pdfs.semanticscholar.org/d533/63bc770f93f876277136a5c6a9fba12a27e6.pdf

Original Problem from Galan and Grossmann, 1998, Example 1

Full GDP model and explanation of model for wastewater treatment network from
Lee and Grossmann (2003)
https://ac-els-cdn-com.proxy.library.cmu.edu/S009813540300098X/1-s2.0-S009813540300098X-main.pdf?_tid=b8fa9ebf-0035-4764-8cbd-aff04fe33274&acdnat=1544402786_0af0a426c9d1438838d007e0cd023399

'''

from pyomo.environ import *
from pyomo.gdp import *

def build_water_treatment_network_model():
    """Build the water treatment network model"""
    m = ConcreteModel(name = "Water Treatment Network")

    """Set declarations"""
    m.streams = RangeSet(3, 19, doc="Process streams")
    m.units = RangeSet(1, 9, doc="Process units")

    """Parameter and initial point declarations"""
    #Initial point information for stream flows
    initFlow = {3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0,
                14:0, 15:0, 16:0, 17:0, 18:0, 19:0}

    initConcA = {9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0, 17:0, 18:0, 
                19:0}

    initConcB = {9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0, 17:0, 18:0, 
                19:0}

    #Flow rate [=] k/h
    #Concentration [=] ppm

    #Process liquid streams 1 and 2 have known flows and pollutant conc
    s1F, s1A, s1B = 40, 20, 20
    s2F, s2A, s2B = 40, 15, 200
    #Streams 3-8 have set concentrations (split off from s1 and s2)
    s3A, s3B = 100, 20
    s4A, s4B = 100, 20
    s5A, s5B = 100, 20
    s6A, s6B = 15, 200
    s7A, s7B = 15, 200
    s8A, s8B = 15, 200
    #Stream 19 (final output) flow and composition is also known
    #Must meet the discharge composition regulations for each pollutant
    #10 ppm limit for each pollutant A and B
    s19F, s19A, s19B = 80, 10, 10 


    """Variable Declarations"""
    m.f = Var(m.streams, domain=NonNegativeReals, initialize=initFlow,
                  bounds=(0,80))

    m.cA = Var(m.streams, domain=NonNegativeReals, initialize=initConcA,
                  bounds=(0,100))

    m.cB = Var(m.streams, domain=NonNegativeReals, initialize=initConcB,
                  bounds=(0,100))


    """Constraint definitions"""

    #Mass Balances for S1 

    m.massbal_S1 = Constraint(expr=s1F = m.f[3] + m.f[4] + m.f[5])
    m.Abal_S1 = Constraint(expr=s1F*s1A = m.f[3]*s3A + m.f[4]*s4A + m.f[5]*s5A)
    m.Bbal_S1 = Constraint(expr=s1F*s1B = m.f[3]*s3B + m.f[4]*s4B + m.f[5]*s5B)
    z13 = m.f[3]/s1F
    z14 = m.f[3]/s1F
    z15 = m.f[3]/s1F
    m.splitA13 = Constraint(expr=s1F*s1A*z13 = m.f[3]*s3A)
    m.splitA14 = Constraint(expr=s1F*s1A*z14 = m.f[4]*s4A)
    m.splitA15 = Constraint(expr=s1F*s1A*z15 = m.f[5]*s5A)
    m.splitB13 = Constraint(expr=s1F*s1B*z13 = m.f[3]*s3B)
    m.splitB14 = Constraint(expr=s1F*s1B*z14 = m.f[4]*s4B)
    m.splitB15 = Constraint(expr=s1F*s1B*z15 = m.f[5]*s5B)

    #Mass Balances for S2
    
    m.massbal_S2 = Constraint(expr=s2F = m.f[6] + m.f[7] + m.f[8])
    m.Abal_S2 = Constraint(expr=s2F*s2A = m.f[6]*s6A + m.f[7]*s7A + m.f[8]*s8A)
    m.Bbal_S2 = Constraint(expr=s2F*s2B = m.f[6]*s6B + m.f[7]*s7B + m.f[8]*s8B)
    z26 = m.f[6]/s2F
    z27 = m.f[7]/s2F
    z28 = m.f[8]/s2F
    m.splitA26 = Constraint(expr=s2F*s2A*z26 = m.f[6]*s6A)
    m.splitA27 = Constraint(expr=s2F*s2A*z27 = m.f[7]*s7A)
    m.splitA28 = Constraint(expr=s2F*s2A*z28 = m.f[8]*s8A)
    m.splitB26 = Constraint(expr=s2F*s2B*z26 = m.f[6]*s6B)
    m.splitB27 = Constraint(expr=s2F*s2B*z27 = m.f[7]*s7B)
    m.splitB28 = Constraint(expr=s2F*s2B*z28 = m.f[8]*s8B)

    #Mass Balances for M1

    m.massbal_M1 = Constraint(expr=m.f[3] + m.f[7] + m.f[15] = m.f[9])
    m.Abal_M1 = Constraint(expr=m.f[3]*s3A + m.f[7]*s7A + m.f[15]*m.cA[15] = m.f[9]*m.cA[9])
    m.Bbal_M1 = Constraint(expr=m.f[3]*s3B + m.f[7]*s7B + m.f[15]*m.cB[15] = m.f[9]*m.cB[9])

    #Mass Balances for M2

    m.massbal_M2 = Constraint(expr=m.f[4] + m.f[6] + m.f[18] = m.f[10])
    m.Abal_M2 = Constraint(expr=m.f[4]*s4A + m.f[6]*s6A + m.f[18]*m.cA[18] = m.f[10]*m.cA[10])
    m.Bbal_M2 = Constraint(expr=m.f[4]*s4B + m.f[6]*s6B + m.f[18]*m.cB[18] = m.f[10]*m.cB[10])

    #Mass Balance for TX

    m.massbal_TX = Constraint(expr=m.f[9] = m.f[11])

    #Mass Balance for TXX

    m.massbal_TXX = Constraint(expr=m.f[10] = m.f[12])

    #Mass Balances for S3

    m.massbal_S3 = Constraint(expr=m.f[11] = m.f[13] + m.f[14] + m.f[15])
    m.Abal_S3 = Constraint(expr=m.f[11]*m.cA[11] = m.f[13]*m.cA[13] + m.f[14]*m.cA[14] + m.f[15]*m.cA[15])
    m.Bbal_S3 = Constraint(expr=m.f[11]*m.cB[11] = m.f[13]*m.cB[13] + m.f[14]*m.cB[14] + m.f[15]*m.cB[15])
    z1113 = m.f[13]/m.f[11]
    z1114 = m.f[14]/m.f[11]
    z1115 = m.f[15]/m.f[11]
    m.splitA1113 = Constraint(expr=m.f[11]*m.cA[11]*z1113 = m.f[13]*m.cA[13])
    m.splitA1114 = Constraint(expr=m.f[11]*m.cA[11]*z1114 = m.f[14]*m.cA[14])
    m.splitA1115 = Constraint(expr=m.f[11]*m.cA[11]*z1115 = m.f[15]*m.cA[15])
    m.splitB1113 = Constraint(expr=m.f[11]*m.cB[11]*z1113 = m.f[13]*m.cB[13])
    m.splitB1114 = Constraint(expr=m.f[11]*m.cB[11]*z1114 = m.f[14]*m.cB[14])
    m.splitB1115 = Constraint(expr=m.f[11]*m.cB[11]*z1115 = m.f[15]*m.cB[15])

    #Mass Balances for S4

    m.massbal_S4 = Constraint(expr=m.f[12] = m.f[16] + m.f[17] + m.f[18])
    m.Abal_S4 = Constraint(expr=m.f[12]*m.cA[12] = m.f[16]*m.cA[16] + m.f[17]*m.cA[17] + m.f[18]*m.cA[18])
    m.Bbal_S4 = Constraint(expr=m.f[12]*m.cB[12] = m.f[16]*m.cB[16] + m.f[17]*m.cB[17] + m.f[18]*m.cB[18])
    z1216 = m.f[16]/m.f[12]
    z1217 = m.f[17]/m.f[12]
    z1218 = m.f[18]/m.f[12]
    m.splitA1216 = Constraint(expr=m.f[12]*m.cA[12]*z1216 = m.f[16]*m.cA[16])
    m.splitA1217 = Constraint(expr=m.f[12]*m.cA[12]*z1217 = m.f[17]*m.cA[17])
    m.splitA1218 = Constraint(expr=m.f[12]*m.cA[12]*z1218 = m.f[18]*m.cA[18])
    m.splitB1216 = Constraint(expr=m.f[12]*m.cB[12]*z1216 = m.f[16]*m.cB[16])
    m.splitB1217 = Constraint(expr=m.f[12]*m.cB[12]*z1217 = m.f[17]*m.cB[17])
    m.splitB1218 = Constraint(expr=m.f[12]*m.cB[12]*z1218 = m.f[18]*m.cB[18])

    #Mass Balances for M3

    m.massbal_M3 = Constraint(expr=m.f[5] + m.f[8] + m.f[13] + m.f[16] = s19F)
    m.Abal_M3 = Constraint(expr=m.f[5]*s5A + m.f[8]*s8A + m.f[13]*m.cA[13] + m.f[16]*m.cA[16] = s19F*s19A)
    m.Bbal_M3 = Constraint(expr=m.f[5]*s5B + m.f[8]*s8B + m.f[13]*m.cB[13] + m.f[16]*m.cB[16] = s19F*s19B)

    """Disjunctions"""

    #Different types of equipment/technologies can be used for each treatment unit
    #Technologies have constant removal ratios for each pollutant

    #Treatment unit TX (Equipment TX1 v Equipment TX2)
    m.use_TX1 = Disjunct()
    m.use_TX2 = Disjunct()
    m.use_TX1.removeA = Constraint(expr=m.cA[11] = 0.1*m.cA[9])
    m.use_TX1.removeB = Constraint(expr=m.cB[11] = m.cB[9])
    m.use_TX2.removeA = Constraint(expr=m.cA[11] = 0.5*m.cA[9])
    m.use_TX2.removeB = Constraint(expr=m.cB[11] = 0.3*m.cB[9])
    m.use_TX1_or_TX2 = Disjunction(expr=[m.use_TX1, m.use_TX2])

    #Treatment unit TXX (Equipment TXX1 v Equipment TXX2)

    m.use_TXX1 = Disjunct()
    m.use_TXX2 = Disjunct()
    m.use_TXX1.removeA = Constraint(expr=m.cA[12] = m.cA[10])
    m.use_TXX1.removeB = Constraint(expr=m.cB[12] = 0.1*m.cB[10])
    m.use_TXX2.removeA = Constraint(expr=m.cA[12] = 0.5*m.cA[10])
    m.use_TXX2.removeB = Constraint(expr=m.cB[12] = 0.01*m.cB[10])
    m.use_TXX1_or_TXX2 = Disjunction(expr=[m.use_TX1, m.use_TX2])

    """Objective function definition"""

    # Objective: Minimize total inflet flow to treatment units (simplification
    # of minimizing total cost by minimizing amount of flow that has to be
    # treated)
    
    model.flow_treated = m.f[9] + m.f[10]
    model.treatment_inlet = Objective(expr=model.flow_treated,
                                   doc="Minimize flowrate")

model = build_water_treatment_network_model()


