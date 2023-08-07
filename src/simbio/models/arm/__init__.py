from simbio.components import EmptyGroup, Override, Parameter, Species
from simbio.models.earm import AlbeckAsMatlab
from simbio.reactions.compound import ReversibleSynthesis
from simbio.reactions.enzymatic import MichaelisMenten


class Sensor(EmptyGroup):
    dimer: Species
    monomer: Species = 0


class ARM(AlbeckAsMatlab):
    KF: Parameter
    KR: Parameter
    KC: Parameter

    Apop: Species

    class Bid(AlbeckAsMatlab.Bid):
        pass

    class C3(AlbeckAsMatlab.C3):
        pass

    class C8(AlbeckAsMatlab.C8):
        pass

    L: Species[Override] = 0  # Set extrinsic stimuli to 0
    IntrinsicStimuli: Species = 0  # Add intrinsic stimuli

    class Smac(AlbeckAsMatlab.Smac):
        M: Species[Override] = 1e5

    class CytoC(AlbeckAsMatlab.CytoC):
        M: Species[Override] = 1e5

    class Apaf(AlbeckAsMatlab.Apaf):
        I: Species[Override] = 1e3  # noqa: E741

    C9: Species[Override] = 0
    XIAP: Species[Override] = 1e4

    sCas3 = Sensor(dimer=7.5e5)
    sCas9 = Sensor(dimer=7.5e5)
    sCas8 = Sensor(dimer=7.5e5)

    # Apoptosome formation
    #   aApaf + pC3 <-->  aApaf:pC3 --> aApaf + C3
    #   C3 + aApaf <-->  C3:aApaf --> C3 + Apop
    Apaf_activates_C3 = MichaelisMenten(
        E=Apaf.A,
        S=C3.pro,
        ES=0,
        P=C3.A,
        forward_rate=5e-9,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    Apaf_and_C3_form_Apop = MichaelisMenten(
        E=C3.A,
        S=Apaf.A,
        ES=0,
        P=Apop,
        forward_rate=1.3e-6,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    # Apoptosome-related inhibitors
    #   Apaf + XIAP <-->  Apaf:XIAP
    XIAP_inhibits_Apaf = ReversibleSynthesis(
        A=Apaf.A,
        B=XIAP,
        AB=0,
        forward_rate=2e-6,
        reverse_rate=KR,
    )

    # Sensors
    C3_cleaves_sCas3 = MichaelisMenten(
        E=C3.A,
        S=sCas3.dimer,
        ES=0,
        P=2 * sCas3.monomer,
        forward_rate=2 * 2.8e-7,
        reverse_rate=1e-2,
        catalytic_rate=KC,
    )
    Apop_cleaves_sCas9 = MichaelisMenten(
        E=Apop,
        S=sCas9.dimer,
        ES=0,
        P=2 * sCas9.monomer,
        forward_rate=2 * 2.8e-7,
        reverse_rate=KR,
        catalytic_rate=KC,
    )
    C8_cleaves_sCas8 = MichaelisMenten(
        E=C8.A,
        S=sCas8.dimer,
        ES=0,
        P=2 * sCas8.monomer,
        forward_rate=2 * 5.4e-8,
        reverse_rate=KR,
        catalytic_rate=KC,
    )

    # Add interaction between single activation
    Apaf_cleaves_sCas9 = MichaelisMenten(
        E=Apaf.A,
        S=sCas9.dimer,
        ES=0,
        P=2 * sCas9.monomer,
        forward_rate=2 * 2e-10,
        reverse_rate=KR,
        catalytic_rate=KC,
    )

    IntrinsicStimuli_activates_Bid = MichaelisMenten(
        E=IntrinsicStimuli,
        S=Bid.U,
        ES=0,
        P=Bid.T,
        forward_rate=KF,
        reverse_rate=KR,
        catalytic_rate=KC,
    )


# Remove some ARM reactions
for name in [
    "Apaf_and_C9_to_Apop",  # aApaf + pC9 <-->  Apop
    "XIAP_inhibits_Apop",  # Apop + XIAP <-->  Apop:XIAP
]:
    del ARM._contents[name]


class ARM_extrinsic(ARM):
    L: Species[Override] = 1e3


class ARM_intrinsic(ARM):
    IntrinsicStimuli: Species[Override] = 1e2
