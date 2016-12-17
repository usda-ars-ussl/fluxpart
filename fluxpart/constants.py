"""Physical constants"""

from collections import namedtuple


class MolecularWeight(namedtuple('MolecularWeight', 'dryair vapor co2')):
    """
    Attributes
    ----------
    dryair, vapor, co2 : float
        kg/mol
    """


class SpecificGasConstant(namedtuple('SpecificGasConstant',
                                     'dryair vapor co2')):
    """
    Attributes
    ----------
    dryair, vapor, co2 : float
        J/K/mol == m^3 Pa/kg/K
    """


class SpecificHeatCapacity(namedtuple('SpecificHeatCapacity', 'dryair')):
    """
    Attributes
    ----------
    dryair : float
        J/K/kg
    """


GRAVITY = 9.81             # m/s^2
VON_KARMAN = 0.40          # unitless
UNIVERSAL_GAS = 8.3144598  # J/K/mol

MOLECULAR_WEIGHT = MolecularWeight(
    dryair=0.0289645,
    vapor=0.018016,
    co2=0.044010)

SPECIFIC_GAS_CONSTANT = SpecificGasConstant(
    dryair=UNIVERSAL_GAS / MOLECULAR_WEIGHT.dryair,
    vapor=UNIVERSAL_GAS / MOLECULAR_WEIGHT.vapor,
    co2=UNIVERSAL_GAS / MOLECULAR_WEIGHT.co2)

SPECIFIC_HEAT_CAPACITY = SpecificHeatCapacity(
    dryair=1004.67)
