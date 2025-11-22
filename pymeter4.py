import fractions, numpy

# Library written by Noah Dye,
# Copyright 2025

# Most Units defined here come from the wikipedia article on SI Units
# https://en.wikipedia.org/wiki/International_System_of_Units
# https://en.wikipedia.org/wiki/SI_derived_unit

baseUnitNames = ["s", "m", "kg", "A", "K", "mol", "cd"]

namedUnits = {
    "Hz":"s^-1",
    "r":"m / m",
    "sr":"m^2 / m^2",
    "N":"kg m / s^2",
    "Pa":"N / m^2",
    "J":"N m",
    "W":"J / s",
    "C":"A s",
    "V":"W / A",
    "F":"C / V",
    "Ω":"V / A",
    "S":"A / V",
    "Wb":"V s",
    "T":"Wb / m^2",
    "H":"Wb / A",
    "lm":"cd sr",
    "lx":"lm / m^2",
    "Bq":"s^-1",
    "Gy":"m^2 / s^2",
    "Sv":"m^2 / s^2",
    "kat":"mol / s",
    "nat":"J / K"
}

quantities = {
    # Base Units
    "Time":"s",
    "Distance":"m",
    "Mass":"kg",
    "ElectricalCurrent":"A",
    "Tempurature":"K",
    "Amount":"mol",
    "LuminousIntensity":"cd",

    # Named Units
    "Frequency":"Hz",
    "Angle":"r",
    "SolidAngle":"sr",
    "Force":"N",
    "Pressure":"Pa",
    "Stress":"Pa",
    "Energy":"J",
    "Power":"W",
    "RadiantFlux":"W",
    "Charge":"C",
    "Voltage":"V",
    "Capacitance":"F",
    "Resistance":"Ω",
    "Reactance":"Ω",
    "Impedance":"Ω",
    "Conductance":"S",
    "MagneticFlux":"Wb",
    "MagneticFluxDensity":"T",
    "Inductance":"H",
    "LuminousFlux":"lm",
    "Illuminance":"lx",
    "Activity":"Bq",
    "Kerma":"Gy",
    "RadiationDose":"Sv",
    "CataliticActivity":"kat",
    "Entropy":"nat",

    #Kinematics
    "Velocity":"m / s",
    "Acceleration":"m / s^2",
    "Jolt":"m / s^3",
    "Yank":"kg m / s^3",
    "Snap":"m / s^4",
    "AngularVelocity":"r / s",
    "AngularAcceleration":"r / s^2",
    "FrequencyDrift":"Hz/s",
    "VolumetricFlow":"m^3 / s",

    #Mechanics
    "Area":"m^2",
    "Volume":"m^3",
    "Momentum":"N s",
    "AngularMomentum":"N m s",
    "Torque":"N m",
    "Yank":"N / s",
    "Wavenumber":"m^-1",
    "OpticalPower":"m^-1",
    "Curvature":"m^-1",
    "SpatialFrequency":"m^-1",
    "AreaDensity":"kg / m^2",
    "Density":"kg / m^3",
    "SpecificVolume":"m^3 / kg",
    "Action":"J s",
    "SpecificEnergy":"J / kg",
    "EnergyDensity":"J / m^3",
    "Stiffness":"N / m",
    "Irradiance":"W / m^2",
    "KinematicViscosity":"m^2 / s",
    "Diffusivity":"m^2 / s",
    "DynamicViscosity":"Pa s",
    "LinearMassDensity":"kg / m",
    "MassFlowRate":"kg / s",
    "Radiance":"W / sr m^2",
    "VolumetricRadiance":"W / sr m^3",
    "SpectralPower":"W / m",
    "AbsorbedDoseRate":"Gy / s",
    "FuelEfficiency":"m / m^3",
    "SpectralIrradiance":"W / m^3",
    "PowerDensity":"W / m^3",
    "EnergyFluxDensity":"J / m^2 s",
    "Compressibility":"Pa^-1",
    "RadiantExposure":"J / m^2",
    "InertialMoment":"kg m^2",
    "SpecificAngularMomentum":"N m s / kg",
    "RadiantIntensity":"W / sr",
    "SpectralIntensity":"W / sr m",

    #Chemistry
    "Molarity":"mol / m^3",
    "MolarVolume":"m^3 / mol",
    "MolarHeatCapacity":"J / K mol",
    "MolarEntropy":"J / K mol",
    "MolarEnergy":"J / mol",
    "MolarConductivity":"S m^2 / mol",
    "Molality":"mol / kg",
    "MolarMass":"kg / mol",
    "CatalyticEfficiency":"m^3 / mol s",

    #Electromagnetics
    "ElectricDisplacement":"C / m^2",
    "PolarizationDensity":"C / m^2",
    "ElectricCurrentDensity":"A / m^2",
    "ElectricalConductivity":"S / m",
    "Permitivity":"F / m",
    "MagneticPermeability":"H / m",
    "ElectricFieldStrength":"V / m",
    "Magnetization":"A / m",
    "Exposure":"C / kg",
    "Resistivity":"Ω m",
    "LinearChargeDensity":"C / m",
    "MagneticDipoleMoment":"J / T",
    "ElectronMobility":"m^2 / V s",
    "Reluctance":"H^-1",
    "MagneticPotential":"Wb / m",
    "MagneticMoment":"Wb m",
    "MagneticRigidity":"T m",
    "Magnetomotivation":"A r",
    "MagneticSusceptibility":"m / H",

    #Photomotry
    "LuminousEnergy":"lm s",
    "LuminousExoisure":"lx s",
    "Luminance":"cd / m^2",
    "LuminousEfficacy":"lm / W",

    #Thermodynamics
    "SpecificHeatCapacity":"J / K kg",
    "ThermalConductivity":"W / m K",
    "ThermalResistance":"K / W",
    "TemperatureGradient":"K / m",

    #Disambiguations
    "Length":"m",
}

aliases = {
    "Ω":["Ohm"],
}

exceptions = {
    "Force Length":"Torque",
    "Distance / Length":"Angle",
    "Area / Length^2":"SolidAngle",
}

prefixes = {
    "q": -30,
    "r": -27,
    "y": -24,
    "z": -21,
    "a": -18,
    "f": -15,
    "p": -12,
    "n": -9,
    "u": -6,
    
    "m": -3,
    "c": -2,
    "d": -1,

    "" : 0,

    "da": 1,
    "h" : 2,
    "k" : 3,
    
    "M" : 6,
    "G" : 9,
    "T" : 12,
    "P" : 15,
    "E" : 18,
    "Z" : 21,
    "Y" : 24,
    "R" : 27,
    "Q" : 30,
}


baseUnits = {}
for baseUnitName, i in zip(baseUnitNames, range(len(baseUnitNames))):
    baseUnits[baseUnitName] = [fractions.Fraction(0),] * len(baseUnitNames)
    baseUnits[baseUnitName][i] = fractions.Fraction(1)

class Unit():
    def __init__(self, rep):
        self.givenRep = rep
        if isinstance(rep, list):
            if len(rep) == 3:
                if isinstance(rep[0], Unit) and isinstance(rep[1], str):
                    self.repType = 'child'
                else:
                    self.repType = 'vector'
                    self.givenRep = numpy.ndarray((len(rep),), fractions.Fraction)
                    for i in range(len(rep)):
                        self.givenRep[i] = fractions.Fraction(rep[i])
            else:
                self.repType = 'vector'
                self.givenRep = numpy.ndarray((len(rep),), fractions.Fraction)
                for i in range(len(rep)):
                    self.givenRep[i] = fractions.Fraction(rep[i])
        elif isinstance(rep, numpy.ndarray):
            self.repType = 'vector'
            if not isinstance(rep.dtype, fractions.Fraction):
                for i in range(len(rep)):
                    self.givenRep[i] = fractions.Fraction(rep[i])
        elif isinstance(rep, str):
            self.repType = 'fDim'
    def __getattr__(self, name):
        match self.repType:
            case 'vector':
                match name:
                    case 'vector':
                        attr = self.givenRep
                    case 'fDim':
                        pass # disassemble
                    case 'scalar':
                        attr = 1
                    case 'quantity':
                        pass # decompile
            case 'fDim':
                match name:
                    case 'vector':
                        numerator = self.givenRep.split(' / ')
                        if len(numerator) == 2:
                            denominator = numerator[1]
                        else:
                            denominator = ''
                        numerator = numerator[0]
                        attr = nullUnit.vector
                        for term in numerator:
                            unit = term.split('^')
                            if len(unit) == 2:
                                exp = fractions.Fraction(unit[1])
                            else:
                                exp = fractions.Fraction(1)
                            unit = unit[0]
                            attr += units[unit].vector * exp
                        for term in denominator:
                            unit = term.split('^')
                            if len(unit) == 2:
                                exp = fractions.Fraction(unit[1])
                            else:
                                exp = fractions.Fraction(1)
                            unit = unit[0]
                            attr -= units[unit].vector * exp
                    case 'fDim':
                        attr = self.givenRep
                    case 'scalar':
                        numerator = self.givenRep.split(' / ')
                        if len(numerator) == 2:
                            denominator = numerator[1]
                        else:
                            denominator = ''
                        numerator = numerator[0]
                        attr = 1
                        for term in numerator:
                            unit = term.split('^')
                            if len(unit) == 2:
                                exp = fractions.Fraction(unit[1])
                            else:
                                exp = fractions.Fraction(1)
                            unit = unit[0]
                            attr *= units[unit].scalar * exp
                        for term in denominator:
                            unit = term.split('^')
                            if len(unit) == 2:
                                exp = fractions.Fraction(unit[1])
                            else:
                                exp = fractions.Fraction(1)
                            unit = unit[0]
                            attr /= units[unit].scalar * exp
                    case 'quantity':
                        pass # decompile
            case 'child':
                mother = self.givenRep[0]
                father = self.givenRep[2]
                op = self.givenRep[1]
                match op:
                    case 'mul':
                        match name:
                            case 'vector':
                                attr = mother.vector + father.vector
                            case 'fDim':
                                pass # disassemble
                            case 'scalar':
                                pass
                            case 'quantity':
                                pass
                    case 'div':
                        match name:
                            case 'vector':
                                attr = mother.vector - father.vector
                            case 'fDim':
                                pass
                            case 'scalar':
                                pass
                            case 'quantity':
                                pass
                    case 'pow':
                        match name:
                            case 'vector':
                                attr = mother.vector * father
                            case 'fDim':
                                pass
                            case 'scalar':
                                attr = mother.scalar * father
                            case 'quantity':
                                pass
        self.__dict__[name] = attr
        return attr
    def __str__(self):
        return self.fDim
    def __repr__(self):
        return f'u<{self}>'
    def __mul__(self, other):
        return Unit([self, 'mul', other])
    def __truediv__(self, other):
        return Unit([self, 'div', other])
    def __pow__(self, other):
        return Unit([self, 'pow', other])

units = {}
for unit in baseUnits.keys():
    units[unit] = Unit(baseUnits[unit])
for unit in namedUnits.keys():
    units[unit] = Unit(namedUnits[unit])
for unit in units:
    units[unit].givenRep = units[unit].givenRep[1]

class UnequalUnitArithmeticError(Exception):
    pass

nullUnit = Unit([0] * len(baseUnitNames))

class Value(object):
    def __init__(self, *args):
        self.args = args
        self.isChild = len(args) == 1
    def __getattr__(self, name):
        #Superposition Collapse
        if self.isChild:
            mother = self.args[0][0]
            father = self.args[0][2]
            match self.args[0][1]:
                case 'sum':
                    match name:
                        case 'value':
                            attr = mother.value + father.value
                        case 'unit':
                            if mother.fDim == father.fDim:
                                attr = mother.unit
                            else:
                                raise UnequalUnitArithmeticError
                case 'sub':
                    match name:
                        case 'value':
                            attr = mother.value - father.value
                        case 'unit':
                            if mother.fDim == father.fDim:
                                attr = mother.unit
                            else:
                                raise UnequalUnitArithmeticError
                case 'mul':
                    match name:
                        case 'value':
                            attr = mother.value * father.value
                        case 'unit':
                            if isinstance(father, Value):
                                attr = mother.unit * father.unit
                            else:
                                attr = mother.unit
                case 'div':
                    match name:
                        case 'value':
                            attr = mother.value / father.value
                        case 'unit':
                            if isinstance(father, Value):
                                attr = mother.unit * father.unit
                            else:
                                attr = mother.unit
                case 'pow':
                    match name:
                        case 'value':
                            attr = mother.value ** father
                        case 'unit':
                            attr = mother.unit ** father
        else:
            match name:
                case 'value':
                    attr = self.args[0]
                case 'unit':
                    if isinstance(self.args[1], Unit):
                        attr = self.args[1]
                    else:
                        attr = Unit(self.args[1])
        self.__dict__[name] = attr
        return attr
    def __sum__(self, other):
        return Value([self, 'sum', other])
    def __sub__(self, other):
        return Value([self, 'sub', other])
    def __mul__(self, other):
        return Value([self, 'mul', other])
    def __truediv__(self, other):
        return Value([self, 'div', other])
    def __pow__(self, other):
        return Value([self, 'pow', other])
    def __str__(self):
        return f'{self.value} {self.unit}'
    def __repr__(self):
        return f'v<{self.value} u<{self.unit}>>'

class Convertion():
    def __init__(self, value, dim):
        self.value = value
        self.dim = dim
    def __mul__(self, other: Value):
        pass