# Library written by Noah Dye,
# Copyright 2025

# Most Units defined here come from the wikipedia article on SI Units
# https://en.wikipedia.org/wiki/International_System_of_Units
# https://en.wikipedia.org/wiki/SI_derived_unit

# No, this library wasn't written by AI, I'm just an idiot

import numpy, fractions

baseUnitNames = ['s', 'm', 'kg', 'A', 'K', 'mol', 'cd']

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

quantities = { # <------ NEEDS REWRITE
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

inferedPriority = {}

def sortUnits(*units):
    pass

baseUnits = {}
for baseUnitName, i in zip(baseUnitNames, range(len(baseUnitNames))):
    baseUnits[baseUnitName] = [fractions.Fraction(0),] * len(baseUnitNames)
    baseUnits[baseUnitName][i] = fractions.Fraction(1)

class Unit():
    def __init__(self, rep):
        self.givenRep = rep
        if isinstance(rep, list):
            if len(rep) == 3:
                self.repType = 'child'
            else:
                self.repType = 'vector'
                for i in range(len(rep)):
                    self.givenRep[i] = fractions.Fraction(rep[i])
        elif isinstance(rep, numpy.ndarray):
            self.repType = 'vector'
            if rep.dtype != fractions.Fraction:
                for i in range(len(rep)):
                    self.givenRep[i] = fractions.Fraction(rep[i])
        elif isinstance(rep, str):
            self.repType = 'fDim'
        elif isinstance(rep, Unit):
            self = rep
    def disassemble(self):
        remainingUnit = self.vector
        terms = {}
        while remainingUnit != nullUnit.vector:
            closestUnit = ''
            closestDistance = 0
            for unit in units.keys():
                distance = sum(abs(remainingUnit - units[unit].vector))
                if closestUnit == '':
                    closestUnit = unit
                    closestDistance = distance
                else:
                    if distance < closestDistance:
                        closestUnit = unit
                        closestDistance = distance
                distance = sum(abs(remainingUnit + units[unit].vector))
                if closestUnit == '':
                    closestUnit = f'-{unit}'
                    closestDistance = distance
                else:
                    if distance < closestDistance:
                        closestUnit = f'-{unit}'
                        closestDistance = distance
            remainingUnit -= units[unit].vector if closestUnit[0] != '-' else -units[unit].vector
            if unit not in terms.keys():
                terms[unit] = 1 if closestUnit[0] != '-' else -1
            else:
                terms[unit] += 1 if closestUnit[0] != '-' else -1
        numer = []
        denom = []
        for term in terms.keys():
            if terms[term] == 1:
                numer += [term]
            elif terms[term] == -1:
                denom += [term]
            elif terms[term] > 1:
                numer += [f'{term}^{terms[term]}']
            elif terms[term] < 1:
                denom += [f'{term}^{-terms[term]}']
        return ' / '.join([sortUnits(numer), sortUnits(denom)])
    def __getattr__(self, name):
        match self.repType:
            case 'vector':
                match name:
                    case 'vector':
                        attr = self.givenRep
                    case 'fDim':
                        attr = self.disassemble()
                    case 'scalar':
                        attr = 1
                    case 'quantity':
                        pass
            case 'fDim':
                match name:
                    case 'vector':
                        attr = nullUnit.vector
                        numerator = self.givenRep.split(' / ')
                        if len(numerator) == 2:
                            denominator = numerator[1].split(' ')
                            numerator = numerator[0].split(' ')
                            for term in denominator:
                                unit = term.split('^')
                                if len(unit) == 2:
                                    exp = fractions.Fraction(unit[1])
                                else:
                                    exp = 1
                                unit = unit[0]
                                attr -= units[unit].vector * exp
                        for term in numerator:
                            unit = term.split('^')
                            if len(unit) == 2:
                                exp = fractions.Fraction(unit[1])
                            else:
                                exp = 1
                            unit = unit[0]
                            attr += units[unit].vector * exp
                    case 'fDim':
                        attr = self.givenRep
                    case 'scalar':
                        attr = 1
                        numerator = self.givenRep.split(' / ')
                        if len(numerator) == 2:
                            denominator = numerator[1].split(' ')
                            numerator = numerator[0].split(' ')
                            for term in denominator:
                                unit = term.split('^')
                                if len(unit) == 2:
                                    exp = fractions.Fraction(unit[1])
                                else:
                                    exp = 1
                                unit = unit[0]
                                attr /= units[unit].scalar * exp
                        for term in numerator:
                            unit = term.split('^')
                            if len(unit) == 2:
                                exp = fractions.Fraction(unit[1])
                            else:
                                exp = 1
                            unit = unit[0]
                            attr *= units[unit].scalar * exp
                    case 'quantity':
                        pass
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
                                attr = self.disassemble()
                            case 'scalar':
                                attr = mother.scalar * father.scalar
                            case 'quantity':
                                pass
                    case 'div':
                        match name:
                            case 'vector':
                                attr = mother.vector - father.vector
                            case 'fDim':
                                attr = self.disassemble()
                            case 'scalar':
                                attr = mother.scalar * father.scalar
                            case 'quantity':
                                pass
                    case 'pow':
                        match name:
                            case 'vector':
                                attr = mother.vector * father
                            case 'fDim': # Could be written a little bit better
                                attr = ''
                                numerator = mother.fDim.split(' / ')
                                numer = []
                                denom = []
                                hasDenom = len(numerator) == 2
                                if hasDenom:
                                    denominator = numerator[1].split(' ')
                                    for term in denominator:
                                        unit = term.split('^')
                                        if len(unit) == 2:
                                            exp = fractions.Fraction(unit[1]) * father
                                        else:
                                            exp = fractions.Fraction(1) * father
                                        unit = unit[0]
                                        if exp > 0:
                                            denom += [f'{unit}^{exp}']
                                        elif exp < 0:
                                            numer += [f'{unit}^{-exp}']
                                for term in numerator:
                                    unit = term.split('^')
                                    if len(unit) == 2:
                                        exp = fractions.Fraction(unit[1]) * father
                                    else:
                                        exp = fractions.Fraction(1) * father
                                    unit = unit[0]
                                    if exp > 0:
                                        numer += [f'{unit}^{exp}']
                                    elif exp < 0:
                                        denom += [f'{unit}^{-exp}']
                                attr = ' / '.join([' '.join(numer), ' '.join(denom)])
                            case 'scalar':
                                attr = mother.scalar ** father
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

nullUnit = Unit([0] * len(baseUnitNames))

units = {}

class Value():
    def __init__(self, *rep):
        self.givenRep = rep
        if len(rep) == 2:
            self.repType = 'parent'
        elif len(rep) == 3:
            self.repType = 'child'
    def __getattr__(self, name):
        match self.repType:
            case 'parent':
                match name:
                    case 'value':
                        attr = self.givenRep[0]
                    case 'unit':
                        attr = Unit(self.givenRep[1])
            case 'child':
                mother = self.givenRep[0]
                father = self.givenRep[2]
                op = self.givenRep[1]
                match op:
                    case 'sum':
                        pass
                    case 'sub':
                        pass
                    case 'mul':
                        pass
                    case 'div':
                        pass
                    case 'pow':
                        pass
        self.__dict__[name] = attr
    def __str__(self):
        return f'{self.value} {self.unit}'
    def __repr__(self):
        return f'v<{self.value} {self.unit.__repr__()}>'
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
    def to(self, target):
        pass

class Tensor():
    def __init__(self, data):
        pass
    def __getattr__(self, name):
        pass
    def __sum__(self, other):
        pass
    def __sub__(self, other):
        pass
    def __mul__(self, other):
        pass
    def __truediv__(self, other):
        pass
    def to(self, target):
        pass

class Expression:
    def __init__(self, expr: function):
        self.expr = expr
        self.stdCall = {}
        self.vars = {}
    def __call__(self, *args, **kwargs):
        for var, val in zip(self.stdCall.keys(), self.stdCall.items()):
            self.expr.__globals__[var] = val
        return self.expr(*self.vars)

