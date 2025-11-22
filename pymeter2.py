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

units = {}

class reverseableDict(dict):
    def __init__(self, Dict: dict):
        self[Dict]
        self.reverseDict = {}
        for key, item in zip(Dict.keys(), Dict.items()):
            self.reverseDict[item] = key
    def __setitem__(self, index, value):
        self[index] = value
        self.reverseDict[value] = index
    def __getitem__(self, key):
        if key in self.keys():
            return self[key]
        elif key in self.items():
            return self.reverseDict[key]
        else:
            raise KeyError()


def combineFDims(fDim1, fDim2, op):
    return None

def combineQuantities(q1, q2, op):
    # Combine the quantities
    return None

def findQuantity(self):
    pass

class Unit():
    def __init__(self, rep):
        self.initDict = {'initDict': True, '__setattribute__': True}
        self.initDict["initDict"] = True
        self.givenRep = rep
        self.initDict["givenRep"] = True
        self.isChild = False
        self.scalar = 1 # Combined scaling of all included units
        self.vector = numpy.ndarray((len(baseUnitNames),), fractions.Fraction) # exponent representation with base units
        self.fDim = '' # The Displayed representation of the unit
        self.quantity = '' # Physical quantity being measured
        # The above is helpful for distinguishing certain units which cancel or decompile incorrectly
        self.simple = False
    def __setattribute__(self, name, value):
        self.__dict__[name] = value
        self.__dict__["initDict"][name] = True
    def __getattribute__(self, name):
        try:
            if not name == "__dict__":
                isInited = self.__dict__["initDict"][name]
            else:
                attr = object.__getattribute__(self, "__dict__")
                isInited = True
                return attr
        except KeyError:
            isInited = False
        if not isInited: # Do if the value has not yet been calculated
            # Collapse "superposition" of the attribute
            if name == "isChild": # Edge case
                # Check if self.givenRep looks like a Child Descriptor
                if isinstance(self.givenRep, list):
                    attr = isinstance(self.givenRep[0], Unit) and isinstance(self.givenRep[1], str)
                else:
                    attr = False
            else:
                if self.isChild:
                    # Infer the value of name from the parents and the operation
                    op = self.givenRep[1]
                    match op:
                        case "mul":
                            father = self.givenRep[0] # Arbitrary names for the two parents
                            mother = self.givenRep[2]
                            match name:
                                case "scalar":
                                    attr = father.scalar * mother.scalar
                                case "vector":
                                    attr = father.vector + mother.vector
                                case "fDim":
                                    # Combine both fDims
                                    attr = combineFDims(father.fDim, mother.fDim)
                                case "quantity":
                                    quantity = combineQuantities(father.quantity, mother.quantity, op)
                                    # Check if the parent's quantities fit an exception
                                    if quantity in exceptions.keys():
                                        # Replace the correct quantity
                                        attr = exceptions[quantity]
                                    else:
                                        # Construct the new quantity based on self.vector
                                        attr = findQuantity(self)
                        case "div":
                            father = self.givenRep[0] # Arbitrary names for the two parents
                            mother = self.givenRep[2]
                            match name:
                                case "scalar":
                                    attr = father.scalar / mother.scalar
                                case "vector":
                                    attr = father.vector - mother.scalar
                                case "fDim":
                                    # Combine fDims
                                    # Recursive method
                                    attr = combineFDims(father.fDim, mother.fDim)
                                case "quantity":
                                    quantity = combineQuantities(father.quantity, mother.quantity, op)
                                    # Check if the parent's quantities fit an exception
                                    if quantity in exceptions.keys():
                                        # Replace the correct quantity
                                        attr = exceptions[quantity]
                                    else:
                                        # Construct the new quantity based on self.vector
                                        attr = findQuantity(self)
                        case "pow":
                            father = self.givenRep[0] # Arbitrary names for the two parents
                            mother = self.givenRep[2]
                            match name:
                                case "scalar":
                                    attr = father.scalar ** mother
                                case "vector":
                                    attr = father.vector * mother
                                case "fDim":
                                    # Combine fDims
                                    attr = combineFDims(father.fDim, mother.fDim, op)
                                case "quantity":
                                    quantity = combineQuantities(father.quantity, mother, op)
                                    # Check if the parent's quantities fit an exception
                                    if quantity in exceptions.keys():
                                        # Replace the correct quantity
                                        attr = exceptions[quantity]
                                    else:
                                        # Construct the new quantity based on self.vector
                                        attr = findQuantity(self)
                else:
                    # Parse self.givenRep
                    if isinstance(self.givenRep, str):
                        if any(self.givenRep.split(' ').remove('/')) == any(quantities.keys()):
                            # Probably a quantity rather than an fDim
                            match name:
                                case "quantity":
                                    attr = self.givenRep
                                case "fDim":
                                    attr = quantities[self.givenRep]
                                case "vector":
                                    fDim = quantities[self.givenRep]
                                    self.__setattribute__("fDim", fDim)
                                    vector = nullUnit.vector
                                    expSign = 1
                                    for term in fDim.split(' '):
                                        if term == '/':
                                            expSign = -1
                                        else:
                                            unit = term.split('^')
                                            if len(unit) - 1:
                                                exp = int(unit[1]) * expSign
                                                unit = unit[0]
                                            vector += units[unit].vector * exp
                                    attr = vector
                                case "scalar":
                                    attr = 1
                        elif any(self.givenRep.split(' ').remove('/')) == any(units.fDim):
                            # Probably an fDim rather than a quantity
                            try:
                                scalar = int(self.givenRep.split(' ')[0])
                            except ValueError():
                                scalar = 1
                            match name:
                                case "fDim":
                                    attr = self.givenRep
                                case "scalar":
                                    try:
                                        attr = int(self.givenRep.split(' ')[0])
                                    except ValueError():
                                        attr = 1
                        else:
                            # Idk what to do with this
                            pass
                    elif isinstance(self.givenRep, numpy.ndarray): # probably vector
                        match name:
                            case "vector": # Parse the vector to make sure its the correct dtype
                                match self.givenRep.dtype:
                                    case fractions.Fraction:
                                        attr = self.givenRep
                                    case type(0):
                                        attr = numpy.ndarray((len(baseUnitNames),), fractions.Fraction)
                                        attr[:] = (fractions.Fraction(i) for i in self.givenRep)
                                    case type(0.0):
                                        attr = numpy.ndarray((len(baseUnitNames),), fractions.Fraction)
                                        attr[:] = (fractions.Fraction(i) for i in self.givenRep)
                            case "fDim":
                                # Find Taxicab distance to closest named unit
                                rep = {}
                                remainder = self.vector
                                numerator = []
                                denominator = []
                                while True:
                                    closestUnit = None
                                    closestDistance = -1
                                    for name, unit in zip(units.keys(), unit.items()):
                                        diff = sum(abs(selfExp + otherExp) for selfExp, otherExp in zip(remainder, unit.vector))
                                        if closestDistance == -1:
                                            closestDistance = diff
                                            closestUnit = name
                                        elif diff < closestDistance:
                                            closestDistance = diff
                                            closestUnit = name
                                        diff = sum(abs(selfExp - otherExp) for selfExp, otherExp in zip(remainder, unit.vector))
                                        if diff < closestDistance:
                                            closestDistance = diff
                                            closestUnit = f'-{name}'
                                    if closestUnit in rep.keys():
                                        rep[closestUnit] += 1
                                    else:
                                        rep[closestUnit] = 1
                                    if closestUnit[0] == '-':
                                        remainder -= units[closestUnit[1:]].vector
                                        denominator += [closestUnit[1:]]
                                    else:
                                        remainder += units[closestUnit].vector
                                        numerator += [closestUnit]
                                    if all(remainder) == all(nullUnit.vector):
                                        break
                                for term, exp in zip(rep.keys, rep.items):
                                    if term[0] == '-':
                                        denominator += [f'{term[1:]}^{exp}']
                                    else:
                                        numerator += [f'{term}^{exp}']
                                # Sort the terms based on examples


                                # then we combine them
                                attr = ' / '.join([' '.join(numerator), ' '.join(denominator)])
                            case "quantity":
                                pass
                    elif isinstance(self.givenRep, list):
                        attr = numpy.ndarray((len(baseUnitNames),), fractions.Fraction)
                        attr[:] = (fractions.Fraction(i) for i in self.givenRep)
                    else:
                        pass
            object.__setattribute__(self, name, attr) # Update the value so we don't recaluclate it on every access
        else: # If the value has already been calculated, just grab it
            attr = self.__dict__[name]
        return attr
    def decompile(self):
        for attr in self.initDict.keys():
            self.__getattribute__(attr)
    def __str__(self):
        return self.fDim
    def __repr__(self):
        return f'u<{self}>'
    def __mul__(self, other):
        newUnit = Unit([self, 'mul', other])
        return newUnit
    def __truediv__(self, other):
        newUnit = Unit([self, 'div', other])
        return newUnit
    def __pow__(self, other):
        return Unit([self, 'pow', other])

# breakpoint()
nullUnit = Unit([0] * len(baseUnitNames))

class UnequalUnitArithmeticError(Exception):
    pass

class Value():
    def __init__(self, value, unit):
        self.value = value
        if not isinstance(unit, Unit):
            self.unit = Unit(unit)
        else:
            self.unit = unit
    def __sum__(self, other):
        if other.unit == self.unit:
            return Value(self.value + other.value, self.unit)
        else:
            raise UnequalUnitArithmeticError
    def __sub__(self, other):
        if other.unit == self.unit:
            return Value(self.value - other.value, self.unit)
        else:
            raise UnequalUnitArithmeticError
    def __mul__(self, other):
        return Value(self.value * other.value, self.unit * other.unit)
    def __truediv__(self, other):
        return Value(self.value / other.value, self.unit / other.unit)
    def __str__(self):
        return f'{self.value} {self.unit}'
    def __repr__(self):
        return f'v<{self.value} u<{self.unit}>>'
    def __pow__(self, other):
        return Value(self.value ** other, self.unit ** other)
    
# class Expression():
#     def __init__(self, exp: str):
#         self.exp = exp
#         self.expectedUnit = None
#     def simplifyUnits(self):
#         pass
#     def eval(self, *args):
#         pass
#     def derivative(self, relativeTerm):
#         terms = self.exp.split(' ')
#     def integrate(self, relativeTerm):
#         pass
# nonSI = {
#     # Time
#     "min":"60 s / min",
#     "hr":"60 min / hr",
#     "d":"24 hr / d",

#     # distance
#     ""

#     # Volume
#     "L":"1 dm^3 / L",
# }