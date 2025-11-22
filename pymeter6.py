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

baseQuantities = {
    'Time':'s',
    'Length':'m',
    'Mass':'kg',
    'Current':'A',
    'Tempurature':'K',
    'Amount':'mol',
    'LuminousIntensity':'cd'
}

quantities = {
    # Quantities with named units
    'Frequency':'Time^-1',
    'Angle':'Length / Radius',
    'SolidAngle':'Area / Radius^2',
    'Force':'Mass Acceleration',
    'Pressure':'Force / Area',
    'Energy':'Force / Length',
    'Power':'Energy / Time',
    'Charge':'Current Time',
    'Voltage':'Power / Current',
    'Capacitance':'Charge / Voltage',
    'Resistance':'Voltage / Current',
    'Impedance':'Voltage / Current',
    'Reactance':'Voltage / Current',
    'Conductance':'Current / Voltage',
    'MagneticFlux':'Voltage Time',
    'MagneticFluxDensity':'MagneticFlux / Area',
    'Inductance':'MagneticFlux / Current',
    'LuminousFlux':'LuminousIntensity / SolidAngle',
    'Illuminance':'LuminousFlux / Area',
    'Radioactivity':'Time^-1',
    'AbsorbedRadiationDose':'Energy / Mass',
    'EquivalentRadiationDose':'Energy / Mass',
    'Entropy':'Energy / Tempurature',

    # Kinematics
    'Velocity':'Length / Time',
    'Acceleration':'Velocity / Time',
    'Jerk':'Acceleration / Time',
    'Snap':'Jerk / Time',
    'Yank':'Mass Jerk',
    'AngularVelocity':'Angle / Time',
    'AngularAcceleration':'AngularVelocity / Time',
    'FrequencyDrift':'Frequency / Time',
    'VolumetricFlow':'Volume / Time',

    # Mechanics
    'Area':'Length^2',
    'Volume':'Length Area',
    'Momentum':'Force Time',
    'AngularMomentum':'Momentum Radius',
    'Torque':'Force Radius',
    'Wavenumber':'Length^-1',
    'OpticalPower':'Length^-1',
    'Curvature':'Radius^-1',
    'SpatialFrequency':'Length^-1',
    'AreaDensity':'Mass / Area',
    'Density':'Mass / Volume',
    'SpecificVolume':'Volume / Mass',
    'Action':'Energy Time',
    'SpecificEnergy':'Energy / Mass',
    'EnergyDensity':'Energy / Volume',
    'SurfaceTension':'Force / Length',
    'Stiffness':'Energy / Area',
    'Irradiance':'Power / Area',
    'KinematicViscosity':'Area / Time',
    'Diffusivity':'Area / Time',
    'DynamicViscosity':'Pressure Time',
    'LinearMassDensity':'Mass / Length',
    'MassFlowRate':'Mass / Time',
    'Radiance':'Power / SolidAngle Area',
    'SpectralPower':'Power / Length',
    'AbsorbedDoseRate':'AbsorbedDose / Time',
    'FuelEfficiency':'Length / Volume',
    'SpectralIrradiance':'Power / Volume',
    'PowerDensity':'Power / Volume',
    'EnergyFluxDensity':'Energy / Area Time',
    'Compressibility':'Pressure^-1',
    'RadiantExposure':'Energy / Area',
    'InertialMoment':'Mass Radius^2',
    'SpecificAngularMomentum':'AngularMoment / Mass',
    'RadiantIntensity':'Power / SolidAngle',
    'SpectralIntensity':'RadiantIntensity / Length',

    # Chemistry
    'Molarity':'Amount / Volume',
    'MolarVolume':'Volume / Amount',
    'MolarHeatCapacity':'Entropy / Amount',
    'MolarEnergy':'Energy / Amount',
    'MolarConductivity':'Conductivity Length / Amount',
    'Molality':'Amount / Mass',
    'MolarMass':'Mass / Amount',
    'CatalyticEfficiency':'MolarVolume / Time',

    # Electromagnetics
    'PolarizationDensity':'Charge / Area',
    'ChargeDensity':'Charge / Volume',
    'CurrentDensity':'Current / Area',
    'Conductivity':'Conductance / Length',
    'Permittivity':'Capacitance / Length',
    'MagneticPermeability':'Inductance / Length',
    'ElectricFieldStrength':'Voltage / Length',
    'Magnetization':'Current / Length',
    'Exposure':'Charge / Mass',
    'Resistivity':'Resistance Length',
    'LinearChargeDensity':'Charge / Length',
    'MagneticDipoleMoment':'Energy / MagneticFluxDensity',
    'ElectronMobility':'Area / MagneticFlux',
    'Reluctance':'Inductance^-1',
    'MagneticVectorPotential':'MagneticFlux / Length',
    'MagneticMoment':'MagneticFlux Length',
    'MagneticRigidity':'MagneticFluxDensity Length',
    'MagnetomotiveForce':'Current Angle',
    'MagneticSusceptibility':'Length / Inductance',

    # Photometry
    'LuminousEnergy':'LuminousFlux Time',
    'LuminousExposure':'Illuminance Time',
    'Luminance':'LuminousIntensity / Area',
    'LuminousEfficacy':'LuminousFlux / Power',

    # Thermodynamics
    'HeatCapacity':'Energy / Tempurature',
    'SpecificHeatCapacity':'HeatCapacity / Mass',
    'ThermalConductance':'Power / Tempurature',
    'ThermalConductivity':'ThermalConductance / Length',
    'ThermalResistance':'Tempurature / Power',
    'ThermalResistivity':'ThermalResistance / Length',
    'ThermalExpansion':'Tempurature^-1',
    'TempuratureGradiant':'Tempurature / Length',
}

exceptions = {

}

aliases = {

}

inferedPriority = {

}

nullVector = numpy.ndarray((len(baseUnitNames),), fractions.Fraction)
nullVector[:] = fractions.Fraction(0)

def formatVector(vector: numpy.ndarray) -> numpy.ndarray:
    newVector = numpy.ndarray((len(baseUnitNames),), fractions.Fraction)
    newVector[:] = [fractions.Fraction(exp) for exp in vector]
    return newVector

def assembleVector(fDim: str) -> numpy.ndarray:
    vector = numpy.ndarray((len(baseUnitNames),), fractions.Fraction)
    vector[:] = fractions.Fraction(0)
    numer = fDim.split(' / ')
    if len(numer) == 2:
        denom = numer[1]
        for term in denom.split(' '):
            unit = term.split('^')
            if len(unit) == 2:
                exp = fractions.Fraction(unit[1])
            else:
                exp = 1
            unit = unit[0]
            vector -= units[unit].vector * exp
    numer = numer[0]
    for term in numer.split(' '):
        unit = term.split('^')
        if len(unit) == 2:
            exp = fractions.Fraction(unit[1])
        else:
            exp = 1
        unit = unit[0]
        vector += units[unit].vector * exp
    return vector

def assembleScalar(fDim: str) -> float:
    scalar = 1
    numer = fDim.split(' / ')
    if len(numer) == 2:
        denom = numer[1]
        for term in denom.split(' '):
            unit = term.split('^')
            if len(unit) == 2:
                exp = fractions.Fraction(unit[1])
            else:
                exp = 1
            unit = unit[0]
            scalar /= units[unit].scalar * exp
    numer = numer[0]
    for term in numer.split(' '):
        unit = term.split('^')
        if len(unit) == 2:
            exp = fractions.Fraction(unit[1])
        else:
            exp = 1
        unit = unit[0]
        scalar *= units[unit].scalar * exp
    return scalar

def getDistance(vector1: numpy.ndarray, vector2: numpy.ndarray) -> fractions.Fraction:
    distance = sum([abs(vector1[i] - vector2[i]) for i in range(len(vector1))])
    return distance

def disassemble(vector: numpy.ndarray) -> str:
    remainder = vector
    numer = {}
    denom = {}
    while getDistance(remainder, nullVector) != fractions.Fraction(0):
        closestDistance = 0
        closestUnit = ''
        for unit in units.keys():
            if closestUnit == '':
                closestUnit = unit
                closestDistance = getDistance(remainder, units[unit].vector)
                distance = getDistance(remainder, -units[unit].vector)
                if closestDistance > distance:
                    closestUnit = f'-{unit}'
                    closestDistance = distance
            else:
                distance = getDistance(remainder, units[unit].vector)
                if closestDistance > distance:
                    closestDistance = distance
                    closestUnit = unit
                distance = getDistance(remainder, -units[unit].vector)
                if closestDistance > distance:
                    closestDistance = distance
                    closestUnit = f'-{unit}'
            if closestDistance == fractions.Fraction(0):
                break
        if closestUnit[0] == '-':
            if closestUnit[1:] in denom.keys():
                denom[closestUnit[1:]] += 1
            else:
                denom[closestUnit[1:]] = 1
            remainder += units[closestUnit[1:]].vector
        else:
            if closestUnit in numer.keys():
                numer[closestUnit] += 1
            else:
                numer[closestUnit] = 1
            remainder -= units[closestUnit].vector
    # ' / '.join([' '.join([term if numer[term] == 1 else f'{term}^{numer[term]}' for term in numer.keys()]), ' '.join([term if denom[term] == 1 else f'{term}^{denom[term]}' for term in denom.keys()])]) if ' '.join([term if denom[term] == 1 else f'{term}^{denom[term]}' for term in denom.keys()]) != {} else ' '.join([term if numer[term] == 1 else f'{term}^{numer[term]}' for term in numer.keys()])
    numerator = ' '.join([term if numer[term] == 1 else f'{term}^{numer[term]}' for term in numer.keys()])
    denominator = ' '.join([term if denom[term] == 1 else f'{term}^{denom[term]}' for term in denom.keys()])
    return ' / '.join([numerator, denominator]) if denom != {} else numerator

def interlace(self) -> str:
    pass

def compile(mother, father, op) -> str:
    pass

def categorize(self) -> str:
    pass

class Unit():
    def __init__(self, rep):
        self.givenRep = rep
        if isinstance(rep, list):
            if isinstance(rep[0], Unit) and isinstance(rep[1], str):
                self.repType = 'child'
            else:
                self.repType = 'vector'
                self.givenRep = formatVector(rep)
        elif isinstance(rep, numpy.ndarray):
            self.repType = 'vector'
            if not rep.dtype == fractions.Fraction:
                self.givenRep = formatVector(rep)
        elif isinstance(rep, str):
            self.repType = 'fDim'
    def __getattr__(self, name):
        match self.repType:
            case 'vector':
                match name:
                    case 'vector':
                        attr = self.givenRep
                    case 'fDim':
                        attr = disassemble(self.vector)
                    case 'scalar':
                        attr = 1
                    case 'quantity':
                        attr = categorize(self)
            case 'fDim':
                match name:
                    case 'vector':
                        attr = assembleVector(self.givenRep)
                    case 'fDim':
                        attr = self.givenRep
                    case 'scalar':
                        attr = assembleScalar(self.givenRep)
                    case 'quantity':
                        attr = categorize(self)
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
                                attr = interlace(self)
                            case 'scalar':
                                attr = mother.scalar * father.scalar
                            case 'quantity':
                                attr = compile(mother, father, op)
                    case 'div':
                        match name:
                            case 'vector':
                                attr = mother.vector - father.vector
                            case 'fDim':
                                attr = interlace(self)
                            case 'scalar':
                                attr = mother.scalar / father.scalar
                            case 'quantity':
                                attr = compile(mother, father, op)
                    case 'pow':
                        match name:
                            case 'vector':
                                attr = mother.vector * father
                            case 'fDim':
                                attr = interlace(self)
                            case 'scalar':
                                attr = mother ** father
                            case 'quantity':
                                attr = compile(mother, father, op)
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

baseUnits = {}
for baseUnitName, i in zip(baseUnitNames, range(len(baseUnitNames))):
    baseUnits[baseUnitName] = [fractions.Fraction(0),] * len(baseUnitNames)
    baseUnits[baseUnitName][i] = fractions.Fraction(1)

units = {}
for unit in baseUnits.keys():
    units[unit] = Unit(baseUnits[unit])
for unit in namedUnits.keys():
    units[unit] = Unit(namedUnits[unit])

# def quantityHandler(parentClass):
#     def wrapper(*rep):
#         self = Value(*rep)
#         return globals()[self.quantity](*rep)

class Value():
    def __init__(self, *rep):
        self.givenRep = rep
        if isinstance(rep[0], list):
            self.repType = 'child'
        elif isinstance(rep[1], Unit):
            self.repType = 'Unit'
        elif isinstance(rep[1], str):
            self.repType = 'fDim'
    def __getattr__(self, name):
        match self.repType:
            case 'Unit':
                match name:
                    case 'value':
                        attr = self.givenRep[0]
                    case 'unit':
                        attr = self.givenRep[1]
                    case 'quantity':
                        attr = self.unit.quantity
            case 'fDim':
                match name:
                    case 'value':
                        attr = self.givenRep[0]
                    case 'unit':
                        attr = Unit(self.givenRep[1])
                    case 'quantity':
                        attr = self.unit.quantity
            case 'child':
                mother = self.givenRep[0][0]
                father = self.givenRep[0][2]
                op = self.givenRep[0][1]
                match op:
                    case 'sum':
                        match name:
                            case 'value':
                                if mother.unit.scalar == father.unit.scalar:
                                    attr = mother.value + father.value
                                elif mother.unit.scalar < father.unit.scalar:
                                    attr = (father.unit.scalar / mother.unit.scalar) * mother.value + father.value
                                elif mother.unit.scalar > father.unit.scalar:
                                    attr = mother.value + (mother.unit.scalar / father.unit.scalar) * father.value
                            case 'unit':
                                if mother.unit.scalar < father.unit.scalar:
                                    attr = father.unit
                                else:
                                    attr = mother.unit
                            case 'quantity':
                                attr = mother.quantity
                    case 'sub':
                        match name:
                            case 'value':
                                if mother.unit.scalar == father.unit.scalar:
                                    attr = mother.value - father.value
                                elif mother.unit.scalar < father.unit.scalar:
                                    attr = (father.unit.scalar / mother.unit.scalar) * mother.value + father.value
                                elif mother.unit.scalar > father.unit.scalar:
                                    attr = mother.value + (mother.unit.scalar / father.unit.scalar) * father.value
                            case 'unit':
                                if mother.unit.scalar < father.unit.scalar:
                                    attr = father.unit
                                else:
                                    attr = mother.unit
                            case 'quantity':
                                attr = mother.quantity
                    case 'mul':
                        if isinstance(father, int):
                            match name:
                                case 'value':
                                    attr = mother.value * father
                                case 'unit':
                                    attr = mother.unit
                                case 'quanitity':
                                    attr = mother.quantity
                        else:
                            match name:
                                case 'value':
                                    attr = mother.value * father.value
                                case 'unit':
                                    attr = mother.unit * father.unit
                                case 'quantity':
                                    attr = self.unit.quantity
                    case 'div':
                        if isinstance(father, int):
                            match name:
                                case 'value':
                                    attr = mother.value / father
                                case 'unit':
                                    attr = mother.unit
                                case 'quantity':
                                    attr = mother.quantity
                        else:
                            match name:
                                case 'value':
                                    attr = mother.value / father.value
                                case 'unit':
                                    attr = mother.unit / father.unit
                                case 'quantity':
                                    attr = self.unit.quantity
                    case 'pow':
                        match name:
                            case 'value':
                                attr = mother.value ** father
                            case 'unit':
                                attr = Unit([mother.unit, 'pow', father])
                            case 'quantity':
                                attr = self.unit.quantity
        self.__dict__[name] = attr
        return attr
    def __str__(self):
        return f'{self.value} {self.unit}'
    def __repr__(self):
        return f'u<{self.value} v<{self.unit}>>'
    def __add__(self, other):
        if self.quantity == other.quantity:
            return globals()[self.quantity]([self, 'sum', other])
        else:
            raise TypeError()
    def __sub__(self, other):
        if self.quantity == other.quantity:
            return globals()[self.quantity]([self, 'sub', other])
        else:
            raise TypeError()
    def __mul__(self, other):
        return Value([self, 'mul', other])
    def __rmul__(self, other):
        return self.__mul__(other)
    def __truediv__(self, other):
        return Value([self, 'div', other])
    def __rtruediv__(self, other):
        return self.__truediv__(other)
    def __pow__(self, other):
        return Value([self, 'pow', other])
    def to(self, target: str):
        pass

quantityClasses = {}

for quantity in quantities.keys():
    globals()[quantity] = type(quantity, (Value,), {"quantity": quantity})

# class Tree(list):
#     def __get

class Variable():
    def __init__(self, root):
        self.root = root
    # def inject(self, other):
    #     if self.name not in self.root
    def __call__(self, *args, **kwargs):
        pass
    def __add__(self, other):
        pass
    def __radd__(self, other):
        pass
    def __sub__(self, other):
        pass
    def __rsub__(self, other):
        pass
    def __mul__(self, other):
        pass
    def __rmul__(self, other):
        pass
    def __truediv__(self, other):
        pass
    def __rtruediv__(self, other):
        pass
    def __pow__(self, other):
        pass
    def __rpow__(self, other):
        pass
    def __eq__(self, other):
        pass
    def __neg__(self):
        pass

class Expression():
    def __init__(self, equation):
        self.equation = equation
        self.code = []
        self.vars = {}
        self.lookup = {}
        vars = self.equation.__code__.co_names
        for i in range(len(vars)):
            self.vars[vars[i]] = Variable(self)
            self.vars[vars[i]].name = vars[i]
            self.equation.__globals__[vars[i]] = self.vars[vars[i]]
        self.equation()
    def __call__(self, *args, **kwargs):
        pass
        # if len(self.vars):
        #     for var, arg in self.vars.keys(), self.vars.items():
        #         self.equation.__globals__[var] = arg
        # if len(kwargs):
        #     for var, arg in kwargs.keys(), kwargs.items():
        #         self.equation.__globals__[var] = arg
        # if len(args):
        #     for var, arg in self.args, args:
        #         self.equation.__globals__[var] = arg
        # return self.equation()

# @Expression
# def foo():
#     # Equivalent to 
#     # y == a * x ** 4 + b * x ** 2 + c * x + d
#     (y - d) / x == x * (a * x ** 2 + b) + c

# [[0, 'sub', 1], 'div', 3]
# [[3, 'mul', [[6, 'mul', [3, 'pow', 4]], 'add', 8]], 'add', 12]
# [[[0, 'sub', 1], 'div', 3], 'eql', [[3, 'mul', [[6, 'mul', [3, 'pow', 4]], 'add', 8]], 'add', 12]]
# (foo - foo.c)

# print(foo())