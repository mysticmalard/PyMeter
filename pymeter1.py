import fractions, numpy



baseUnitNames = ["s", "m", "kg", "A", "K", "mol", "cd", "b"]

baseUnits = {}
for baseUnitName, i in zip(baseUnitNames, range(len(baseUnitNames))):
    baseUnits[baseUnitName] = [fractions.Fraction(0),] * len(baseUnitNames)
    baseUnits[baseUnitName][i] = fractions.Fraction(1)

prefixes = {
    "q": -30,
    "r": -27,
    "y": -27,
    "z": -27,
    "a": -27,
    "f": -27,
    "p": -27,
    "n": -27,
    "u": -27,
    
    "m": -27,
    "c": -27,
    "d": -27,

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

units = {}

nextUnitID = 0



# kg A b cd m K mol s

class Unit(object):
    def getBaseRepresentation(self):
        repList = []
        for exp, unit in zip(self.baseUnitVector, baseUnits.keys()):
            try:
                match int(exp / abs(exp)):
                    case 0:
                        pass
                    case 1:
                        if exp == 1:
                            repList = [f'{unit}'] + repList
                        else:
                            repList = [f'{unit}^{exp}'] + repList
                    case -1:
                        if ' / ' not in repList:
                            repList += ['/']
                        if exp == -1:
                            repList += [f'{unit}']
                        else:
                            repList += [f'{unit}^{-exp}']
            except ZeroDivisionError:
                pass
        return ' '.join(repList)
    def __neg__(self):
        newUnit = Unit(baseUnitVector=-self.baseUnitVector)
        # newUnit.decompile()
        return newUnit
    def getDistance(self, other):
        diff = self / other
        distance = 0
        for exp in diff.baseUnitVector:
            distance += abs(exp)
        return distance
    def decompile(self):
        rep = {}
        remainingUnit = self
        while True:
            distances = {}
            # breakpoint()
            for unit in units.keys():
                distances[unit] = self.getDistance(units[unit])
                distances[f'-{unit}'] = self.getDistance(-units[unit])
            smallestDistance = 0
            closestUnit = ''
            for unit, distance in zip(distances.keys(), distances.items()):
                distance = distance[-1]
                if closestUnit == '':
                    closestUnit = unit
                    smallestDistance = distance
                elif abs(distance) < abs(smallestDistance):
                    closestUnit = unit
                    smallestDistance = distance
            if closestUnit in rep.keys():
                rep[closestUnit] += 1
            else:
                rep[closestUnit] = 1
            # breakpoint()
            if closestUnit[0] == '-':
                remainingUnit = remainingUnit * units[closestUnit[1:]]
            else:
                remainingUnit = remainingUnit / units[closestUnit]
            if remainingUnit.baseUnitVector.all() == nullUnit.baseUnitVector.all():
                break
        numerator = []
        denominator = []
        for unit, exp in zip(rep.keys(), rep.items()):
            exp = exp[-1]
            try:
                match (exp / abs(exp)):
                    case 1:
                        if exp == 1:
                            numerator += [f'{unit}']
                        else:
                            numerator += [f'{unit}^{exp}']
                    case -1:
                        if exp == -1:
                            denominator += [f'{unit}']
                        else:
                            denominator += [f'{unit}^{-exp}']
            except ZeroDivisionError:
                pass
        notation = ' '.join(numerator)
        if len(denominator) > 0:
            notation += f' / {' '.join(denominator)}'
        self.displayNotation = notation
    def __init__(self, displayNotation = None, baseUnitVector=None, name='', ID=nextUnitID, scalar = 1):
        global units, nextUnitID
        nextUnitID += 1
        self.scalar = scalar
        self.ID = ID
        self.preference = None
        self.simple = False
        self.displayNotation = displayNotation
        self.baseUnitVector = numpy.ndarray((len(baseUnitNames),), fractions.Fraction)
        self.baseUnitVector[:] = [fractions.Fraction(0)] * len(baseUnitNames)
        try:
            self.simple = bool(_.split(' ') for _ in self.displayNotation.split(' / '))
        except AttributeError:
            pass
        if isinstance(baseUnitVector, type(None)):
            if not isinstance(displayNotation, (type(None), Unit,)):
                # breakpoint()
                frac = self.displayNotation.split(' / ')
                numerator = frac[0]
                for dim in numerator.split(' '):
                    if '^' in dim:
                        dim = dim.split('^')
                        exp = int(dim[1])
                        dim = dim[0]
                    guesses = []
                    for unit in units.keys():
                        if unit in dim:
                            guesses += [unit]
                    try:
                        if len(guesses) == 1:
                            unit = guesses[0]
                        else:
                            unit = guesses[-1]
                    except IndexError:
                        pass
                    if unit == dim:
                        pass
                    else:
                        prefix = dim[:-len(unit)]
                        scale = prefixes[prefix]
                        if unit == 'b':
                            scale = scale * 10 / 3
                            scale *= exp
                            scalar *= (1 << scale)
                        else:
                            scale *= exp
                            scalar *= 10 ** scale
                    try:
                        self.baseUnitVector += units[unit].baseUnitVector * exp
                    except UnboundLocalError:
                        self.baseUnitVector += units[unit].baseUnitVector
                try:
                    denominator = frac[1]
                    for dim in denominator.split(' '):
                        if '^' in dim:
                            dim = dim.split('^')
                            exp = dim[1]
                            dim = dim[0]
                        guesses = []
                        for unit in units.keys():
                            if unit in dim:
                                guesses += [unit]
                        if len(guesses) == 1:
                            unit = guesses[0]
                        else:
                            unit = guesses[len(guesses)-1]
                        if unit == dim:
                            pass
                        else:
                            prefix = dim[:-len(unit)]
                            scale = prefixes[prefix]
                            if unit == 'b':
                                scale = scale * 10 / 3
                                scale *= exp
                                scalar /= (1 << scale)
                            else:
                                scale *= exp
                                scalar /= 10 ** scale
                        try:
                            # self.baseUnitVector -= units[unit].baseUnitVector * exp
                            for i in range(abs(int(exp))):
                                self.baseUnitVector -= units[unit].baseUnitVector
                        except UnboundLocalError:
                            self.baseUnitVector -= units[unit].baseUnitVector
                except IndexError:
                    pass
            elif isinstance(displayNotation, Unit):
                self = displayNotation
        elif isinstance(baseUnitVector, numpy.ndarray):
            if baseUnitVector.dtype == fractions.Fraction:
                self.baseUnitVector = baseUnitVector
            elif baseUnitVector.dtype == int:
                for unit, i in zip(baseUnitVector, range(len(baseUnitVector))):
                    self.baseUnitVector[i] = unit
            self.displayNotation = self.getBaseRepresentation()
        else:
            self.baseUnitVector = numpy.ndarray((len(baseUnitNames),), fractions.Fraction)
            for unit, i in zip(baseUnitVector, range(len(baseUnitVector))):
                self.baseUnitVector[i] = fractions.Fraction(unit)
            self.displayNotation = self.getBaseRepresentation()

    def name(self, displayNotation, preference = None):
        try:
            units.pop(self.displayNotation)
        except KeyError:
            pass
        units[displayNotation] = self
        self.displayNotation = displayNotation
        self.simple = len(self.displayNotation.split(' ')) == 1
        self.preference = preference
        if self.simple:
            units[self.displayNotation] = self
        return self
    def getTex(self):
        pass
    def getBaseRepresentation(self):
        repList = []
        for exp, unit in zip(self.baseUnitVector, baseUnits.keys()):
            try:
                match int(exp / abs(exp)):
                    case 0:
                        pass
                    case 1:
                        if exp == 1:
                            repList = [f'{unit}'] + repList
                        else:
                            repList = [f'{unit}^{exp}'] + repList
                    case -1:
                        if ' / ' not in repList:
                            repList += ['/']
                        if exp == -1:
                            repList += [f'{unit}']
                        else:
                            repList += [f'{unit}^{-exp}']
            except ZeroDivisionError:
                pass
        return ' '.join(repList)
    def __eq__(self, other):
        return self.baseUnitVector == other.baseUnitVector
    def __str__(self):
        if isinstance(self.displayNotation, (type(None), Unit)) or self.displayNotation == '':
            pass
        self.decompile()
        return self.displayNotation
    def __repr__(self):
        return f'u<{self}>'
    def __mul__(self, other):
        newUnit = Unit(baseUnitVector = self.baseUnitVector + other.baseUnitVector)
        # newUnit.decompile()
        return newUnit
    def __truediv__(self, other):
        newUnit = Unit(baseUnitVector = self.baseUnitVector - other.baseUnitVector)
        # newUnit.decompile()
        return newUnit
    def __pow__(self, other):
        newUnit = Unit(baseUnitVector = self.baseUnitVector * other)
        # newUnit.decompile()
        return newUnit

class Value(object):
    def __init__(self, scalar, unit):
        self.scalar = scalar
        if isinstance(unit, Unit):
            self.unit = unit
        else:
            self.unit = Unit(unit)
    def __str__(self):
        return f'{self.scalar} {self.unit}'
    def __repr__(self):
        return f'{self.scalar} {self.unit}'
    def __neg__(self):
        return Value(-self.scaler, self.unit)
    def __add__(self, other):
        if self.unit == other.unit:
            return Value(self.scalar + other.scalar, self.unit)
    def __sub__(self, other):
        return self + (-other)
    def __mul__(self, other):
        return Value(self.scalar * other.scalar, self.unit * other.unit)
    def __truediv__(self, other):
        return Value(self.scalar / other.scalar, self.unit / other.unit)
    def __pow__(self, other):
        return Value(self.scalar ** other, self.unit ** other)

nullUnit = Unit('null', [0] * len(baseUnitNames))

for baseUnit, i in zip(baseUnits.keys(), range(len(baseUnits.keys()))):
    baseUnits[baseUnit] = Unit(baseUnitVector=baseUnits[baseUnit]).name(baseUnit)

derivedUnits = {
    "Hz":"s^-1",
    "N":"kg m / s^2",
    "Pa":"N / m^2",
    "J":"N m",
    "W":"J / s",
    "C":"s A",
    "V":"W / A",
    "F":"C / V",
    "Ω":"V / A",
    "S":"A / V",
    "Wb":"V s",
    "T":"Wb / m^2",
    "H":"Wb / A",
    "lm":"cd",
    "lx":"lm / m^2",
    "Bq":"s^-1",
    "Gy":"m^2 / s^2",
    "Sv":"m^2 / s^2",
    "kat":"mol / s",
}

for derivedUnitSymbol, derivation in derivedUnits.items():
    Unit(derivation).name(derivedUnitSymbol)
