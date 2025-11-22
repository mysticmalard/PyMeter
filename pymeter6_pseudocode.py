format_vector = '''
for exp
    exp = exp -> fraction
'''

Unit_init = '''
match type
    case list
        if rep length == 3
            if rep[0] is Unit and rep[1] is str
                repType = child
            else
                repType = vector
                - Format vector
        else
            repType = vector
            - Format vector
    case array
        repType = vector
        - Format vector
    case str
        repType = fDim
'''

disassemble = '''
loop
    until remainder = 0
    do
        # Find closest unit
        for unit in namedUnits
            - Check standard unit
            - Check inverted unit
        - Pull closest unit out
- Combine terms
    - Inference
'''

interlace = '''

'''

distribute = '''
if exp > 0
    if denom?
        if numer?
            for term in numer
                if uexp?
                    uexp = uexp -> fraction
                else
                    uexp = 1
                term = unit ^ uexp * exp
            for term in denom
                if uexp?
                    uexp = uexp -> fraction
                else
                    uexp = 1
                term = unit ^ uexp * exp
            fDim = numer + ' / ' + denom
        else
            pass
    else
        if numer?
            for term in numer
                if uexp?
                    uexp = uexp -> fraction
                else
                    uexp = 1
                term = unit ^ uexp * exp
            fDim = numer
if exp < 0
    if denom?
        if numer?
            pass
        else
            pass
    else
        if numer?
            pass
        else
            pass
if exp == 0
    fDim = ''
return fDim
'''

categorize = '''
'''

assemble = '''
if denom?
    for term in denom
        if exp?
            exp = exp -> fraction
        else
            exp = 1
        vector -= unit * exp
        scalar /= unit * exp
for term in numer
    if exp?
        exp = exp -> fraction
    else
        exp = 1
    vector += unit * exp
    scalar *= unit * exp
'''

Unit_get_attr = '''
match repType
    case vector
        match name
            case vector
                - Trivial
            case fDim
                - Disassemble
            case scalar
                - 1
            case quantity
                - Categorize
    case fDim
        match name
            case vector
                - Assemble
            case fDim
                - Trivial
            case scalar
                - Assemble
            case quantity
                - Categorize
    case child
        match op
            case mul
                match name
                    case vector
                        - Add vectors
                    case fDim
                        - Disassemble
                    case scalar
                        - Multiply scalars
                    case quantity
                        - Compile
            case div
                match name
                    case vector
                        - Subtract vectors
                    case fDim
                        - Disassemble
                    case scalar
                        - Divide scalars
                    case quantity
                        - Compile
            case pow
                match name
                    case vector
                        - Multiply vector
                    case fDim
                        - Disassemble
                    case scalar
                        - Exp scalar
                    case quantity
                        - Compile
'''

Value_init = '''

'''