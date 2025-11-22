compileInference = '''
for example in examples
    for term in numer
        if term not in ordering
            - find adjacent terms
            - set priority to one greater than the previous term
            if previous term and next term adjacent priority
                - increment all units' priorities for next term and all following it
            else
                pass
        else

'''

Unit_init = '''
match type
    case list
        if rep length == 3
            case True
                if rep[0] is Unit and rep[1] is str
                    case True
                        repType = child
                    case False
                        repType = vector
                        - Format vector
            case False
                repType = vector
                - Format vector
    case array
        repType = vector
        if rep.dtype is fraction
            case True
            case False
                - Format vector
    case str
        repType = fDim
'''

inference = '''
loop until done
    - find lowest priority term in numerator and add it to the list
loop until done
    - find lowest priority term in denominator and add it to the list
- combine terms
'''

disassemble = '''
loop
    until
        remainder = 0
    do
        # Find closest unit
        for term in namedUnits
            # Check standard unit
            # Check inverted unit
        - Pull closest unit out of the record
- Combine terms
    - inference
'''

assemble = '''
if denom?
    case True
        for term in denom
            if exp?
                case True
                case False
    case False
for term in numer
    if exp?
        case True
        case False
'''

decompile = '''

'''

categorize = '''
Check fDim
    
'''

Unit_getattr = '''
match repType
    case vector
        match name
            case vector
                - Trivial
            case fDim
                - disassemble
            case scalar
                - 1
            case quantity
                - decompile
    case fDim
        match name
            case vector
                - assemble
            case fDim
                - Trivial
            case scalar
                - assemble
            case quantity
                - decompile
    case child
        match op
            case mul
                match name
                    case vector
                        - Add parent's vectors
                    case fDim
                        - disassemble
                    case scalar
                        - Multiply mother's scalar by father's scalar
                    case quantity
                        - categorize
            case div
                match name
                    case vector
                        - Subtract parent's vectors
                    case fDim
                        - disassemble
                    case scalar
                        - Divide mother's scalar by father's scalar
                    case quantity
                        - categorize
            case pow
                match name
                    case vector
                        - Multiply mother's vector by father
                    case fDim
                        - disassemble
                    case scalar
                        - Raise mother's scalar to the power of the father
                    case quantity
                        - categorize
'''

Value_init = '''

'''

Value_getattr = '''

'''

Value_convert = '''

'''

Epression_init = '''

'''

Tensor_init = '''

'''