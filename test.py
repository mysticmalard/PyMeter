def expression(equation):
    print(equation.__annotations__['a'])
    def wrapper(*args, **kwargs):
        i = 0
        while True:
            try:
                return equation()
            except NameError as e:
                equation.__globals__[e.name] = args[i]
                i+=1
    return wrapper

# def expression(equation):
#     def wrapper(*args):
#         return equation(*args)
#     return wrapper

@expression
def foo(a: int, x: int) -> int:
    return a**2 * x
# equation.__annotations__ => type hints
# args => positional arguments in a tuple
# kwargs nonpositional arguments in a dict


# print(foo(a = 10))

class A:
    def __init__(self):
        pass
    def __getattribute__(self, name):
        pass
    def a(self):
        self = B()
        return self

class B(A):
    pass

C = A().a()
